!#########################################################################
!     The Doubly-Periodic Single-Layer Magnetised Vertically-Averaged 
!               Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!  Code developed in late 2019 by D G Dritschel & M R Jalali @ St Andrews
!  Adapted from va and msw directoris in late August 2022 by DGD

!       This code simulates the unforced Vertically-Averaged (VA) equations
!       for a magnetised fluid in variables (q_l,delta,gamma_l,b_x,b_y),
!       where q_l is the linearised potential vorticity (PV), delta is the
!       velocity divergence, gamma_l is the linearised acceleration
!       divergence f*zeta-g*lap(h), while b_x and b_y are the horizontal
!       components of the magnetic field. See Dritschel & Jalali, JFM (2020)
!       and Dritschel & Tobias, JFM (2023).

!       Contour advection and generation are done internally now.  For
!       details of the method, see Dritschel & Fontane, J. Comput. Phys.
!       229, pp. 5408--5417 (2010).

!       The full algorithm consists of the following modules:

!       caps.f90      : This source - main program loop, repeats successive 
!                       calls to evolve fields and recontour;
!       parameters.f90: User defined parameters for a simulation;
!       constants.f90 : Fixed constants used throughout the other modules;
!       common.f90    : Common data preserved throughout simulation 
!                       (through recontouring--evolution cycle);
!       spectral.f90  : Fourier transform common storage and routines;
!       contours.f90  : Contour advection common storage and routines;
!       congen.f90    : Source code for contour-to-grid conversion;
!       evolution.f90 : Main time evolution module - advects gridded 
!                       fields using a PS method along with contours.
!----------------------------------------------------------------------------
program casl

use common

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Obtain new PV contours:
  call recont

   !Advect PV and other fields until next recontouring or end:
  call evolve

enddo

 !End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

call finalise

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none
double precision:: qmin,ql1,ql2
integer:: ix,iy

!----------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

!----------------------------------------------------------------------
 !Read in gridded (linearised) PV and convert to spectral space as qs:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,zz
close(11)
qr=zz !Preserve physical space field in qr for contouring below

call ptospc(ng,ng,zz,qs,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
qs(1,1)=zero

!----------------------------------------------------------------------
 !Read in gridded divergence and convert to spectral space as ds:
open(11,file='dd_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,ds,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
ds(1,1)=zero

!----------------------------------------------------------------------
 !Read in gridded acceleration divergence and convert to spectral space
 !as gs:
open(11,file='gg_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,gs,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
gs(1,1)=zero

!----------------------------------------------------------------------
 !Read in x component of B and convert to spectral space as bxs:
open(11,file='bx_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,bxs,xfactors,yfactors,xtrig,ytrig)
 !Domain average can be non-zero.

!----------------------------------------------------------------------
 !Read in y component of B and convert to spectral space as bys:
open(11,file='by_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,bys,xfactors,yfactors,xtrig,ytrig)
 !Domain average can be non-zero.

!----------------------------------------------------------------------
 !Spectrally-truncate all fields for use in de-aliasing:
qs=filt*qs
ds=filt*ds
gs=filt*gs
bxs=filt*bxs
bys=filt*bys

 !Obtain initial dimensionless height anomaly (hh), velocity (uu,vv),
 !relative vorticity (zz) by inversion:
call main_invert(qs,ds,gs,hh,uu,vv,zz)
 !Note: qs, ds & gs are in spectral space while 
 !      hh, uu, vv and zz are in physical space.

 !Choose PV contour interval based on <q^2>/<|q|> for |q| > q_rms:
qmin=sqrt(dsumi*sum(qr**2))
ql1=zero
ql2=zero
do ix=1,ng
  do iy=1,ng
    if (abs(qr(iy,ix)) .gt. qmin) then
      ql1=ql1+abs(qr(iy,ix))
      ql2=ql2+qr(iy,ix)**2
    endif
  enddo
enddo
qjump=ql2/(ql1*dble(ncont))

!Initially there are no contours (they are built from the gridded PV):
nq=0
nptq=0

!--------------------------------------
 !Open all plain text diagnostic files:
open(14,file='contours/complexity.asc',status='replace')
open(15,file='evolution/ecomp.asc',status='replace')
open(16,file='evolution/ro-fr-hm.asc',status='replace')
open(17,file='evolution/norms.asc',status='replace')

 !Open files for coarse grid saves:
open(31,file='evolution/qq.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(32,file='evolution/dd.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(33,file='evolution/gg.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(34,file='evolution/hh.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(35,file='evolution/zz.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(36,file='evolution/pn.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(37,file='evolution/jz.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(38,file='evolution/tt.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)

 !Open files for 1d spectra:
open(51,file='spectra/zspec.asc',status='replace')
open(52,file='spectra/dspec.asc',status='replace')
open(53,file='spectra/gspec.asc',status='replace')
open(54,file='spectra/hspec.asc',status='replace')
open(55,file='spectra/pspec.asc',status='replace')
open(56,file='spectra/jspec.asc',status='replace')

 !Open files for contour writes:
open(80,file='contours/qqsynopsis.asc',status='replace')
open(83,file='contours/qqresi.r8',form='unformatted',access='direct', &
                                status='replace',recl=nbytes)

return
end subroutine initialise

!=======================================================================

subroutine evolve

use evolution

implicit none

 !Advect PV until next recontouring or end:
write(*,*) 'Evolving contours and fields ...'
call advect

return 
end subroutine evolve

!=======================================================================

subroutine recont

use congen

implicit none

 !Obtain new PV contours:
write(*,*) 'Recontouring PV ...'
call recontour(qr)
write(*,'(a,i8,a,i9)') '   nq = ',nq,'   nptq = ',nptq

return 
end subroutine recont

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(14)
close(15)
close(16)
close(17)
close(31)
close(32)
close(33)
close(34)
close(35)
close(36)
close(37)
close(38)
close(51)
close(52)
close(53)
close(54)
close(55)
close(56)
close(80)
close(83)

return
end subroutine finalise

 !End main program
end program casl
!=======================================================================
