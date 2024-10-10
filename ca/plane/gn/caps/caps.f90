!#########################################################################
!                   The Doubly-Periodic Single-Layer 
!       Green-Naghdi Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!  Code developed in late 2017/early 2018 by D G Dritschel @ St Andrews

!       This code simulates the unforced Green-Naghdi Equations (GNE) 
!       in variables (q,delta,gamma), where q is the potential vorticity,
!       delta is the velocity divergence, and gamma is the acceleration 
!       divergence (called ageostrophic vorticity).

!       Contour advection and generation are done internally now.  For
!       details of the method, see Dritschel & Fontane, J. Comput. Phys.
!       229, pp. 5408--5417 (2010).

!       The full algorithm consists of the following modules:

!       caps.f90      : This source - main program loop, repeats successive 
!                       calls to evolve fields and recontour;
!       parameters.f90: User defined parameters for a simulation;
!       constants.f90 : Fixed constants used throughout the other modules;
!       variables.f90 : Global quantities that may change in time;
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

!----------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

!----------------------------------------------------------------------
 !Read in gridded PV anomaly and convert to spectral space as qs:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qr
close(11)
 !Note: qr typically has zero domain average, whereas the actual
 !      PV anomaly may not since this is determined by the 
 !      requirement that the mean relative vorticity is zero;
 !      qr is corrected upon calling main_invert in spectral.f90

 !Convert to spectral space (qr is overwritten; it is recovered below):
call ptospc(ng,ng,qr,qs,xfactors,yfactors,xtrig,ytrig)

 !Ensure domain average qs is zero (this does not matter):
qs(1,1)=zero

!----------------------------------------------------------------------
 !Read in gridded divergence and convert to spectral space as ds:
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,ds,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
ds(1,1)=zero

!----------------------------------------------------------------------
 !Read in gridded acceleration divergence and convert to spectral space
 !as gs:
open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,gs,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
gs(1,1)=zero

!----------------------------------------------------------------------
 !Spectrally-truncate all fields for use in de-aliasing:
qs=filt*qs
ds=filt*ds
gs=filt*gs

 !Obtain initial dimensionless height anomaly (hh), velocity (uu,vv),
 !relative vorticity (zz) and adjusted PV anomaly consistent with 
 !zero domain averaged zz:
call main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
 !Note: qs, ds & gs are in spectral space while 
 !      hh, uu, vv, qq and zz are in physical space.

 !Define qr = qq to build initial PV contours in recontour below:
qr=qq

 !Determine PV contour interval:
qjump=(maxval(qq)-minval(qq))/dble(ncont)

!Initially there are no contours (they are built from the gridded PV):
nq=0
nptq=0

!--------------------------------------
 !Open all plain text diagnostic files:
open(14,file='complexity.asc',status='replace')
open(15,file='ecomp.asc',status='replace')
open(16,file='ubar.asc',status='replace')
open(17,file='monitor.asc',status='replace')

 !Open file for 1d vorticity & divergence spectra:
open(51,file='spectra.asc',status='replace')

 !Open files for coarse grid saves:
open(31,file='qq.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(32,file='dd.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(33,file='gg.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(34,file='hh.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(35,file='zz.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)

 !Open files for contour writes:
open(80,file='cont/qqsynopsis.asc',status='replace')
open(83,file='cont/qqresi.r4',form='unformatted',access='direct', &
                          & status='replace',recl=nbytes)

 !Define number of time steps between grid and contour saves:
ngsave=nint(tgsave/dt)
ncsave=nint(tcsave/dt)
 !*** WARNING: tgsave and tcsave should be an integer multiple of dt

return
end subroutine

!=======================================================================

subroutine evolve

use evolution

implicit none

 !Advect PV until next recontouring or end:
write(*,*) 'Evolving contours and fields ...'
call advect

return 
end subroutine

!=======================================================================

subroutine recont

use congen

implicit none

 !Obtain new PV contours:
write(*,*) 'Recontouring PV ...'
call recontour(qr)
write(*,'(a,i8,a,i9)') '   nq = ',nq,'   nptq = ',nptq

return 
end subroutine

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
close(51)
close(80)
close(83)

return
end subroutine

 !End main program
end program
!=======================================================================
