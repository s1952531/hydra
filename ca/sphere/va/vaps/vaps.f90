!#########################################################################
!              The Spherical Single-Layer Vertically-Averaged
!               Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!       Code developed from SW one in June 2020 by D G Dritschel @ St Andrews.

!       This code simulates the unforced Vertically-Averaged or Green-
!       Naghdi Equations (GNE) in variables (q,delta,gamma), where q is
!       the potential vorticity, delta is the velocity divergence, and
!       gamma is the linear part of the acceleration divergence.

!       Contour advection and generation are done internally now.  For
!       details of the method, see Dritschel & Fontane, J. Comput. Phys.
!       229, pp. 5408--5417 (2010).

!       The full algorithm consists of the following modules:

!       vaps.f90      : This source - main program loop, repeats successive 
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
program clam

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

! Local variables:
double precision:: qt(ng,nt)
integer:: i  

!----------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

!----------------------------------------------------------------------
 !Read in gridded PV, qr = (zeta+f)/(1+h)+(H^2/3)*J(delta,h), where
 !zeta is the relative vorticity, and h is the dimensionless height
 !anomaly and delta is the (velocity) divergence:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qr
close(11)

 !Copy into qs for proper start (see subroutine init in evolution.f90):
qs=qr

!----------------------------------------------------------------------
 !Read in gridded divergence, ds:
open(11,file='dd_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,ds
close(11)

 !Ensure domain average is zero:
call zeroavg(ds)

 !De-alias and convert to semi-spectral space as ds:
call dealias(ds)

!----------------------------------------------------------------------
 !Read in gridded acceleration divergence, gs:
open(11,file='gg_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,gs
close(11)

 !Ensure domain average is zero:
call zeroavg(gs)

 !De-alias and convert to semi-spectral space as gs:
call dealias(gs)

!----------------------------------------------------------------------
 !Obtain initial dimensionless height anomaly (hh), velocity (uu,vv),
 !relative vorticity (zz) and adjusted PV anomaly consistent with 
 !zero domain averaged zz:

 !Define PV anomaly (qt) needed for inversion below:
do i=1,nt
  qt(:,i)=qs(:,i)-cof
enddo

 !Convert qt to semi-spectral space:
call forfft(ng,nt,qt,trig,factors) 

call main_invert(qt,ds,gs,hh,uu,vv,qq,zz)
 !Note: qt, ds & gs are in semi-spectral space while 
 !      hh, uu, vv, qq and zz are in physical space.

!----------------------------------------------------------------------
 !Determine PV contour interval from range of Coriolis frequency:
qjump=two*fpole/dble(ncont)

!Initially there are no contours (they are built from the gridded PV):
n=0
npt=0

!--------------------------------------
 !Open all plain text diagnostic files:
open(14,file='complexity.asc',status='replace')
open(15,file='ecomp.asc',status='replace')
open(17,file='monitor.asc',status='replace')

 !Open file for longitudinal vorticity & divergence spectra
 !(averaged over z = sin(latitude)):
open(51,file='spectra.asc',status='replace')

 !Open files for coarse grid saves:
open(31,file='qq.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(32,file='dd.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(33,file='gg.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(34,file='hh.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(35,file='zz.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(36,file='pp.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)

 !Open files for contour writes:
open(80,file='cont/qqsynopsis.asc',status='replace')
open(83,file='cont/qqresi.r4',form='unformatted',access='direct', &
                            status='replace',recl=nbytes)

 !Define number of time steps between grid and contour saves:
ngsave=nint(tgsave/dt)
ncsave=nint(tcsave/dt)
 !*** WARNING: tgsave and tcsave should be an integer multiple of dt

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
call recontour
write(*,'(a,i8,a,i9)') '   n = ',n,'   npt = ',npt

return 
end subroutine recont

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(14)
close(15)
close(17)
close(31)
close(32)
close(33)
close(34)
close(35)
close(36)
close(51)
close(80)
close(83)

return
end subroutine finalise

 !End main program
end program clam
!=======================================================================
