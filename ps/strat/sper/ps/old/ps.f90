!#########################################################################
!                      The pseudo-spectral method for 
!      2D non-rotating Boussinesq flow in a periodic channel geometry
!#########################################################################

!        Code adapted from strat/caps by David Dritschel @ St Andrews
!                ***Version 1.0 completed 11 January 2022***

!     This code solves: 
!            zeta_t + (u,v)*grad(zeta) = db/dx + D[zeta]
!               b_t + (u,v)*grad(b) = D[b]
!               u_x + v_y = 0
!     in the domain xmin < x < xmax ; ymin < y < ymax
!     (free slip boundary conditions in y, periodic in x).
!     Here D is a hyperviscous dissipation operator.

!     Incompressibility is handled via the introduction of a stream-
!     function such that:
!             Lap(psi) = zeta

!     The full algorithm consists of the following modules:
!        ps.f90        : This source - main program loop
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        common.f90    : Common data preserved throughout simulation 
!        spectral.f90  : Fourier transform common storage and routines;
!        evolution.f90 : Main time evolution module - advects gridded fields 
!----------------------------------------------------------------------------
program ps

use common
use evolution

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

 !Advect buoyancy & vorticity until end of simulation:
call advect

 !Close all files and finish:
call finalise

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none

!-----------------------------------------------------------------
 !Read in initial vorticity (zz) and buoyancy (bb):
open(11,file='zz_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)
open(12,file='bb_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(12,rec=1) t,bb
close(12)

!--------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!--------------------------------------
 !Open all plain text diagnostic files:
open(12,file='monitor.asc',status='unknown')
open(13,file='norms.asc',status='unknown')
open(15,file='ene.asc',status='unknown')
 !Open file for 1d vorticity spectrum:
open(50,file='zspec.asc',status='unknown')
 !Open files for coarse grid saves:
open(31,file='zz.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(32,file='bb.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

!------------------------------------------------------------
 !Set flag for computing PE:
iene=0

 !Time to save next data:
tgrid=zero

 !Initialise counter for writing direct files to the correct counter:
igrids=0

return
end subroutine initialise

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(12) 
close(13)
close(15)
close(31)
close(32)
close(50)

return
end subroutine finalise

 !End main program
end program ps
!=======================================================================
