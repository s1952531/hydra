!#########################################################################
!              The Doubly-Periodic Single-Layer Quasi-Geostrophic
!                  Contour-Advective Semi-Lagrangian (CASL) algorithm
!#########################################################################

!        Code written by Stuart King & David Dritschel @ St Andrews

!        Principally adapted from vclam2d.F and ancillary codes:
!                     dgd, 11 November 2010, St Andrews

!          This code simulates the following system of equations:

!                                  Dq/Dt = 0                           (1)
!          L(psi) = (d^2/dx^2 + d^2/dy^2 - 1/L_D^2)psi = q - beta*y    (2)
!                        u = -dpsi/dy ; v = dpsi/dx                    (3)

!          where:
!             L_D    is the Rossby deformation length 
!             beta   is the planetary vorticity gradient

!     The full algorithm consists of the following modules:
!        ca.f90        : This source - main program loop, repeats successive 
!                        calls to evolve contours and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation 
!                        (through recontouring--evolution cycle);
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        generic.f90   : Generic service routines for CASL;
!        congen.f90    : Source code for contour-to-grid conversion;
!        evolution.f90 : Main time evolution module - advects contours.
!----------------------------------------------------------------------------
program ca

use common

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tfin)

   !Obtain new PV contours:
  call recont

   !Advect PV until next recontouring or end:
  call evolve

enddo

!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

call finalise

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit double precision(a-h,o-z)
implicit integer(i-n)

!-----------------------------------------------------------------
 !Read in initial gridded PV:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qq
close(11)

 !Compute contour interval:
qqmax=zero
qqmin=zero
do ix=1,nx
  do iy=1,ny
    qqmax=max(qqmax,qq(iy,ix))
    qqmin=min(qqmin,qq(iy,ix))
  enddo
enddo
qjump=(qqmax-qqmin)/dble(ncontq)

 !Write information to log file:
write(*,*)
write(*,'(a,3(1x,f9.5))') ' q_min, q_max, qjump = ',qqmin,qqmax,qjump

!------------------------------------------------------------
 !Initially there are no contours:
nq=0
nptq=0

 !Initialise time step so that subroutine adapt chooses a suitable one:
dt=zero

 !Set final time for simulation end:
itime=int((t+small)/tgsave)
tgrid=tgsave*dble(itime)
tfin=tgrid+tsim

!--------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

!--------------------------------------
 !Open all plain text diagnostic files:
open(14,file='complexity.asc',status='unknown')
open(15,file='ene.asc',status='unknown')
open(16,file='monitor.asc',status='unknown')
 !Open file for 1d PV spectrum:
open(50,file='qspec.asc',status='unknown')
 !Open files for coarse grid saves:
open(31,file='qq.r4',form='unformatted', &
    access='direct',status='replace',recl=nbytes)
 !Open files for contour writes:
open(80,file='cont/qqsynopsis.asc',status='unknown')
open(83,file='cont/qqresi.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

 !Initialise counter for writing direct files to the correct counter:
igrids=0

return
end subroutine

!=======================================================================

subroutine evolve

use evolution

implicit none

!Advect PV until next recontouring or end:
write(*,*) 'Evolving contours ...'
call advect

return 
end subroutine

!=======================================================================

subroutine recont

use congen

implicit none

!Obtain new PV contours:
write(*,*) 'Recontouring PV ...'
call recontour(qq)
write(*,'(a,i8,a,i9,a,f9.5)') '   nq = ',nq,'   nptq = ',nptq,'   dq = ',qjump

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
close(31)
close(50)
close(80)
close(83)

return
end subroutine


 !End main program
end program
!=======================================================================
