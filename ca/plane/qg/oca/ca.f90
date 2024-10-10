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
!        ca.f90        : This source - main program loop
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation 
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        generic.f90   : Generic service routines for CASL;
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

   !Advect PV until end time reached:
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
 !Read in initial PV contours:
open(11,file='pvcont.dat',status='old')
read(11,*) nq,nptq,t
do j=1,nq
  read(11,*) npq(j),i1q(j),indq(j)
enddo
do i=1,nptq
  read(11,*) xq(i),yq(i)
enddo
close(11)

 !Set ending contour indices:
do j=1,nq
  i2q(j)=i1q(j)+npq(j)-1
enddo

 !Define nextq(i), which gives the node following node i:
do i=1,nptq-1
  nextq(i)=i+1
enddo
do j=1,nq
  nextq(i2q(j))=i1q(j)
enddo

 !Ensure points all lie inside the domain:
do i=1,npt
  xq(i)=(xq(i)-ellx*int(xq(i)*hlxi))*oms
  yq(i)=(yq(i)-elly*int(yq(i)*hlyi))*oms
enddo

!------------------------------------------------------------
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
open(17,file='steady.asc',status='unknown')
 !Open file for 1d PV spectrum:
open(50,file='qspec.asc',status='unknown')
 !Open files for coarse grid saves:
open(31,file='qq.r4',form='unformatted', &
    access='direct',status='replace',recl=nbytes)
 !Open files for contour writes:
open(80,file='cont/qqsynopsis.asc',status='unknown')

 !Initialise counter for writing direct files to the correct counter:
igrids=0

return
end subroutine

!=======================================================================

subroutine evolve

use evolution

implicit none

!Advect PV until end:
write(*,*) 'Evolving contours ...'
call advect

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
close(50)
close(80)

return
end subroutine

 !End main program
end program
!=======================================================================
