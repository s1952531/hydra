!#########################################################################
!              The Doubly-Periodic Single-Layer Quasi-Geostrophic
!                  Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!        Code written by Stuart King & David Dritschel @ St Andrews

!        Principally adapted from vclam2d.F and ancillary codes:
!                     dgd, 11 November 2010, St Andrews

!          This code simulates the following system of equations:

!             Dq/Dt = S = -r*zeta + (psi-psi_eq)/(tau*L_D^2) + F       (1)
!          L(psi) = (d^2/dx^2 + d^2/dy^2 - 1/L_D^2)psi = q - beta*y    (2)
!                        u = -dpsi/dy ; v = dpsi/dx                    (3)

!          where:
!             r      is the Ekman damping rate
!             tau    is the thermal damping rate
!             L_D    is the Rossby deformation length
!             psi_eq is the thermal equilibrium flow, obtained from a
!                    prescribed PV q_eq by inverting L(psi_eq) = q_eq
!             F      is stochastic forcing (by vortex injection)

!          We split (1) into *three* equations:
!             Dq_a/Dt = 0                          :  contour advection
!             Dq_s/Dt = 0                          :  pseudo-spectral
!             Dq_d/Dt = S + (weak hyperdiffusion)  :  pseudo-spectral
!          such that, in spectral space, q is a weighted sum of these:
!                         q = w*q_s + (1-w)*q_a + q_d
!          where the Butterworth 2nd order filter is used:
!                      w(k) = 1 / (1 + (k/kc)^(2*n))
!              with    n is the order of the filter (default n=2),
!              and     kc is the cut-off wavenumber (default kc=kmax/3).
!          Hence advection at scales k <= kc is controlled by the spectal
!          method (where it is most accurate), whereas advection at scales
!          k >= kc is controlled by contour advection (where it is most
!          accurate).  The source term S is handled entirely by the
!          spectral method, at all k (up to the maximum wavenumber
!          allowed by the grid).

!          At the beginning of each time step, q_s is replaced by q,
!          while q_d is replaced by (1-w)*(q-q_a).  q_a remains as
!          contours for a period of time determined by twistmax below.
!          After this period, q is obtained on an ultra-fine grid and
!          re-contoured so that the accumulated forcing in q_d is
!          given to the contours in q_a (to the extent possible).

!     The full algorithm consists of the following modules:
!        casl.f90      : This source - main program loop, repeats successive 
!                        calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation 
!                        (through recontouring--evolution cycle);
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        generic.f90   : Generic service routines for CASL;
!        congen.f90    : Source code for contour-to-grid conversion;
!        evolution.f90 : Main time evolution module - advects gridded fields 
!                        using SL method along with contours.
!----------------------------------------------------------------------------
program casl

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
 !Offsets used to determine trajectory location in the grid during
 !semi-lagrangian advection:
do ix=1,nx
  xig(ix)=dble(ix-1)
enddo

do iy=1,ny
  yig(iy)=dble(iy-1)
enddo
 
!-----------------------------------------------------------------
 !Read in initial gridded PV:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qs
close(11)

 !Write PV contour interval to log file:
write(*,*)
if (beffect) then
   !Choose an integral number of PV jumps (required in recontour):
  qjump=beta*elly/dble(ncontq)
  write(*,'(a,f14.8,a,f12.8)') 'beta = ',beta,'   qjump = ',qjump
   !Subtract beta*y to define the PV anomaly:
  do ix=1,nx
    do iy=1,ny
      qs(iy,ix)=qs(iy,ix)-beta*(ymin+gly*dble(iy-1))
    enddo
  enddo
else
   !Choose contour interval based on range of PV values:
  qqmax=zero
  qqmin=zero
  do ix=1,nx
    do iy=1,ny
      qqmax=max(qqmax,qs(iy,ix))
      qqmin=min(qqmin,qs(iy,ix))
    enddo
  enddo
  qjump=(qqmax-qqmin)/dble(ncontq)
  write(*,'(a,2(1x,f9.5),1x,f12.8)') ' q_min, q_max, qjump = ', &
                                     & qqmin, qqmax, qjump
endif

 !Copy qs into qd for recontouring:
do ix=1,nx
  do iy=1,ny
    qd(iy,ix)=qs(iy,ix)
  enddo
enddo

!------------------------------------------------------------
 !Initially there are no contours:
nq=0
nptq=0

 !Initialise stochastic vortex forcing variables if used:
if (stoch) then
  if (ivor .eq. 1) then
     !Point vortices of mean (grid) vorticity of +/-vorvor are
     !added at an average rate of dnvor vortices per unit time:
    dnvor=two*esr*(three*pi/vorvor)**2/garea
  else
     !Point vortex dipoles (concentrated to points) are added
     !at an average rate of dnvor dipoles per unit time, with 
     !a maximum absolute grid vorticity of vorvor:
    dnvor=six*esr*(pi/vorvor)**2/garea
  endif
!Above, esr is the enstrophy input rate.

   !Initialize random # generator on first call:
  do i=1,iseed
    uni=rand(0)
  enddo

   !Initialise total added vortices:
  totnvor=zero
endif

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
call recontour(qd)
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
