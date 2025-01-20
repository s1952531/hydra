!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!            The Combined Lagrangian Advection Method for
!            Multi-Layer Quasi-Geostrophic Flow in a Basin

!    Code written in early-mid 2024 by David Dritschel @ St Andrews
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!        This code solves Dq_l/Dt = F_l, for layers l = 1, ..., nz,
!     where F_l may be wind-stress forcing (l = 1) or Ekman damping
!     (l = nz). Note, l = 1 is the uppermost layer, and l = nz is the
!     lowest layer. Additionally, bathymetry (at z = -H + eta_b) is
!     included by prescribing the scaled bathymetry f*eta_b/H_nz,
!     referred to as qb in the code (H_l = the mean depth of layer l).

!     The rectangular horizontal domain, xmin < x < xmax and
!     ymin < y < ymax, has free-slip boundaries in each layer.

!     Contour advection + a pseudo-spectral scheme is used for the
!     time evolution (see Dritschel & Fontane, J Comput. Phys. 229,
!     5408-5417, 2010).

!     The vertical layer structure is read in from vertical.asc and
!     used to invert PV to find the layerwise-2D velocity field (u,v).
!     Details may be found in spectral.f90.

!     The full algorithm consists of the following modules:
!        caps.f90      : This source - main program loop, repeats successive 
!                        calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        congen.f90    : Source code for contour-to-grid conversion;
!        evolution.f90 : Main time evolution module - advects gridded fields 
!                        using a SL method along with contour advection.
!----------------------------------------------------------------------------
program caps

 !Import common areas:
use common

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t <= tsim)

   !Obtain PV contours:
   call recont

   !Advect PV until next recontouring or end:
   call evolve

enddo

!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Close open files and terminate:
call finalise

!===============================================================

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none

 !Local variables:
double precision:: wkp(0:ny,0:nx),wks(0:nx,0:ny)
integer:: ix,iz

!------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral

 !Initialise constants and arrays for contour advection:
call init_contours

!------------------------------------------------------------------
 !Read in initial PV (qq) as a double-precision field:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qq
close(11)

 !Compute horizontal average in each layer (to be preserved):
do iz=1,nz
   qavg(iz)=sum(qq(:,:,iz)*danorm)
enddo
 !Here, danorm = dx * dy / (L_x * L_y) essentially.

 !Transform qq to spectral space as qs for time stepping:
do iz=1,nz
   wkp=qq(:,:,iz)
   call ptospc_cc(nx,ny,wkp,wks,xfactors,yfactors,xtrig,ytrig)
   qs(:,:,iz)=wks
enddo

!------------------------------------------------------------------
 !Initially there are no contours (they are built from qq):
nq=0
nptq=0

 !Initialise counters used for data saves:
igrids=0
iconts=0
 !These are used in module evolution

!------------------------------------------------------------------
 !Open all diagnostic plain text files:     
open(14,file='evolution/complexity.asc',status='unknown')
open(15,file='evolution/energy.asc',status='unknown')
open(16,file='evolution/monitor.asc',status='unknown')

 !Open files for saving gridded PV, streamfunction & vorticity:
open(31,file='evolution/qq.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(32,file='evolution/pp.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(33,file='evolution/zz.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)

 !Open files for writing PV contour data:
open(80,file='cont/synopsis.asc',status='unknown')
open(83,file='cont/resi.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)

return
end subroutine initialise

!=======================================================================

subroutine evolve

use evolution

implicit none

!------------------------------------------------------------------
 !Advect PV until next recontouring or end:
call advect

return 
end subroutine evolve

!=======================================================================

subroutine recont

use congen

implicit none

!------------------------------------------------------------------
 !Obtain new PV contours:
if (t < small) write(*,*) ' Contouring initial PV field ...'

call recontour(qq)

if (t < small) write(*,'(a,i8,a,i9)') ' Number of contours = ',nq, &
                                      '   Number of nodes = ',nptq

return 
end subroutine recont

!=======================================================================

subroutine finalise

implicit none

!------------------------------------------------------------------
write(*,*)
write(*,*) ' @@@    Code completed normally    @@@'
write(*,*)

 !Close output files (opened in subroutine initialise):
close(14)
close(15)
close(16) 
close(31)
close(32)
close(33)
close(80)
close(83)

return
end subroutine finalise

 !End main program
end program caps
!=======================================================================
