!#########################################################################
!            The Balanced Spherical Single-Layer Shallow-Water
!               Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!   Adapted from sw code in February 2021 by D G Dritschel @ St Andrews

!    This code simulates the unforced balanced Shallow-Water Equations
!    enforcing first-order delta-gamma balance, in which the first time
!    derivatives of delta and gamma are set to zero.  Here, delta is
!    the velocity divergence and gamma is the acceleration divergence.

!    To initialise, simply create a PV field in qq_init.r8 using one
!    of the data generation routines in the init subdirectory, or use
!    the flow-setup script in the scripts subdirectory (called by the
!    master "hydra" script).

!    For details of CLAM, see Dritschel & Fontane, J. Comput. Phys. 229,
!    pp. 5408--5417 (2010).

!    The full algorithm consists of the following modules:

!    caps.f90      : This source - main program loop, repeats successive 
!                    calls to evolve fields and recontour;
!    parameters.f90: User defined parameters for a simulation;
!    constants.f90 : Fixed constants used throughout the other modules;
!    common.f90    : Common data preserved throughout simulation 
!                    (through recontouring--evolution cycle);
!    spectral.f90  : Fourier transform common storage and routines;
!    contours.f90  : Contour advection common storage and routines;
!    congen.f90    : Source code for contour-to-grid conversion;
!    evolution.f90 : Main time evolution module - advects gridded PV
!                    fields using a PS method along with PV contours.
!----------------------------------------------------------------------------
program caps

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

   !Advect PV until next recontouring or end:
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

! Routine initialises fixed constants and arrays, reads in input data
! files, and opens output data files ready for writing to. 

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
 !Read in gridded PV, qr = (zeta+f)/(1+h), where zeta is the relative
 !vorticity and h is the dimensionless height anomaly:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qr
close(11)

 !Copy into qs for proper start (see subroutine init in evolution.f90):
qs=qr

 !Initialise other fields to zero values until inversion is performed:
ds=zero
gs=zero
hh=zero
uu=zero
vv=zero
 !Inversion (and balance) defines all of these fields from the PV field.

!----------------------------------------------------------------------
 !Initially there are no contours (they are built from the gridded PV):
n=0
npt=0

!--------------------------------------
 !Open all plain text diagnostic files:
open(14,file='contours/complexity.asc',status='replace')
open(15,file='evolution/ecomp.asc',status='replace')
open(16,file='evolution/ro-fr-hm.asc',status='replace')

 !Open files for coarse grid saves:
open(31,file='evolution/qq.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(32,file='evolution/dd.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(33,file='evolution/gg.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(34,file='evolution/hh.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(35,file='evolution/zz.r4',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)

 !Open files for 1d longitudinal spectra (averaged over cos(latitude)):
open(51,file='spectra/zspec.asc',status='replace')
open(52,file='spectra/dspec.asc',status='replace')
open(53,file='spectra/gspec.asc',status='replace')
open(54,file='spectra/hspec.asc',status='replace')

 !Open files for contour writes:
open(80,file='contours/qqsynopsis.asc',status='replace')
open(83,file='contours/qqresi.r4',form='unformatted',access='direct', &
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
write(*,*) 'Evolving PV contours and fields ...'
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
close(16)
close(31)
close(32)
close(33)
close(34)
close(35)
close(51)
close(52)
close(53)
close(54)
close(80)
close(83)

return
end subroutine finalise

 !End main program
end program caps
!=======================================================================
