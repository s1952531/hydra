!#########################################################################
!                 The Doubly-Periodic Single-Layer Surface
!                 Quasi-Geostrophic Pseudo-Spectral Method
!#########################################################################

!   Code adapted by dgd from the CAPS code on 23 July 2025 @ Wildwood Crest

!          This code simulates the following system of equations:

!             Dq/Dt = 0

!          where q = -b_0/N; here b_0 is the surface buoyancy and
!          N is the uniform buoyancy frequency.

!          The velocity field (u,v) is found by 

!             u = -dpsi/dy ; v = dpsi/dx                   

!          where in spectral space psi_hat = -q_hat*cosh(K*D)/(K*sinh(KD))
!          where K is the wavenumber magnitude and D = NH/f is the scaled
!          depth.  Note, 1/D is specified in parameters.f90 so that the
!          conventional limit D -> infinity can be studied.

!          The fields of b_0/N and zeta are output to bb.r4 & zz.r4.

!     The full algorithm consists of the following modules:
!        ps.f90        : This source - main program loop;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation ;
!        spectral.f90  : Fourier transform common storage and routines;
!        evolution.f90 : Main time evolution module - advects gridded 
!                        fields using the pseudo-spectral method.
!----------------------------------------------------------------------------
program ps

use common

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

 !Evolve flow until end:
call evolve

 !Close all files and stop:
call finalise

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none

 !Local variables:
double precision:: ff(ny,nx)
integer:: itime

!--------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral

!-----------------------------------------------------------------
 !Read in full buoyancy and convert to spectral space as qs:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,ff
close(11)

call ptospc(nx,ny,ff,qs,xfactors,yfactors,xtrig,ytrig)

!------------------------------------------------------------
 !Initialise time step so that subroutine adapt chooses a suitable one:
dt=zero

 !Set final time for simulation end:
itime=int((t+small)/tgsave)
tgrid=tgsave*dble(itime)
tfin=tgrid+tsim

!--------------------------------------
 !Open all plain text diagnostic files:
open(15,file='ene-ens.asc',status='unknown')
open(17,file='monitor.asc',status='unknown')

 !Open files for 1d spectra:
open(51,file='spectra/bspec.asc',status='replace')
open(52,file='spectra/zspec.asc',status='replace')

 !Open files for coarse grid saves of b_0/N and zeta:
open(31,file='bb.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(32,file='zz.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)

 !Initialise counter for writing direct files to the correct counter:
igrids=0

return
end subroutine initialise

!=======================================================================

subroutine evolve
!Advects buoyancy until end:

use evolution

implicit none

write(*,*) 'Evolving buoyancy ...'
call advect

return 
end subroutine evolve

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(15)
close(17)
close(31)
close(32)
close(51)
close(52)

return
end subroutine finalise

 !End main program
end program ps
!=======================================================================
