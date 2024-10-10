!#########################################################################
!                    The Doubly-Periodic Single-Layer 
!        Quasi-Geostrophic Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!        Code adapted by D G Dritschel from the imhd code on 17/1/2014
!        Revised by dgd on 15 Jan 2024 to include -beta*v in PV anomaly
!        source, rather than use contours for the total PV (difficult
!        in a y periodic domain, especially for small beta).

!          This code evolves the QG PV anomaly q from

!             Dq/Dt = -r*zeta + (psi-psi_eq)/(tau*L_D^2) - beta*v = S

!          where
!             r      is the Ekman damping rate,
!             tau    is the thermal damping rate,
!             L_D    is the Rossby deformation length,
!             zeta   is the relative (vertical) vorticity,
!             psi    is the streamfunction,
!             psi_eq is the thermal equilibrium streamfunction, and
!             beta   is the planetary vorticity gradient.

!          The 2D velocity field (u,v) is found by inverting 

!             Lap{psi} - psi/L_D^2 = q ; u = -dpsi/dy ; v = dpsi/dx.

!          We split the PV evolution equation into *three* equations,

!             Dq_c/Dt = 0  :  contour advection
!             Dq_s/Dt = 0  :  pseudo-spectral
!             Dq_d/Dt = S  :  pseudo-spectral

!          such that q is a weighted sum of these,

!             q = F(q_s) + (1-F)*(q_c) + q_d

!          where F is a low-pass filter defined in spectral.f90 and 1-F is
!          a complementary high pass filter (see Dritschel & Fontane, JCP,
!          2010).

!          Hence advection at large to intermediate scales is controlled 
!          by the pseudo-spectral method, whereas advection at intermediate
!          to small scales is controlled by contour advection (where it 
!          is most accurate).  The source term S is handled entirely by
!          the pseudo-spectral method.

!          At the beginning of each time step, q_s is replaced by q,
!          while q_d is replaced by (1-F)(q-q_c).  q_c remains as
!          contours for a period of time determined by twistmax below.
!          After this period, q is obtained on an ultra-fine grid and
!          re-contoured so that the accumulated forcing in q_d is
!          given to the contours in q_c (to the extent possible).

!     The full algorithm consists of the following modules:
!        caps.f90      : This source - main program loop, repeats successive 
!                        calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation 
!                        (through recontouring--evolution cycle);
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        generic.f90   : Generic service routines for CAPS;
!        congen.f90    : Source code for contour-to-grid conversion;
!        evolution.f90 : Main time evolution module - advects gridded 
!                        fields using a PS method along with contours.
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

   !Advect PV anomaly until next recontouring or end:
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

implicit none

 !Local variables:
double precision:: ff(ny,nx)

!--------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

!--------------------------------------------------------------------
 !Read in thermal equilibrium streamfunction (if present) and
 !convert to spectral coefficients (ppeq):
if (heating) then
  open(12,file='psieq.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(12,rec=1) t,ff
  close(12)
   !Allocate memory for thermal equilibrium streamfunction if heating
  allocate(ppeq(nx,ny))
  call ptospc(nx,ny,ff,ppeq,xfactors,yfactors,xtrig,ytrig)
endif

!--------------------------------------------------------------------
 !Read in PV anomaly and FFT to initialise qs:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qr
close(11)

 !Choose contour interval based on range of PV values:
call contint(qr,ncontq,qjump)
write(*,*)
write(*,'(a,1x,f13.8)') ' qjump = ',qjump

 !Copy PV anomaly qr into ff before taking Fourier transform 
 !(ff is overwritten):
ff=qr
 !Convert PV anomaly to spectral space as qs:
call ptospc(nx,ny,ff,qs,xfactors,yfactors,xtrig,ytrig)

 !Initially there are no contours:
nq=0
nptq=0

!--------------------------------------------------------------------
 !Open all plain text diagnostic files:
open(14,file='complexity.asc',status='replace')
open(15,file='ene.asc',status='replace')
open(16,file='norms.asc',status='replace')
open(17,file='monitor.asc',status='replace')

 !Open file for 1d PV spectra:
open(51,file='spectra.asc',status='replace')

 !Open file for coarse grid saves:
open(31,file='qq.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)

 !Open files for contour writes:
open(80,file='cont/qqsynopsis.asc',status='replace')
open(83,file='cont/qqresi.r4',form='unformatted',access='direct', &
                            status='replace',recl=nbytes)

 !Initialise counter for writing direct access grid files:
igrec=0
 !Initialise counter for writing direct access contour files:
icrec=1

return
end subroutine initialise

!=======================================================================

subroutine evolve

use evolution

implicit none

 !Advect PV anomaly until next recontouring or end:
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

if (forcing) call contint(qr,ncontq,qjump)

call recontour(qr)
write(*,'(a,i9,a,i10,a,f9.5)') '   nq = ',nq,'   nptq = ',nptq,'   dq = ',qjump

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
close(51)
close(80)
close(83)

return
end subroutine finalise

 !End main program
end program caps
!=======================================================================
