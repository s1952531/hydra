!#########################################################################
!           The Doubly-Periodic Two-Layer Vertically-Averaged 
!              Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!     This code simulates the unforced Vertically-Averaged (VA) equations
!     in variables (q,delta,gamma_l), where q is the potential vorticity,
!     delta is the velocity divergence, and gamma_l is the linearised 
!     acceleration divergence f*zeta-lap(p_h) with p_h the hydrostatic
!     pressure in the given layer (see Dritschel & Jalali, JFM, 2021).

!     Contour advection and generation are done internally now.  For
!     details of the method, see Dritschel & Fontane, J. Comput. Phys.
!     229, pp. 5408--5417 (2010).

!     The full algorithm consists of the following modules:

!     caps.f90      : This source - main program loop, repeats successive 
!                     calls to evolve fields and recontour;
!     parameters.f90: User defined parameters for a simulation;
!     constants.f90 : Fixed constants used throughout the other modules;
!     common.f90    : Common data preserved throughout simulation 
!                     (through recontouring--evolution cycle);
!     spectral.f90  : Fourier transform common storage and routines;
!     contours.f90  : Contour advection common storage and routines;
!     congen.f90    : Source code for contour-to-grid conversion;
!     evolution.f90 : Main time evolution module - advects gridded 
!                     fields using a PS method along with contours.

!     Code completed in December 2020 by D G Dritschel @ St Andrews
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
 !Read in gridded PV anomaly and convert to spectral space as qs1,2:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,z1
read(11,rec=2) t,z2
close(11)
 !Note: The input fields typically have zero domain average, whereas
 !      the actual PV anomalies may not since this is determined by the 
 !      requirement that the mean relative vorticities are zero; the
 !      correct means are found upon calling main_invert in spectral.f90

 !Convert to spectral space (qr is overwritten; it is recovered below):
call ptospc(ng,ng,z1,qs1,xfactors,yfactors,xtrig,ytrig)
call ptospc(ng,ng,z2,qs2,xfactors,yfactors,xtrig,ytrig)

 !Ensure domain average qsj is zero (this does not matter):
qs1(1,1)=zero
qs2(1,1)=zero

!----------------------------------------------------------------------
 !Read in gridded divergence and convert to spectral space as dsj:
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,z1
read(11,rec=2) t,z2
close(11)

call ptospc(ng,ng,z1,ds1,xfactors,yfactors,xtrig,ytrig)
call ptospc(ng,ng,z2,ds2,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
ds1(1,1)=zero
ds2(1,1)=zero

!----------------------------------------------------------------------
 !Read in gridded linearised acceleration divergence and convert to 
 !spectral space as gsj:
open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,z1
read(11,rec=2) t,z2
close(11)

call ptospc(ng,ng,z1,gs1,xfactors,yfactors,xtrig,ytrig)
call ptospc(ng,ng,z2,gs2,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
gs1(1,1)=zero
gs2(1,1)=zero

!----------------------------------------------------------------------
 !Spectrally-truncate all fields for use in de-aliasing:
qs1=filt*qs1
qs2=filt*qs2
ds1=filt*ds1
ds2=filt*ds2
gs1=filt*gs1
gs2=filt*gs2

 !Obtain initial dimensionless height anomaly (hj), velocity (uj,vj),
 !relative vorticity (zj) and adjusted PV anomaly consistent with 
 !zero domain averaged zj (j = 1 & 2):
h1=zero
h2=zero
call main_invert(qs1,qs2,ds1,ds2,gs1,gs2,h1,h2,u1,u2,v1,v2,q1,q2,z1,z2)
 !Note: qsj, dsj & gsj (j = 1 & 2) are in spectral space while 
 !      hj, uj, vj, qj and zj are in physical space.

 !Define qr to build initial PV contours in recontour below:
qr(:,:,1)=q1
qr(:,:,2)=q2

 !Determine PV contour interval:
qjump=(max(maxval(q1),maxval(q2))-min(minval(q1),minval(q2)))/dble(ncont)

 !Initially there are no contours (they are built from the gridded PV):
nq=0
nptq=0

!--------------------------------------
 !Open all plain text diagnostic files:
open(14,file='complexity.asc',status='replace')
open(15,file='ecomp.asc',status='replace')
open(16,file='ubar.asc',status='replace')
open(17,file='monitor.asc',status='replace')

 !Open files for 1d spectra projected onto each vertical mode:
open(51,file='zspec.asc',status='replace')
open(52,file='dspec.asc',status='replace')
open(53,file='gspec.asc',status='replace')
open(54,file='hspec.asc',status='replace')
open(55,file='pspec.asc',status='replace')

 !Open files for coarse grid saves:
open(31,file='q1.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(41,file='q2.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(32,file='d1.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(42,file='d2.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(33,file='g1.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(43,file='g2.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(34,file='h1.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(44,file='h2.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(35,file='z1.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(45,file='z2.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(36,file='p1.r4',form='unformatted',access='direct', &
                   status='replace',recl=nbytes)
open(46,file='p2.r4',form='unformatted',access='direct', &
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
close(36)
close(41)
close(42)
close(43)
close(44)
close(45)
close(46)
close(51)
close(52)
close(53)
close(54)
close(55)
close(80)
close(83)

return
end subroutine

 !End main program
end program
!=======================================================================
