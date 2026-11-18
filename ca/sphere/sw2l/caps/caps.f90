!#########################################################################
!                The Spherical Two-Layer Shallow-Water
!              Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!      Code developed from plane/sw2l and sphere/sw in November 2026
!      by D G Dritschel @ Kyoto.

!      This code simulates the Shallow-Water Equations (SWE) in variables
!      (q,delta,gamma), where q is the potential vorticity,
!      delta is the velocity divergence, and gamma is the acceleration 
!      divergence (called ageostrophic vorticity).

!      The baroclinic PV can be forced randomly using a Markovian process.

!      Contour advection and generation are done internally now.  For
!      details of the method, see Dritschel & Fontane, J. Comput. Phys.
!      229, pp. 5408--5417 (2010) i.e. CLAM paper.

!      The full algorithm consists of the following modules:

!      caps.f90      : This source - main program loop, repeats successive 
!                      calls to evolve fields and recontour;
!      parameters.f90: User defined parameters for a simulation;
!      constants.f90 : Fixed constants used throughout the other modules;
!      common.f90    : Common data preserved throughout simulation 
!                      (through recontouring--evolution cycle);
!      spectral.f90  : Fourier transform common storage and routines;
!      contours.f90  : Contour advection common storage and routines;
!      congen.f90    : Source code for contour-to-grid conversion;
!      evolution.f90 : Main time evolution module - advects gridded 
!                      fields using a PS method along with contours.
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
integer:: i, iz

!----------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral
 !Initialise constants and arrays for contour advection:
call init_contours

 !If baroclinic PV forcing is present, allocate memory and generate
 !its initial form:
if (forcing) then
  allocate(fq(ng,nt))
  call generate_forcing(fq, fqrms)
endif

 !If fixed topography is present, allocate memory and read in its form:
if (topogr) then
  allocate(bb(ng,nt), tsds(ng,nt), tsgs(ng,nt))   
  open(11,file = 'bb_init.r8',form = 'unformatted', &
        access = 'direct',status = 'old',recl = 2*nbytes)
  read(11,rec=1) t, bb
  close(11)
  !Store its contribution to the gamma and delta tendencies:
  tsgs = csq1*bb
  call forfft(ng, nt, tsgs, trig, factors)
  call laplace(tsgs, tsds)
  call revfft(ng, nt, tsds, trig, factors)
endif

!----------------------------------------------------------------------
 !Read in the total PV, (zeta+f)/(1+h), where zeta is the relative
 !vorticity and h is the dimensionless height anomaly in each layer:
open(11, file = 'qq_init.r8', form = 'unformatted', &
      access = 'direct', status = 'old', recl = 2*nbytes)
read(11,rec=1) t, qs1
read(11,rec=2) t, qs2
close(11)

 !Subtract planetary vorticity (cof), de-alias and convert to
 !semi-spectral space as qs1 & qs2:
do i = 1,nt
  qs1(:,i) = qs1(:,i)-cof
  qs2(:,i) = qs2(:,i)-cof
enddo
call dealias(qs1)
call dealias(qs2)

!----------------------------------------------------------------------
 !Read in gridded divergence, ds:
open(11, file = 'dd_init.r8', form = 'unformatted', &
      access = 'direct', status = 'old', recl = 2*nbytes)
read(11,rec=1) t, ds1
read(11,rec=2) t, ds2
close(11)

 !Ensure domain average is zero:
call zeroavg(ds1)
call zeroavg(ds2)

 !De-alias and convert to semi-spectral space as ds:
call dealias(ds1)
call dealias(ds2)

!----------------------------------------------------------------------
 !Read in gridded acceleration divergence, gs:
open(11, file = 'gg_init.r8', form = 'unformatted', &
      access = 'direct', status = 'old', recl = 2*nbytes)
read(11,rec=1) t, gs1
read(11,rec=2) t, gs2
close(11)

 !Ensure domain average is zero:
call zeroavg(gs1)
call zeroavg(gs2)

 !De-alias and convert to semi-spectral space as gs:
call dealias(gs1)
call dealias(gs2)

!----------------------------------------------------------------------
 !Obtain initial dimensionless height anomaly (h), velocity (u,v),
 !relative vorticity (z) and adjusted PV anomaly (q) consistent
 !with zero domain averaged relative vorticity in each layer:
h1 = zero
h2 = zero
u1 = zero
u2 = zero
v1 = zero
v2 = zero
 !These quantities need initialisation before the first call in order
 !to measure convergence in subroutine main_invert:
call main_invert(qs1,qs2,ds1,ds2,gs1,gs2,h1,h2,u1,u2,v1,v2,q1,q2,z1,z2)
 !Note: qsj, dsj & gsj (j = 1 & 2) are in semi-spectral space while 
 !      hj, uj, vj, qj and zj are in physical space.

 !Define qr to build initial PV contours in recontour below:
do i = 1,nt
  qr(:,i,1) = q1(:,i)+cof
  qr(:,i,2) = q2(:,i)+cof
enddo

 !Obtain unique polar values of qr for use in recontouring:
do iz = 1,2
  qsp(iz) = sum(qr( 1,:,iz))/dble(nt)
  qnp(iz) = sum(qr(ng,:,iz))/dble(nt)
enddo
 !These are averages of points half a latitude spacing from the poles

 !Initially there are no contours (they are built from the total PV qr):
n = 0
npt = 0
jl2 = 0 !Used in congen initially; must be zero to start
 !Note, the PV jump dq is fixed in constants.f90

!--------------------------------------------------------------------
 !Open all plain text diagnostic files:
open(14, file = 'evolution/complexity.asc', status = 'replace')
open(15, file = 'evolution/ecomp.asc', status = 'replace')
open(16, file = 'evolution/rossby-froude.asc', status = 'replace')

 !Open files for 1d longitudinal power spectra for various fields,
 !averaged over cos(latitude) and projected onto each vertical mode:
open(51, file = 'spectra/zspec.asc', status = 'replace')
open(52, file = 'spectra/dspec.asc', status = 'replace')
open(53, file = 'spectra/gspec.asc', status = 'replace')
open(54, file = 'spectra/hspec.asc', status = 'replace')

 !Open files for coarse (inversion) grid saves:
open(31, file = 'evolution/q1.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)
open(41, file = 'evolution/q2.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)
open(32, file = 'evolution/d1.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)
open(42, file = 'evolution/d2.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)
open(33, file = 'evolution/g1.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)
open(43, file = 'evolution/g2.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)
open(34, file = 'evolution/h1.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)
open(44, file = 'evolution/h2.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)
open(35, file = 'evolution/z1.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)
open(45, file = 'evolution/z2.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)

 !Open files for contour writes:
open(80, file = 'contours/qqsynopsis.asc', status = 'replace')
open(83, file = 'contours/qqresi.r4', form = 'unformatted', &
       access = 'direct', status = 'replace', recl = nbytes)

 !Define number of time steps between grid and contour saves:
ngsave = nint(tgsave/dt)
ncsave = nint(tcsave/dt)
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
close(16)
close(31)
close(32)
close(33)
close(34)
close(35)
close(41)
close(42)
close(43)
close(44)
close(45)
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
