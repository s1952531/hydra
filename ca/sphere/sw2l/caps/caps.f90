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
double precision:: qTot(numLatGridPts,numLongGridPts) 
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
  allocate(qForcing(numLatGridPts,numLongGridPts))
  call generate_forcing(qForcing, rmsPvForcing)
endif

 !If fixed topography is present, allocate memory and read in its form:
if (topoFlag) then
  allocate(botTopo(numLatGridPts,numLongGridPts), specDeltaTopoSrc(numLatGridPts,numLongGridPts), specGammaTopoSrc(numLatGridPts,numLongGridPts))   
  open(11,file = 'bb_init.r8',form = 'unformatted', &
        access = 'direct',status = 'old',recl = 2*nbytes)
  read(11,rec=1) t, botTopo
  close(11)
  !Store its contribution to the gamma and delta tendencies:
  specGammaTopoSrc = csq1*botTopo
  call forfft(numLatGridPts, numLongGridPts, specGammaTopoSrc, trig, factors)
  call laplace(specGammaTopoSrc, specDeltaTopoSrc)
  call revfft(numLatGridPts, numLongGridPts, specDeltaTopoSrc, trig, factors)
endif

!----------------------------------------------------------------------
 !Read in the total PV, (zeta+f)/(1+h), where zeta is the relative
 !vorticity and h is the dimensionless height anomaly in each layer:
open(11, file = 'qq_init.r8', form = 'unformatted', &
      access = 'direct', status = 'old', recl = 2*nbytes)
read(11,rec=1) t, qSpec1
read(11,rec=2) t, qSpec2
close(11)

 !Subtract planetary vorticity (cof), de-alias and convert to
 !semi-spectral space as qSpec1 & qSpec2:
do i = 1,numLongGridPts
  qSpec1(:,i) = qSpec1(:,i)-cof
  qSpec2(:,i) = qSpec2(:,i)-cof
enddo
call dealias(qSpec1)
call dealias(qSpec2)

!----------------------------------------------------------------------
 !Read in gridded divergence, ds:
open(11, file = 'dd_init.r8', form = 'unformatted', &
      access = 'direct', status = 'old', recl = 2*nbytes)
read(11,rec=1) t, deltaSpec1
read(11,rec=2) t, deltaSpec2
close(11)

 !Ensure domain average is zero:
call zeroavg(deltaSpec1)
call zeroavg(deltaSpec2)

 !De-alias and convert to semi-spectral space as ds:
call dealias(deltaSpec1)
call dealias(deltaSpec2)

!----------------------------------------------------------------------
 !Read in gridded acceleration divergence, gs:
open(11, file = 'gg_init.r8', form = 'unformatted', &
      access = 'direct', status = 'old', recl = 2*nbytes)
read(11,rec=1) t, gammaSpec1
read(11,rec=2) t, gammaSpec2
close(11)

 !Ensure domain average is zero:
call zeroavg(gammaSpec1)
call zeroavg(gammaSpec2)

 !De-alias and convert to semi-spectral space as gs:
call dealias(gammaSpec1)
call dealias(gammaSpec2)

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
call main_invert(qSpec1,qSpec2,deltaSpec1,deltaSpec2,gammaSpec1,gammaSpec2,h1,h2,u1,u2,v1,v2,qAnom1,qAnom2,z1,z2)
 !Note: qSpecj, deltaSpecj & gammaSpecj (j = 1 & 2) are in semi-spectral space while 
 !      hj, uj, vj, qAnomj and zj are in physical space.

 !Define qResid to build initial PV contours in recontour below:
do i = 1,numLongGridPts
  qResid(:,i,1) = qAnom1(:,i)+cof
  qResid(:,i,2) = qAnom2(:,i)+cof
enddo

 !Obtain unique polar values of qResid for use in recontouring:
do iz = 1,2
  qSPole(iz) = sum(qResid( 1,:,iz))/dble(numLongGridPts)
  qNPole(iz) = sum(qResid(numLatGridPts,:,iz))/dble(numLongGridPts)
enddo
 !These are averages of points half a latitude spacing from the poles

 !Initially there are no contours (they are built from the total PV qResid):
numContoursTot = 0
numNodesTot = 0
layerLastContIdx = 0 !Used in congen initially; must be zero to start
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
open(31, file = 'evolution/qAnom1.r4', form = 'unformatted', access = 'direct', &
      status = 'replace', recl = nbytes)
open(41, file = 'evolution/qAnom2.r4', form = 'unformatted', access = 'direct', &
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
gridSaveInterval = nint(tGridSaveInterval/dt)
contourSaveInterval = nint(tContourSaveInterval/dt)
 !*** WARNING: tGridSaveInterval and tContourSaveInterval should be an integer multiple of dt

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
write(*,'(a,i8,a,i9)') '   numContoursTot = ',numContoursTot,'   numNodesTot = ',numNodesTot

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
