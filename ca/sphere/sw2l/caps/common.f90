module common

 !Module containing all global common areas

 !Import constants, parameters and common arrays:
use constants
use contours
use spectral
use force

!-----------------------------------------------------------------------
 !Define quantities to be preserved between recontouring and evolution:
!-----------------------------------------------------------------------

 !PV anomaly, relative vorticity and PV residual for recontouring:
double precision:: qAnom1(numLatGridPts,numLongGridPts), z1(numLatGridPts,numLongGridPts)
double precision:: qAnom2(numLatGridPts,numLongGridPts), z2(numLatGridPts,numLongGridPts)
double precision:: qResid(numLatGridPts,numLongGridPts,2), qSPole(2), qNPole(2)

 !Velocity field and dimensionless height anomaly in physical space:
double precision:: u1(numLatGridPts,numLongGridPts), v1(numLatGridPts,numLongGridPts), h1(numLatGridPts,numLongGridPts)
double precision:: u2(numLatGridPts,numLongGridPts), v2(numLatGridPts,numLongGridPts), h2(numLatGridPts,numLongGridPts)

 !Spectral prognostic fields:
double precision:: qSpec1(numLatGridPts,numLongGridPts), deltaSpec1(numLatGridPts,numLongGridPts), gammaSpec1(numLatGridPts,numLongGridPts)
double precision:: qSpec2(numLatGridPts,numLongGridPts), deltaSpec2(numLatGridPts,numLongGridPts), gammaSpec2(numLatGridPts,numLongGridPts)

 !PV forcing (if used):
double precision,allocatable,dimension(:,:):: qForcing

 !Bottom topography and spectral fields needed in evolution (if used):
double precision,allocatable,dimension(:,:):: botTopo, specDeltaTopoSrc, specGammaTopoSrc

 !Time and twist variable:
double precision:: t, twist

 !Number of time steps between grid and contour saves:
double precision:: gridSaveInterval, contourSaveInterval

end module common
