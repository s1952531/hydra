module common

 !Module containing all global common areas

 !Import contants, parameters and common arrays:
use constants
use contours
use spectral
use force

!-----------------------------------------------------------------------
 !Define quantities to be preserved between recontouring and evolution:
!-----------------------------------------------------------------------

 !PV anomaly (full and residual for recontouring) and relative vorticity:
double precision:: qAnomFull(nLatGridPts,nLongGridPts),qAnomResid(nLatGridPts,nLongGridPts),relVort(nLatGridPts,nLongGridPts)

 !Velocity field and dimensionless height anomaly in physical space:
double precision:: uVel(nLatGridPts,nLongGridPts),vVel(nLatGridPts,nLongGridPts),heightAnom(nLatGridPts,nLongGridPts)

 !Spectral prognostic fields:
double precision:: qSpec(nLatGridPts,nLongGridPts),velDiv(nLatGridPts,nLongGridPts),accelDiv(nLatGridPts,nLongGridPts)

 !Bottom topography at current time (if used):
double precision:: bTopog(nLatGridPts,nLongGridPts)

 !Time and twist variable:
double precision:: t,twist

 !Number of time steps between grid and contour saves:
double precision:: nDTGridSave,nDTContSave

end module
