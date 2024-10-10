module common

 !Module containing all global common areas

 !Import contants, parameters and common arrays:
use constants
use contours
use spectral
use generic

!-----------------------------------------------------------------------
 !Define quantities to be preserved between recontouring and evolution:
!-----------------------------------------------------------------------

 !PV fields used in CLAM and recontouring:
double precision:: qq(0:ny,0:nxm1),qs(0:ny,0:nxm1),qr(0:ny,0:nxm1)

 !Momentum divergence div((1+h)*(u,v)) and acceleration divergence:
double precision:: cc(0:ny,0:nxm1),gg(0:ny,0:nxm1)

 !Velocity field and dimensionless height anomaly:
double precision:: uu(0:ny,0:nxm1),vv(0:ny,0:nxm1),hh(0:ny,0:nxm1)

 !Zonal velocities at y = y_min and y_max:
double precision:: uum(0:nxm1),uup(0:nxm1)

 !Mean values of the boundary zonal velocities (remain constant):
double precision:: uumbar,uupbar,qoff

 !Time and twist parameter (used for contour surgery):
double precision:: t,twist

 !Number of time steps between grid and contour saves:
double precision:: ngsave,ncsave

end module common
