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

 !Prognostic fields for b = Dv/Dt and q (semi-spectral):
double precision:: bs(0:ny,0:nxm1),qs(0:ny,0:nxm1)

 !PV field and residue used for recontouring (physical):
double precision:: qq(0:ny,0:nxm1),qr(0:ny,0:nxm1)

 !Dimensionless height anomaly and relative vorticity (physical):
double precision:: hh(0:ny,0:nxm1),zz(0:ny,0:nxm1)

 !Velocity field (physical):
double precision:: uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)

 !Linearised PV anomaly and balanced meridional velocity (physical):
double precision:: qt(0:ny,0:nxm1),vb(0:ny,0:nxm1)

 !Zonal velocities at y = y_min and y_max:
double precision:: uum(0:nxm1),uup(0:nxm1)

 !Mean values of the boundary zonal velocities (remain constant):
double precision:: uumbar,uupbar,qoff

 !Time and twist parameter (used for contour surgery):
double precision:: t,twist

 !Number of time steps between grid and contour saves:
double precision:: ngsave,ncsave

end module common
