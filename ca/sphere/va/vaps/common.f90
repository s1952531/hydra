module common

 !Module containing all global common areas

 !Import contants, parameters and common arrays:
use constants
use variables
use contours
use spectral

!-----------------------------------------------------------------------
 !Define quantities to be preserved between recontouring and evolution:
!-----------------------------------------------------------------------

 !PV anomaly (full and residual for recontouring) and relative vorticity:
double precision:: qq(ng,nt),qr(ng,nt),zz(ng,nt)

 !Velocity field and dimensionless height anomaly in physical space:
double precision:: uu(ng,nt),vv(ng,nt),hh(ng,nt)

 !Non-hydrostatic pressure (\bar{P}_n/H in notes):
double precision:: ppn(ng,nt)

 !Spectral prognostic fields:
double precision:: qs(ng,nt),ds(ng,nt),gs(ng,nt)

 !PV jump between contours:
double precision:: qjump

 !Number of time steps between grid and contour saves:
double precision:: ngsave,ncsave

end module
