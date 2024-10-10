module common

 !Module containing all global common areas

 !Import contants, parameters and common arrays:
use constants
use contours
use spectral

!-----------------------------------------------------------------------
 !Define quantities to be preserved between recontouring and evolution:
!-----------------------------------------------------------------------

 !PV anomaly, relative vorticity and PV residual for recontouring:
double precision:: q1(ng,ng),z1(ng,ng)
double precision:: q2(ng,ng),z2(ng,ng)
double precision:: qr(ng,ng,nz)

 !Velocity field and dimensionless height anomaly in physical space:
double precision:: u1(ng,ng),v1(ng,ng),h1(ng,ng)
double precision:: u2(ng,ng),v2(ng,ng),h2(ng,ng)

 !Non-hydrostatic pressure:
double precision:: pn1(ng,ng),pn2(ng,ng)

 !Spectral prognostic fields:
double precision:: qs1(ng,ng),ds1(ng,ng),gs1(ng,ng)
double precision:: qs2(ng,ng),ds2(ng,ng),gs2(ng,ng)

 !Time and twist parameter (used for contour surgery):
double precision:: t,twist

 !Number of time steps between grid and contour saves:
double precision:: ngsave,ncsave

end module
