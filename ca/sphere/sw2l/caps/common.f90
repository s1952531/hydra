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
double precision:: q1(ng,nt), z1(ng,nt)
double precision:: q2(ng,nt), z2(ng,nt)
double precision:: qr(ng,nt,2), qsp(2), qnp(2)

 !Velocity field and dimensionless height anomaly in physical space:
double precision:: u1(ng,nt), v1(ng,nt), h1(ng,nt)
double precision:: u2(ng,nt), v2(ng,nt), h2(ng,nt)

 !Spectral prognostic fields:
double precision:: qs1(ng,nt), ds1(ng,nt), gs1(ng,nt)
double precision:: qs2(ng,nt), ds2(ng,nt), gs2(ng,nt)

 !PV forcing (if used):
double precision,allocatable,dimension(:,:):: fq

 !Bottom topography and spectral fields needed in evolution (if used):
double precision,allocatable,dimension(:,:):: bb, tsds, tsgs

 !Time and twist variable:
double precision:: t, twist

 !Number of time steps between grid and contour saves:
double precision:: ngsave, ncsave

end module common
