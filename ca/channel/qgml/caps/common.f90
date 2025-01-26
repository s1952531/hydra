module common

 !Module containing all global common areas

 !Import contour advection and spectral modules:
use constants
use contours
use spectral

!-----------------------------------------------------------------------
 !Define quantities to be preserved between recontouring and evolution:
!-----------------------------------------------------------------------

 !Gridded PV field used for recontouring:
double precision:: qq(0:ny,0:nxm1,nz)

 !Spectral PV fields used in time stepping (note array order):
double precision:: qs(0:nxm1,0:ny,nz),qd(0:nxm1,0:ny,nz)

 !Current model time since beginning of simulation:
double precision:: t

 !Counters for saving gridded and contour data periodically:
integer:: igrids,iconts

end module common
