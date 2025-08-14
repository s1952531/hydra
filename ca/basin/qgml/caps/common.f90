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
double precision:: qq(0:ny,0:nx,nz)

 !Spectral PV fields used in time stepping (note array order):
double precision:: qs(0:nx,0:ny,nz),qd(0:nx,0:ny,nz)

 !Current model time since beginning of simulation:
double precision:: t

 !Optionally used for checking energy equilibration:
double precision:: rme,rmep,dtavg_rme

 !Counters for saving gridded and contour data periodically:
integer:: igrids,iconts

end module common
