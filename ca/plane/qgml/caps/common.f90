module common

 !Module containing all global common areas

 !Import contants, parameters and common arrays:
use constants
use contours
use spectral

!-----------------------------------------------------------------------
 !Define quantities to be preserved between recontouring and evolution:
!-----------------------------------------------------------------------

 !Gridded PV field used for recontouring:
double precision:: qq(ng,ng,nz)

 !Spectral PV fields used in time stepping:
double precision:: qs(ng,ng,nz),qd(ng,ng,nz)

 !Time and twist parameter (used for contour surgery):
double precision:: t

 !Counters for saving gridded and contour data periodically:
integer:: igrids,iconts

end module common
