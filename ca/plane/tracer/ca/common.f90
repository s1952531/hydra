module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral
use generic

 !Define quantities which need to be preserved between recontouring and evolution:

 !Gridded PV fields:
double precision:: qq(ny,nx)
 !Time stepping parameters:
double precision:: tgrid
 !Velocity parameters:
double precision:: theta1,theta2
integer:: nsteps,nper,iperkeep,jtkeep,itkeep

end module
