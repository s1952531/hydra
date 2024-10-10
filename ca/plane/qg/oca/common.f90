module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral
use generic

 !Main PV field:
double precision:: qq(ny,nx)

 !Time stepping parameters:
double precision:: dtmax,tfin,tgrid

end module
