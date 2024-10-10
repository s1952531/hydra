module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral
use generic

 !Define quantities to be preserved between recontouring and evolution:

 !PV residue:
double precision:: qr(ny,nx)

 !Spectral prognostic fields (note array order):
double precision:: qs(nx,ny)

 !Exponential integrating factors:
double precision:: emq(nx,ny),epq(nx,ny)

 !Time stepping parameters:
double precision:: dtmax,tfin,tgrid


end module
