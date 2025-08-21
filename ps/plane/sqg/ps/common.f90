module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use spectral

 !Spectral prognostic field (note array order):
double precision:: qs(nx,ny),qspre(nx,ny)

 !Time stepping parameters:
double precision:: dtmax,tfin,tgrid

end module common
