module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use spectral

 !Define quantities which need to be preserved between recontouring and evolution:

 !Variables:
double precision:: t

 !Physical fields:
double precision:: hh(nx,ny),shh(nx,ny),dd(nx,ny),sdd(nx,ny)
double precision:: uu(ny,nx),vv(ny,nx)
double precision:: hhp(ny,nx),ddp(ny,nx)
double precision:: hhb(ny,nx)

end module
