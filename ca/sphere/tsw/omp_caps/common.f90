module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral

 !Quantities which need to be preserved between recontouring and evolution:

 !Height & divergence (semi-spectral fields):
double precision:: hh(ng,nt),dd(ng,nt)

 !Height & divergence (physical space):
double precision:: hhp(ng,nt),ddp(ng,nt)

 !Equilibrium height & topographic source:
double precision:: hhe(ng,nt),ttdd(ng,nt) 

end module
