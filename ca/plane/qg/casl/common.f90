module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral
use generic

 !Define quantities which need to be preserved between recontouring and evolution:

 !Gridded PV fields:
double precision:: qs(ny,nx),qd(ny,nx)

 !For semi-Lagrangian advection: 
double precision:: xig(nx),yig(ny)
double precision,parameter:: xigmax=dble(nx),yigmax=dble(ny) 

 !Time stepping parameters:
double precision:: dtmax,tfin,tgrid

 !Parameter for vortex forcing:
double precision:: dnvor

end module
