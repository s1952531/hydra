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

 !PV anomaly evolved in spectral space (note array order):
double precision:: qs(nx,ny)

 !Thermal equilibrium streamfunction (spectral):
double precision,allocatable,dimension(:,:):: ppeq

 !Time stepping parameters:
double precision:: dtmax,tgrid

end module common
