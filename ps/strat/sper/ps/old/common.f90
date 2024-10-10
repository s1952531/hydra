module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use spectral

 !Buoyancy & vorticity fields:
double precision:: bb(0:ny,0:nxm1),zz(0:ny,0:nxm1)

 !Time stepping parameters:
double precision:: dt,hfdt,qudt
double precision:: t,tgrid
integer:: igrids

 !Background potential energy calculation flag:
integer:: iene

end module
