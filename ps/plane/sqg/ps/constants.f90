module constants

 !Include all modifiable parameters for use below:
use parameters

! Contains all the non-modifiable parameters as well as all 
! quantities which never change throughout a simulation
! for the suite of pseudo-spectral f90 codes.

integer,parameter:: nxm1=nx-1,nym1=ny-1

 !Grid dimensions used in write statements:
integer,parameter:: ngridp=nx*ny,nbytes=4*(ngridp+1)

 !Generic double precision numerical constants: 
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0,three=3.d0
double precision,parameter:: four=4.d0,six=6.d0,f12=one/two,f13=one/three
double precision,parameter:: f14=one/four,f16=one/six
double precision,parameter:: twopi=two*pi
double precision,parameter:: small=1.d-12
double precision,parameter:: oms=one-small

 !Domain limits:
double precision,parameter:: hlx=f12*ellx,xmin=-hlx,xmax=hlx
double precision,parameter:: hly=f12*elly,ymin=-hly,ymax=hly

 !Basic constants:
double precision,parameter:: domarea=ellx*elly,aspect=ellx/elly
double precision,parameter:: glx=ellx/dble(nx),glxi=dble(nx)/ellx
double precision,parameter:: gly=elly/dble(ny),glyi=dble(ny)/elly
double precision,parameter:: garea=glx*gly,dsumi=one/dble(nx*ny)

 !Logical for using fixed or varying hyperviscosity:
logical,parameter:: fixed_viscosity=cdamp < zero

end module constants
