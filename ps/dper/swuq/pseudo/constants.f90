module constants

 !Include all modifiable parameters for use below:
use parameters

double precision,parameter:: zero=0.d0,  one=1.d0, two=2.d0,  three=3.d0
double precision,parameter:: four=4.d0, five=5.d0, six=6.d0, twelve=12.d0
double precision,parameter:: f12=one/two,  f13=one/three,  f23=two/three
double precision,parameter:: f14=one/four, f32=three/two,  f43=four/three
double precision,parameter:: f53=five/three, f56=five/six, f74=7.d0/four
double precision,parameter:: f16=one/six, f112=one/twelve, f1112=11.d0/twelve
double precision,parameter:: hpi=pi/two,twopi=two*pi

double precision,parameter:: hlx=ellx/two,hly=elly/two
double precision,parameter:: glx=ellx/dble(nx),gly=elly/dble(ny)
double precision,parameter:: dsumi=one/dble(nx*ny)
! The x range is -L_x/2 <= x <= L_x/2   
! The y range is -L_y/2 <= y <= L_y/2   

 !Basic constants:
double precision,parameter:: cof=two*omega,fsq=cof**2,csq=cgw**2
double precision,parameter:: cfl=0.2d0,dt=cfl*min(glx,gly)/cgw
double precision,parameter:: dt2=dt/two,alp=one/dt2,alpsq=alp**2
end module
