module constants

 !Include all modifiable parameters for use below:
use parameters

 !Grid dimensions +/- 1:
integer,parameter:: nxm1=nx-1,nym1=ny-1
integer,parameter:: nxp1=nx+1,nyp1=ny+1

 !Grid dimensions used in write statements:
integer,parameter:: nbytes=4*(nx+1)*(ny+1)

 !Generic double precision numerical constants: 
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0,four=4.d0
double precision,parameter:: f12=one/two,f14=one/four

 !Domain widths:
double precision,parameter:: ellx=xmax-xmin,elly=ymax-ymin
double precision,parameter:: hlx=ellx/two,hly=elly/two

 !Basic constants:
double precision,parameter:: dsumi=one/dble(nx*ny)

end module
