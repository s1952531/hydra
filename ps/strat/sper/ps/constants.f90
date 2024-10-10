module constants

 !Include all modifiable parameters for use below:
use parameters

! Contains all the non-modifiable parameters as well as all 
! quantities which never change throughout a simulation
! for the suite of ps f90 codes.

 !Grid dimensions - 1:
integer,parameter:: nxm1=nx-1,nym1=ny-1

 !Grid dimensions used in write statements:
integer,parameter:: ngridp=nx*(ny+1),nbytes=4*(ngridp+1)

 !Generic double precision numerical constants: 
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0
double precision,parameter:: three=3.d0,four=4.d0
double precision,parameter:: f12=one/two,f14=one/four
double precision,parameter:: f13=one/three,f23=two/three
double precision,parameter:: pi=3.14159265358979d0,twopi=two*pi
double precision,parameter:: small=1.d-12

 !Domain lengths and inverses:
double precision,parameter:: xmax=ellx/two,xmin=-xmax
double precision,parameter:: elly=ymax-ymin

 !Basic constants:
double precision,parameter:: glx=ellx/dble(nx),gly=elly/dble(ny)
double precision,parameter:: glmin=min(glx,gly)
double precision,parameter:: garea=glx*gly,dsumi=one/dble(nx*ny)

end module
