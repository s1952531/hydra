module constants

use parameters

integer,parameter:: nwx=nx/2
integer,parameter:: nyf=mgf*ny
integer,parameter:: nym1=ny-1

 !Generic double precision numerical constants:
double precision,parameter:: zero=0.d0,  one=1.d0, two=2.d0,  three=3.d0
double precision,parameter:: four=4.d0, five=5.d0, six=6.d0, twelve=12.d0
double precision,parameter:: f12=one/two,  f13=one/three,  f23=two/three
double precision,parameter:: f14=one/four, f32=three/two,  f43=four/three
double precision,parameter:: f53=five/three, f56=five/six, f74=7.d0/four
double precision,parameter:: f16=one/six, f112=one/twelve, f1112=11.d0/twelve
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0
double precision,parameter:: hpi=pi/two,twopi=two*pi
double precision,parameter:: small=1.d-12,small3=small*small*small

 !Domain constants:
double precision,parameter:: hlx=ellx/two,hly=elly/two
double precision,parameter:: y1=y1nd*elly,y2=y2nd*elly
double precision,parameter:: gly=elly/dble(ny),glx=ellx/dble(nx)
double precision,parameter:: hgly=gly/two,hglx=glx/two
double precision,parameter:: glyi=one/gly,hglyi=one/(two*gly)
double precision,parameter:: glyf=elly/dble(nyf),hglyf=glyf/two,glyfi=one/glyf

end module
