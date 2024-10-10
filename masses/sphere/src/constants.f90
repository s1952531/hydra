module constants

use parameters

! Contains all the non-modifiable parameters as well as all 
! quantities which never change throughout a simulation

 !Point masses / 4*pi:
double precision:: s(n),stot

 !Number of bytes written when saving direct access data:
integer,parameter:: nbytes=4*(3*n+1)

 !Generic double precision numerical constants: 
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,six=6.d0
double precision,parameter:: half=one/two,third=one/three,sixth=one/six
double precision,parameter:: big=1.d14

double precision,parameter:: pi=3.14159265358979d0,twopi=two*pi

end module
