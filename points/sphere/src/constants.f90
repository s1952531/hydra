module constants

use parameters

! Contains all the non-modifiable parameters as well as all 
! quantities which never change throughout a simulation

 !Point vortex strengths:
double precision:: s(n)

 !Generic double precision numerical constants: 
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0,three=3.d0,six=6.d0
double precision,parameter:: half=one/two,third=one/three,sixth=one/six
double precision,parameter:: big=1.d14

end module
