module parameters

! This module contains all the modifiable parameters (except pi) for 
! the suite of caps f90 files.

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

 !Domain grid dimensions:
integer,parameter:: nx=400,ny=200

 !Conformal domain limits:
double precision,parameter:: xmin=1e-08,xmax=11.84333444237903
double precision,parameter:: ymin=-1.5707963171193782,ymax=1.5707963171193782

end module
