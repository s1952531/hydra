module parameters

! This module contains all the modifiable parameters for 
! the suite of pms f90 files.

 !Total number of point masses:
integer,parameter:: n=2

 !Approximate time between data saves & simulation duration:
double precision,parameter:: tsave=0.1d0,tsim=10.0d0

 !For display purposes, the number of latitudes used to convert points
 !to small circles for imaging (used in p2g.f90 and image.f90):
integer,parameter:: ng=512

end module
