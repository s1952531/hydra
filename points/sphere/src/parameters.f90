module parameters

! This module contains all the modifiable parameters for 
! the suite of pvs f90 files.

 !Total number of point vortices:
integer,parameter:: n=N_POINTS

 !Initial frame to start from (use 0 to start from t = 0):
integer,parameter:: iniframe=0

 !Approximate time between data saves & simulation duration:
double precision,parameter:: tsave=T_SAV,tsim=T_MAX

 !For display purposes, the number of latitudes used to convert points
 !to small circles for imaging (used in p2g.f90 and image.f90):
integer,parameter:: ng=N_G

end module
