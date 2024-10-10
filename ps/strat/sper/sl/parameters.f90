module parameters

! This module contains all the modifiable parameters for 
! the suite of casl f90 files.

 !Domain grid dimensions:
integer,parameter:: nx=N_X,ny=N_Y

 !Domain width in x and limits in y:
double precision,parameter:: ellx=L_X
double precision,parameter:: ymin=Y_MIN,ymax=Y_MAX

 !Simulation time length etc..
double precision,parameter:: tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE,tcsave=T_CSAVE
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)
 
 !Reference translational velocity (added to u):
double precision,parameter:: uref=U_REF

end module
