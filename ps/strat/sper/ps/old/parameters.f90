module parameters

! This module contains all the modifiable parameters for 
! the suite of ps f90 files.

 !Domain grid dimensions:
integer,parameter:: nx=N_X,ny=N_Y

 !Domain width in x and limits in y:
double precision,parameter:: ellx=L_X
double precision,parameter:: ymin=Y_MIN,ymax=Y_MAX

 !Simulation duration and data save interval:
double precision,parameter:: tsim=T_SIM,tgsave=T_GSAVE
 
 !Reference translational velocity (added to u):
double precision,parameter:: uref=U_REF

 !Hyperviscosity parameters:
integer,parameter:: nnu=N_NU
double precision,parameter:: prediss=PRE_DISS
! nnu      : power of hyperviscous damping term (1 -> Laplacian)
! prediss  : nu = prediss/kmax^(2*nnu) where kmax is the maximum
!          : wavenumber magnitude.
!          : (Note for nnu=3, prediss=80 is recommended)

end module
