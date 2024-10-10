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
 
 !(Hyper)viscosity parameters:
integer,parameter:: nnu=N_NU
double precision,parameter:: prediss=PRE_DISS
! If nnu = 1, this is the molecular viscosity case.  Then, we 
! choose the viscosity nu = prediss*((b_max-b_min)/k_{x,max}^3)
! where k_{x_max} is the maximum x wavenumber.
! Note: prediss = 2 is recommended.
! ----------------------------------------------------------------
! If nnu > 1, this is the hyperviscosity case.  Then, the damping
! rate is prediss*zeta_char*(k/k_max)^(2*nnu) on wavenumber k
! where k_max is the maximum x or y wavenumber and zeta_char is
! a characteristic vorticity (see subroutine adapt of strat.f90).
! Note: nnu = 3 and prediss = 10 are recommended.

end module
