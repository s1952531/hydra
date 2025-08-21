module parameters

! This module contains all the modifiable parameters (except pi below)
! for the suite of pseudo-spectral f90 files.

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!***Numerical parameters:***
integer,parameter:: nx=N_X,ny=N_Y
! nx, ny: number of grid intervals in x and y
double precision,parameter:: tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 

!***Physical parameters:***
double precision,parameter:: ellx=L_X,elly=L_Y
double precision,parameter:: depth=L_Z
double precision,parameter:: cdamp=C_DAMP,nnu=POW_HYPER
! ellx   : domain width in x (periodic, centred at 0)
! elly   : domain width in y (periodic, centred at 0)
! depth  : N*H/f, the scaled depth
! cdamp  : hyperviscosity prefactor
! nnu    : hyperviscosity power
! ***note: if cdamp < 0, the damping rate is -cdamp*(k/k_max)^(2*nnu)
!          otherwise, it varies as cdamp*|zeta|_rms*(k/k_max)^(2*nnu)
!----------------------------------------------------------------

end module parameters
