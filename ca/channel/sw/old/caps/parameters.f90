module parameters

! This module contains all modifiable parameters (except pi below)

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! Number of grid boxes in the x & y directions (inversion grid):
integer,parameter:: nx=N_X
integer,parameter:: ny=N_Y

! Number of contours used for representing PV:
integer,parameter:: ncontq=N_CONTQ
! The PV contour interval is found from dq = (qq_max-qq_min)/ncontq
! where qq_min & qq_max are the minimum & maximum PV values initially.

! Simulation time length etc..
double precision,parameter:: dt=T_STEP,tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE,tcsave=T_CSAVE
! dt     : time step (fixed) - usually gravity-wave resolving
! tsim   : total duration of the simulation
! tgsave : grid data save time increment
! tcsave : contour data save time increment
! ***NOTE*** tcsave should always be an integer multiple of tgsave

!***Physical parameters:***
double precision,parameter:: ellx=L_X
double precision,parameter:: elly=L_Y
double precision,parameter:: cof=COR_FREQ,cgw=C_GW
double precision,parameter:: cdamp=C_DAMP,nnu=POW_HYPER
! ellx   : domain length in x
! elly   : domain  width in y
! cof    : Constant Coriolis frequency f
! cgw    : Short-scale gravity wave speed c
! cdamp  : This times f is the damping rate on the maximum x wavenumber
! nnu    : power of hyperviscosity (used only in the periodic x direction)
!----------------------------------------------------------------

end module parameters
