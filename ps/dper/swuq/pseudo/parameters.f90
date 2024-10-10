module parameters

! This module contains all the modifiable parameters for 
! the suite of pseudosp f90 files.

 !Define pi for use below: ***don't alter***
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!***Numerical parameters:***
integer,parameter:: nx=N_X, ny=N_Y
double precision,parameter:: tdsave=T_SAV,tsim=T_SIM
double precision,parameter:: ellx=L_X,elly=L_Y
double precision,parameter:: drate=D_RATE,hyppow=H_POW
! nx,ny    : number of grid boxes in the x & y directions
! tdsave   : approximate time interval between data saves (of fields)
! tsim     : simulation duration
! ellx,elly: domain lengths L_x & L_y in x & y directions resp.
! drate    : divergence damping has the form: drate*omega/K_max^2p where
!            K_max^2 equals the product of the maximum x & y wavenumbers
! hyppow   : power of hyper-diffusion (p in above definition)
!----------------------------------------------------------------

!***Physical parameters:***
double precision,parameter:: omega=OMEGA,cgw=CGW
logical,parameter:: topo=.false.
! omega  : the planetary rotation rate
! cgw    : gravity wave speed 
! topo   : logical specifying existence of topography
!----------------------------------------------------------------

end module
