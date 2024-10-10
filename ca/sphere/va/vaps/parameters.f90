module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: ng=N_G, nt=2*ng
double precision,parameter:: ncont=N_CONTQ
 !Simulation time length etc..
double precision,parameter:: dt=T_STEP,tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE,tcsave=T_CSAVE
! ng     : number of  latitudinal grid points
! nt     : number of longitudinal grid points
! ncont  : Number of PV jumps (contours) spanning q_min to q_max
! dt     : time step (fixed)
! tsim   : total duration of the simulation
! tgsave : grid data save time increment
! tcsave : contour data save time increment

! ==> Physical parameters <==
double precision,parameter:: fpole=COR_FREQ,cgw=C_GW,hbar=H_BAR
double precision,parameter:: drate=D_RATE,nnu=POW_HYPER
! fpole  : Coriolis frequency f at the north pole
! cgw    : Short-scale gravity wave speed c
! hbar   : Mean fluid depth
! drate  : This times fpole is the damping rate on wavenumber ng
! nnu    : power of hyperviscosity
!----------------------------------------------------------------

end module parameters
