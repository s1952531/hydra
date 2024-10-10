module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: ng=N_G, nz=2
double precision,parameter:: ncont=N_CONTQ
 !Simulation time length etc..
double precision,parameter:: dt=T_STEP, tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE, tcsave=T_CSAVE
! ng     : inversion grid resolution in both x and y
!          (Note: the domain is a 2*pi periodic box.)
! nz     : number of layers (here 2, but the contour modules are general)
! ncont  : Number of PV jumps (contours) spanning q_min to q_max
! dt     : time step (fixed)
! tsim   : total duration of the simulation
! tgsave : grid data save time increment
! tcsave : contour data save time increment

! ==> Physical parameters <==
double precision,parameter:: hbar1=H_BAR1, hbar2=H_BAR2, alpha=DEN_RAT
double precision,parameter:: cof=COR_FREQ, cgw=C_GW
double precision,parameter:: cdamp=C_DAMP, nnu=POW_HYPER
! hbar1  : Mean fluid depth of layer 1, H_1
! hbar2  : Mean fluid depth of layer 2, H_2
! alpha  : density ratio, rho_2/rho_1
! cof    : Constant Coriolis frequency
! cgw    : Short-scale gravity wave speed, sqrt(g*(H_1+H_2))
! cdamp  : This times f is the damping rate on wavenumber ng/2
!----------------------------------------------------------------

end module parameters
