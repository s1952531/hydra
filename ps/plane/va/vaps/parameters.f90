module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: ng=N_G
double precision,parameter:: dt=T_STEP,tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE
! ng     : inversion grid resolution in both x and y
!          (Note: the domain is a 2*pi periodic box.)
! dt     : time step (fixed)
! tsim   : total duration of the simulation
! tgsave : grid data save time increment

! ==> Physical parameters <==
double precision,parameter:: cof=COR_FREQ,cgw=C_GW,hbar=H_BAR
double precision,parameter:: cdamp=C_DAMP,nnu=POW_HYPER
! cof    : Constant Coriolis frequency f
! cgw    : Short-scale gravity wave speed c
! hbar   : Mean fluid height H
! cdamp  : This times f is the damping rate on wavenumber ng/2
!----------------------------------------------------------------

end module parameters
