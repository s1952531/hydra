module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: ng=256
double precision,parameter:: ncont=80
 !Simulation time length etc..
double precision,parameter:: dt=0.0025d0,tsim=25.d0
double precision,parameter:: tgsave=0.25d0,tcsave=2.5d0
! ng     : inversion grid resolution in both x and y
!          (Note: the domain is a 2*pi periodic box.)
! ncont  : Number of PV jumps (contours) spanning q_min to q_max
! dt     : time step (fixed)
! tsim   : total duration of the simulation
! tgsave : grid data save time increment
! tcsave : contour data save time increment

! ==> Physical parameters <==
double precision,parameter:: cof=4.d0*pi,cgw=2.d0*pi,hbar=0.2d0
double precision,parameter:: cdamp=10.d0,nnu=3
! cof    : Constant Coriolis frequency
! cgw    : Short-scale gravity wave speed
! hbar   : Mean fluid height H
! cdamp  : This times f is the damping rate on wavenumber ng/2
!----------------------------------------------------------------

end module parameters
