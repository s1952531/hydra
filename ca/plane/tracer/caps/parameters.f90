module parameters

! This module contains all the modifiable parameters (except nz & pi below)
! for the suite of casl f90 files.

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!***Numerical parameters:***
integer,parameter:: nx=256,ny=256

 !Number of contours used for representing the PV variation:
integer,parameter:: ncontq=50
! ncontq : used to compute the PV contour interval from 
!          dq = (qq_max-qq_min)/ncontq

 !Simulation time length etc..
double precision,parameter:: tsim=250.d0
double precision,parameter:: tgsave=1.0d0,tcsave=20.d0
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)

!***Physical parameters:***
double precision,parameter:: ellx=2.d0*pi,elly=2.d0*pi
double precision,parameter:: kd=0.d0,beta=0.d0
double precision,parameter:: rtherm=0.d0,rekman=0.d0
double precision,parameter:: cdamp=10.d0,nnu=3
double precision,parameter:: esr=0.d0,vorvor=1.d0
double precision,parameter:: srcalp=0.044d0,uscale=pi
double precision,parameter:: ku=1.d0,kf=1.d0
double precision,parameter:: tper=1.d0
integer,parameter:: ivor=0,iseed=764
! ellx   : domain width in x (periodic, centred at 0)
! elly   : domain width in y (periodic, centred at 0)
! kd     : Rossby deformation wavenumber associated with the baroclinic mode
! beta   : planetary vorticity gradient
! rtherm : thermal damping rate
! rekman : Ekman   damping rate
! cdamp  : this times |zeta|_rms*(k/k_max)^(2*nnu) is the hyperdiffusivity 
!          coefficient, where nnu is the power specified above
! nnu    : hyperviscous power, e.g. Lap^nnu.
! esr    : enstrophy input rate (via pairs of point vortices which
!          are converted to gridded vorticity and added to qd)
! vorvor : the mean vorticity associated with the point vortex 
!          upon placement on the grid
! ivor   : use 1 for arbitrarily spaced pairs (effectively monopoles)
!          use 2 for dipoles concentrated at a point, and 0 otherwise
! iseed  : seed for initialising point vortex forcing
!----------------------------------------------------------------

end module
