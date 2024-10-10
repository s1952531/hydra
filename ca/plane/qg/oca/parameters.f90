module parameters

! This module contains all the modifiable parameters (except nz & pi below)
! for the suite of casl f90 files.

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!***Numerical parameters:***
integer,parameter:: nx=256,ny=256

 !Simulation time length etc..
double precision,parameter:: tsim=25.d0
double precision,parameter:: tgsave=0.25d0,tcsave=1.d0
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)

 !PV jump across all contours:
integer,parameter:: qjump=4.d0*pi

!***Physical parameters:***
double precision,parameter:: ellx=2.d0*pi,elly=2.d0*pi
double precision,parameter:: kd=2.d0,beta=0.d0
! ellx   : domain width in x (periodic, centred at 0)
! elly   : domain width in y (periodic, centred at 0)
! kd     : Rossby deformation wavenumber associated with the baroclinic mode
! beta   : planetary vorticity gradient

end module
