module parameters

! This module contains all the modifiable parameters (except nz & pi below)
! for the suite of casl f90 files.

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!***Numerical parameters:***
integer,parameter:: nx=512,ny=512

 !Number of contours used for representing the PV variation:
integer,parameter:: ncontq=100
! ncontq : used to compute the PV contour interval from 
!          dq = (qq_max-qq_min)/ncontq

 !Simulation time length etc..
double precision,parameter:: tper=1.d0
integer,parameter:: totper=200
double precision,parameter:: tgsave=0.2d0,tcsave=25.0d0
! tper   : time period of the velocity field
! totper : total number of periods to simulate
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)

!***Other physical parameters:***
double precision,parameter:: ellx=2.d0*pi,elly=2.d0*pi
double precision,parameter:: uscale=pi
double precision,parameter:: ku=1.d0
integer,parameter:: iseed=7438
! ellx   : domain width in x (periodic, centred at 0)
! elly   : domain width in y (periodic, centred at 0)
!...

!----------------------------------------------------------------

end module
