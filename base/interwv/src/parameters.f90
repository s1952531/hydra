module parameters

! This module contains all the modifiable parameters 
! for the epbeq quite of f90 files.

integer,parameter:: nx=1024, ny=128
integer,parameter:: mgf=16
integer,parameter:: nitmax=10000
double precision,parameter:: ellx=16.0d0,elly=1.0d0
double precision,parameter:: y1nd=0.8d0,y2nd=0.9d0
double precision,parameter:: bf1=0.0d0,bf2=1.d0,bf3=1.0d0
double precision,parameter:: bfsmoond=2.0d0
double precision,parameter:: outsmoond=2.0d0
double precision,parameter:: orelax=0.2d0
double precision,parameter:: tol=1.d-7,derrmax=4.d0
! nx,ny      : Basic grid dimensions in x & y
! mgf        : Ratio of the fine grid to coarse used in y direction
! nitmax     : Maximum number of iterations allowed for convergence 
!              to a given amplitude solution
! ellx,elly  : The domain is the rectangle -L_x/2 < x < L_x/2
!              and 0 < y < L_y
! bf1,bf2,bf3: Buoyancy frequencies in lower, middle and upper layers
!              respectively (bf2=1 can always be guaranteed by 
!              suitable scaling)
! bfsmoond   : Smoothing distance between layers within the domain 
!              given in dimensionless units of number of grid cells 
!              (typically =2)
! outsmoond  : Smoothing distance to buoyancy nominally outside the domain
!              given in dimensionless units of number of grid cells 
!              (typically =2)
! orelax     : Relaxation coefficient used in the iterative solver
!              (typically a value around 0.2 helps convergence)
! tol,derrmax: Tolerance of rms error in converged solution for a 
!              given amplitude, and maximum error allowed before 
!              iterations are halted
end module
