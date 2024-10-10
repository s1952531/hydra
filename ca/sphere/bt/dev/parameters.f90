# 1 "parameters.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 13 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4

# 17 "/usr/include/stdc-predef.h" 3 4















































# 13 "<command-line>" 2
# 1 "parameters.f90"
module parameters

! This module contains all the modifiable parameters for 
! the suite of spe f90 files.

!***Numerical parameters:***
integer,parameter:: ng=256, nt=2*ng
integer,parameter:: nperiod=40
double precision,parameter:: dtmax=.05000000d0,tsave=0.5d0,tsim=5.0d0
double precision,parameter:: dq=0.2d0,cdamp=2.0d0
! ng, nt : number of grid boxes in latitude & longitude (inversion grid)
! nperiod: total number of periods to run (each of length tsim)
! dtmax  : maximum allowed timestep
! tsave  : approximate time interval between data saves (of fields)
! tsim   : simulation duration (one "period")
! dq     : PV jump across all contours
! cdamp  : this times zz_rms is the damping rate on wavenumber ng;
!          cdamp*zz_rms*(m/(ng*cos(lat)))^4 is the damping applied to
!          wavenumber m; latitude damping at the same rate is also done

!***Physical parameters:***
double precision,parameter:: asp=0.5d0,omega=0.0d0
double precision,parameter:: rekman=0.0d0,esr=0.0d0
integer,parameter:: ksr=1,iseed=646483
! asp    : the height:width aspect ratio of the ellipsoidal surface
! omega  : the planetary rotation rate
! rekman : Ekman friction rate (the relative vorticity zeta
!          is relaxed back to zero at the rate rekman)
! esr    : Enstrophy injection rate per unit of time
! ksr    : Wavenumber centroid of enstrophy injection; 
!          the enstrophy spectrum is proportional to 
!          k^5*exp(-2k^2/k0^2) --- see ranspec
! iseed  : seed for initialising stochastic forcing
!----------------------------------------------------------------

end module
