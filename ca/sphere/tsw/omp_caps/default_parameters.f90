module parameters

! This module contains all the modifiable parameters for 
! the suite of spe f90 files.

!***Numerical parameters:***
integer,parameter:: ng=256, nt=2*ng, nlay=1
integer,parameter:: nperiod=20
double precision,parameter:: dtmax=0.0017469281074166689d0
double precision,parameter:: tsave=0.25d0,tsim=1.0d0,tramp=0.d0
double precision,parameter:: dq=0.196349540850d0,drate=8.d0
! nt,ng  : number of grid boxes in the x & y directions (inversion  grid)
! nlay   : number of layers 
! nperiod: total number of periods to run (each of length tsim)
! dtmax  : maximum allowed timestep (note this should be < or = pi/(ng*cgw) to
!          resolve the short scale gravity waves) 
! tsave  : approximate time interval between data saves (of fields)
! tsim   : simulation duration (one "period")
! tramp  : ramp period for PV anomaly to reach final form specified
!          on initialisation (use tramp = 0 for no ramping)
! dq     : PV jump across all contours
! drate  : residual PV damping rate per day on wavenumber ng;
!          drate*omega*(m/(ng*cos(lat)))^4 is the damping
!          applied to wavenumber m of the residual PV;
!          latitude damping at the same rate is also done

!***Physical parameters:***
double precision,parameter:: omega=6.2831853072d0,cgw=7.024814731061d0
double precision,parameter:: rtherm=0.d0,atherm=0.d0
double precision,parameter:: rekman=0.d0,esr=0.d0
integer,parameter:: ksr=1,iseed=0
double precision,parameter:: atopo=0.d0,btopo=0.d0,ftopo=0.d0
! omega  : the planetary rotation rate
! cgw    : short-scale gravity wave speed, sqrt{gH}
! rtherm, atherm : The thermal relaxation rate 
!          r_th = rtherm*(atherm+rtherm*t)/(1+rtherm*t);
!          hh is relaxed back to a specified field hhe at the
!          rate r_th)
! rekman : Ekman friction rate (the relative vorticity zeta
!          is relaxed back to zero at the rate rekman)
! esr    : Enstrophy injection rate per unit of time
! ksr    : Wavenumber centroid of enstrophy injection; 
!          the enstrophy spectrum is proportional to 
!          k^5*exp(-2k^2/k0^2) --- see ranspec
! iseed  : seed for initialising stochastic forcing
! atopo  : time mean topographic amplitude
! btopo  : amplitude of time variation of topography
! ftopo  : frequency of time variation of topography
!          Topography varies as: atopo+btopo*sin(ftopo*t)
!----------------------------------------------------------------


end module
