module parameters

! This module contains all the modifiable parameters for 
! the suite of tsw f90 files.

!***Numerical parameters:***
integer,parameter:: ng=N_G, nt=2*ng, nlay=1, nperiod=N_PER
double precision,parameter:: dtmax=DT_MAX, tsave=T_SAV, tsim=T_SIM
double precision,parameter:: dq=DQ, drate=D_RAT
! nt,ng  : number of grid boxes in the x & y directions (inversion  grid)
! nlay   : number of layers 
! nperiod: total number of periods to run (each of length tsim)
! dtmax  : maximum allowed timestep (note this should be < or = pi/(ng*cgw) to
!          resolve the short scale gravity waves) 
! tsave  : approximate time interval between data saves (of fields)
! tsim   : simulation duration (one "period")
! dq     : PV jump across all contours
! drate  : height & divergence damping rate per day on wavenumber ng;
!          drate*omega*(m/(ng*cos(lat)))^4 is the damping applied to
!          wavenumber m; latitude damping at the same rate is also done

!***Physical parameters:***
double precision,parameter:: omega=OMEGA, cgw=C_GW
double precision,parameter:: kappa=2.d0/7.d0, aa0=1.d0-kappa
double precision,parameter:: alpha=ALPHA, esr=ESR
integer,parameter:: ksr=KSR, iseed=ISEED
! omega  : the planetary rotation rate
! cgw    : short-scale gravity wave speed, sqrt{g*H_char}
! kappa  : ratio of specific heats, normally 2/7
! aa0    : associated with vertical temperature structure, normally 1 - kappa
! alpha  : thermal damping rate
! esr    : Enstrophy injection rate per unit of time
! ksr    : Wavenumber centroid of enstrophy injection; 
!          the enstrophy spectrum is proportional to 
!          k^5*exp(-2k^2/k0^2) --- see ranspec
! iseed  : seed for initialising stochastic forcing
!----------------------------------------------------------------

end module
