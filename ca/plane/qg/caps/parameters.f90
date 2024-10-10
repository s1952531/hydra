module parameters

! This module contains all the modifiable parameters (except pi below)
! for the suite of casl f90 files.

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!***Numerical parameters:***
integer,parameter:: nx=N_X, ny=N_Y
! nx, ny: number of grid intervals in x and in y

 !Number of contours used for representing the PV variation:
integer,parameter:: ncontq=N_CONTQ
! ncontq : used to compute the PV contour interval from 
!          dq = (qq_max-qq_min)/ncontq

 !Simulation time length etc..
double precision,parameter:: tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE, tcsave=T_CSAVE
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)

!***Physical parameters:***
double precision,parameter:: ellx=L_X, elly=L_Y
double precision,parameter:: kd=K_D, beta=PV_GRAD, ubar=U_BAR
double precision,parameter:: rtherm=R_THERM, rekman=R_EKMAN
double precision,parameter:: cdamp=C_DAMP, nnu=POW_HYPER
double precision, parameter:: fbase=F_BASE, fwid=F_WID
integer, parameter:: kfmid=KF_MID, nfamp=NF_AMP, iseed=I_SEED
! ellx   : domain width in x (periodic, centred at 0)
! elly   : domain width in y (periodic, centred at 0)
! kd     : Rossby deformation wavenumber associated with the baroclinic mode
! beta   : planetary vorticity gradient
! ubar   : mean velocity in the x direction
! rtherm : thermal damping rate
! rekman : Ekman   damping rate
! cdamp  : this times |zeta|_rms*(k/k_max)^(2*nnu) is the hyperdiffusivity 
!          coefficient, where nnu is the power specified above
! fbase  : prefactor for energy input rate (narrow band spectral forcing)
! nfamp  : exponent  for energy input rate
!          *** Note, the energy input rate epsilon = fbase*2**nfamp 
!              (Suggested values for nfamp: -5, -7, -10, -12, -14).
! kfmid  : centre of band
! fwid   : width  of band
! iseed  : seed for initialising forcing
!----------------------------------------------------------------

end module parameters
