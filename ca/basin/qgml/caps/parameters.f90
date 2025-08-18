module parameters

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! This module contains all the modifiable parameters for the suite of
! contour advection (CLAM) codes in a closed rectangular basin.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 !Number of vertical layers:
integer,parameter:: nz=N_Z

 !Number of grid intervals in x:
integer,parameter:: nx=N_X

 !Number of grid intervals in y:
integer,parameter:: ny=N_Y

 !Horizontal domain limits:
double precision,parameter:: xmin=X_MIN
double precision,parameter:: xmax=X_MAX
double precision,parameter:: ymin=Y_MIN
double precision,parameter:: ymax=Y_MAX

 !Number of contours used for representing the total PV contrast:
integer,parameter:: ncontq=N_CONTQ
! ncontq : used to compute the PV contour interval from 
!          dq = (qq_max-qq_min)/ncontq

! Duration of simulation and data save intervals:
double precision,parameter:: tsim=T_SIM
double precision,parameter:: tgsave=T_GSAVE
double precision,parameter:: tcsave=T_CSAVE
! tsim   : total duration of the simulation
! tgsave : time interval between gridded (field) data saves
! tcsave : time interval between contour data saves
!          *** should be a multiple of tgsave ***
double precision,parameter:: tavg_rme=T_AVG
! tavg_rme : if > 0, this time interval is used to accumulate
!            the running-mean total energy E_rm, which, if
!            within a fraction frme of its previous value,
!            signals the simulation to stop and save the
!            current PV field in qq_final.r8.
double precision,parameter:: frme=F_RME
! frme     : see above. Use 0 to avoid running-mean calculation

! Hyperdiffusion damping parameters:
double precision,parameter:: cdamp=C_DAMP,nnu=POW_HYPER
! cdamp  : this times |PVA|_rms*(k/k_max)^(2*nnu) is the hyperdiffusivity
!          coefficient for qd, where nnu is the power specified above and
!          PVA is the PV anomaly (q - beta*y).

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

!***Physical parameters:***
double precision,parameter:: beta=PV_GRAD
double precision,parameter:: rekman=R_EKMAN
double precision,parameter:: fwind=F_WIND
logical,parameter:: bath=BATH_FLAG
! beta   : planetary vorticity gradient
! rekman : Ekman damping rate (1/tau_E)
! fwind  : wind stress forcing applied to upper layer PV is
!          fwind*sin(pi*(y-y_min)/(y_max-y_min))
! bath   : logical variable to indicate presence of bathymetry
!          (the form of bathymetry is read in from a file)
!----------------------------------------------------------------

 !Random number seed:
integer,parameter:: iseed=2073600

end module parameters
