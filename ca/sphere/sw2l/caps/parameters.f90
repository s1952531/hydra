module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: numLayers=2, numLatGridPts=N_G, numLongGridPts=2*numLatGridPts
integer,parameter:: numContoursTot=N_CONTQ
double precision,parameter:: dt=T_STEP,tsim=T_SIM
double precision,parameter:: tGridSaveInterval=T_GSAVE,tContourSaveInterval=T_CSAVE
! numLayers     : number of layers (here 2, but the contour modules are general)
! numLatGridPts     : number of  latitudinal grid points (using a half grid)
! numLongGridPts     : number of longitudinal grid points (always equal to 2*numLatGridPts)
! numContoursTot  : number of PV jumps used to represent f = 2*Omega*sin(lat)
! dt     : time step (*fixed*)
! tsim   : total duration of the simulation
! tGridSaveInterval : grid data save time increment
! tContourSaveInterval : contour data save time increment

! ==> Physical parameters <==
double precision,parameter:: vm11=V_11, vm12=V_12, vm21=V_21, vm22=V_22
double precision,parameter:: hAvgLayer1=H_BAR1, hAvgLayer2=H_BAR2, densityRatio=DEN_RAT
double precision,parameter:: gravity=ACC_GRAV, ctsq1=C_MODE1, ctsq2=C_MODE2
double precision,parameter:: thermDampRate=R_THERM, hyperviscDampRate=D_RATE, hyperviscPower=POW_HYPER
double precision,parameter:: rmsPvForcing=FQ_RMS, tPvForcingDecorrelation=T_FQ
integer,parameter:: minSpherHarmDegree=N_BEG, maxSpherHarmDegree=N_END
! vmij       : Coefficients to convert layer quantities to mode quantities;
!              for a field F_j in layers j = 1 & 2, its projection onto
!              vertical mode 1 is F_1*vm11 + F_2*vm12, while for mode 2,
!              its projection is  F_1*vm21 + F_2*vm22
! hAvgLayer1      : Mean fluid depth of layer 1, H_1
! hAvgLayer2      : Mean fluid depth of layer 2, H_2
!              *** Note, wlog we take hAvgLayer1 + hAvgLayer2 = 1 throughout ***
! densityRatio     : density ratio, rho_2/rho_1
!              *** Note, densityRatio must be strictly less than 1 ***
! gravity    : The acceleration due to gravity (see hinvert script)
! ctsq1      : Short-scale gravity-wave speed associated with 1st vertical mode
! ctsq2      : Short-scale gravity-wave speed associated with 2nd vertical mode
! thermDampRate        : thermal damping rate
! hyperviscDampRate      : hyperviscDampRate*fpole is the hyperviscous damping rate on wavenumber numLatGridPts
! hyperviscPower        : power of hyperviscosity
! rmsPvForcing      : r.m.s. PV forcing amplitude
! tPvForcingDecorrelation        : Markovian de-correlation time of the forcing
! minSpherHarmDegree, maxSpherHarmDegree : min/max degree n of spherical harmonics Y_n^m used

!----------------------------------------------------------------------------
! For implementing topography if desired:
logical,parameter:: topoFlag=TOP_FLAG
 !Random number seed:
integer,parameter:: rngSeed=I_SEED

end module parameters
