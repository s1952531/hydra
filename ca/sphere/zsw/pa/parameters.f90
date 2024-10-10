module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision, parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

double precision, parameter:: fpole=4.d0*pi, rrad=0.05d0
! fpole  : Coriolis frequency f at the north pole
! rrad   : Rossby deformation length, sqrt(g*H)/fpole (H = mean depth)

integer, parameter:: n=500
! n      : number of latitudinal divisions *** MUST BE DIVISIBLE BY 10 ***

double precision, parameter:: cfl0=0.1d0
! cfl0   : CFL number = c*dt/dphi initially (c = GW speed)

integer, parameter:: nsteps=4*n/10, nloops=5000
! nsteps : number of time steps between field data saves
! nloops : total number of field data saves (apart from first)

end module parameters
