# 1 "parameters.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 13 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4

# 17 "/usr/include/stdc-predef.h" 3 4























# 13 "<command-line>" 2
# 1 "parameters.f90"
module parameters

! This module contains all the modifiable parameters (except pi) for 
! the suite of caps f90 files.

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

 !Domain grid dimensions:
integer,parameter:: nx=400,ny=200

 !Conformal domain limits:
double precision,parameter:: xmin=1e-08,xmax=11.84333444237903
double precision,parameter:: ymin=-1.5707963171193782,ymax=1.5707963171193782

 !Physical domain x & y coordinates specifying weir, river & sea states:
double precision:: xsea,xwsea,xwriv,xriv
double precision:: ysea,ybsea,ybriv,yriv

 !River and sea uniform flow speeds:
double precision:: usea,vsea,uriv,vriv

 !Number of contours used for representing buoyancy and vorticity:
integer,parameter:: ncontb=100,ncontz=50
! ncontb : used to compute the buoyancy contour interval from 
!          dbb = (bb_max-bb_min)/ncontb 
! ncontz : used to compute the vorticity contour interval from 
!          dzz = (zz_2/zz_1)/ncontz  where zz_n is the L_n norm

 !Simulation time length etc..
double precision,parameter:: tsim=100.0d0
double precision,parameter:: tgsave=0.25d0,tcsave=5.0d0
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)

 !(Hyper-)Viscosity parameters (for residual vorticity only):
integer,parameter:: nnu=3
double precision,parameter:: prediss=100.0d0
! nnu      : power of damping term (1 -> Laplacian)
! prediss  : if nnu=1 --> nu = prediss * (bbmax-bbmin)/kmax^(3/2) 
!          : where kmax is the max wavenumber magnitude 
!          : if nnu>1 --> nu = prediss/kmax^(2*nnu)   
!          : (Note for nnu=3, prediss=80 is recommended)

end module
