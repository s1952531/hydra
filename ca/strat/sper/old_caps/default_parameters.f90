module parameters

! This module contains all the modifiable parameters for 
! the suite of casl f90 files.

 !Domain grid dimensions:
integer,parameter:: nx=1024,ny=128

 !Domain width in x and limits in y:
double precision,parameter:: ellx=8.d0
double precision,parameter:: ymin=0.d0,ymax=1.d0

 !Number of contours used for representing buoyancy and vorticity:
integer,parameter:: ncontb=100,ncontz=20
! ncontb : used to compute the buoyancy contour interval from 
!          dbb = (bb_max-bb_min)/ncontb 
! ncontz : used to compute the vorticity contour interval from 
!          dzz = (zz_2/zz_1)/ncontz  where zz_n is the L_n norm

 !Simulation time length etc..
double precision,parameter:: tsim=20.d0
double precision,parameter:: tgsave=0.2d0,tcsave=2.0d0
! tsim   : total duration of the simulation
! tgsave : grid data save time increment 
! tcsave : contour data save time increment (approximate)
 
 !Reference translational velocity (added to u):
double precision,parameter:: uref=0.d0

 !(Hyper-)Viscosity parameters (for residual vorticity only):
integer,parameter:: nnu=3
double precision,parameter:: prediss=100.d0
! nnu      : power of damping term (1 -> Laplacian)
! prediss  : if nnu=1 --> nu = prediss * (bbmax-bbmin)/kmax^(3/2) 
!          : where kmax is the max wavenumber magnitude 
!          : if nnu>1 --> nu = prediss/kmax^(2*nnu)   
!          : (Note for nnu=3, prediss=80 is recommended)

end module
