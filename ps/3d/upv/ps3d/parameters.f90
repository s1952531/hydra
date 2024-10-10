module parameters

! Module containing all the modifiable parameters (except pi below).

! The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

! ==> Numerical parameters <==
integer,parameter:: ng=1024,nz=32
double precision,parameter:: dt=0.01d0,tsim=10.d0
double precision,parameter:: tgsave=1.d0
! nz     : number of vertical layers
! ng     : inversion grid resolution in both x and y
!          (Note: the domain is a 2*pi periodic box horizontally)
! dt     : time step (fixed)
! tsim   : total duration of the simulation
! tgsave : grid data save time increment

! ==> Physical parameters <==
double precision,parameter:: cof=1.d0,bvf=8.d0
double precision,parameter:: depth=pi/2.d0
double precision,parameter:: cdamp=10.d0,nnu=3
! cof    : Coriolis frequency f
! bvf    : Buoyancy frequency N of the background linear stratification
! depth  : NH/f, total fluid depth H stretched by N/f; note the domain
!          width is L = 2*pi. To have an isotropic grid in unscaled
!          coordinates, ensure L/ng = H/nz, i.e.
!              2*pi/ng = cof*depth/(bvf*nz)
! cdamp  : This times f is the damping rate on wavenumber ng/2
! nnu    : Power of hyperviscosity
!----------------------------------------------------------------

end module parameters
