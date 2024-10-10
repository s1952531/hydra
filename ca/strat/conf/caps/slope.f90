program conform

use stafft

! =========================================================================
! This routine sets up the conformal transformation of an aperiodic domain
! given the bottom topography h(x) in the function hp below.

! Use precomp to pass the parameter values for nx, ny, ellx & ymax.
! Then compile using
! gfortran -O3 -o conform conform.f90 ~/hydra/lib/stafft/stafft.f90

! This is usually done through the "flow_setup" script.

! Adapted on 23 Aug 2015 by dgd from the python code topo_final.py 
! originally written by sek.

! Extended to aperiodic case on 19-20 Sep 2017 by dgd
! ===================================================================

implicit double precision (a-h,o-z)

! Grid dimensions:
integer,parameter:: nx=N_X,ny=N_Y

! Physical domain width (0 < x < ellx) & upper surface location (y = ymax):
double precision,parameter:: ellx=L_X,ymax=L_Y
! The conformal domain has the same width, without loss of generality
double precision,parameter:: tol=1.d-14
! tol: error tolerance in the conformal domain length

 !The number pi (*** non-modifiable ***):
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0

 !Local work arrays:
double precision:: x(0:nx),xp(0:nx),xpo(0:nx),h(0:nx)
double precision:: gh(nx),d(nx)
double precision:: rkx(nx),rky(ny)
double precision:: xo(0:ny,0:nx),yo(0:ny,0:nx)
double precision:: yox(0:ny,0:nx),yoy(0:ny,0:nx)
double precision:: yh0(0:ny),y(0:ny)
double precision:: dco(0:ny),dci(0:ny)

double precision:: trig(2*nx)
integer:: factors(5)

!-------------------------------------------------------------------
write(*,'(2(a,f6.2))') ' We consider a domain of width L_x = ',ellx, &
                           ' and height L_y = ',ymax
write(*,*)
write(*,*) ' Enter the ramp width,  l:'
read(*,*) wramp
write(*,*) ' Enter the ramp height, h:'
read(*,*) hramp
slope=hramp/wramp
write(*,*) ' Enter a smoothing length epsilon:'
read(*,*) eps
x2=ellx
x1=x2-wramp
hm=slope/2.d0
epsi=1.d0/eps

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! Note, quantities relevant to the conformally mapped domain will NOT
! have a "p" appended to the end of their name.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! Initialise regular rectangular transformed grid (x,y):
glx=ellx/dble(nx)
do ix=0,nx
  x(ix)=glx*dble(ix)
enddo

div=1.d0/dble(ny)
do iy=0,ny
  yh0(iy)=1.d0-div*dble(iy)
enddo

! Set up FFTs:
call initfft(nx,factors,trig)

! Define x wavenumbers (these are redefined below as ellx changes):
scx=pi/ellx
do kx=1,nx
  rkx(kx)=scx*dble(kx)
enddo

! Start by choosing the conformal domain height elly so that the area
! of the conformal domain is the same as that of the physical domain:
do ix=0,nx
  h(ix)=hp(xp(ix))
enddo
elly=ymax-(0.5d0*(h(0)+h(nx))+sum(h(1:nx-1)))/dble(nx)
write(*,*)
write(*,'(a,f10.6)') ' Starting with a conformal domain height L_y = ',elly
write(*,*)

!------------------------------------------------------------------
! Iterate to find conformal transform:
xperr=0.5d0
xp=x
do while (xperr .gt. tol .and. xperr .lt. 1.d0)
  !Assign h to hp(xp) and compute average (havg):
  xpo=xp
  do ix=0,nx
    h(ix)=hp(xp(ix))
  enddo
  havg=(0.5d0*(h(0)+h(nx))+sum(h(1:nx-1)))/dble(nx)

  ! Perform cosine Fourier transform of h (overwrites h by the coeffs):
  call dct(1,nx,h,trig,factors)

  ! Multiply by spectral gh operator:
  d=h(1:nx)/tanh(rkx*elly)

  ! Perform (inverse) sine Fourier transform of d (overwrites d):
  call dst(1,nx,d,trig,factors)

  ! Obtain new estimate for xp:
  xp(1:nx)=x(1:nx)-d

  ! Obtain new estimate for L_y:
  elly=ymax-havg

  ! Compute error:
  xperr=sum(abs(xp-xpo))/dble(nx)

  write(*,*) ' Error = ',xperr
enddo

!------------------------------------------------------------------
if (abs(xerr) .le. tol) then

  ! Finalise values of h(x):
  do ix=0,nx
    h(ix)=hp(xp(ix))
  enddo

  ! Write the conformal domain dimensions to domdim.asc:
  open(12,file='domdim.asc',status='replace')
  write(12,*) 0.d0,ellx
  write(12,*) 0.d0,elly
  close(12)

  !------------------------------------------------------------
  ! Compute X and Y throughout the domain and write to a file:

   !Perform cosine Fourier transform of h (overwrites h by the coeffs):
  call dct(1,nx,h,trig,factors)

   !Define the interior semi-spectral fields of X and Y and well
   !as the derivatives Y_x and Y_y:
  xo(:,0)=0.d0
  yo(:,0)=h(0)
  yox(:,0)=0.d0
  yoy(:,0)=0.d0
  do kx=1,nx
    fac=rkx(kx)*elly
    div=h(kx)/sinh(fac)
    dco=div*cosh(fac*yh0)
    dci=div*sinh(fac*yh0)
    xo(:,kx)=-dco
    yo(:,kx)= dci
    yox(:,kx)=-rkx(kx)*dci
    yoy(:,kx)=-rkx(kx)*dco
  enddo
   !Invert using sine/cosine transforms in x:
  call dst(ny+1,nx,xo(:,1:nx),trig,factors)
  call dct(ny+1,nx,yo,trig,factors)
  call dst(ny+1,nx,yox(:,1:nx),trig,factors)
  call dct(ny+1,nx,yoy,trig,factors)
   !Enforce boundary conditions on xo and yox:
  xo(:,0)=0.d0
  xo(:,nx)=0.d0
  yox(:,0)=0.d0
  yox(:,nx)=0.d0

   !Add on the transformed x & y coordinates to the deviations:
  y=elly*(1.0-yh0)
  do ix=0,nx
    xo(:,ix)=xo(:,ix)+x(ix)
    yo(:,ix)=yo(:,ix)+y
    yoy(:,ix)=yoy(:,ix)+1.d0
  enddo

   !Write coordinates X and Y:
  open(20,file='coords.r8',form='unformatted',access='stream',status='replace')
  write(20) xo
  write(20) yo
  close(20)

   !Also write derivatives Y_x and Y_y:
  open(20,file='derivs.r8',form='unformatted',access='stream',status='replace')
  write(20) yox
  write(20) yoy
  close(20)

  !Write basic grid dimensions to params.asc:
  open(10,file='params.asc',status='replace')
  write(10,'(2(1x,i5),1x,f17.13)') nx,ny,elly
  close(10)

else

  write(*,*) ' *** Not converging! ***'

endif

contains 

!=======================================================================

double precision function hp(s)
 !Returns the bottom topography height hp at position s:

double precision:: s

hp=hm*(eps*log(0.5d0*cosh(epsi*(s-x1))/cosh(epsi*(s-x2))**2)-s+2.d0*x2-x1)

return
end function

!=======================================================================
      
end program
