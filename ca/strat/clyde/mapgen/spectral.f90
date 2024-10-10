module spectral

! Module contains routines for inversion and other spectral operations.

use constants
use sta2dfft

 !Common arrays, constants:
double precision:: cecx(0:nx,nym1),cecy(0:ny,nxm1)
double precision:: decx(0:nx,nym1),decy(0:ny,nxm1)
double precision:: xh0(0:nx),xh1(0:nx)
double precision:: yh0(0:ny),yh1(0:ny)

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

 !Arrays related to the conformal map:
double precision:: xori(0:ny,0:nx),yori(0:ny,0:nx)
double precision:: dyoridx(0:ny,0:nx),dyoridy(0:ny,0:nx)

 !For calculating derivatives and integrating:
double precision:: rkx(nx),rky(ny)
double precision:: rkxi(nx),rkyi(ny)

contains

!================================================================

subroutine init_spectral

! Initialises module

implicit none

 !Local variables:
double precision:: fac,scx,scy
double precision:: div,wdi,argm,argp
integer:: ix,iy,kx,ky

!-----------------------------------------------------------------
 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,rkx,rky)

 !Weights for removing corner values of the streamfunction:
fac=one/dble(nx)
do ix=0,nx
  xh1(ix)=fac*dble(ix)
  xh0(ix)=one-xh1(ix)
enddo

fac=one/dble(ny)
do iy=0,ny
  yh1(iy)=fac*dble(iy)
  yh0(iy)=one-yh1(iy)
enddo

 !Define x wavenumbers:
scx=pi/ellx
do kx=1,nx
  rkx(kx)=scx*dble(kx)
enddo

 !Define y wavenumbers:
scy=pi/elly
do ky=1,ny
  rky(ky)=scy*dble(ky)
enddo

 !Define inverses:
rkxi=one/rkx
rkyi=one/rky

!---------------------------------------------------------------
 !Hyperbolic functions used for solutions of Laplace's equation:
do kx=1,nxm1
  fac=rkx(kx)*elly
  div=one/(one-exp(-two*fac))
  wdi=rkx(kx)*div
  do iy=0,ny
    argm=fac*(one-yh1(iy))
    argp=fac*(one+yh1(iy))
    cecy(iy,kx)=(exp(-argm)-exp(-argp))*div
    decy(iy,kx)=(exp(-argm)+exp(-argp))*wdi
  enddo
enddo

do ky=1,nym1
  fac=rky(ky)*ellx
  div=one/(one-exp(-two*fac))
  wdi=rky(ky)*div
  do ix=0,nx
    argm=fac*(one-xh1(ix))
    argp=fac*(one+xh1(ix))
    cecx(ix,ky)=(exp(-argm)-exp(-argp))*div
    decx(ix,ky)=(exp(-argm)+exp(-argp))*wdi
  enddo
enddo

!-------------------------------------------------------------
 !Read in conformal map:
open(11,file='coords.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) xori
read(11,rec=2) yori
close(11)

open(11,file='derivs.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) dyoridx
read(11,rec=2) dyoridy
close(11)

return 
end subroutine

!================================================================

subroutine conform(ybot,ytop)

! Given Y at the bottom and top of the domain (in arrays ybot and 
! ytop), this routine generates the conformal map everywhere in 
! the domain, both Y and X, the latter by using the Cauchy 
! relations X_x = Y_y and X_y = -Y_x.

implicit none

 !Passed arrays and variables:
double precision:: ybot(0:nx),ytop(0:nx)

 !Local arrays and variables:
double precision:: abot(nx),atop(nx)
double precision:: yh(0:ny,nx),yy(0:ny,nx)
double precision:: xc(0:ny),yc(0:nx)
double precision:: ylb,yrb,dycdx,xcor,yinc
integer:: ix,iy,kx,ky

!------------------------------------------------------------------
 ! Use corner values of Y in order to construct residual boundary 
 ! values which are zero at the corners:
ylb=ybot(0)
yrb=ybot(nx)
 ! Note, in doing this we define a solution to Laplace's equation
 !                 Y_c = ylb*(1-xi)+yrb*xi
 ! where  xi = (x-x_min)/L_x = xh1 in the code & 1-xi = xh0.
 ! Note also, ytop(0)=ybot(0) and ytop(nx)=ybot(nx).

 ! Define Y_c:
yc=ylb*xh0+yrb*xh1

 ! Obtain dY_c/dx (constant):
dycdx=(yrb-ylb)/ellx

!------------------------------------------------------------------
 ! Define Y - Y_c on lower & upper boundaries (use abot & atop):
do ix=1,nxm1
  abot(ix)=ybot(ix)-yc(ix)
  atop(ix)=ytop(ix)-yc(ix)
enddo
 ! By construction, each of these functions vanish at ix = 0 & nx.
 ! These arrays only need to be defined where they are non-zero.

!------------------------------------------------------------------
 ! Do a sine transform of abot at y = ymin and atop at y = ymax
 ! to construct half of the homogeneous solution Y_h which, when
 ! added to Y_c above, has the correct boundary conditions in y
 ! (Y = ybot & ytop at ymin & ymax):
call dst(1,nx,abot,xtrig,xfactors)
call dst(1,nx,atop,xtrig,xfactors)

 ! Define the interior semi-spectral field and its spatial derivatives:
do kx=1,nxm1
  do iy=0,ny
    yh(iy,kx)=atop(kx)*cecy(iy,kx)+abot(kx)*cecy(ny-iy,kx)
    yy(iy,kx)=atop(kx)*decy(iy,kx)-abot(kx)*decy(ny-iy,kx)
    dyoridx(iy,kx)=rkx(kx)*yh(iy,kx)
    xori(iy,kx)= -rkxi(kx)*yy(iy,kx)
  enddo
enddo
dyoridx(:,0)=zero
dyoridx(:,nx)=zero
xori(:,0)=zero
xori(:,nx)=zero

 ! Invert yh & yy using a sine transform:
call dst(nyp1,nx,yh,xtrig,xfactors)
call dst(nyp1,nx,yy,xtrig,xfactors)
 ! Invert dyoridx & xori using a cosine transform:
call dct(nyp1,nx,dyoridx,xtrig,xfactors)
call dct(nyp1,nx,xori,xtrig,xfactors)

!------------------------------------------------------------------
 ! Add Y_c and Y_h to obtain the final map yori:
yori(:, 0)=ylb
yori(:,nx)=yrb
yori(0 ,:)=ybot
yori(ny,:)=ytop
do ix=1,nxm1
  do iy=1,nym1
    yori(iy,ix)=yc(ix)+yh(iy,ix)
  enddo
enddo

 ! Add X_c and X_h to obtain the final map xori after ensuring
 ! xori = 0 in lower left corner:
xcor=-xori(0,0)
yinc=elly*dycdx
xc=xcor-yinc*yh1
do ix=0,nx
  xori(:,ix)=xori(:,ix)+xc
enddo

 ! Also obtain the x and y derivatives of Y:
dyoridx=dyoridx+dycdx
dyoridy(:, 0)=zero
dyoridy(:,nx)=zero
do ix=1,nxm1
  dyoridy(:,ix)=yy(:,ix)
enddo

return
end subroutine

!============================================================
end module spectral
