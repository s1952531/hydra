module spectral

use constants
use sta2dfft

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Common arrays, constants:
double precision:: pbar(0:ny,0:nx)
double precision:: decx(nxm1,nym1),decy(nym1,nxm1)
double precision:: xh0(0:nx),xh1(0:nx)
double precision:: yh0(0:ny),yh1(0:ny)
double precision:: green(0:nx,0:ny),diss(0:nx,0:ny)
double precision:: dafx(0:nx),dafy(0:ny) 

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

 !Arrays related to the conformal map:
double precision:: xo(0:ny,0:nx),yo(0:ny,0:nx)
double precision:: yox(0:ny,0:nx),yoy(0:ny,0:nx)
double precision:: confac(0:ny,0:nx),confaci(0:ny,0:nx)

double precision:: rkx(nx),rky(ny)

contains

!================================================================

subroutine init_spectral

implicit none

 !Local variables:
double precision:: fac,scx,scy,rkxmax,rkymax
double precision:: delk,delki,snorm,div,argm,argp,visc
integer:: ix,iy,kx,ky,k

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

 !Define part of streamfunction proportional to the mean vorticity:
do ix=0,nx
  do iy=0,ny
    pbar(iy,ix)=-f14*(ellx**2*xh0(ix)*xh1(ix)+elly**2*yh0(iy)*yh1(iy))
  enddo
enddo

 !Define x wavenumbers:
scx=pi/ellx
rkxmax=scx*dble(nx)
do kx=1,nx
  rkx(kx)=scx*dble(kx)
enddo

 !Define y wavenumbers:
scy=pi/elly
rkymax=scy*dble(ny)
do ky=1,ny
  rky(ky)=scy*dble(ky)
enddo

!--------------------------------------------------------
 !Define de-aliasing filter (2/3 rule):
dafx(0)=one
do kx=1,nx
  if (rkx(kx) .lt. f23*rkxmax) then
    dafx(kx)=one
  else
    dafx(kx)=zero
  endif
enddo

dafy(0)=one
do ky=1,ny
  if (rky(ky) .lt. f23*rkymax) then
    dafy(ky)=one
  else
    dafy(ky)=zero
  endif
enddo
 
!-----------------------------------------------------------------
 !Define Green function for solving Poisson's equation spectrally:
green(0,0)=zero
do kx=1,nx
  green(kx,0)=-one/rkx(kx)**2
enddo
do ky=1,ny
  green(0,ky)=-one/rky(ky)**2
enddo
do ky=1,ny
  do kx=1,nx
    green(kx,ky)=-one/(rkx(kx)**2+rky(ky)**2)
  enddo
enddo

!---------------------------------------------------------------
 !Hyperbolic functions used for solutions of Laplace's equation:
do kx=1,nxm1
  fac=rkx(kx)*elly
  div=one/(one-exp(-two*fac))
  do iy=1,nym1
    argm=fac*(one-yh1(iy))
    argp=fac*(one+yh1(iy))
    decy(iy,kx)=(exp(-argm)-exp(-argp))*div
  enddo
enddo

do ky=1,nym1
  fac=rky(ky)*ellx
  div=one/(one-exp(-two*fac))
  do ix=1,nxm1
    argm=fac*(one-xh1(ix))
    argp=fac*(one+xh1(ix))
    decx(ix,ky)=(exp(-argm)-exp(-argp))*div
  enddo
enddo

!---------------------------------------------------------
 !Read in conformal map:
open(11,file='coords.r8',form='unformatted',access='stream',status='old')
read(11) xo
read(11) yo
close(11)

open(11,file='derivs.r8',form='unformatted',access='stream',status='old')
read(11) yox
read(11) yoy
close(11)

 !Define conformal factor for Poisson problem, and rescale derivatives 
 !of transform for use in calculating grad(b) in conformal space:
confac=yox**2+yoy**2
confaci=one/confac
yox=yox*confaci
yoy=yoy*confaci

!---------------------------------------------------------
 !Define hyperdiffusion:
visc=prediss/(rkxmax**2+rkymax**2)**nnu

 !Write viscosity to log file:
write(*,*)
write(*,'(a,1x,1p,e14.7)') ' hyperviscosity = ',visc
write(*,*)

diss(0,0)=zero
do kx=1,nx
  diss(kx,0)=f12*visc*rkx(kx)**(2*nnu)
enddo
do ky=1,ny
  diss(0,ky)=f12*visc*rky(ky)**(2*nnu)
enddo
do ky=1,ny
  do kx=1,nx
    diss(kx,ky)=f12*visc*(rkx(kx)**2+rky(ky)**2)**nnu
  enddo
enddo

return 
end subroutine

!================================================================

subroutine main_invert(zz,uu,vv)
! Main routine for computing the velocity field (uu,vv) from the
! vorticity field (zz).

implicit none

 !Passed arrays:
double precision:: zz(0:ny,0:nx),uu(0:ny,0:nx),vv(0:ny,0:nx)

 !Local arrays and variables:
double precision:: ss(0:nx,0:ny),pp(0:ny,0:nx)
double precision:: pbot(nx),ptop(nx),cppy(nym1,nx)
double precision:: plft(ny),prgt(ny),cppx(nxm1,ny)
double precision:: zbar,sw00,sw10,sw01,sw11
integer:: ix,iy,kx,ky

!-----------------------------------------------------------------
 !Apply conformal factor to vorticity (store in pp momentarily):
pp=zz*confac

 !Solve for psi (pp):

 !(1) compute mean "vorticity" (zbar) over conformal domain:
zbar=zero
do ix=1,nxm1
  zbar=zbar+pp(0,ix)+pp(ny,ix)
enddo
do iy=1,nym1
  zbar=zbar+pp(iy,0)+pp(iy,nx)
enddo
zbar=(f12*zbar+f14*(pp(0,0)+pp(ny,0)+pp(0,nx)+pp(ny,nx))+ &
      sum(pp(1:nym1,1:nxm1)))/dble(nx*ny)

 !(2) Remove mean vorticity from pp:
pp=pp-zbar
 
 !(3) FFT pp and invert to get uncorrected streamfunction pp:
call ptospc_cc(nx,ny,pp,ss,xfactors,yfactors,xtrig,ytrig)
ss=green*ss
call spctop_cc(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)

 !(4) Add part of pp due to mean vorticity:
pp=pp+zbar*pbar

 !(5) Remove a bi-linear function so that pp is zero at the corners:
sw00=pp(0,0)
sw10=pp(ny,0)
sw01=pp(0,nx)
sw11=pp(ny,nx)
do ix=0,nx
  do iy=0,ny
    pp(iy,ix)=pp(iy,ix)-(sw00*xh0(ix)+sw01*xh1(ix))*yh0(iy) &
                       -(sw10*xh0(ix)+sw11*xh1(ix))*yh1(iy)
  enddo
enddo
 !Note:  xh0 = (xmax - x)/ellx, xh1 = (x - xmin)/ellx etc.

 !(6) Do a sine transform of pp at y = ymin and ymax and obtain the
 !    interior field (cppy) that must be subtracted to give pp = 0
 !    at y = ymin and ymax:
do ix=1,nxm1
  pbot(ix)=pp(0,ix)
  ptop(ix)=pp(ny,ix)
enddo
call dst(1,nx,pbot,xtrig,xfactors)
call dst(1,nx,ptop,xtrig,xfactors)

 !Define the interior semi-spectral field:
do kx=1,nxm1
  do iy=1,nym1
    cppy(iy,kx)=pbot(kx)*decy(ny-iy,kx)+ptop(kx)*decy(iy,kx)
  enddo
enddo
 !Invert using a sine transform:
call dst(nym1,nx,cppy,xtrig,xfactors)

 !(7) Do a sine transform of pp at x = xmin and xmax and obtain the
 !    interior field (cppx) that must be subtracted to give pp = 0
 !    at x = xmin and xmax:
do iy=1,nym1
  plft(iy)=pp(iy,0)
  prgt(iy)=pp(iy,nx)
enddo
call dst(1,ny,plft,ytrig,yfactors)
call dst(1,ny,prgt,ytrig,yfactors)

 !Define the interior semi-spectral field:
do ky=1,nym1
  do ix=1,nxm1
    cppx(ix,ky)=plft(ky)*decx(nx-ix,ky)+prgt(ky)*decx(ix,ky)
  enddo
enddo
 !Invert using a sine transform:
call dst(nxm1,ny,cppx,ytrig,yfactors)

 !(8) Remove cppx and cppy to obtain the final streamfunction pp:
pp(:, 0)=zero
pp(:,nx)=zero
pp(0 ,:)=zero
pp(ny,:)=zero
do ix=1,nxm1
  do iy=1,nym1
    pp(iy,ix)=pp(iy,ix)-cppx(ix,iy)-cppy(iy,ix)
  enddo
enddo

 !(9) Obtain the velocity field by differentiation:
call getvel(pp,uu,vv)

return
end subroutine

!================================================================

subroutine getvel(pp,uu,vv)
! Computes the velocity components uu & vv from the streamfunction
! pp via uu = -lambda^{-1}d(pp)/dy and vv = lambda^{-1}d(pp)/dx,
! where lambda is the conformal factor (confac).
! *** pp, uu & vv are all in physical space
! *** and include the domain edges.

implicit none

 !Passed arrays:
double precision:: pp(0:ny,0:nx),uu(0:ny,0:nx),vv(0:ny,0:nx)

 !Local arrays and variables:
double precision:: ppi(ny,nx),pps(nx,ny)
double precision:: ppx(0:nx,ny),vvi(ny,0:nx)
double precision:: ppy(nx,0:ny),uui(0:ny,nx)
integer:: ix,iy,kx,ky

!-----------------------------------------------------------------
 !Copy non-zero interior values of pp to ppi:
do ix=1,nxm1
  do iy=1,nym1
    ppi(iy,ix)=pp(iy,ix)
  enddo
enddo

 !Transform ppi to spectral space:
call ptospc_ss(nx,ny,ppi,pps,xfactors,yfactors,xtrig,ytrig)

 !Apply de-aliasing filter:
do ky=1,ny
  do kx=1,nx
    pps(kx,ky)=pps(kx,ky)*dafx(kx)*dafy(ky)
  enddo
enddo

 !Compute d(ppi)/dx = ppx spectrally:
call xderiv_ss(nx,ny,rkx,pps,ppx)

 !Transform ppx back to physical space as vvi:
call spctop_cs(nx,ny,ppx,vvi,xfactors,yfactors,xtrig,ytrig)

 !Copy vvi into vv after applying the conformal factor and add on 
 !zero edge values at iy = 0 & ny:
do ix=0,nx
  vv(0,ix)=zero
  do iy=1,nym1
    vv(iy,ix)=confaci(iy,ix)*vvi(iy,ix)
  enddo
  vv(ny,ix)=zero
enddo

 !Compute d(ppi)/dy = ppy spectrally:
call yderiv_ss(nx,ny,rky,pps,ppy)

 !Transform ppy back to physical space as uui:
call spctop_sc(nx,ny,ppy,uui,xfactors,yfactors,xtrig,ytrig)

 !Copy -uui into uu after applying the conformal factor and add on 
 !zero edge values at ix = 0 & nx:
do ix=1,nxm1
  do iy=0,ny
    uu(iy,ix)=-confaci(iy,ix)*uui(iy,ix)
  enddo
enddo
do iy=0,ny
  uu(iy, 0)=zero
  uu(iy,nx)=zero
enddo

return
end subroutine

!============================================================
end module
