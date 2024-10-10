module spectral

! Module contains routines for inversion and other spectral operations.

use constants
use sta2dfft

 !Common arrays, constants:
double precision:: p0(0:ny,0:nx),u0(0:ny),v0(0:nx)
double precision:: cecx(0:nx,nym1),cecy(0:ny,nxm1)
double precision:: decx(0:nx,nym1),decy(0:ny,nxm1)
double precision:: xh0(0:nx),xh1(0:nx)
double precision:: yh0(0:ny),yh1(0:ny)
double precision:: green(0:nx,0:ny),diss(0:nx,0:ny)
double precision:: dafx(0:nx),dafy(0:ny) 

double precision:: xtrig(2*nx),ytrig(2*ny)
integer:: xfactors(5),yfactors(5)

 !Arrays related to the conformal map:
double precision:: xori(0:ny,0:nx),yori(0:ny,0:nx)
double precision:: dyoridx(0:ny,0:nx),dyoridy(0:ny,0:nx)
double precision:: confac(0:ny,0:nx),confaci(0:ny,0:nx)
double precision:: xsea,ysea,ybsea,xriv,yriv,ybriv,xwsea,xwriv
integer:: isea,iriv

 !For calculating derivatives:
double precision:: rkx(nx),rky(ny)

contains

!================================================================

subroutine init_spectral(bbdif)

! Initialises module; bbdif = bb_max - bb_min (fixed in time).

implicit none

 !Passed variable
double precision:: bbdif

 !Local variables:
double precision:: fac,scx,scy,rkxmax,rkymax
double precision:: delk,delki,snorm,div,wdi,argm,argp,visc
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
    p0(iy,ix)=-f14*(ellx**2*xh0(ix)*xh1(ix)+elly**2*yh0(iy)*yh1(iy))
  enddo
enddo

 !Do same for the x and y derivatives:
do ix=0,nx
  v0(ix)=hlx*(xh1(ix)-f12)
enddo
do iy=0,ny
  u0(iy)=hly*(yh1(iy)-f12)
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

!---------------------------------------------------------
 !Define dissipation operator:
if (nnu .eq. 1) then 
  visc=prediss*sqrt(bbdif)/(rkxmax**2+rkymax**2)**0.75d0
   !Note: bbdif = bb_max - bb_min
else
  visc=prediss/(rkxmax**2+rkymax**2)**nnu
endif

 !Write viscosity to log file:
write(*,*)
write(*,'(a,1x,1p,e14.7)') ' viscosity = ',visc
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

!-------------------------------------------------------------
 !Read in conformal map:
open(11,file='coords.r8',form='unformatted', &
    & access='direct',status='old',recl=2*(nbytes-4))
read(11,rec=1) xori
read(11,rec=2) yori
close(11)

open(11,file='derivs.r8',form='unformatted', &
    & access='direct',status='old',recl=2*(nbytes-4))
read(11,rec=1) dyoridx
read(11,rec=2) dyoridy
close(11)

 !Define conformal factor and rescale derivatives Y_x and Y_y
 !for use in calculating grad(b) in conformal space:
confac=dyoridx**2+dyoridy**2
confaci=one/confac      !lambda in the notes
dyoridx=dyoridx*confaci
dyoridy=dyoridy*confaci

!-------------------------------------------------------------
 !Determine grid points lying along the left (sea) edge:
xsea=xori(0,0)
ysea=yori(0,0)
isea=1
do
  if (xori(0,isea) .gt. xsea+1.d-7) exit
  isea=isea+1
enddo
isea=isea-1
ybsea=yori(0,isea)

 !Determine grid points lying along the right (river) edge:
xriv=xori(0,nx)
yriv=yori(0,nx)
iriv=nx-1
do
  if (xori(0,iriv) .lt. xriv-1.d-7) exit
  iriv=iriv-1
enddo
iriv=iriv+1
ybriv=yori(0,iriv)

 !Determine x coordinates on either edge of the weir:
xwsea=xori(ny,0)
xwriv=xori(ny,nx)

return 
end subroutine

!================================================================

subroutine main_invert(zz,uu,vv,vsea,uriv,vriv)

! Main routine for computing the velocity field (uu,vv) from the
! vorticity field (zz) as well as the river and sea flow speeds.

implicit none

 !Passed arrays and variables:
double precision:: zz(0:ny,0:nx),uu(0:ny,0:nx),vv(0:ny,0:nx)
double precision:: vsea,uriv,vriv

 !Local arrays and variables:
double precision:: ss(0:nx,0:ny),pp(0:ny,0:nx)
double precision:: ppx(nx,0:ny),ppy(0:nx,ny)
double precision:: vvi(0:ny,nx),uui(ny,0:nx)
double precision:: pbot(0:nx),ptop(0:nx),plft(0:ny),prgt(0:ny)
double precision:: abot(nx),atop(nx),alft(ny),argt(ny)
double precision:: phbt(0:ny,nx),pybt(0:ny,nx),pxbt(0:ny,0:nx)
double precision:: phlr(0:nx,ny),pxlr(0:nx,ny),pylr(0:nx,0:ny)
double precision:: uc(0:nx),vc(0:ny)
double precision:: usea,cw,z0,plb,prb,plt,prt
double precision:: dprgt,dplft,dptop,dpbot
integer:: ix,iy,kx,ky

!-----------------------------------------------------------------------
 !Determine usea by mass (volume) conservation:
usea=(uriv*(yriv-ybriv)+vriv*(xriv-xwriv)+vsea*(xwsea-xsea))/(ysea-ybsea)

 !Set streamfunction on boundaries (called "psi_bar" in comments below):
cw=-uriv*(yriv-ybriv)-vriv*(xriv-xwriv)

 !The "bottom" of the conformal domain consists of the sea edge
 !from ix = 0 to isea, the bottom from isea+1 to iriv-1, and the
 !river edge from iriv to nx.  The streamfunction psi is zero along
 !the bottom:
do ix=0,isea
  pbot(ix)=usea*(ybsea-yori(0,ix))
enddo
do ix=isea+1,iriv-1
  pbot(ix)=zero
enddo
do ix=iriv,nx
  pbot(ix)=uriv*(ybriv-yori(0,ix))
enddo

 !The "top" of the conformal domain is the weir where psi = constant:
do ix=0,nx
  ptop(ix)=cw
enddo

 !The "left edge" of the conformal domain is the sea surface:
do iy=0,ny
  plft(iy)=vsea*(xori(iy,0)-xwsea)+cw
enddo

 !The "right edge" of the conformal domain is the river surface:
do iy=0,ny
  prgt(iy)=vriv*(xori(iy,nx)-xwriv)+cw
enddo

!-----------------------------------------------------------------
 !Apply conformal factor to vorticity (store in pp momentarily):
pp=zz*confac

 !Solve for psi_p, the particular solution (store result in pp):

 !(1) compute mean "vorticity" (z0) over conformal domain:
z0=(f14*(pp(0,0)+pp(0,nx)+pp(ny,0)+pp(ny,nx)) &
   +f12*(sum(pp(1:nym1,0)+pp(1:nym1,nx))+sum(pp(0,1:nxm1)+pp(ny,1:nxm1))) &
   +sum(pp(1:nym1,1:nxm1)))*dsumi

 !(2) Remove mean vorticity (pp then contains zz - z0):
pp=pp-z0
 !Note: this is necessary as the mean part cannot be inverted spectrally.
 
 !(3) FFT pp and invert to get psi_p - psi_0:
call ptospc_cc(nx,ny,pp,ss,xfactors,yfactors,xtrig,ytrig)
ss=green*ss ! Solves Poisson's equation
 !Get x and y derivatives for use in constructing u_p & v_p below:
call xderiv_cc(nx,ny,rkx,ss,ppx)
call yderiv_cc(nx,ny,rky,ss,ppy)
 !Transform fields back to physical space:
call spctop_cc(nx,ny,ss,pp,xfactors,yfactors,xtrig,ytrig)
call spctop_sc(nx,ny,ppx,vvi,xfactors,yfactors,xtrig,ytrig)
call spctop_cs(nx,ny,ppy,uui,xfactors,yfactors,xtrig,ytrig)

 !(4) Add back psi_0 (due to mean vorticity) to define psi_p:
pp=pp+z0*p0
 !Note: p0 is the streamfunction for unit mean vorticity.

 ! Do same for the x and y derivatives for velocity components:
do iy=0,ny
  vv(iy, 0)=z0*v0(0)
  vv(iy,nx)=z0*v0(nx)
enddo
do ix=1,nxm1
  do iy=0,ny
    vv(iy,ix)=vvi(iy,ix)+z0*v0(ix)
  enddo
enddo
 ! vv now contains v_p.

do ix=0,nx
  uu(0,ix)=z0*u0(0)
  do iy=1,nym1
    uu(iy,ix)=uui(iy,ix)+z0*u0(iy)
  enddo
  uu(ny,ix)=z0*u0(ny)
enddo
 ! uu now contains -u_p; the sign of this component is changed at the end.

 !(5) Determine corner values of psi' = psi_bar - psi_p in order to
 !    construct residual boundary values which are zero at the corners:
plb=pbot(0) -pp( 0, 0)
prb=pbot(nx)-pp( 0,nx)
plt=ptop(0) -pp(ny, 0)
prt=ptop(nx)-pp(ny,nx)
 ! Note, in doing this we define another solution to Laplace's equation
 ! psi_c = plb*(1-xi)*(1-eta)+prb*xi*(1-eta)+plt*(1-xi)*eta+prt*xi*eta
 ! where  xi = (x-x_min)/L_x = xh1 in the code, 1- xi = xh0 in the code,
 ! while eta = (y-y_min)/L_y = yh1 in the code, 1-eta = yh0 in the code.

 ! Obtain dpsi_c/dx & dpsi_c/dy (vc and uc below):
dprgt=(prt-prb)/elly
dplft=(plt-plb)/elly
dptop=(prt-plt)/ellx
dpbot=(prb-plb)/ellx
do iy=0,ny
  vc(iy)=dpbot*yh0(iy)+dptop*yh1(iy)
enddo
do ix=0,nx
  uc(ix)=dplft*xh0(ix)+dprgt*xh1(ix)
enddo

 ! Define psi_tilde = psi' - psi_c on boundaries (use abot, atop ...):
do ix=1,nxm1
  abot(ix)=pbot(ix)-pp( 0,ix)-plb*xh0(ix)-prb*xh1(ix)
  atop(ix)=ptop(ix)-pp(ny,ix)-plt*xh0(ix)-prt*xh1(ix)
enddo
do iy=1,nym1
  alft(iy)=plft(iy)-pp(iy, 0)-plb*yh0(iy)-plt*yh1(iy)
  argt(iy)=prgt(iy)-pp(iy,nx)-prb*yh0(iy)-prt*yh1(iy)
enddo
 ! By construction, each of these functions vanish at the corners
 ! (ix,iy) = (0,0), (nx,0), (0,ny) and (nx,ny).  These arrays only
 ! need to be defined where they are non-zero.

 !(6) Do a sine transform of abot at y = ymin and atop at y = ymax
 !    to construct half of the homogeneous solution psi_h which, when
 !    added to psi_c and psi_bar above, has the correct boundary 
 !    conditions in y (psi = psi_bar = pbot & ptop at ymin & ymax):
call dst(1,nx,abot,xtrig,xfactors)
call dst(1,nx,atop,xtrig,xfactors)

 !Define the interior semi-spectral field and its spatial derivatives:
do kx=1,nxm1
  do iy=0,ny
    phbt(iy,kx)=atop(kx)*cecy(iy,kx)+abot(kx)*cecy(ny-iy,kx)
    pybt(iy,kx)=atop(kx)*decy(iy,kx)-abot(kx)*decy(ny-iy,kx)
    pxbt(iy,kx)=rkx(kx)*phbt(iy,kx)
  enddo
enddo
pxbt(:,0)=zero

 !Invert phbt & pybt using a sine transform:
call dst(nyp1,nx,phbt,xtrig,xfactors)
call dst(nyp1,nx,pybt,xtrig,xfactors)
 !Invert pxbt using a cosine transform:
call dct(nyp1,nx,pxbt,xtrig,xfactors)

 !(7) Do a sine transform of alft at x = xmin and argt at x = xmax
 !    to construct half of the homogeneous solution psi_h which, when
 !    added to psi_c and psi_bar above, has the correct boundary 
 !    conditions in x (psi = psi_bar = plft & prgt at xmin & xmax):
call dst(1,ny,alft,ytrig,yfactors)
call dst(1,ny,argt,ytrig,yfactors)

 !Define the interior semi-spectral field and its spatial derivatives:
do ky=1,nym1
  do ix=0,nx
    phlr(ix,ky)=argt(ky)*cecx(ix,ky)+alft(ky)*cecx(nx-ix,ky)
    pxlr(ix,ky)=argt(ky)*decx(ix,ky)-alft(ky)*decx(nx-ix,ky)
    pylr(ix,ky)=rky(ky)*phlr(ix,ky)
  enddo
enddo
pylr(:,0)=zero

 !Invert phlr & pxlr using a sine transform:
call dst(nxp1,ny,phlr,ytrig,yfactors)
call dst(nxp1,ny,pxlr,ytrig,yfactors)
 !Invert pylr using a cosine transform:
call dct(nxp1,ny,pylr,ytrig,yfactors)

 !(8) Add phlr and phbt to obtain the final streamfunction pp:
pp(:, 0)=plft
pp(:,nx)=prgt
pp(0 ,:)=pbot
pp(ny,:)=ptop
do ix=1,nxm1
  do iy=1,nym1
    pp(iy,ix)=pp(iy,ix)+phlr(ix,iy)+phbt(iy,ix)+ &  !psi_p + psi_h
              (plb*xh0(ix)+prb*xh1(ix))*yh0(iy)+ &  !psi_c
              (plt*xh0(ix)+prt*xh1(ix))*yh1(iy)
  enddo
enddo

 !(9) Finally, obtain the velocity field:

 ! d(psi_p + psi_h + psi_c)/dx (use pxlr = 0 when iy = 0 & ny):
do ix=0,nx
  vv(0, ix)=vv(0 ,ix)+pxbt(0 ,ix)+vc(0)
  do iy=1,nym1
    vv(iy,ix)=vv(iy,ix)+pxlr(ix,iy)+pxbt(iy,ix)+vc(iy)
  enddo
  vv(ny,ix)=vv(ny,ix)+pxbt(ny,ix)+vc(ny)
enddo

 ! d(psi_p + psi_h + psi_c)/dy (use pybt = 0 when ix = 0 & nx):
do iy=0,ny
  uu(iy, 0)=uu(iy, 0)+pylr(0 ,iy)+uc(0)
  uu(iy,nx)=uu(iy,nx)+pylr(nx,iy)+uc(nx)
enddo
do ix=1,nxm1
  do iy=0,ny
    uu(iy,ix)=uu(iy,ix)+pylr(ix,iy)+pybt(iy,ix)+uc(ix)
  enddo
enddo

 ! Apply conformal factor lambda (= confaci):
uu=-confaci*uu
vv= confaci*vv

! Insert prescribed normal flow at edges:
do ix=0,isea
  vv(0,ix)=-usea*dyoridx(0,ix)
enddo
do ix=isea+1,iriv-1
  vv(0,ix)=zero
enddo
do ix=iriv,nx
  vv(0,ix)=-uriv*dyoridx(0,ix)
enddo

do ix=0,nx
  vv(ny,ix)=zero
enddo

do iy=0,ny
  uu(iy, 0)=vsea*dyoridx(iy, 0)
  uu(iy,nx)=vriv*dyoridx(iy,nx)
enddo

return
end subroutine

!============================================================
end module spectral
