module spectral

! Module containing subroutines for spectral operations, inversion, etc.

use constants
use sta2dfft

 !Common arrays, constants:
double precision:: opak(ng,ng),helm(ng,ng),rlap(ng,ng)
double precision:: c2g2(ng,ng),simp(ng,ng),rdis(ng,ng)
double precision:: filt(ng,ng),diss(ng,ng)

 !For 2D FFTs:
double precision:: hrkx(ng),hrky(ng),rk(ng)
double precision:: xtrig(2*ng),ytrig(2*ng)
integer:: xfactors(5),yfactors(5)

double precision:: spmf(0:ng),alk(ng)
integer:: kmag(ng,ng),kmax,kmaxred


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================
subroutine init_spectral
! Initialises this module

implicit none

!Local variables:
double precision:: rkmax,rks,snorm
double precision:: anu,rkfsq,fsq
integer:: kx,ky,k

!----------------------------------------------------------------------
 !Set up 2D FFTs:
call init2dfft(ng,ng,twopi,twopi,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define wavenumbers and filtered wavenumbers:
rk(1)=zero
do k=1,ng/2-1
  rk(k+1)   =hrkx(2*k)
  rk(ng+1-k)=hrkx(2*k)
enddo
rk(ng/2+1)=hrkx(ng)

!-----------------------------------------------------------------------
 !Initialise arrays for computing the spectrum of any field:
rkmax=dble(ng/2)
kmax=nint(rkmax*sqrt(two))
do k=0,kmax
  spmf(k)=zero
enddo
do ky=1,ng
  do kx=1,ng
    k=nint(sqrt(rk(kx)**2+rk(ky)**2))
    kmag(kx,ky)=k
    spmf(k)=spmf(k)+one
  enddo
enddo
 !Compute spectrum multiplication factor (spmf) to account for unevenly
 !sampled shells and normalise spectra by 8/(ng*ng) so that the sum
 !of the spectrum is equal to the L2 norm of the original field:
snorm=four*pi/dble(ng*ng)
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(dble(k))
enddo
 !Only output shells which are fully occupied (k <= kmaxred):
kmaxred=ng/2

!-----------------------------------------------------------------------
 !Define a variety of spectral operators:

 !Hyperviscosity coefficient (Dritschel, Gottwald & Oliver, JFM (2017)):
anu=cdamp*cof/rkmax**(2*nnu)
 !Assumes Burger number = 1.

 !Used for de-aliasing filter below:
rkfsq=(dble(ng)/3.d0)**2

 !Squared Coriolis frequency:
fsq=cof**2

do ky=1,ng
  do kx=1,ng
    rks=rk(kx)**2+rk(ky)**2
     !Spectral c^2*grad^2 - f^2 operator:
    opak(kx,ky)=-(fsq+csq*rks)
     !Hyperviscous operator:
    diss(kx,ky)=anu*rks**nnu
     !De-aliasing filter:
    if (rks .gt. rkfsq) then
      filt(kx,ky)=zero
      c2g2(kx,ky)=zero
      rlap(kx,ky)=zero
      helm(kx,ky)=zero
      rdis(kx,ky)=zero
    else
      filt(kx,ky)=one
       !c^2*grad^2:
      c2g2(kx,ky)=-csq*rks
       !grad^{-2}:
      rlap(kx,ky)=-one/(rks+1.d-20)
       !(c^2*grad^2 - f^2)^{-1}:
      helm(kx,ky)=one/opak(kx,ky)
       !R operator in paper:
      rdis(kx,ky)=dt2i+diss(kx,ky)
    endif
     !Semi-implicit operator for delta and gamma_l:
    simp(kx,ky)=one/((dt2i+diss(kx,ky))**2-opak(kx,ky))
     !Re-define damping operator for use in q_l evolution:
    diss(kx,ky)=two/(one+dt2*diss(kx,ky))
  enddo
enddo

 !Ensure mean potentials and height anomaly remain zero:
rlap(1,1)=zero

return 
end subroutine

!======================================================================
subroutine main_invert(qs,ds,gs,hh,uu,vv,zz)
! Given the PV anomaly qs, divergence ds and acceleration divergence gs
! (all in spectral space), this routine computes the dimensionless depth 
! anomaly hh and the velocity field (uu,vv) in physical space.  It also 
! returns the relative vorticity (zz) in physical space.

implicit none

 !Passed variables:
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng)           !Spectral
double precision:: hh(ng,ng),uu(ng,ng),vv(ng,ng),zz(ng,ng) !Physical

 !Local variables:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng)
double precision:: wke(ng,ng),wkf(ng,ng),wkg(ng,ng)
double precision:: uio,vio

!--------------------------------------------------------------
 !Invert linearised PV for the depth anomaly in spectral space:
wka=helm*(cof*qs-gs)

 !Define relative vorticity, zeta:
wkb=qs+cof*wka

 !Invert Laplace operator on zeta & delta to define velocity:
wkc=rlap*wkb
wkd=rlap*ds

 !Calculate derivatives spectrally:
call xderiv(ng,ng,hrkx,wkd,wke)
call yderiv(ng,ng,hrky,wkd,wkf)
call xderiv(ng,ng,hrkx,wkc,wkd)
call yderiv(ng,ng,hrkx,wkc,wkg)

 !Define velocity components:
wke=wke-wkg
wkf=wkf+wkd

 !Bring quantities back to physical space:
call spctop(ng,ng,wka,hh,xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,wkb,zz,xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,wke,uu,xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,wkf,vv,xfactors,yfactors,xtrig,ytrig)

 !Add mean flow (uio,vio):
uio=-sum(hh*uu)*dsumi
vio=-sum(hh*vv)*dsumi
uu=uu+uio
vv=vv+vio

return
end subroutine

!=================================================================
subroutine gradient(ff,ffx,ffy)
! Computes the gradient ffx = dF/dx & ffy = dF/dy of a field F.
! *** ff is in spectral space whereas (ffx,ffy) are in physical space

implicit none

 !Passed arrays:
double precision:: ff(ng,ng)             !Spectral
double precision:: ffx(ng,ng),ffy(ng,ng) !Physical

 !Local array:
double precision:: vtmp(ng,ng)           !Spectral

 !Get derivatives of F:
call xderiv(ng,ng,hrkx,ff,vtmp)
call spctop(ng,ng,vtmp,ffx,xfactors,yfactors,xtrig,ytrig)

call yderiv(ng,ng,hrky,ff,vtmp)
call spctop(ng,ng,vtmp,ffy,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!=================================================================
subroutine jacob(aa,bb,cc)
! Computes the Jacobian of aa and bb and returns it in cc.
! All passed variables are in physical space.

implicit none

 !Passed arrays:
double precision:: aa(ng,ng),bb(ng,ng),cc(ng,ng)           !Physical

 !Work arrays:
double precision:: ax(ng,ng),ay(ng,ng),bx(ng,ng),by(ng,ng) !Physical
double precision:: wka(ng,ng),wkb(ng,ng)                   !Spectral

!---------------------------------------------------------
cc=aa
call ptospc(ng,ng,cc,wka,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
wka=filt*wka
 !Get derivatives of aa:
call xderiv(ng,ng,hrkx,wka,wkb)
call spctop(ng,ng,wkb,ax,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,wkb)
call spctop(ng,ng,wkb,ay,xfactors,yfactors,xtrig,ytrig)

cc=bb
call ptospc(ng,ng,cc,wka,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
wka=filt*wka
 !Get derivatives of bb:
call xderiv(ng,ng,hrkx,wka,wkb)
call spctop(ng,ng,wkb,bx,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,wkb)
call spctop(ng,ng,wkb,by,xfactors,yfactors,xtrig,ytrig)

cc=ax*by-ay*bx

return
end subroutine

!=================================================================
subroutine divs(aa,bb,cs)
! Computes the divergence of (aa,bb) and returns it in cs.
! Both aa and bb in physical space but cs is in spectral space.

implicit none

 !Passed arrays:
double precision:: aa(ng,ng),bb(ng,ng)   !Physical
double precision:: cs(ng,ng)             !Spectral

 !Work arrays:
double precision:: wkp(ng,ng)            !Physical
double precision:: wka(ng,ng),wkb(ng,ng) !Spectral

!---------------------------------------------------------
wkp=aa
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call xderiv(ng,ng,hrkx,wka,wkb)

wkp=bb
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,cs)

cs=wkb+cs

return
end subroutine

!===================================================================

subroutine spec1d(ss,spec)
! Computes the 1d spectrum of a spectral field ss and returns the
! result in spec.

implicit none

 !Passed variables:
double precision:: ss(ng,ng),spec(0:ng)

 !Local variables:
integer:: kx,ky,k

!--------------------------------------------------------
do k=0,kmax
  spec(k)=zero
enddo

 !x and y-independent mode:
k=kmag(1,1)
spec(k)=spec(k)+f14*ss(1,1)**2

 !y-independent mode:
do kx=2,ng
  k=kmag(kx,1)
  spec(k)=spec(k)+f12*ss(kx,1)**2
enddo

 !x-independent mode:
do ky=2,ng
  k=kmag(1,ky)
  spec(k)=spec(k)+f12*ss(1,ky)**2
enddo

 !All other modes:
do ky=2,ng
  do kx=2,ng
    k=kmag(kx,ky)
    spec(k)=spec(k)+ss(kx,ky)**2
  enddo
enddo

return
end subroutine

!===================================================================

end module     
