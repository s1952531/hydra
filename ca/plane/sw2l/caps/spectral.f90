module spectral

! Module containing subroutines for spectral operations, inversion, etc.

use constants
use sta2dfft

 !Common arrays, constants:
double precision:: gw11(ng,ng),gw12(ng,ng),gw21(ng,ng),gw22(ng,ng)
double precision:: ho11(ng,ng),ho12(ng,ng),ho21(ng,ng),ho22(ng,ng)
double precision:: bb11(ng,ng),bb12(ng,ng),bb21(ng,ng),bb22(ng,ng)
double precision:: rksq(ng,ng),rlap(ng,ng), filt(ng,ng)
double precision:: diss(ng,ng),rdis(ng,ng),rdisi(ng,ng)
double precision:: bflo(ng,ng),bfhi(ng,ng)
double precision:: vm11,vm12,vm21,vm22

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
double precision:: fac,rkmax,rks,snorm
double precision:: anu,rkfsq,fsq,deni
double precision:: beta1,beta2
integer:: kx,ky,k,kc

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
 !Define vertical mode projection coefficients:
beta1=f12*(hrat-one)+sqrt(f14*(hrat+one)**2-(one-alpha)*hrat)
vm11=one/sqrt(one+beta1**2)
vm12=beta1*vm11
beta2=beta1/mubar
vm22=one/sqrt(one+beta2**2)
vm21=-beta2*vm22
 !For a field F_j in layers j = 1 & 2, its projection on vertical mode 1
 !is F_1*vm11 + F_2*vm12, while for mode 2, it is F_1*vm21 + F_2*vm22.

!-----------------------------------------------------------------------
 !Define a variety of spectral operators:

 !Hyperviscosity coefficient (Dritschel, Gottwald & Oliver, JFM (2017)):
anu=cdamp*cof/rkmax**(2*nnu)
 !Assumes Burger number = 1.

 !Used for de-aliasing filter below:
rkfsq=(dble(ng)/3.d0)**2

 !Used for Butterworth filter below:
kc=ng/6
fac=one/dble(kc**2)

 !Squared Coriolis frequency:
fsq=cof**2

 !Define all spectrally-truncated operators:
do ky=1,ng
  do kx=1,ng
    rks=rk(kx)**2+rk(ky)**2
    if (rks .gt. rkfsq) then
      filt(kx,ky)=zero
      diss(kx,ky)=zero
      rdis(kx,ky)=zero
      rdisi(kx,ky)=zero
      bflo(kx,ky)=zero
      bfhi(kx,ky)=zero
      rksq(kx,ky)=zero
      rlap(kx,ky)=zero
      gw11(kx,ky)=zero
      gw12(kx,ky)=zero
      gw21(kx,ky)=zero
      gw22(kx,ky)=zero
      ho11(kx,ky)=zero
      ho12(kx,ky)=zero
      ho21(kx,ky)=zero
      ho22(kx,ky)=zero
      bb11(kx,ky)=zero
      bb12(kx,ky)=zero
      bb21(kx,ky)=zero
      bb22(kx,ky)=zero
    else
       !De-aliasing filter:
      filt(kx,ky)=one
       !Hyperviscous operator:
      diss(kx,ky)=anu*rks**nnu
       !R operator that appears in time integration:
      rdis(kx,ky)=dt2i+diss(kx,ky)
       !R^{-1} needed in the gamma_l evolution:
      rdisi(kx,ky)=one/rdis(kx,ky)
       !Butterworth low-pass (F) & high-pass (1-F) filters:
      bflo(kx,ky)=one/(one+(fac*rks)**2)
      bfhi(kx,ky)=one-bflo(kx,ky)
       !-grad^2:
      rksq(kx,ky)=rks
       !grad^{-2}:
      rlap(kx,ky)=-one/(rks+1.d-12)
       !Spectral GW operators (G_11 ... G_22):
      gw11(kx,ky)=fsq+csq1*rks
      gw12(kx,ky)=alpha*csq2*rks
      gw21(kx,ky)=csq1*rks
      gw22(kx,ky)=fsq+csq2*rks
       !Operators for inverting the height field:
      deni=one/(gw11(kx,ky)*gw22(kx,ky)-gw21(kx,ky)*gw12(kx,ky))
      ho11(kx,ky)=gw22(kx,ky)*deni
      ho12(kx,ky)=-gw12(kx,ky)*deni
      ho21(kx,ky)=-gw21(kx,ky)*deni
      ho22(kx,ky)=gw11(kx,ky)*deni
       !Semi-implicit inversion operators for delta and gamma_l:
      bb11(kx,ky)=gw11(kx,ky)+rdis(kx,ky)**2
      bb12(kx,ky)=gw12(kx,ky)
      bb21(kx,ky)=gw21(kx,ky)
      bb22(kx,ky)=gw22(kx,ky)+rdis(kx,ky)**2
      deni=one/(bb11(kx,ky)*bb22(kx,ky)-bb12(kx,ky)*bb21(kx,ky))
      bb11(kx,ky)=bb11(kx,ky)*deni
      bb12(kx,ky)=bb12(kx,ky)*deni
      bb21(kx,ky)=bb21(kx,ky)*deni
      bb22(kx,ky)=bb22(kx,ky)*deni
       !Re-define damping operator for use in q_d evolution:
      diss(kx,ky)=two/(one+dt2*diss(kx,ky))
    endif
  enddo
enddo

 !Ensure mean potentials and height anomalies remain zero:
rlap(1,1)=zero
ho11(1,1)=zero
ho12(1,1)=zero
ho21(1,1)=zero
ho22(1,1)=zero
bflo(1,1)=zero
diss(1,1)=zero

return 
end subroutine

!============================================================================
subroutine main_invert(qs1,qs2,ds1,ds2,gs1,gs2,h1,h2,u1,u2,v1,v2,q1,q2,z1,z2)
! Given the PV anomaly qs, divergence ds and acceleration divergence gs
! (all in spectral space and in both layers), this routine computes the
! dimensionless depth anomaly h and the velocity field (u,v) in physical
! space.  It also returns the corrected PV anomaly (q) and relative
! vorticity (z) in physical space (in both layers).

! Note: we assume zero momentum,
! H_1<(1+h1)(u1,v1)> + H_2<(1+h2)(u2,v2)> = 0.  
! The mean flow is determined from this condition.

implicit none

 !Passed variables:
double precision:: qs1(ng,ng),qs2(ng,ng)
double precision:: ds1(ng,ng),ds2(ng,ng)
double precision:: gs1(ng,ng),gs2(ng,ng)
double precision:: h1(ng,ng),h2(ng,ng)
double precision:: u1(ng,ng),u2(ng,ng)
double precision:: v1(ng,ng),v2(ng,ng)
double precision:: q1(ng,ng),q2(ng,ng)
double precision:: z1(ng,ng),z2(ng,ng)

 !Local variables:
double precision,parameter:: tole=1.d-10
 !tole: relative energy norm error in successive iterates when finding
 !      hj, uj & vj from qj, dj & gj (for j = 1 & 2).  The energy norm is
 !      <du1^2+dv1^2+c^2*dh1^2> + alpha*(H_2/H_1)*<du2^2+dv2^2+c^2*dh2^2>
 !      where <:> means a horizontal domain average and (duj,dvj,dhj)
 !      is either the current guess for (uj,vj,hj) or the difference
 !      from the previous guess.

 !Local work arrays:
double precision:: uds1(ng,ng),uds2(ng,ng),vds1(ng,ng),vds2(ng,ng)
double precision:: un1(ng,ng),un2(ng,ng),vn1(ng,ng),vn2(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng)

 !Other constants:
double precision:: qadd,dhrms,durms,enorm
double precision:: uio,vio

!-----------------------------------------------------------------
 !Define spectral divergent velocity (never changes in iteration):
wkc=rlap*ds1
 !This solves Lap(wkc) = delta_1 in spectral space
call xderiv(ng,ng,hrkx,wkc,uds1)
call yderiv(ng,ng,hrky,wkc,vds1)

wkc=rlap*ds2
 !This solves Lap(wkc) = delta_2 in spectral space
call xderiv(ng,ng,hrkx,wkc,uds2)
call yderiv(ng,ng,hrky,wkc,vds2)

 !Obtain physical space copies of qsj and dsj:
wka=qs1
call spctop(ng,ng,wka,q1,xfactors,yfactors,xtrig,ytrig)
wka=qs2
call spctop(ng,ng,wka,q2,xfactors,yfactors,xtrig,ytrig)

!---------------------------------------------------------------
 !Iteratively solve for hj, uj & vj in the two layers j = 1 & 2:

 !Energy norm error (must be > tole to start):
enorm=f12
do while (enorm .gt. tole)
   !Get average PV from the requirement of zero average vorticity:
  qadd=-dsumi*sum(q1*(one+h1))
  q1=q1+qadd
  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd

   !Compute "sources" (wka,wkb) to invert height in spectral space:
  wkc=(one+h1)*q1
  call ptospc(ng,ng,wkc,wka,xfactors,yfactors,xtrig,ytrig)
  wka=gs1-cof*wka
  wkc=(one+h2)*q2
  call ptospc(ng,ng,wkc,wkb,xfactors,yfactors,xtrig,ytrig)
  wkb=gs2-cof*wkb

   !Solve for height anomalies (hold in zj temporarily):
  wkc=ho11*wka+ho12*wkb
  call spctop(ng,ng,wkc,z1,xfactors,yfactors,xtrig,ytrig)
  wkc=ho21*wka+ho22*wkb
  call spctop(ng,ng,wkc,z2,xfactors,yfactors,xtrig,ytrig)
   !See init_spectral for definitions of ho11 ... ho22.

   !Compute rms error in height fields (mubar = layer mass ratio):
  wkc=csq1*(z1-h1)**2+mubar*csq2*(z2-h2)**2
  dhrms=sum(wkc)

   !Re-assign fields:
  h1=z1
  h2=z2

   !Compute relative vorticities (z1 & z2):
  qadd=-dsumi*sum(q1*(one+h1))
  q1=q1+qadd
  z1=(one+h1)*(q1+cof)-cof

  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd
  z2=(one+h2)*(q2+cof)-cof

   !Find the non-divergent part of the velocity:

   !===> Layer 1: Solve Lap(wka) = z1 spectrally:
  wka=z1
  call ptospc(ng,ng,wka,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkc

   !Compute derivatives in spectral space:
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)

   !New velocity components in spectral space, written in (wkb,wkd):
  wkb=uds1-wkb  !uds1 is the fixed divergent part of u1
  wkd=vds1+wkd  !vds1 is the fixed divergent part of v1

   !Convert "new" velocity to physical space as (un1,vn1):
  call spctop(ng,ng,wkb,un1,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,vn1,xfactors,yfactors,xtrig,ytrig)

   !===> Layer 2: Solve Lap(wka) = z2 spectrally:
  wka=z2
  call ptospc(ng,ng,wka,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkc

   !Compute derivatives in spectral space:
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)

   !New velocity components in spectral space, written in (wkb,wkd):
  wkb=uds2-wkb  !uds2 is the fixed divergent part of u2
  wkd=vds2+wkd  !vds2 is the fixed divergent part of v2

   !Convert "new" velocity to physical space as (un2,vn2):
  call spctop(ng,ng,wkb,un2,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,vn2,xfactors,yfactors,xtrig,ytrig)

   !Compute and add mean flow (uio,vio):
  uio=cio*sum(h1*un1+mubar*h2*un2)
  vio=cio*sum(h1*vn1+mubar*h2*vn2)
   !cio=1/(ng^2*(1+mubar)); mubar = (rho_2*H_2)/(rho_1*H_1)
  un1=un1+uio
  un2=un2+uio
  vn1=vn1+vio
  vn2=vn2+vio

   !Compute rms error in uj & vj (mubar = layer mass ratio):
  wkc=(u1-un1)**2+(v1-vn1)**2+mubar*((u2-un2)**2+(v2-vn2)**2)
  durms=sum(wkc)

   !Re-assign velocity components:
  u1=un1
  u2=un2
  v1=vn1
  v2=vn2

   !Compute overall error:
  wkc=u1**2+v1**2+csq1*h1**2+mubar*(u2**2+v2**2+csq2*h2**2)
  enorm=sqrt((durms+dhrms)/sum(wkc))
enddo
 !Passing this, we have converged.
!------------------------------------------------------------------

 !Filter relative vorticity:
call dealias(z1)
call dealias(z2)

return
end subroutine

!=================================================================
subroutine dealias(aa)
! De-aliases the array physical space array aa.  Returns aa in
! physical space.

implicit none

 !Passed array:
double precision:: aa(ng,ng)   !Physical

 !Work array:
double precision:: wka(ng,ng)  !Spectral

!---------------------------------------------------------
call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka
call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)

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
