!############################################################################
!  Initialises a flow with balanced fields obtained from the conditions 
!  delta_t=gamma_t=0 using the PV anomaly field q previously set up with
!  a data generation routine.

!  Adapted from swbalinit.f90 on 13/8/2020 by D G Dritschel @ St Andrews
!############################################################################

program dgbalini

 !Import contants, parameters and common arrays:
use spectral

implicit none

 !Global variables:

 !PV anomaly and relative vorticity (q & zeta):
double precision:: q1(ng,ng),z1(ng,ng)
double precision:: q2(ng,ng),z2(ng,ng)

 !Velocity field and dimensionless height anomaly (u, v & h):
double precision:: u1(ng,ng),v1(ng,ng),h1(ng,ng)
double precision:: u2(ng,ng),v2(ng,ng),h2(ng,ng)

 !Vertically-averaged non-hydrostatic pressure (p_n):
double precision:: pn1(ng,ng),pn2(ng,ng)

 !Spectral fields of delta = div(u,v) and gamma = f*zeta - c^2*Lap(h):
double precision:: ds1(ng,ng),gs1(ng,ng)
double precision:: ds2(ng,ng),gs2(ng,ng)

double precision:: t,deni
integer:: kx,ky

!----------------------------------------------------------------------
! Initialise inversion constants and arrays:
call init_spectral

! Define operators needed to find gamma below:
do ky=1,ng
  do kx=1,ng
    if (abs(pm11(kx,ky)) .gt. small) then
      deni=one/(pm11(kx,ky)*pm22(kx,ky)-pm21(kx,ky)*pm12(kx,ky))
      bb11(kx,ky)= pm22(kx,ky)*deni
      bb12(kx,ky)=-pm12(kx,ky)*deni
      bb21(kx,ky)=-pm21(kx,ky)*deni
      bb22(kx,ky)= pm11(kx,ky)*deni
    else
      bb11(kx,ky)=zero
      bb12(kx,ky)=zero
      bb21(kx,ky)=zero
      bb22(kx,ky)=zero
    endif
  enddo
enddo

!----------------------------------------------------------------------
! Read in gridded ***shallow-water*** PV anomaly in each layer:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,q1
read(11,rec=2) t,q2
close(11)
! This is overwritten at the end with the VA PV anomaly.

! De-alias:
call dealias(q1)
call dealias(q2)

! Compute balanced flow and write data:
call genbal


! Subroutine definitions follow:
contains

!=======================================================================

subroutine genbal
    
double precision,parameter:: tole=1.d-10
 !tole: relative energy norm error in successive iterates when finding
 !      hj, uj & vj from qj, dj & gj (for j = 1 & 2).  The energy norm is
 !      <du1^2+dv1^2+c^2*dh1^2> + alpha*(H_2/H_1)*<du2^2+dv2^2+c^2*dh2^2>
 !      where <:> means a horizontal domain average and (duj,dvj,dhj)
 !      is either the current guess for (uj,vj,hj) or the difference
 !      from the previous guess.  Note: dj stands for delta in layer j
 !      and gj stands for gamma in layer j.

 !Fields (PV anomalies etc):
double precision:: d1(ng,ng),d2(ng,ng)
double precision:: htot1(ng,ng),htot2(ng,ng)
double precision:: u1pre(ng,ng),v1pre(ng,ng),h1pre(ng,ng)
double precision:: u2pre(ng,ng),v2pre(ng,ng),h2pre(ng,ng)
double precision:: cc1(ng,ng),cc2(ng,ng)
double precision:: hd1(ng,ng),hd2(ng,ng)

double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng)
double precision:: wkf(ng,ng),wkg(ng,ng),wkp(ng,ng),wkq(ng,ng)
double precision:: wks(ng,ng)

 !Other scalar variables:
double precision:: qadd,dhrms,durms,enorm
double precision:: uio,vio,t
integer:: iopt

!-------------------------------------------------------------------
! Iteratively solve for hj, uj & vj (initial guess for dj = gj = 0):

! Energy norm error (must be > tole to start):
enorm=f12
do while (enorm .gt. tole)
  ! Get average PV from the requirement of zero average vorticity:
  qadd=-dsumi*sum(q1*(one+h1))
  q1=q1+qadd
  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd

   !Solve for height anomalies (hold in zj temporarily):
  wkc=(one+h1)*q1   !d1 = 0 at this stage: no Jacobian terms
  call ptospc(ng,ng,wkc,wka,xfactors,yfactors,xtrig,ytrig)
  wka=-ho1*cof*wka  !g1 = 0 at this stage: simpler inversion
  call spctop(ng,ng,wka,z1,xfactors,yfactors,xtrig,ytrig)

  wkc=(one+h2)*q2   !d2 = 0 at this stage: no Jacobian terms
  call ptospc(ng,ng,wkc,wka,xfactors,yfactors,xtrig,ytrig)
  wka=-ho2*cof*wka  !g2 = 0 at this stage simpler inversion
  call spctop(ng,ng,wka,z2,xfactors,yfactors,xtrig,ytrig)

  ! Compute rms error in height fields (mubar = layer mass ratio):
  wkc=csq1*(z1-h1)**2+mubar*csq2*(z2-h2)**2
  dhrms=sum(wkc)

  ! Re-assign updated height anomalies hj:
  h1=z1
  h2=z2

  ! Compute relative vorticities zj:
  qadd=-dsumi*sum(q1*(one+h1))
  q1=q1+qadd
  z1=(one+h1)*(q1+cof)-cof   !ds1 = 0 at this stage: no Jacobian terms

  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd
  z2=(one+h2)*(q2+cof)-cof   !ds2 = 0 at this stage: no Jacobian terms

  ! Find the non-divergent velocity (here the only part as dj = 0):

  ! ===> Layer 1: Solve Lap(wka) = z1 spectrally:
  wka=z1
  call ptospc(ng,ng,wka,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkc

  ! Compute derivatives in spectral space:
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)

  ! Flip sign for x component:
  wkb=-wkb

  ! Convert "new" velocity to physical space as (wkf,wkg):
  call spctop(ng,ng,wkb,wkf,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,wkg,xfactors,yfactors,xtrig,ytrig)

  ! ===> Layer 2: Solve Lap(wka) = z2 spectrally:
  wka=z2
  call ptospc(ng,ng,wka,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkc

  ! Compute derivatives in spectral space:
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)

  ! Flip sign for x component:
  wkb=-wkb

  ! Convert "new" velocity to physical space as (wkp,wkq):
  call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,wkq,xfactors,yfactors,xtrig,ytrig)

  ! Compute and add mean flow (uio,vio):
  uio=cio*sum(h1*wkf+mubar*h2*wkp)
  vio=cio*sum(h1*wkg+mubar*h2*wkq)
  ! cio=1/(ng^2*(1+mubar)); mubar = (rho_2*H_2)/(rho_1*H_1)
  wkf=wkf+uio
  wkp=wkp+uio
  wkg=wkg+vio
  wkq=wkq+vio

  ! Compute rms error in uj & vj (mubar = layer mass ratio):
  wkc=(u1-wkf)**2+(v1-wkg)**2+mubar*((u2-wkp)**2+(v2-wkq)**2)
  durms=sum(wkc)

  ! Re-assign updated velocity components uj & vj:
  u1=wkf
  u2=wkp
  v1=wkg
  v2=wkq

  ! Compute overall relative error:
  wkc=u1**2+v1**2+csq1*h1**2+mubar*(u2**2+v2**2+csq2*h2**2)
  enorm=sqrt((durms+dhrms)/sum(wkc))
enddo
! Passing this, we have converged (but only for when dj = gj = 0).

!-----------------------------------------------------------------------
! Iterate to find the balanced fields of dj & gj as well as hj, uj & vj:
h1pre=h1
h2pre=h2
u1pre=u1
u2pre=u2
v1pre=v1
v2pre=v2

! Initialise dsj = gsj = 0:
ds1=0
ds2=0
gs1=0
gs2=0

! Flag needed in first call to nhpsolve below to initialise pnj:
iopt=0

! Energy norm error (must be > tole to start):
enorm=f12
do while (enorm .gt. tole)
  !------------------------------------------------------------
  ! Solve for the non-hydrostatic pressure pnj in each layer j:
  call nhpsolve(d1,d2,hd1,hd2,cc1,cc2,iopt)
  ! On return, dj = divergence in layer j in physical space,
  ! hdj = -xi_j where xi_j = D(delta_j)/Dt-delta_j^2, while
  ! ccj = J_j (non-hydrostatic part of zeta tendency).
  iopt=1

  !---------------------------------------------------------------
  ! Obtain balanced estimates for gsj (spectral gj) from dj_t = 0:
  wkc=d1*u1
  wkd=d1*v1
  call divs(wkc,wkd,wka)
  ! wka contains div(delta_1*u_1,delta_1*v_1) in spectral space
  hd1=hd1-two*d1**2
  ! hd1 now contains -xi_1 - 2*delta_1^2
  call ptospc(ng,ng,hd1,wkb,xfactors,yfactors,xtrig,ytrig)
  ! Add linear terms pm11*gs1+pm12*gs2 to cancel those in wka+wkb:
  hd1=wka+wkb+pm11*gs1+pm12*gs2
  ! hd1 is now (minus) the entire nonlinear delta_1 source (spectral)

  wkc=d2*u2
  wkd=d2*v2
  call divs(wkc,wkd,wka)
  ! wka contains div(delta_2*u_2,delta_2*v_2) in spectral space
  hd2=hd2-two*d2**2
  ! hd2 now contains -xi_2 - 2*delta_2^2
  call ptospc(ng,ng,hd2,wkb,xfactors,yfactors,xtrig,ytrig)
  ! Add linear terms pm21*gs1+pm22*gs2 to cancel those in wka+wkb:
  hd2=wka+wkb+pm21*gs1+pm22*gs2
  ! hd2 is now (minus) the entire nonlinear delta_2 source (spectral)

  ! Solve 2x2 problem to update gsj:
  gs1=bb11*hd1+bb12*hd2
  gs2=bb21*hd1+bb22*hd2
  ! Note: bbjk are re-defined at the beginning of this code.

  !---------------------------------------------------------------
  ! Obtain balanced estimates for dsj (spectral dj) from gj_t = 0:
  wkf=csq1*h1 ! Hydrostatic pressure in layer 1
  wkg=csq2*h2 ! Hydrostatic pressure in layer 2

  ! Compute div(zeta_1*u_1,zeta_1*v_1) and store in wkc (spectral):
  wkp=z1*u1
  wkq=z1*v1
  call divs(wkp,wkq,wkc)

  ! Combine with J_1 (in cc1) from nhpsolve above:
  call ptospc(ng,ng,cc1,wkq,xfactors,yfactors,xtrig,ytrig)
  wkc=cof*(wkq-wkc)

  ! Form term involving the layer depth:
  wkp=wkf*u1       !c_1^2*h_1*u_1
  wkq=wkf*v1       !c_1^2*h_1*v_1
  call divs(wkp,wkq,wks)
  ! wks = div(c_1^2*h_1*(u_1,v_1)) in spectral space

  ! Update balanced delta in layer 1:
  ds1=ho1*(wkc-rksq*wks) !ho1 = -H_1^{-1} = 1/G_1 in the notes
  wkd=ds1
  call spctop(ng,ng,wkd,d1,xfactors,yfactors,xtrig,ytrig)

  ! Compute div(zeta_2*u_2,zeta_2*v_2) and store in wkc (spectral):
  wkp=z2*u2
  wkq=z2*v2
  call divs(wkp,wkq,wkc)

  ! Combine with J_2 (in cc2) from nhpsolve above:
  call ptospc(ng,ng,cc2,wkq,xfactors,yfactors,xtrig,ytrig)
  wkc=cof*(wkq-wkc)

  ! Form terms involving the layer depth:
  wkp=wkg*u2       !c_2^2*h_2*u_2
  wkq=wkg*v2       !c_2^2*h_2*v_2
  call divs(wkp,wkq,wks)
  ! wks = div(c_2^2*h_2*(u_2,v_2)) in spectral space

  ! Update balanced delta in layer 2:
  ds2=ho2*(wkc-rksq*wks) !ho2 = -H_2^{-1} = 1/G_2 in the notes
  wkd=ds2
  call spctop(ng,ng,wkd,d2,xfactors,yfactors,xtrig,ytrig)

  !--------------------------------------
  ! Find height anomaly fields (h1 & h2):
  htot1=hbar1*(one+h1)
  htot2=hbar2*(one+h2)

  ! Form products htot_j*delta_j for Jacobian terms below:
  hd1=htot1*d1
  call dealias(hd1)
  hd2=htot2*d2
  call dealias(hd2)

   !Re-compute cc1 = (1/3)*J(htot_1,htot_1*delta_1):
  call jacob(htot1,hd1,cc1)
  cc1=f13*cc1

   !Re-compute cc2 = J(htot_1,htot_1*delta_1+(1/2)*htot_2*delta_2) +
   !                 J(htot_2,(1/2)*htot_1*delta_1+(1/3)*htot_2*delta_2):
  wkc=hd1+f12*hd2
  call jacob(htot1,wkc,wka)
  wkc=f12*hd1+f13*hd2
  call jacob(htot2,wkc,wkb)
  cc2=wka+wkb

  ! Get average PV from the requirement of zero average vorticity:
  qadd=-dsumi*sum(q1*(one+h1))
  q1=q1+qadd
  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd

   !Solve for height anomalies hj:
  wkc=(one+h1)*q1-cc1
  call ptospc(ng,ng,wkc,wka,xfactors,yfactors,xtrig,ytrig)
  wkc=ho1*(gs1-cof*wka) ! -H_1^{-1}(gamma_1 - eta_1) in spectral space
  call spctop(ng,ng,wkc,h1,xfactors,yfactors,xtrig,ytrig)

  wkc=(one+h2)*q2-cc2
  call ptospc(ng,ng,wkc,wkb,xfactors,yfactors,xtrig,ytrig)
  wkc=ho2*(gs2-cof*wkb) ! -H_2^{-1}(gamma_2 - eta_2) in spectral space
  call spctop(ng,ng,wkc,h2,xfactors,yfactors,xtrig,ytrig)

  ! Compute rms error in height fields (mubar = layer mass ratio):
  wkc=csq1*(h1-h1pre)**2+mubar*csq2*(h2-h2pre)**2
  dhrms=sum(wkc)

  !-------------------------------------------------------
  ! Find the velocity (uj,vj) in the two layers j = 1 & 2:

  wkc=rlap*ds1
  ! This solves Lap(wkc) = delta_1 in spectral space
  call xderiv(ng,ng,hrkx,wkc,wkf)
  call yderiv(ng,ng,hrky,wkc,wkg)

  wkc=rlap*ds2
  ! This solves Lap(wkc) = delta_2 in spectral space
  call xderiv(ng,ng,hrkx,wkc,wkp)
  call yderiv(ng,ng,hrky,wkc,wkq)

  ! Define total height in each layer:
  htot1=hbar1*(one+h1)
  htot2=hbar2*(one+h2)

  ! Re-form products htot_j*delta_j for Jacobian terms below:
  hd1=htot1*d1
  call dealias(hd1)
  hd2=htot2*d2
  call dealias(hd2)

  ! Re-compute cc1 = (1/3)*J(htot_1,htot_1*delta_1):
  call jacob(htot1,hd1,cc1)
  cc1=f13*cc1

  ! Re-compute cc2 = J(htot_1,htot_1*delta_1+(1/2)*htot_2*delta_2) +
  !                  J(htot_2,(1/2)*htot_1*delta_1+(1/3)*htot_2*delta_2):
  wkc=hd1+f12*hd2
  call jacob(htot1,wkc,wka)
  wkc=f12*hd1+f13*hd2
  call jacob(htot2,wkc,wkb)
  cc2=wka+wkb

   ! Compute relative vorticities zj:
  qadd=-dsumi*sum(q1*(one+h1))
  q1=q1+qadd
  z1=(one+h1)*(q1+cof)-cof-cc1

  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd
  z2=(one+h2)*(q2+cof)-cof-cc2

  wka=z1
  call ptospc(ng,ng,wka,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkc
  ! This solves Lap(wka) = zeta_1 in spectral space

  ! Compute derivatives of streamfunction wka in spectral space:
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)

  ! New velocity components in spectral space, written in (wkb,wkd):
  wkb=wkf-wkb
  wkd=wkg+wkd

  ! Convert to physical space as (u1,v1):
  call spctop(ng,ng,wkb,u1,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,v1,xfactors,yfactors,xtrig,ytrig)

  wka=z2
  call ptospc(ng,ng,wka,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkc
  ! This solves Lap(wka) = zeta_1 in spectral space

  ! Compute derivatives of streamfunction wka in spectral space:
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)

  ! New velocity components in spectral space, written in (wkb,wkd):
  wkb=wkp-wkb
  wkd=wkq+wkd

  ! Convert to physical space as (u2,v2):
  call spctop(ng,ng,wkb,u2,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,v2,xfactors,yfactors,xtrig,ytrig)

  ! Compute and add mean flow (uio,vio):
  uio=cio*sum(h1*u1+mubar*h2*u2)
  vio=cio*sum(h1*v1+mubar*h2*v2)
  ! cio=1/(ng^2*(1+mubar)); mubar = (rho_2*H_2)/(rho_1*H_1)
  u1=u1+uio
  u2=u2+uio
  v1=v1+vio
  v2=v2+vio

  ! Compute rms error in uj & vj (mubar = layer mass ratio):
  wkc=(u1-u1pre)**2+(v1-v1pre)**2+mubar*((u2-u2pre)**2+(v2-v2pre)**2)
  durms=sum(wkc)

   !Compute overall error:
  wkc=u1pre**2+v1pre**2+csq1*h1pre**2+mubar*(u2pre**2+v2pre**2+csq2*h2pre**2)
  enorm=sqrt((durms+dhrms)/(sum(wkc)+1.d-20))

  write(*,*) ' Relative energy error = ',enorm

  !Otherwise continue with another iteration:
  h1pre=h1
  h2pre=h2
  u1pre=u1
  u2pre=u2
  v1pre=v1
  v2pre=v2
enddo

!-----------------------------------------------------------------
 !Write delta:
open(11,file='dd_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
call spctop(ng,ng,ds1,d1,xfactors,yfactors,xtrig,ytrig)
write(11,rec=1) zero,d1
call spctop(ng,ng,ds2,d2,xfactors,yfactors,xtrig,ytrig)
write(11,rec=2) zero,d2
close(11)

 !Write gamma:
open(11,file='gg_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
call spctop(ng,ng,gs1,wka,xfactors,yfactors,xtrig,ytrig)
write(11,rec=1) zero,wka
call spctop(ng,ng,gs2,wka,xfactors,yfactors,xtrig,ytrig)
write(11,rec=2) zero,wka
close(11)

 !Write height anomaly (only for interest):
open(11,file='hh_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,h1
write(11,rec=2) zero,h2
close(11)

 !Write relative vorticity (only for interest):
open(11,file='zz_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,z1
write(11,rec=2) zero,z2
close(11)

 !Write vertically-averaged non-hydrostatic pressure (only for interest):
open(11,file='pp_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,pn1
write(11,rec=2) zero,pn2
close(11)

 !Write PV anomaly:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
 !Define total height in each layer:
htot1=hbar1*(one+h1)
htot2=hbar2*(one+h2)
 !Form products H_j*(1+h_j)*delta_j for Jacobian terms below:
d1=htot1*d1
call dealias(d1)
d2=htot2*d2
call dealias(d2)
 !Compute wkf = (1/3)*J(htot_1,htot_1*delta_1):
call jacob(htot1,d1,wkf)
call dealias(wkf)
wkf=f13*wkf/(one+h1)
call dealias(wkf)
 !Compute wkg = J(htot_1,htot_1*delta_1+(1/2)*htot_2*delta_2) +
 !              J(htot_2,(1/2)*htot_1*delta_1+(1/3)*htot_2*delta_2):
wkc=d1+f12*d2
call jacob(htot1,wkc,wka)
wkc=f12*d1+f13*d2
call jacob(htot2,wkc,wkb)
wkg=wka+wkb
call dealias(wkg)
wkg=wkg/(one+h2)
call dealias(wkg)
q1=q1+wkf
q2=q2+wkg
write(11,rec=1) zero,q1
write(11,rec=2) zero,q2
close(11)

write(*,*)
write(*,*) &
     ' Balanced initial fields balanced and re-written to qq, dd & gg_init.r8'

return
end subroutine genbal

!=======================================================================

subroutine nhpsolve(d1,d2,s1,s2,wka,wkb,iopt)

! Finds the scaled vertically-averaged non-hydrostatic pressure
! pnj = H_j^{-2}*bar{P}_{nj}/h_j, in layers j = 1 & 2.
! On return, dj = delta in layer j, sj = -xi_j where x_j = 
! [D(delta_j)/Dt - delta_j^2, while (wka,wkb) = (J_1,J_2) =
! non-hydrostatic part of zeta_j tendency.

! iopt = 0 on initialisation when a previous estimate for pnj
! is unavailable.  

! *** All passed variables are in physical space ***
  
implicit none

! Passed variables:
double precision:: d1(ng,ng),d2(ng,ng),s1(ng,ng),s2(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng)
integer:: iopt

! Local variables:
double precision,parameter:: ptol=1.d-7
! ptol: maximum relative rms NH pressure error

double precision:: htot1(ng,ng),htot2(ng,ng),hinv1(ng,ng),hinv2(ng,ng)
double precision:: h1x(ng,ng),h1y(ng,ng),tt(ng,ng)
double precision:: r1x(ng,ng),r1y(ng,ng),t1x(ng,ng),t1y(ng,ng)
double precision:: r2x(ng,ng),r2y(ng,ng),sgx(ng,ng),sgy(ng,ng)
double precision:: ht1x(ng,ng),ht1y(ng,ng)
double precision:: ch1(ng,ng),ch2(ng,ng),ch3(ng,ng)
double precision:: wkc(ng,ng),wkd(ng,ng),wkp(ng,ng)
double precision:: perr

!---------------------------------------------------------------
! Get total dimensionless layer thicknesses (1 + h_j):
htot1=one+h1
htot2=one+h2

! Get their inverses:
hinv1=one/htot1
call dealias(hinv1)
hinv2=one/htot2
call dealias(hinv2)

! Define T = 6/(4 + 3*mu) where mu = alpha*(H_2/H_1)*(1+h_2)/(1+h_1):
wkp=mubar*htot2*hinv1
call dealias(wkp)
tt=six/(four+three*wkp)
call dealias(tt)

! Define c_hat_1,2,3:
wkp=three*tt*hinv1
call dealias(wkp)
ch1=wkp*hinv1-cona1
call dealias(ch1)
ch2=hrati*wkp*hinv2-cona2
call dealias(ch2)
wkd=tt*hinv2          !wkd = T/(1+h_2): *** preserve to define sgx & sgy
call dealias(wkd)
ch3=two*wkd*(hinv2+mubar3*hinv1)-cona3   !mubar3 = 3*mubar
call dealias(ch3)

! Calculate (h1x,h1y) = grad{h_1} and (r1x,r1y) = H_1^2*grad{h_1}/(1+h_1):
wkp=h1
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,h1x,h1y)
r1x=hbsq1*hinv1*h1x
call dealias(r1x)
r1y=hbsq1*hinv1*h1y
call dealias(r1y)

! Calculate (t1x,t1y) = tau = T*(r1x,r1y) = H_1^2*T*grad{h_1}/(1+h_1):
t1x=tt*r1x
call dealias(t1x)
t1y=tt*r1y
call dealias(t1y)

! Define tau/2 for efficiency in pressure iteration below:
ht1x=f12*t1x
ht1y=f12*t1y

! Calculate (r2x,r2y) = H_2^2*grad{h_2}/(1+h_2):
wkp=hbsq2*h2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,r2x,r2y)
r2x=hinv2*r2x
call dealias(r2x)
r2y=hinv2*r2y
call dealias(r2y)

! Calculate (sgx,sgy) = (r2x,r2y) + H_1*H_2*T*grad{h_1}/(1+h_2):
wka=wkd*h1x
call dealias(wka)
sgx=r2x+hb1hb2*wka
wka=wkd*h1y
call dealias(wka)
sgy=r2y+hb1hb2*wka
! Note: here wkd contains T/(1+h_2) from above; wkd can now be re-used

!--------------------------------------------------------------------
! Form the fixed parts (s1,s2) of the sources needed to calculate the
! non-hydrostatic pressure:

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute gamma_tilde = gamma_l + 2*(J(u,v) - delta^2) in layer 1:
wka=ds1
call spctop(ng,ng,wka,d1,xfactors,yfactors,xtrig,ytrig)
! d1 contains delta_1 in physical space
call jacob(u1,v1,wkp)
! wkp contains J(u_1,v_1) in physical space
wkp=wkp-d1**2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
! Form eta_2 = f*(zeta_2 - f*h_2) -> wkb
wkb=cof*(z2-cof*h2)
call ptospc(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
wka=gs1+ee2*(gs2-wkc)+two*filt*wka
! Above gs1 + ee2*(gs2-wkc) = gamma_1l (spectral); ee2 = E_2 operator
call spctop(ng,ng,wka,s1,xfactors,yfactors,xtrig,ytrig)
! s1 now contains gamma_tilde_1 in physical space (de-aliased)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute gamma_tilde = gamma_l + 2*(J(u,v) - delta^2) in layer 2:
wka=ds2
call spctop(ng,ng,wka,d2,xfactors,yfactors,xtrig,ytrig)
! d2 contains delta_2 in physical space
call jacob(u2,v2,wkp)
! wkp contains J(u_2,v_2) in physical space
wkp=wkp-d2**2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
! Form eta_1 = f*(zeta_1 - f*h_1) -> wkb
wkb=cof*(z1-cof*h1)
call ptospc(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
wka=gs2+ee1*(gs1-wkc)+two*filt*wka
! Above gs2 + ee1*(gs1-wkc) = gamma_2l (spectral); ee1 = E_1 operator
call spctop(ng,ng,wka,s2,xfactors,yfactors,xtrig,ytrig)
! s2 now contains gamma_tilde_2 in physical space (de-aliased)

!======================================================================
! Iterate to find approximate solution:

if (iopt .eq. 0) then
  ! This is done at t = 0 only when pn1 & pn2 are not yet defined:
  wkc=htot1*s1
  call dealias(wkc)
  wkd=htot2*s2
  call dealias(wkd)

  ! Use the exact solution valid in the long-wave limit (k*H_j)^2 -> 0:
  pn2=-htot2*(hhrati*wkc+f13*wkd)
  pn1=1.5d0*ahrsq*pn2-(f13*htot1+f14*mubar*htot2)*wkc
endif
  
perr=one
do while (perr .gt. ptol)

  ! Find full S_1 (rhs of 1st pressure equation) in spectral space (wka):
  wkp=f23*pn1-ahrsq*pn2     !P = (2/3)*pn1 - alpha*(H_2/H_1)^2*pn2
  wkc=wkp*t1x
  wkd=wkp*t1y
  call divs(wkc,wkd,wka)    !wka = div{P*tau} (spectral)
  wkp=s1+ch1*wkp            !wkp = gamma_tilde_1 + c_hat_1*P
  call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=wkc-wka               !wka = gamma_tilde_1 + c_hat_1*P - div{P*tau}

  ! Find full S_2 (rhs of 2nd pressure equation) in spectral space (wkb):
  wkc=pn1*ht1x+pn2*sgx
  wkd=pn1*ht1y+pn2*sgy
  call divs(wkc,wkd,wkb)    !wkb = div{(1/2)*pn1*tau+pn2*sigma} (spectral)
  wkp=s2-ch2*pn1+ch3*pn2    !wkp = gamma_tilde_2 - c_hat_2*pn1 + c_hat_3*pn2
  call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
  wkb=wkc-wkb               !wkb = gamma_tilde_2 - c_hat_2*pn1 + c_hat_3*pn2
                            !     -div{(1/2)*pn1*tau+pn2*sigma}

  ! Obtain next approximation for pn1 & pn2 (hold in wkc & wkd temporarily):
  wkp=pi11*wka+pi12*wkb
  call spctop(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)

  wkp=pi21*wka+pi22*wkb
  call spctop(ng,ng,wkp,wkd,xfactors,yfactors,xtrig,ytrig)

  ! Compute relative rms difference error:
  perr=sqrt(sum((wkc-pn1)**2+(wkd-pn2)**2)/sum(pn1**2+pn2**2))

  ! Copy wkc & wkd into pn1 & pn2 with relaxation (essential!):
  pn1=f12*(pn1+wkc)
  pn2=f12*(pn2+wkd)
  ! From tests, this 50/50 blend of old and new solutions leads to
  ! exponential convergence (50% error reduction each iteration).
enddo

!----------------------------------------------------------
! Compute -xi_1 = -[D(delta_1)/Dt - delta_1^2] -> s1
!     and -xi_2 = -[D(delta_2)/Dt - delta_2^2] -> s2:
wkp=f23*pn1-ahrsq*pn2
s1=(ch1+cona1)*wkp
s2=(ch3+cona3)*pn2-(ch2+cona2)*pn1

! Compute non-hydrostatic parts of zeta_1 tendency (-> wka):
wkp=tt*wkp
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call gradient(wkb,wkc,wkd)
wka=r1x*wkd-r1y*wkc                     ! wka = J_1

! Compute non-hydrostatic parts of zeta_2 tendency (-> wkb):
wkp=pn2
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
call gradient(wkb,sgx,sgy)

wkp=pn2*hinv2+hhrati*pn1*hinv1          ! pn2/(1+h_2)+(H_1/2*H_2)*pn1/(1+h_1)
call dealias(wkp)
wkp=tt*wkp                              ! wkp = kappa in the notes
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call gradient(wkb,wkc,wkd)              ! (wkc,wkd) = grad{kappa}
wkb=r2x*sgy-r2y*sgx+hb1hb2*(h1x*wkd-h1y*wkc)   ! wkb = J_2

!------------------------------------------------------------------
wkp=pn1*(one+h1)

return
end subroutine

 !End main program
end program dgbalini
!=======================================================================
