!###############################################===##########################
!  Initialises a flow with balanced fields obtained from the conditions 
!  delta_t=gamma_t=0 using the PV anomaly field q previously set up with
!  a data generation routine.

!  The balance conditions are those relevant to a two-layer shallow-water
!  flow (hydrostatic).  Here, we use the same balance conditions to ensure
!  identical fields of h_j, u_j and v_j when comparing sw2l & va2l results.

! Adapted from analogous sw2l code on 3/8/2020 by D G Dritschel @ St Andrews
!############################################===#############################

program dgbalini

 !Import contants, parameters and common arrays:
use spectral

implicit none

double precision,parameter:: tole=1.d-10
 !tole: relative energy norm error in successive iterates when finding
 !      hj, uj & vj from qj, dj & gj (for j = 1 & 2).  The energy norm is
 !      <du1^2+dv1^2+c^2*dh1^2> + alpha*(H_2/H_1)*<du2^2+dv2^2+c^2*dh2^2>
 !      where <:> means a horizontal domain average and (duj,dvj,dhj)
 !      is either the current guess for (uj,vj,hj) or the difference
 !      from the previous guess.

 !Fields (PV anomaly etc):
double precision:: q1(ng,ng),z1(ng,ng),d1(ng,ng)
double precision:: q2(ng,ng),z2(ng,ng),d2(ng,ng)
double precision:: u1(ng,ng),v1(ng,ng),h1(ng,ng),htot1(ng,ng)
double precision:: u2(ng,ng),v2(ng,ng),h2(ng,ng),htot2(ng,ng)
double precision:: u1pre(ng,ng),v1pre(ng,ng),h1pre(ng,ng)
double precision:: u2pre(ng,ng),v2pre(ng,ng),h2pre(ng,ng)
double precision:: ds1(ng,ng),gs1(ng,ng),sgs1(ng,ng)
double precision:: ds2(ng,ng),gs2(ng,ng),sgs2(ng,ng)

double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng)
double precision:: wkf(ng,ng),wkg(ng,ng),wkp(ng,ng),wkq(ng,ng)
double precision:: wks(ng,ng)

 !Other constants:
double precision:: qadd,dhrms,durms,enorm
double precision:: uio,vio,t

!----------------------------------------------------------------------
! Initialise inversion constants and arrays:
call init_spectral

!----------------------------------------------------------------------
! Read in gridded ***shallow-water*** PV anomaly in each layer:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,q1
read(11,rec=2) t,q2
close(11)

! De-alias:
call dealias(q1)
call dealias(q2)

!----------------------------------------------------------------------
! Start with zero divergence and acceleration divergence:
ds1=zero
ds2=zero
gs1=zero
gs2=zero

! Initialise previous guess for hj, uj & vj:
h1=zero
h2=zero
u1=zero
u2=zero
v1=zero
v2=zero

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

  ! Compute "sources" (wka,wkb) to invert height in spectral space:
  wkc=-cof*(one+h1)*q1
  call ptospc(ng,ng,wkc,wka,xfactors,yfactors,xtrig,ytrig)
  wkc=-cof*(one+h2)*q2
  call ptospc(ng,ng,wkc,wkb,xfactors,yfactors,xtrig,ytrig)

  ! Solve for height anomalies (hold in zj temporarily):
  wkc=ho11*wka+ho12*wkb
  call spctop(ng,ng,wkc,z1,xfactors,yfactors,xtrig,ytrig)
  wkc=ho21*wka+ho22*wkb
  call spctop(ng,ng,wkc,z2,xfactors,yfactors,xtrig,ytrig)
  ! See init_spectral for definitions of ho11 ... ho22.

  ! Compute rms error in height fields (mubar = layer mass ratio):
  wkc=csq1*(z1-h1)**2+mubar*csq2*(z2-h2)**2
  dhrms=sum(wkc)

  ! Re-assign fields:
  h1=z1
  h2=z2

  ! Compute relative vorticities (z1 & z2):
  qadd=-dsumi*sum(q1*(one+h1))
  q1=q1+qadd
  z1=(one+h1)*(q1+cof)-cof

  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd
  z2=(one+h2)*(q2+cof)-cof

  ! Find the non-divergent velocity (here the only part):

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

  ! Compute rms error in uu & vv (mubar = layer mass ratio):
  wkc=(u1-wkf)**2+(v1-wkg)**2+mubar*((u2-wkp)**2+(v2-wkq)**2)
  durms=sum(wkc)

  ! Re-assign velocity components:
  u1=wkf
  u2=wkp
  v1=wkg
  v2=wkq

  ! Compute overall error:
  wkc=u1**2+v1**2+csq1*h1**2+mubar*(u2**2+v2**2+csq2*h2**2)
  enorm=sqrt((durms+dhrms)/sum(wkc))
enddo
! Passing this, we have converged.

!-----------------------------------------------------------------
! Iterate to find the balanced fields:
h1pre=h1
h2pre=h2
u1pre=u1
u2pre=u2
v1pre=v1
v2pre=v2

! Energy norm error (must be > tole to start):
enorm=f12
do while (enorm .gt. tole)
  !----------------------------------
  ! Obtain balanced estimate for gs1:
  wka=ds1
  call spctop(ng,ng,wka,wkb,xfactors,yfactors,xtrig,ytrig)
  ! wkb contains delta_1 in physical space
  wkc=wkb*u1
  wkd=wkb*v1
  call divs(wkc,wkd,gs1)
  ! gs1 contains div(delta_1*u_1,delta_1*v_1) in spectral space
  call jacob(u1,v1,wka)
  call ptospc(ng,ng,wka,wkb,xfactors,yfactors,xtrig,ytrig)
  ! wkb contains J(u_1,v_1) in spectral space
  gs1=filt*(gs1-two*wkb)
  ! gs1 is now the de-aliased estimate for gamma_1 (spectral)
  
  !----------------------------------
  ! Obtain balanced estimate for gs2:
  wka=ds2
  call spctop(ng,ng,wka,wkb,xfactors,yfactors,xtrig,ytrig)
  ! wkb contains delta_2 in physical space
  wkc=wkb*u2
  wkd=wkb*v2
  call divs(wkc,wkd,gs2)
  ! gs2 contains div(delta_2*u_2,delta_2*v_2) in spectral space
  call jacob(u2,v2,wka)
  call ptospc(ng,ng,wka,wkb,xfactors,yfactors,xtrig,ytrig)
  ! wkb contains J(u_2,v_2) in spectral space
  gs2=filt*(gs2-two*wkb)
  ! gs2 is now the de-aliased estimate for gamma_2 (spectral)

  !-------------------------------------------------------------
  ! Obtain balanced estimates for ds1 & ds2 (solve 2x2 problem):
  wkf=csq1*h1 ! Hydrostatic pressure in layer 1
  wkg=csq2*h2 ! Hydrostatic pressure in layer 2

  wka=wkf*u1  ! c_1^2*h_1*u_1
  wkb=wkg*u2  ! c_2^2*h_2*u_2
  wkf=wkf*v1  ! c_1^2*h_1*v_1
  wkg=wkg*v2  ! c_2^2*h_2*v_2

  ! Compute div(zeta_1*u_1,zeta_1*v_1) and store in sgs1 (spectral):
  wkp=z1*u1
  wkq=z1*v1
  call divs(wkp,wkq,sgs1)

  ! Form terms involving the layer depths:
  wkp=wka+alpha*wkb
  wkq=wkf+alpha*wkg
  call divs(wkp,wkq,wks)
  ! wks = div(c_1^2*h_1*(u_1,v_1)+alpha*c_2^2*h_2*(u_2,v_2)) in spectral space

  ! Finalise gs source in layer 1:
  sgs1=-cof*sgs1-rksq*wks
  ! Here -rksq = -k^2 = filtered Laplacian in spectral space

  ! Compute div(zeta_2*u_2,zeta_2*v_2) and store in sgs2 (spectral):
  wkp=z2*u2
  wkq=z2*v2
  call divs(wkp,wkq,sgs2)

  ! Form terms involving the layer depths:
  wkp=wka+wkb
  wkq=wkf+wkg
  call divs(wkp,wkq,wks)
  ! wks = div(c_1^2*h_1*(u_1,v_1)+c_2^2*h_2*(u_2,v_2)) in spectral space

  ! Finalise gamma_l source in layer 2:
  sgs2=-cof*sgs2-rksq*wks
  ! Here -rksq = -k^2 = filtered Laplacian in spectral space

  ! Solve 2x2 problem for ds1 & ds2:
  ds1=ho11*sgs1+ho12*sgs2
  ds2=ho21*sgs1+ho22*sgs2

  !--------------------------------------
  ! Find height anomaly fields (h1 & h2):

  ! Get average PV from the requirement of zero average vorticity:
  qadd=-dsumi*sum(q1*(one+h1))
  q1=q1+qadd
  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd

  ! Compute "sources" (wka,wkb) to invert height in spectral space:
  wkc=(one+h1)*q1
  call ptospc(ng,ng,wkc,wka,xfactors,yfactors,xtrig,ytrig)
  wka=gs1-cof*wka
  wkc=(one+h2)*q2
  call ptospc(ng,ng,wkc,wkb,xfactors,yfactors,xtrig,ytrig)
  wkb=gs2-cof*wkb

  ! Solve for height anomalies (h1 & h2):
  wkc=ho11*wka+ho12*wkb
  call spctop(ng,ng,wkc,h1,xfactors,yfactors,xtrig,ytrig)
  wkc=ho21*wka+ho22*wkb
  call spctop(ng,ng,wkc,h2,xfactors,yfactors,xtrig,ytrig)
  ! See init_spectral for definitions of ho11 ... ho22.

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

  ! Compute relative vorticities (z1 & z2):
  qadd=-dsumi*sum(q1*(one+h1))
  q1=q1+qadd
  z1=(one+h1)*(q1+cof)-cof

  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd
  z2=(one+h2)*(q2+cof)-cof

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
 !Write velocity divergence:
open(11,file='dd_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
call spctop(ng,ng,ds1,d1,xfactors,yfactors,xtrig,ytrig)
write(11,rec=1) zero,d1
call spctop(ng,ng,ds2,d2,xfactors,yfactors,xtrig,ytrig)
write(11,rec=2) zero,d2
close(11)

 !Write acceleration divergence:
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

 !End main program
end program dgbalini
!=======================================================================
