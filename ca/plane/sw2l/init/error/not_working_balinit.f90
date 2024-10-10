!############################################################################
!  Initialises a flow with balanced fields obtained from the conditions 
!  delta_t=gamma_t=0 using the PV anomaly field q previously set up with
!  a data generation routine.

!  Adapted from va2l/init/balinit.f90 on 27/11/2020 by D G Dritschel @ St A
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
    if (abs(ee1(kx,ky)) .gt. small) then
      deni=one/(one-ee1(kx,ky)*ee2(kx,ky))
      bb11(kx,ky)=deni
      bb12(kx,ky)=-ee2(kx,ky)*deni
      bb21(kx,ky)=-ee1(kx,ky)*deni
      bb22(kx,ky)=deni
    else
      bb11(kx,ky)=zero
      bb12(kx,ky)=zero
      bb21(kx,ky)=zero
      bb22(kx,ky)=zero
    endif
  enddo
enddo
! Ensure zero domain average:
bb11(1,1)=zero
bb12(1,1)=zero
bb21(1,1)=zero
bb22(1,1)=zero

!----------------------------------------------------------------------
! Read in gridded PV anomaly in each layer:
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
double precision:: u1pre(ng,ng),v1pre(ng,ng),h1pre(ng,ng)
double precision:: u2pre(ng,ng),v2pre(ng,ng),h2pre(ng,ng)

double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng)
double precision:: wkf(ng,ng),wkg(ng,ng),wkp(ng,ng),wkq(ng,ng)
double precision:: wks(ng,ng)

 !Other scalar variables:
double precision:: qadd,dhrms,durms,enorm
double precision:: uio,vio,t

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
  wkc=(one+h1)*q1
  call ptospc(ng,ng,wkc,wka,xfactors,yfactors,xtrig,ytrig)
  wka=-ho1*cof*wka  !g1 = 0 at this stage: simpler inversion
  call spctop(ng,ng,wka,z1,xfactors,yfactors,xtrig,ytrig)

  wkc=(one+h2)*q2
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
  z1=(one+h1)*(q1+cof)-cof

  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd
  z2=(one+h2)*(q2+cof)-cof

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

! Energy norm error (must be > tole to start):
enorm=f12
do while (enorm .gt. tole)

  !---------------------------------------------------------------
  ! Obtain balanced estimates for gsj (spectral gj) from dj_t = 0:
  wka=ds1
  call spctop(ng,ng,wka,wkf,xfactors,yfactors,xtrig,ytrig)
  ! wkf contains delta_1 in physical space
  wkp=wkf*u1
  wkq=wkf*v1
  call divs(wkp,wkq,wkf)
  ! wkf contains div(delta_1*u_1,delta_1*v_1) in spectral space
  call jacob(u1,v1,wkp)
  call ptospc(ng,ng,wkp,wkq,xfactors,yfactors,xtrig,ytrig)
  ! wkq contains J(u_1,v_1) in spectral space
  wkb=cof*(z2-cof*h2)
  ! Here wkb = eta_2 = f*(zeta_2 - f*h_2)
  call ptospc(ng,ng,wkb,wks,xfactors,yfactors,xtrig,ytrig)
  wkf=wkf+ee2*wks-two*wkq
  ! wkf is now (minus) the full nonlinear ds source in layer 1

  wka=ds2
  call spctop(ng,ng,wka,wkg,xfactors,yfactors,xtrig,ytrig)
  ! wkg contains delta_2 in physical space
  wkp=wkg*u2
  wkq=wkg*v2
  call divs(wkp,wkq,wkg)
  ! wkg contains div(delta_2*u_2,delta_2*v_2) in spectral space
  call jacob(u2,v2,wkp)
  call ptospc(ng,ng,wkp,wkq,xfactors,yfactors,xtrig,ytrig)
  ! wkq contains J(u_2,v_2) in spectral space
  wkb=cof*(z1-cof*h1)
  ! Here wkb = eta_1 = f*(zeta_1 - f*h_1)
  call ptospc(ng,ng,wkb,wks,xfactors,yfactors,xtrig,ytrig)
  wkg=wkg+ee1*wks-two*wkq
  ! wkg is now (minus) the full nonlinear ds source in layer 2

  ! Solve 2x2 problem to update gsj:
  gs1=bb11*wkf+bb12*wkg
  gs2=bb21*wkf+bb22*wkg
  ! Note: bbjk are re-defined at the beginning of this code.

  !---------------------------------------------------------------
  ! Obtain balanced estimates for dsj (spectral dj) from gj_t = 0:
  wkf=csq1*h1 ! Hydrostatic pressure in layer 1
  wkg=csq2*h2 ! Hydrostatic pressure in layer 2

  ! Compute div(zeta_1*u_1,zeta_1*v_1) and store in wkc (spectral):
  wkp=z1*u1
  wkq=z1*v1
  call divs(wkp,wkq,wkc)

  ! Form term involving the layer depth:
  wkp=wkf*u1       !c_1^2*h_1*u_1
  wkq=wkf*v1       !c_1^2*h_1*v_1
  call divs(wkp,wkq,wks)
  ! wks = div(c_1^2*h_1*(u_1,v_1)) in spectral space

  ! Update balanced delta in layer 1:
  ds1=-ho1*(cof*wkc+rksq*wks) !ho1 = -H_1^{-1} = 1/G_1 in the notes
  wkd=ds1
  call spctop(ng,ng,wkd,d1,xfactors,yfactors,xtrig,ytrig)

  ! Compute div(zeta_2*u_2,zeta_2*v_2) and store in wkc (spectral):
  wkp=z2*u2
  wkq=z2*v2
  call divs(wkp,wkq,wkc)

  ! Form terms involving the layer depth:
  wkp=wkg*u2       !c_2^2*h_2*u_2
  wkq=wkg*v2       !c_2^2*h_2*v_2
  call divs(wkp,wkq,wks)
  ! wks = div(c_2^2*h_2*(u_2,v_2)) in spectral space

  ! Update balanced delta in layer 2:
  ds2=-ho2*(cof*wkc+rksq*wks) !ho2 = -H_2^{-1} = 1/G_2 in the notes
  wkd=ds2
  call spctop(ng,ng,wkd,d2,xfactors,yfactors,xtrig,ytrig)

  !--------------------------------------
  ! Get average PV from the requirement of zero average vorticity:
  qadd=-dsumi*sum(q1*(one+h1))
  q1=q1+qadd
  qadd=-dsumi*sum(q2*(one+h2))
  q2=q2+qadd

   !Solve for height anomalies hj:
  wkc=(one+h1)*q1
  call ptospc(ng,ng,wkc,wka,xfactors,yfactors,xtrig,ytrig)
  wkc=ho1*(gs1-cof*wka) ! -H_1^{-1}(gamma_1 - eta_1) in spectral space
  call spctop(ng,ng,wkc,h1,xfactors,yfactors,xtrig,ytrig)

  wkc=(one+h2)*q2
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

   ! Compute relative vorticities zj:
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

 !Write PV anomaly:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,q1
write(11,rec=2) zero,q2
close(11)

write(*,*)
write(*,*) &
     ' Balanced initial fields balanced and re-written to qq, dd & gg_init.r8'

return
end subroutine genbal

 !End main program
end program dgbalini
!=======================================================================
