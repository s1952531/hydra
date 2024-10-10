module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

! Various PVs
double precision:: qt1(ng,ng),qc1(ng,ng),qd1(ng,ng)
double precision:: qt2(ng,ng),qc2(ng,ng),qd2(ng,ng)

! Logical needed for calling PV inversion routine:
logical:: ggen

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================
subroutine advect

! Main subroutine for advecting fields and contours

implicit none

! Local variables:
double precision,parameter:: twistmax=2.5d0
!      twistmax: the maximum value of the time integral of |zeta|_max
!                between regularisations of the contours.
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, the code stops
!      to rebuild the contours in a separate memory space.
integer:: ireg,itime,jtime

!-----------------------------------------------------------------------
! Define fixed arrays and constants and read initial data:
call init

! Used for regularising contours:
twist=zero

! Counter used for counting number of contour regularisations done:
ireg=0

! Logical for deciding when to call inversion in subroutine advance:
ggen=.true.

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Start the time loop:
do while (t .le. tsim)
  ! Save data periodically:
  itime=nint(t/dt)

  jtime=itime/ngsave
  if (ngsave*jtime .eq. itime) call savegrid(jtime+1)

  jtime=itime/ncsave
  if (ncsave*jtime .eq. itime) call savecont(jtime+1)

  ! Perform contour surgery or recontouring when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1

    ! Don't continue if maximum number of regularisations reached:
    if (ireg .eq. nregmax) then
      ! Prepare PV residual qr for recontouring (preserve qs1 & qs2):
      call prepare
      ! Exit module and go to recontouring:
      return
    endif

    ! Regularise the PV contours (surgery + node redistribution):
    call surgery
    ! Record active contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i9,1x,i10)') t,nq,nptq

    twist=twist-twistmax
  endif

  ! Advect flow from time t to t + dt:
  call advance
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! Possibly save final data:
itime=nint(t/dt)

jtime=itime/ngsave
if (ngsave*jtime .eq. itime) call savegrid(jtime+1)

jtime=itime/ncsave
if (ncsave*jtime .eq. itime) call savecont(jtime+1)

return
end subroutine

!=======================================================================

subroutine init

! Initialises residual PV for normal time integration following 
! contour regeneration

implicit none

! Local variables:
double precision:: qa(ng,ng,nz)

!------------------------------------------------------------------
! Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

!-------------------------------------------------------------
! Convert PV contours (xq,yq) to gridded values as qa:
call con2grid(qa)

! Convert each component of qa to spectral space as qc1,2:
call ptospc(ng,ng,qa(:,:,1),qc1,xfactors,yfactors,xtrig,ytrig)
call ptospc(ng,ng,qa(:,:,2),qc2,xfactors,yfactors,xtrig,ytrig)

! Define (spectral) residual PV qdj = (1-F)[qsj-qcj] (j = 1,2):
qd1=bfhi*(qs1-qc1)
qd2=bfhi*(qs2-qc2)
! Here bfhi = 1-F is a high-pass spectral filter

return
end subroutine

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current (spectral) PV anomaly field is stored in qs1,2, and the 
! residual PV needed in congen.f90 is stored in qr.

implicit none

! Local variables:
double precision:: qa(ng,ng,nz)

!-----------------------------------------------------------------
! Convert PV contours (xq,yq) to gridded values as qa:
call con2grid(qa)

! Convert each component of qa to spectral space as qc1,2:
call ptospc(ng,ng,qa(:,:,1),qc1,xfactors,yfactors,xtrig,ytrig)
call ptospc(ng,ng,qa(:,:,2),qc2,xfactors,yfactors,xtrig,ytrig)

! Define spectral PV anomaly (qs1,2) and PV residual (qd1,2):
qs1=bflo*qs1+bfhi*qc1+qd1
qd1=qs1-qc1
qs2=bflo*qs2+bfhi*qc2+qd2
qd2=qs2-qc2

! Convert qd1,2 to physical space as qr (used in recontouring):
call spctop(ng,ng,qd1,qr(:,:,1),xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,qd2,qr(:,:,2),xfactors,yfactors,xtrig,ytrig)
! Note: qd1,2 are overwritten, but we are leaving this module next
!       and qd1,2 will be redefined upon re-entry in subroutine init.

return
end subroutine

!=======================================================================

subroutine advance

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method 
! (for all spectral fields, namely qs, qd, ds & gs) in each layer.

! Uses an iterative implicit method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}-F^n)/2] + N[(F^{n+1}-F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms, and N refers to the nonlinear source
! terms.  We start with a guess for F^{n+1} in N and iterate 
! niter times (see parameter statement below).

implicit none

! Local variables:
integer,parameter:: niter=2

! Spectral fields needed in time stepping:
double precision:: qs1i(ng,ng),sqs1(ng,ng)
double precision:: qs2i(ng,ng),sqs2(ng,ng)
double precision:: qd1i(ng,ng),qd1m(ng,ng),sqd1(ng,ng)
double precision:: qd2i(ng,ng),qd2m(ng,ng),sqd2(ng,ng)
double precision:: ds1i(ng,ng),sds1(ng,ng),nds1(ng,ng)
double precision:: ds2i(ng,ng),sds2(ng,ng),nds2(ng,ng)
double precision:: gs1i(ng,ng),sgs1(ng,ng),ngs1(ng,ng)
double precision:: gs2i(ng,ng),sgs2(ng,ng),ngs2(ng,ng)
double precision:: ttt1(ng,ng),ttt2(ng,ng)
! Velocity arrays needed for calling velint:
double precision:: uu(ng,ng,nz),vv(ng,ng,nz)
! Contour positions needed in time stepping:
double precision:: xqi(nptq),yqi(nptq)
! Contour velocities:
double precision:: uq(nptq),vq(nptq)
! Other local quantities:
double precision:: xx,yy
integer:: i,iter

!-----------------------------------------------------------------------
! Re-initialise qs & qd in each layer at the beginning of the time step:
!           Reset qs = F*(qs-qc) + qc + qd
!             and qd = (1-F)*(qs-qc)
! Here F is a low pass filter (see spectral.f90)
qs1=bflo*qs1+bfhi*qc1+qd1
qd1=bfhi*(qs1-qc1)

qs2=bflo*qs2+bfhi*qc2+qd2
qd2=bfhi*(qs2-qc2)

!--------------------------------------------------------------------
! Invert PV and compute hj, uj & vj at current time level, say t=t^n:
if (ggen) call inversion
! If ggen is false, inversion was called previously at this time level.
ggen=.true.

! Compute twist parameter and save various diagnostics each time step:
call diagnose

!------------------------------------------------------------------
! Start with a guess for F^{n+1} for all contours and fields:

! Contours:
uu(:,:,1)=u1
uu(:,:,2)=u2
vv(:,:,1)=v1
vv(:,:,2)=v2
call velint(uu,vv,uq,vq)
! Here (uq,vq) stands for the velocity at time level n, i.e. u(x^n,t^n)
do i=1,nptq
  ! Store x^n+0.5*dt*u(x^n,t^n) for efficiency in the iteration below:
  xx=xq(i)+dt2*uq(i)
  yy=yq(i)+dt2*vq(i)
  xqi(i)=oms*(xx-twopi*dble(int(xx*pinv)))
  yqi(i)=oms*(yy-twopi*dble(int(yy*pinv)))
  ! Preliminary guess for x^{n+1}:
  xx=xq(i)+dt*uq(i)
  yy=yq(i)+dt*vq(i)
  xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
  yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
enddo

! Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
! (sds,sgs) for divergence and acceleration divergence (ds,gs):
call source(sqs1,sqs2,sqd1,sqd2,sds1,sds2,sgs1,sgs2)

! Update PV fields:
qs1i=qs1+dt2*sqs1
qs1=qs1+dt*sqs1
qd1i=qd1
qd1m=qd1+dt4*sqd1
qd1=diss*(qd1m+dt4*sqd1)-qd1i

qs2i=qs2+dt2*sqs2
qs2=qs2+dt*sqs2
qd2i=qd2
qd2m=qd2+dt4*sqd2
qd2=diss*(qd2m+dt4*sqd2)-qd2i

! Update divergence and acceleration divergence:
ds1i=ds1
gs1i=gs1
nds1=sds1+dt4i*ds1i
ngs1=sgs1+dt4i*gs1i
sds1=nds1+sds1             !2*N_tilde_delta_1
sgs1=ngs1+sgs1             !2*N_tilde_gamma_1

ds2i=ds2
gs2i=gs2
nds2=sds2+dt4i*ds2i
ngs2=sgs2+dt4i*gs2i
sds2=nds2+sds2             !2*N_tilde_delta_2
sgs2=ngs2+sgs2             !2*N_tilde_gamma_2

ttt1=rdis*sds1+pm11*sgs1+pm12*sgs2
ttt2=rdis*sds2+pm21*sgs1+pm22*sgs2
sds1=bb22*ttt1-bb12*ttt2   !2*bar_delta_1
sds2=bb11*ttt2-bb21*ttt1   !2*bar_delta_2

ds1=sds1-ds1i
ds2=sds2-ds2i
gs1=rdisi*(sgs1-gw11*sds1-gw12*sds2)-gs1i
gs2=rdisi*(sgs2-gw21*sds1-gw22*sds2)-gs2i

!------------------------------------------------------------------
! Iterate to improve estimates of F^{n+1}:
do iter=1,niter
  ! Perform inversion at t^{n+1} from estimated quantities:
  call inversion

  ! Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
  ! (sds,sgs) for divergence and acceleration divergence (ds,gs):
  call source(sqs1,sqs2,sqd1,sqd2,sds1,sds2,sgs1,sgs2)

  ! Interpolate gridded velocity in each layer, (u1,v1) & (u2,v2),
  ! at contour nodes as (uq,vq):
  uu(:,:,1)=u1
  uu(:,:,2)=u2
  vv(:,:,1)=v1
  vv(:,:,2)=v2
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqi(i)+dt2*uq(i)
    yy=yqi(i)+dt2*vq(i)
    xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
    yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
  enddo
  ! Now (xq,yq) contain a new guess for x^{n+1}.

  ! Update PV fields:
  qs1=qs1i+dt2*sqs1
  qd1=diss*(qd1m+dt4*sqd1)-qd1i

  qs2=qs2i+dt2*sqs2
  qd2=diss*(qd2m+dt4*sqd2)-qd2i
  
  ! Update divergence and acceleration divergence:
  sds1=nds1+sds1
  sgs1=ngs1+sgs1

  sds2=nds2+sds2
  sgs2=ngs2+sgs2

  ttt1=rdis*sds1+pm11*sgs1+pm12*sgs2
  ttt2=rdis*sds2+pm21*sgs1+pm22*sgs2
  sds1=bb22*ttt1-bb12*ttt2
  sds2=bb11*ttt2-bb21*ttt1

  ds1=sds1-ds1i
  ds2=sds2-ds2i
  gs1=rdisi*(sgs1-gw11*sds1-gw12*sds2)-gs1i
  gs2=rdisi*(sgs2-gw21*sds1-gw22*sds2)-gs2i
enddo

! Advance time:
t=t+dt

return
end subroutine

!=======================================================================

subroutine source(sqs1,sqs2,sqd1,sqd2,sds1,sds2,sgs1,sgs2)

! Gets the source terms (sqs,sqd) for the PV (qs,qd) as well as those
! (sds,sgs) for divergence and acceleration divergence (ds,gs) ---
! --- all in spectral space and in both layers (1 & 2).

! Note that (sds,sgs) only include the nonlinear terms for a 
! semi-implicit treatment, closely analogous to that described in 
! the appendix of Dritschel & Jalali (JFM, 2020).

! The spectral fields ds, gs, qd and qs are all spectrally truncated.
! Note, h, u, v, q & z obtained by main_invert before calling this 
! routine are all spectrally truncated.

! The sources are *not* spectrally truncated as that occurs when they
! are used in subroutine advance.
  
implicit none

! Passed variables:
double precision:: sqs1(ng,ng),sqd1(ng,ng),sds1(ng,ng),sgs1(ng,ng)
double precision:: sqs2(ng,ng),sqd2(ng,ng),sds2(ng,ng),sgs2(ng,ng)

! Local variables:
double precision:: d1(ng,ng),d2(ng,ng)
double precision:: htot1(ng,ng),htot2(ng,ng)
double precision:: hinv1(ng,ng),hinv2(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wkf(ng,ng),wkg(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng)

!---------------------------------------------------------------
! qd source --- only NL advection term is needed:
call gradient(qd1,wka,wkb)
wkp=-u1*wka-v1*wkb
call ptospc(ng,ng,wkp,sqd1,xfactors,yfactors,xtrig,ytrig)
sqd1(1,1)=zero ! ensures zero domain average

call gradient(qd2,wka,wkb)
wkp=-u2*wka-v2*wkb
call ptospc(ng,ng,wkp,sqd2,xfactors,yfactors,xtrig,ytrig)
sqd2(1,1)=zero ! ensures zero domain average

!---------------------------------------------------------------
! qs source --- only NL advection term is needed:
call gradient(qs1,wka,wkb)
wkp=-u1*wka-v1*wkb
call ptospc(ng,ng,wkp,sqs1,xfactors,yfactors,xtrig,ytrig)
sqs1(1,1)=zero ! ensures zero domain average

call gradient(qs2,wka,wkb)
wkp=-u2*wka-v2*wkb
call ptospc(ng,ng,wkp,sqs2,xfactors,yfactors,xtrig,ytrig)
sqs2(1,1)=zero ! ensures zero domain average

!----------------------------------------------------------------
! Find the scaled vertically-averaged non-hydrostatic pressure
! (pn1,pn2) in each layer (note: pnj = p_tilde_j/H_j^2):
call nhpsolve(d1,d2,wkf,wkg,wka,wkb,1)
! On return, d1,d2 = divergence in layers 1,2 in physical space,
! wkf,wkg = -xi_1,-xi_2 where xi_j = D(delta_j)/Dt-delta_j^2,
! while wka,wkb = J_1,J_2 (non-hydrostatic part of zeta tendency).

!---------------------------------------------------------------
! Compute next the nonlinear part of delta source in each layer.

! ==>  Layer 1  <==
! Compute div(delta_1*u_1,delta_1*v_1) and store in sds1 (spectral):
wkp=d1*u1
wkq=d1*v1
call divs(wkp,wkq,sds1)

! Compute 2*delta_1^2 + xi_1 and store in spectral space as wkq:
wkp=two*d1**2-wkf
call ptospc(ng,ng,wkp,wkq,xfactors,yfactors,xtrig,ytrig)
! *** wkf is now free to re-use

! Finalise delta_1 source, N_delta_1 (de-aliased):
sds1=wkq-sds1-pm11*gs1-pm12*gs2
! Note, the linear terms are subtracted here but they largely cancel
! similar terms in wkq-sds1 on rhs.
sds1(1,1)=zero ! ensures zero domain average

! ==>  Layer 2  <==
! Compute div(delta_2*u_2,delta_2*v_2) and store in sds2 (spectral):
wkp=d2*u2
wkq=d2*v2
call divs(wkp,wkq,sds2)

! Compute 2*delta_2^2 + xi_2 and store in spectral space as wkq:
wkp=two*d2**2-wkg
call ptospc(ng,ng,wkp,wkq,xfactors,yfactors,xtrig,ytrig)
! *** wkg is now free to re-use

! Finalise delta_2 source, N_delta_2 (de-aliased):
sds2=wkq-sds2-pm21*gs1-pm22*gs2
! Note, the linear terms are subtracted here but they largely cancel
! similar terms in wkq-sds2 on rhs.
sds2(1,1)=zero ! ensures zero domain average

!-----------------------------------------------------------------
! Compute next the nonlinear part of gamma_l source in each layer.

wkf=csq1*h1 ! Hydrostatic pressure in layer 1
wkg=csq2*h2 ! Hydrostatic pressure in layer 2

! Re-use the arrays htotj & hinvj for the fluxes needed below:
htot1=wkf*u1
htot2=wkg*u2
hinv1=wkf*v1
hinv2=wkg*v2

! ==>  Layer 1  <==
! Compute div(zeta_1*u_1,zeta_1*v_1) and store in sgs1 (spectral):
wkp=z1*u1
wkq=z1*v1
call divs(wkp,wkq,sgs1)

! Combine with J_1 (in wka) from nhpsolve above:
call ptospc(ng,ng,wka,wkq,xfactors,yfactors,xtrig,ytrig)
sgs1=cof*(wkq-sgs1)
! *** wka is now free to re-use

! Form terms involving the layer depths:
wkp=htot1+alpha*htot2
wkq=hinv1+alpha*hinv2
call divs(wkp,wkq,wka)
! wka = div(c_1^2*h_1*(u_1,v_1)+alpha*c_2^2*h_2*(u_2,v_2)) in spectral space

! Finalise gamma_l source in layer 1:
sgs1=sgs1-rksq*wka
! Here -rksq = -k^2 = filtered Laplacian in spectral space
sgs1(1,1)=zero ! ensures zero domain average

! ==>  Layer 2  <==
! Compute div(zeta_2*u_2,zeta_2*v_2) and store in sgs2 (spectral):
wkp=z2*u2
wkq=z2*v2
call divs(wkp,wkq,sgs2)

! Combine with J_2 (in wkb) from nhpsolve above:
call ptospc(ng,ng,wkb,wkq,xfactors,yfactors,xtrig,ytrig)
sgs2=cof*(wkq-sgs2)
! *** wkb is now free to re-use

! Form terms involving the layer depths:
wkp=htot1+htot2
wkq=hinv1+hinv2
call divs(wkp,wkq,wkb)
! wkb = div(c_1^2*h_1*(u_1,v_1)+c_2^2*h_2*(u_2,v_2)) in spectral space

! Finalise gamma_l source in layer 2:
sgs2=sgs2-rksq*wkb
! Here -rksq = -k^2 = filtered Laplacian in spectral space
sgs2(1,1)=zero ! ensures zero domain average

return
end subroutine

!=======================================================================

subroutine inversion

! Finds the gridded dimensionless height anomaly (hj) and velocity
! field (uj,vj) from the PV contours and the (spectral) divergence (dsj)
! and acceleration divergence (gsj) in both layers j = 1 & 2.

implicit none

! Local variables:
double precision:: qa(ng,ng,nz)

!-----------------------------------------------------------------
! Convert PV contours (xq,yq) to gridded values as qa:
call con2grid(qa)

! Convert each component of qa to spectral space as qc1,2:
call ptospc(ng,ng,qa(:,:,1),qc1,xfactors,yfactors,xtrig,ytrig)
call ptospc(ng,ng,qa(:,:,2),qc2,xfactors,yfactors,xtrig,ytrig)

! Combine fields to update qt1,2 with full (spectral) field,
! qt1,2 = F[qs1,2-qc1,2]+qc1,2+qd1,2, where F is a low pass filter:
qt1=bflo*qs1+bfhi*qc1+qd1
qt2=bflo*qs2+bfhi*qc2+qd2

! Invert PV, divergence and acceleration divergence to obtain the
! dimensionless height anomaly and velocity field, as well as the
! gridded PV anomaly and relative vorticity (see spectral.f90):
call main_invert(qt1,qt2,ds1,ds2,gs1,gs2,h1,h2,u1,u2,v1,v2,q1,q2,z1,z2)
! Note: qtj, dsj & gsj are in spectral space while 
!       hj, uj, vj, qj and zj are in physical space (j = 1 & 2).

return
end subroutine

!=======================================================================

subroutine nhpsolve(d1,d2,s1,s2,wka,wkb,iopt)

! Finds the scaled vertically-averaged non-hydrostatic pressure
! pnj = H_j^{-2}*bar{P}_{nj}/h_j, in layers j = 1 & 2.
! On return, dj = divergence in layer j, and if ggen is true, 
! sj = -xi_j = -[D(delta_j)/Dt - delta_j^2], while 
! (wka,wkb) = (J_1,J_2) = non-hydrostatic part of zeta_j tendency.

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
wka=gs1+two*filt*wka
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
wka=gs2+two*filt*wka
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

  ! Use the exact solution valid in the SW limit (k*H_j)^2 -> 0:
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

! If this was called from savegrid, only pn1 & pn2 are needed:
if (.not. ggen) return

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

!=======================================================================

subroutine diagnose

! Computes the twist parameter, the time integral of |zeta|_max, and
! various quantities at every time step to monitor the flow evolution.

implicit none

! Local variables:
double precision:: umax,zrms,zmax
double precision:: uio,vio

!----------------------------------------------------------------------
! Compute diagnostics:
umax=sqrt(max(maxval(u1**2+v1**2),maxval(u2**2+v2**2)))
zmax=max(maxval(abs(z1)),maxval(abs(z2)))
zrms=sqrt(cio*sum(z1**2+mubar*z2**2))
! above, mubar = (rho_2*H_2)/(rho_1*H_1) and cio = 1/(ng^2*(1+mubar))

! Increment the integral of |zeta|_max:
twist=twist+dt*zmax

! Record various diagnostics to monitor.asc:
write(17,'(1x,f12.5,4(1x,f12.6))') t,f12*zrms**2,zrms,zmax,umax

! Compute and write mean flow:
uio=cio*sum(h1*u1+mubar*h2*u2)
vio=cio*sum(h1*v1+mubar*h2*v2)
write(16,'(1x,f12.5,2(1x,f14.10))') t,uio,vio

return
end subroutine

!=======================================================================

subroutine savegrid(igrids)

! Saves PV, energy and various spectra at the desired save time

implicit none

! Passed variable:
integer:: igrids

! Local variables:
double precision:: d1(ng,ng),d2(ng,ng)
double precision:: htot1(ng,ng),htot2(ng,ng)
double precision:: wkp(ng,ng),wka(ng,ng),wkb(ng,ng)
double precision:: spec1(0:ng),spec2(0:ng)
double precision:: ekih,ekiv,epot,etot,fac
integer:: k

!-----------------------------------------------------------------
! Invert PV and compute hj, uj, vj etc at the current time:
call inversion

ggen=.false.
! This means there is no need to call inversion again immediately
! after calling this routine in advance.

! Solve for the non-hydrostatic pressure (pn1,pn2) in each layer:
call nhpsolve(d1,d2,htot1,htot2,wka,wkb,igrids-1)
! On return, d1,d2 = divergence in layers 1,2 in physical space.
! The remaining variables are not used here.

!-----------------------------------------------------------------
! Compute energy components, divided by rho_1*H where H = H_1+H_2:
htot1=hbar1*(one+h1)
wka=d1*htot1

htot2=hbar2*(one+h2)
wkb=d2*htot2

! Horizontal and vertical parts of the kinetic energy:
fac=f12*garea/hbar
ekih=fac*sum(htot1*(u1**2+v1**2)+alpha*htot2*(u2**2+v2**2))
ekiv=fac*sum(htot1*f13*wka**2+alpha*htot2*(wka**2+wka*wkb+f13*wkb**2))

! Potential energy:
fac=fac*gravity
epot=fac*sum(htot1**2+alpha*htot2*(two*htot1+htot2)-dape)

! Total energy:
etot=ekih+epot+ekiv

! Write energies to ecomp.asc:
write(15,'(f12.5,5(1x,f16.9))') t,ekiv,ekih,ekih+ekiv,epot,etot
write(*,'(a,f12.5,a,f13.6)') ' t = ',t,'  E_tot = ',etot

!-------------------------------------------------------------------
! Compute 1d vorticity, divergence, acceleration divergence, height
! and non-hydrostatic pressure spectra for each vertical mode:

! Relative vorticity:
wkp=z1
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
wkp=z2
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkp=wka*vm11+wkb*vm12    !Vertical mode 1
call spec1d(wkp,spec1)
spec1=log10(spmf*spec1+1.d-32)
wkp=wka*vm21+wkb*vm22    !Vertical mode 2
call spec1d(wkp,spec2)
spec2=log10(spmf*spec2+1.d-32)
write(51,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(51,'(4(1x,f12.8))') alk(k),spec1(k),spec2(k)
enddo

! Velocity divergence:
wkp=ds1*vm11+ds2*vm12    !Vertical mode 1
call spec1d(wkp,spec1)
spec1=log10(spmf*spec1+1.d-32)
wkp=ds1*vm21+ds2*vm22    !Vertical mode 2
call spec1d(wkp,spec2)
spec2=log10(spmf*spec2+1.d-32)
write(52,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(52,'(4(1x,f12.8))') alk(k),spec1(k),spec2(k)
enddo

! Linearised acceleration divergence:
wkp=gs1*vm11+gs2*vm12    !Vertical mode 1
call spec1d(wkp,spec1)
spec1=log10(spmf*spec1+1.d-32)
wkp=gs1*vm21+gs2*vm22    !Vertical mode 2
call spec1d(wkp,spec2)
spec2=log10(spmf*spec2+1.d-32)
write(53,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(53,'(4(1x,f12.8))') alk(k),spec1(k),spec2(k)
enddo

! Dimensionless height anomaly:
wkp=h1
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
wkp=h2
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkp=wka*vm11+wkb*vm12    !Vertical mode 1
call spec1d(wkp,spec1)
spec1=log10(spmf*spec1+1.d-32)
wkp=wka*vm21+wkb*vm22    !Vertical mode 2
call spec1d(wkp,spec2)
spec2=log10(spmf*spec2+1.d-32)
write(54,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(54,'(4(1x,f12.8))') alk(k),spec1(k),spec2(k)
enddo

! Vertically-averaged non-hydrostatic pressure:
wkp=pn1
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
wkp=pn2
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkp=wka*vm11+wkb*vm12    !Vertical mode 1
call spec1d(wkp,spec1)
spec1=log10(spmf*spec1+1.d-32)
wkp=wka*vm21+wkb*vm22    !Vertical mode 2
call spec1d(wkp,spec2)
spec2=log10(spmf*spec2+1.d-32)
write(55,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(55,'(4(1x,f12.8))') alk(k),spec1(k),spec2(k)
enddo

! spmf: normalisation factor to take into account uneven sampling
!       of wavenumbers in each shell [k-1/2,k+1/2];
! kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
!           kx,ky plane which are not fully populated;
! alk(k) = log_10(k).

!-------------------------------------------------------------------
! Write various gridded fields to direct access files:

! PV anomaly:
write(31,rec=igrids) real(t),real(q1)
write(41,rec=igrids) real(t),real(q2)

! Divergence:
write(32,rec=igrids) real(t),real(d1)
write(42,rec=igrids) real(t),real(d2)

! Linearised acceleration divergence:
wka=gs1
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
write(33,rec=igrids) real(t),real(wkp)
wka=gs2
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
write(43,rec=igrids) real(t),real(wkp)

! Dimensionless depth anomalies:
write(34,rec=igrids) real(t),real(h1)
write(44,rec=igrids) real(t),real(h2)

! Relative vorticities:
write(35,rec=igrids) real(t),real(z1)
write(45,rec=igrids) real(t),real(z2)

! Vertically-averaged non-hydrostatic pressures:
write(36,rec=igrids) real(t),real(pn1)
write(46,rec=igrids) real(t),real(pn2)

return
end subroutine

!=======================================================================
      
subroutine savecont(irec)

! Saves PV contours for post-processing and imaging

implicit none

! Passed variable:
integer:: irec

! Local variables:
double precision:: ss(ng,ng)
double precision:: qa(ng,ng)
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') irec

! Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

! Save residual needed to build ultra-fine-grid PV for plotting purposes:
ss=qt1-qc1
call spctop(ng,ng,ss,qa,xfactors,yfactors,xtrig,ytrig)
write(83,rec=2*irec-1) real(t),real(qa)
ss=qt2-qc2
call spctop(ng,ng,ss,qa,xfactors,yfactors,xtrig,ytrig)
write(83,rec=2*irec) real(t),real(qa)
! Note: both layers are saved above

! Save PV contours:
open(81,file='cont/qqindex'//pind,form='unformatted',status='replace')
write(81) npq(1:nq),i1q(1:nq),indq(1:nq),layq(1:nq)
close(81)

open(82,file='cont/qqnodes'//pind,form='unformatted',status='replace')
write(82) xq(1:nptq),yq(1:nptq)
close(82)

return
end subroutine

!=======================================================================

! Main end module
end module
