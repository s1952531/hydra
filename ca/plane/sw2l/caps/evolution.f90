module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

! Various PVs
double precision:: qt1(ng,ng),qc1(ng,ng),qd1(ng,ng)
double precision:: qt2(ng,ng),qc2(ng,ng),qd2(ng,ng)

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
logical:: ggen

!-----------------------------------------------------------------------
! Define fixed arrays and constants and read initial data:
call init

! Used for regularising contours:
twist=zero

! Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Start the time loop:
do while (t .le. tsim)

  ! Save data periodically:
  itime=nint(t/dt)
  jtime=itime/ngsave
  if (ngsave*jtime .eq. itime) then
    call inversion
    call savegrid(jtime+1)
    ggen=.false.
  else
    ggen=.true.
  endif
  ! ggen is used to indicate if calling inversion is needed in advance below
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
  call advance(ggen)
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! Possibly save final data:
itime=nint(t/dt)
jtime=itime/ngsave
if (ngsave*jtime .eq. itime) then
  call inversion
  call savegrid(jtime+1)
endif
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

subroutine advance(ggen)

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

! Passed variable:
logical:: ggen

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

!------------------------------------------------------------------
! Invert PV and compute velocity at current time level, say t=t^n:
if (ggen) call inversion
! If ggen is false, inversion was called previously at this time level.

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

ttt1=rdis*sds1+sgs1
ttt2=rdis*sds2+sgs2
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

  ttt1=rdis*sds1+sgs1
  ttt2=rdis*sds2+sgs2
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
double precision:: wka(ng,ng),wkb(ng,ng),wkf(ng,ng),wkg(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng),wks(ng,ng)

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

!------------------------------------------------------------------
! ds source:
wka=ds1
call spctop(ng,ng,wka,wkf,xfactors,yfactors,xtrig,ytrig)
! wkf contains delta_1 in physical space
wkp=wkf*u1
wkq=wkf*v1
call divs(wkp,wkq,sds1)
! sds1 contains div(delta_1*u_1,delta_1*v_1) in spectral space
call jacob(u1,v1,wkp)
call ptospc(ng,ng,wkp,wkq,xfactors,yfactors,xtrig,ytrig)
! wkq contains J(u_1,v_1) in spectral space
sds1=two*wkq-sds1
! sds1 is now the full nonlinear ds source in layer 1
sds1(1,1)=zero ! ensures zero domain average

wka=ds2
call spctop(ng,ng,wka,wkg,xfactors,yfactors,xtrig,ytrig)
! wkg contains delta_2 in physical space
wkp=wkg*u2
wkq=wkg*v2
call divs(wkp,wkq,sds2)
! sds2 contains div(delta_2*u_2,delta_2*v_2) in spectral space
call jacob(u2,v2,wkp)
call ptospc(ng,ng,wkp,wkq,xfactors,yfactors,xtrig,ytrig)
! wkq contains J(u_2,v_2) in spectral space
sds2=two*wkq-sds2
! sds2 is now the full nonlinear ds source in layer 2
sds2(1,1)=zero ! ensures zero domain average

!-----------------------------------------------------------------
! gs source:
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
sgs1(1,1)=zero ! ensures zero domain average

! Compute div(zeta_2*u_2,zeta_2*v_2) and store in sgs2 (spectral):
wkp=z2*u2
wkq=z2*v2
call divs(wkp,wkq,sgs2)

! Form terms involving the layer depths:
wkp=wka+wkb
wkq=wkf+wkg
call divs(wkp,wkq,wks)
! wks = div(c_1^2*h_1*(u_1,v_1)+c_2^2*h_2*(u_2,v_2)) in spectral space

! Finalise gs source in layer 2:
sgs2=-cof*sgs2-rksq*wks
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
double precision:: wkp(ng,ng),wka(ng,ng),wkb(ng,ng)
double precision:: spec1(0:ng),spec2(0:ng)
double precision:: ekih,ekiv,epot,etot
integer:: k

!-----------------------------------------------------------------
! Compute energy components, divided by rho_1*H where H = H_1+H_2:
wka=hbar1*(one+h1)
wkb=hbar2*(one+h2)

! Horizontal and vertical parts of the kinetic energy:
ekih=ekmf*sum(wka*(u1**2+v1**2)+alpha*wkb*(u2**2+v2**2))
! ekmf = 0.5*glx*gly/(hbar_1+hbar_2)
ekiv=zero

! Potential energy:
epot=epmf*sum(wka**2+alpha*wkb*(two*wka+wkb)-dape)
! epmf = g*ekmf (see constants.f90)

! Total energy:
etot=ekih+epot+ekiv

! Write energies to ecomp.asc:
write(15,'(f9.2,5(1x,f16.9))') t,ekiv,ekih,ekih+ekiv,epot,etot
write(*,'(a,f9.2,a,f13.6)') ' t = ',t,'  E_tot = ',etot

!-------------------------------------------------------------------
! Compute 1d vorticity, divergence, acceleration divergence and
! height spectra for each vertical mode:

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
write(51,'(f9.2,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(51,'(4(1x,f12.8))') alk(k),spec1(k),spec2(k)
enddo

! Divergence:
wkp=ds1*vm11+ds2*vm12    !Vertical mode 1
call spec1d(wkp,spec1)
spec1=log10(spmf*spec1+1.d-32)
wkp=ds1*vm21+ds2*vm22    !Vertical mode 2
call spec1d(wkp,spec2)
spec2=log10(spmf*spec2+1.d-32)
write(52,'(f9.2,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(52,'(4(1x,f12.8))') alk(k),spec1(k),spec2(k)
enddo

! Acceleration divergence:
wkp=gs1*vm11+gs2*vm12    !Vertical mode 1
call spec1d(wkp,spec1)
spec1=log10(spmf*spec1+1.d-32)
wkp=gs1*vm21+gs2*vm22    !Vertical mode 2
call spec1d(wkp,spec2)
spec2=log10(spmf*spec2+1.d-32)
write(53,'(f9.2,1x,i5)') t,kmaxred
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
write(54,'(f9.2,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(54,'(4(1x,f12.8))') alk(k),spec1(k),spec2(k)
enddo

! spmf: normalisation factor to take into account uneven sampling
!       of wavenumbers in each shell [k-1/2,k+1/2];
! kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
!           kx,ky plane which are not fully populated;
! alk(k) = log_10(k).

!---------------------------------------------------------------
! Write various gridded fields to direct access files:

! PV anomaly:
write(31,rec=igrids) real(t),real(q1)
write(41,rec=igrids) real(t),real(q2)

! Divergence:
wka=ds1
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
write(32,rec=igrids) real(t),real(wkp)
wka=ds2
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
write(42,rec=igrids) real(t),real(wkp)

! Acceleration divergence:
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
