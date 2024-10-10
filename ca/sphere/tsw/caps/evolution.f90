module evolution

! Module contains subroutines to evolve fields as detailed in spe.f90.

use common

implicit none

 !Fields at beginning of a time step (spectral):
double precision:: eeo(ng,nt),tto(ng,nt),ddo(ng,nt)
double precision:: qso(ng,nt),qdo(ng,nt)

 !Velocity fields and contour velocities (physical):
double precision:: uu(0:ng+1,nt),vv(0:ng+1,nt)
double precision:: u(npm),v(npm),w(npm)

 !Widely used work arrays:
double precision:: wke(ng,nt),wku(ng,nt),wkv(ng,nt)

 !Twist parameter used to decide when to regularise PV contours:
double precision:: twist

contains 

!=============================================================
subroutine advect

! Main subroutine for advecting fields and contours

implicit none

 !Local variables:
double precision,parameter:: twistmax=2.5
 !twistmax: the maximum value of the time integral of |zeta|_max
 !          between regularisations of the contours.
integer,parameter:: nregmax=20
 !nregmax: every nregmax regularisations, the code stops to rebuild
 !         the contours in a separate module (congen.f90).
integer:: i,j,nreg

!--------------------------------------------------------------
 !Initialise PV fields after recontouring:
call init

 !Initialise counter for regularisations:
nreg=0

 !Used for regularising contours:
twist=zero

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (nint(tmax/tsim) .le. nperiod) 
  do while (t .lt. tmax) 

     !Evolve contours and fields from t to t + dt:
    call timestep

    if (twist .gt. twistmax) then
       !Check whether to regularise contours or rebuild them
      nreg=nreg+1

      if (nreg .eq. nregmax) then
         !Rebuild the contours using module congen:
        call prepare
        return
      endif
 
       !Regularise the contours (surgery + node redistribution):
      call surgery
      write(19,'(1x,f12.5,1x,i9,1x,i10)') t,n,npt
      twist=twist-twistmax
    endif
  enddo  

   !Write data for later post-processing:
  call writedata
   !Continue with another:
  tmax=t+tsim
enddo
 !End time loop 
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !Code reaching this point has finished entire run

 !If t is a multiple of tsim, write data for later post-processing:
if (abs(t/tsim-dble(nint(t/tsim))) .lt. 1.d-6) call writedata

 !Set t to final time to ensure main code loop stops:
t=tfin

return
end subroutine

!=======================================================================

subroutine init
!  Initialises pv fields after recontouring and before timestepping

implicit none

integer:: i,j

!------------------------------------------------------------------------
 !Update gridded version of contour PV:
call con2grid(qc)

 !Take away coriolis frequency from qs & qc and prepare 
 !the residual PV field qd:
do i=1,nt
  do j=1,ng
    qq(j,i)=qc(j,i)-cof(j)
    qs(j,i)=qs(j,i)-cof(j)
  enddo
enddo
 !FFT qq & qs in longitude (semi-spectral):
call forfft(ng,nt,qq,trig,factors) 
call forfft(ng,nt,qs,trig,factors) 
 !qq contains the PV anomaly arising from the contours in 
 !semi-spectral space (called "qc" in the comments below)
 !qs contains the entire PV anomaly in semi-spectral space.

 !Calculate the residual PV from qd = (1-F)*(q-qc):
qd=qs-qq
qq=qs

 !Apply high-pass Butterworth filter 1-F to (q-qc) in array qd:
call hifilter(qd)

 !Return qq (the entire PV anomaly) to physical space:
call revfft(ng,nt,qq,trig,factors) 

!--------------------------------------------------------------------
 !Time which we are going to compute until:
tmax=tsim*(one+dble(int(t/tsim)))
 !Here, t is the initial time, read in above.

dt=small
 !The time step is worked out in subroutine adapt.

 !Write initial conditions for later post-processing:
if (t .eq. zero) then
  call invert
  call writedata
endif

return
end subroutine

!=======================================================================

subroutine prepare
! Prepares PV fields for recontouring before exiting this module

implicit none

integer:: i,j

!-------------------------------------------
 !Re-initialise qs & qd before recontouring:
 !    Reset qs = q = F*qs + (1-F)*qc + qd
 !      and qd = (1-F)*(q-qc)
call reset
 !qq contains the full gridded PV anomaly after this call.
 !qd contains the *spectral* residual PV.

 !Correct qq so that relative vorticity (zeta, or zz) has zero mean:
call getzet

 !Define total PV and residual needed for recontouring:
do i=1,nt
  do j=1,ng
    qs(j,i)=qq(j,i)+cof(j)
    qd(j,i)=qs(j,i)-qc(j,i)
  enddo
enddo
 !Note: qc is the PV associated with the contours.

return
end subroutine

!=======================================================================

subroutine timestep
!   Integrates the equations of motion from time t to time t + dt,
!   where dt is dynamically adapted in subroutine adapt.
!   *** Uses implicit trapezoidal integration ***

implicit none

 !Define local parameters and arrays:
integer,parameter:: niter=2
! niter+1: number of interations of the implicit scheme

double precision:: qsi(ng,nt),reeo(ng,nt),rddo(ng,nt) !spectral
double precision:: xo(npm),yo(npm),zo(npm)
double precision:: tcrit,xnew,ynew,znew,fac,dti,dzrms
integer:: i,j,m,iter,itime

!------------------------------------------------------------------------
 !Re-initialise qs & qd at the beginning of the time step:
 !    Reset qs = q = F*qs + (1-F)*qc + qd
 !      and qd = (1-F)*(q-qc)
call reset
 !qq contains the full spectral PV anomaly after this call.

 !Find gridded velocity and interpolate it at each contour node:
call invert

 !Save various fields approximately every tsave time units:
itime=int(t/tsave+small)
tcrit=tsave*dble(itime)
if ((t-dt-tcrit)*(t+small*tsave-tcrit) .lt. zero) then
  if (itime .gt. iref) then
    call dump
    iref=itime
  endif
endif

 !Possibly adapt the time step and adjust damping rates:
call adapt

 !Find tendencies:
call tends

!--------------------------------------------------------------------
 !Store quantities at the beginning of a time step needed in the
 !iteration below and take first forward time step:

 !eta and delta (for semi-implicit time stepping):
do m=1,nt
  do j=1,ng
    eeo(j,m)=ee(j,m)
    fac=f12*see(j,m)
    reeo(j,m)=fac+hdti*ee(j,m)
    wkc(j,m)=fac+reeo(j,m)
     !wkc = R_eta in sistep below

    ddo(j,m)=dd(j,m)
    fac=f12*sdd(j,m)
    rddo(j,m)=fac+hdti*dd(j,m)
    wkd(j,m)=fac+rddo(j,m)
     !wkd = R_delta in sistep below
  enddo
enddo

 !Get approximate eta & delta at t + dt:
call sistep

!--------------------------------------------------------
 !Iterate now to solve the implicit problem in time:
do iter=1,niter

   !Generate the PV anomaly at t+dt from the contours & PV fields:
  call genpv

   !Update velocity at this time and interpolate it to the contours:
  call invert

   !Update tendencies at this time:
  call tends

   !Update contours and PV fields:
  do i=1,npt
    xnew=xo(i)+dt2*u(i)
    ynew=yo(i)+dt2*v(i)
    znew=zo(i)+dt2*w(i)
    fac=one/sqrt(xnew**2+ynew**2+znew**2)
    x(i)=fac*xnew
    y(i)=fac*ynew
    z(i)=fac*znew
  enddo
  qs=qso+dt2*sqs
  qd=qdo+dt2*sqd

   !If there is thermal forcing, update log potential temperature:
  if (thermal) tt=tto+dt2*stt

   !Define R_eta & R_delta (use wkc & wkd):
  wkc=f12*see+reeo
  wkd=f12*sdd+rddo

   !Update eta & delta at t + dt:
  call sistep
enddo

 !Apply hyperviscous damping in latitude:
call latdamp(ee)
call latdamp(dd)
if (thermal) call latdamp(tt)

 !Update eta, h & delta in physical space:
eep=ee
call revfft(ng,nt,eep,trig,factors) 
hhp=(one+eep)**kappa
ddp=dd
call revfft(ng,nt,ddp,trig,factors)
 !If there is thermal forcing, also update log potential temperature:
if (thermal) then
  ttp=tt
  call revfft(ng,nt,ttp,trig,factors)
endif

 !Increment time:
t=t+dt

!-----------------------------------------------------------------------
 !Estimate q_t at high wavenumbers (purely diagnostic):
dti=one/dt
qso=(qs-qsi)*dti
call hifilter(qso)
call revfft(ng,nt,qso,trig,factors)
call getabs(qso,dqdthi)
 !Here, dqdthi is the L1 norm of the high-pass part of q_t

 !Add any stochastic forcing to qd here:
if (stoch) then
   !rms value of added vorticity:
  dzrms=sqrt(two*esr*dt)
   !generate random field wka with this rms value:
  call ranspec(wkc,dzrms)

   !Convert to PV forcing:
  wkc=wkc/(one+eep)

   !FFT wkc and add to qd in semi-spectral space:
  call forfft(ng,nt,wkc,trig,factors) 
  qd=qd+wkc
endif

return
end subroutine

!==========================================================================

subroutine invert
! Calculates the velocity field (uu,vv) given eta (ee), delta (dd)
! and the PV anomaly field (qq).  The velocity field is then interpolated
! at all nodes as (u,v).

implicit none

integer:: j,m

!---------------------------------------------------------
 !First compute divergent part of the velocity field:

 !Invert Laplace's operator on the delta:
call laplinv(dd,pp,wkv)
 !Here the divergence potential is pp while d(pp)/dlat = wkv,
 !the divergent meridional velocity.

 !Compute and store divergent zonal velocity (times cos(lat)) in wku:
call deriv(ng,nt,rk,pp,wku) 
 !wku = d(pp)/dlon

 !Compute relative vorticity zz = (1+ee)*(qq+f) - f:
call getzet
 !Note qq+f = q, the full PV.

 !Convert the relative vorticity to semi-spectral space as wke:
wke=zz
call forfft(ng,nt,wke,trig,factors)  

 !Invert Laplace's operator:
call laplinv(wke,pp,wkb)
 !Here the streamfunction is pp while d(pp)/dlat = wkb.
 !Compute d(pp)/dlon = wkc:
call deriv(ng,nt,rk,pp,wkc)

 !Complete calculation of zonal and meridional velocity, wku & wkv:
do m=1,nt
  do j=1,ng
    wku(j,m)=clati(j)*wku(j,m)-wkb(j,m)
    wkv(j,m)=clati(j)*wkc(j,m)+wkv(j,m)
  enddo
enddo

 !Get physical space velocity:
call revfft(ng,nt,wku,trig,factors)
call revfft(ng,nt,wkv,trig,factors)
uu(1:ng,:)=wku
vv(1:ng,:)=wkv
 !uu & vv have extended ranges in latitude for ease of interpolation

 !Spectrally-truncated relative vorticity is in wke:
call revfft(ng,nt,wke,trig,factors)  

!----------------------------------------------------
 !Interpolate the velocity at the contour nodes:
call velint(uu,vv,u,v,w)

return
end subroutine

!==========================================================================

subroutine getzet
! Gets the relative vorticity (zz) given the PV anomaly (qq) and 
! eta = h^{1/kappa} -1 (ee).  Ensures that <zz> = 0.
! This requires qq to be adjusted by a constant.

implicit none

double precision:: hsum,zsum,qadd
integer:: i,j

!----------------------------------------------------------------
 !Store work arrays:
do i=1,nt
  do j=1,ng
    wka(j,i)=one+eep(j,i)
    wkb(j,i)=clat(j)*wka(j,i)*(qq(j,i)+cof(j))
  enddo
enddo

 !Carry out 4th-order average:
hsum=zero
zsum=zero
do i=1,nt
  hsum=hsum+f1112*(wka(1,i)+wka(ng,i))
  zsum=zsum+f1112*(wkb(1,i)+wkb(ng,i))
  do j=2,ngm1
    hsum=hsum+wka(j,i)
    zsum=zsum+wkb(j,i)
  enddo
enddo
qadd=-zsum/hsum

 !Adjust PV anomaly (qq) and define relative vorticity (zz):
do i=1,nt
  do j=1,ng
    qq(j,i)=qq(j,i)+qadd
    zz(j,i)=wka(j,i)*(qq(j,i)+cof(j))-cof(j)
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine tends
! Computes the tendencies for qs & qd (sqs & sqd) as well as for eta
! & delta (see & sdd, omitting the linear part used in the
! semi-implicit time stepping), and (when thermal forcing is present)
! for log potential temperature (stt).

! Thermal forcing, if used, is added to all tendencies at the *end*
! of this subroutine, where stt is defined.

! *** Note: wke must contain the spectrally-truncated relative vorticity
!           coming into this subroutine.

implicit none

double precision:: dtat(ng,nt),dton(ng,nt),tflx(ng,nt)
double precision:: eet(ng,nt)
integer:: i,j

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !Compute A = zU + dV/dlon & B = zV - dU/dlon, where
 !U = u/r & V = v/r, while z = sin(lat) and r = cos(lat):
do i=1,nt
  do j=1,ng
    wku(j,i)=clati(j)*uu(j,i)
    wkv(j,i)=clati(j)*vv(j,i)
    wka(j,i)=wkv(j,i)
    wkb(j,i)=wku(j,i)
  enddo
enddo
 !*** Do not re-use wku ***

 !Get semi-spectral coefficients of wka & wkb:
call forfft(ng,nt,wka,trig,factors)  
call forfft(ng,nt,wkb,trig,factors)  

 !Compute longitudinal derivatives spectrally:
call deriv(ng,nt,rk,wka,wkc)
call deriv(ng,nt,rk,wkb,wkd)

 !Recover physical fields:
call revfft(ng,nt,wkc,trig,factors)  
call revfft(ng,nt,wkd,trig,factors)  

 !Complete definition of A & B:
do i=1,nt
  do j=1,ng
    wka(j,i)=slat(j)*wku(j,i)+wkc(j,i)
    wkb(j,i)=slat(j)*wkv(j,i)-wkd(j,i)
  enddo
enddo

 !Calculate latitudinal and longitudinal derivatives of delta:
call latder(ddp,wkc)        !ddp is delta in physical space
 !wkc = d(dd)/dlat
call deriv(ng,nt,rk,dd,wkv) !dd  is delta in spectral space
call revfft(ng,nt,wkv,trig,factors)  
 !wkv = d(dd)/dlon

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 !Compute the delta tendency (sdd), omitting (c^2*Lap - f^2)eta and forcing:

 !First calculate c^2*Lap(h-kappa*eta) --- needed because h, not eta, 
 !appears in the pressure gradient:
sdd=csq*(hhp-kappa*eep)
call laplace(sdd)

 !Next add up all the terms:
do i=1,nt
  do j=1,ng
    sdd(j,i)=cof(j)*(wke(j,i)-cof(j)*eep(j,i)) - bet(j)*uu(j,i) &
            -uu(j,i)**2-vv(j,i)**2 + sdd(j,i) &
            +two*(wka(j,i)*(wke(j,i)-wka(j,i)) &
                 -wkb(j,i)*(ddp(j,i)+wkb(j,i))) &
            -ddp(j,i)**2-wku(j,i)*wkv(j,i)-vv(j,i)*wkc(j,i)
!    sdd =    f*(zeta-f*eta) - beta*u
!           - u^2 - v^2 + c^2*Lap(h-kappa*eta)
!           + 2*[A*(zeta-A) - B*(delta+B)]
!           - div((u,v)*delta)
  enddo
enddo
 !wku = u/r on rhs above (and throughout this subroutine)

 !If bottom topography is present, add on the topographic contribution:
if (topogr) sdd=sdd-topo
 !See spe.f90 (subroutine initialise) for the form of topo.

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 !Compute the eta tendency (see), omitting delta and forcing:

 !Compute -div((u,v)*eta) and store in see:
call latder(eep,wkc)        !eep is eta in physical space
 !wkc = deta/dlat
call deriv(ng,nt,rk,ee,wkv) !ee  is eta in spectral space (needed in deriv)
call revfft(ng,nt,wkv,trig,factors) 
 !wkv = deta/dlon
see=-ddp*eep-wku*wkv-vv(1:ng,:)*wkc

!- - - - - - - - - - - - - - - - - - - - - - - - - -
 !Get source term (sqs) for spectral PV field qs_t:

 !Calculate dq_s/dlon in spectral space:
call deriv(ng,nt,rk,qs,wka)
 !Return to physical space as wka:
call revfft(ng,nt,wka,trig,factors) 
 !Calculate dq_s/dlat in physical space as wkb:
wkc=qs !Need to copy qs here to avoid overwritting the array in revfft
call revfft(ng,nt,wkc,trig,factors) 
call latder(wkc,wkb)
 !Add up source terms (and include beta = df/dlat):
do i=1,nt
  do j=1,ng
    sqs(j,i)=-wku(j,i)*wka(j,i)-vv(j,i)*(wkb(j,i)+bet(j))
  enddo
enddo
 !wku = u/r throughout this subroutine.

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 !Get source term (sqd) for spectral residual PV field qd_t:

 !Calculate dq_d/dlon in spectral space:
call deriv(ng,nt,rk,qd,wka)
 !Return to physical space as wka:
call revfft(ng,nt,wka,trig,factors) 
 !Calculate dq_d/dlat in physical space as wkb:
wkc=qd !Need to copy qd here to avoid overwritting the array in revfft
call revfft(ng,nt,wkc,trig,factors) 
call latder(wkc,wkb)
 !Add up source terms:
sqd=-wku*wka-vv(1:ng,:)*wkb

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !If thermal forcing is present, compute the log potential temperature
 !tendency (stt) and add appropriate forcing to all other tendencies:

if (thermal) then
  eet=one+eep
   !eet = 1 + eta (denominator in definition of PV)
  wkd=alpha*(one+(hhb-hhe)/hhp)
   !wkd = F = alpha*[1 + (h_b - h_e)/h], the thermal forcing

   !Compute tau (log potential temperature) tendency (stt):
  call latder(ttp,dtat)        !ttp is tau in physical space
   !dtat = dtau/dlat
  call deriv(ng,nt,rk,tt,dton) !tt  is tau in spectral space (needed in deriv)
  call revfft(ng,nt,dton,trig,factors) 
   !dton = dtau/dlon
  tflx=wku*dton+vv(1:ng,:)*dtat
   !tflx = (u,v)*grad(tau)
  stt=-aa0*tflx-wkd
   !stt = S_tau = -A_0*(u,v)*grad(tau) - F

   !*** Do not re-use dtat, dton & tflx ***

   !Compute contribution to other tendencies:
  wka=tflx-kappai*wkd
   !wka = (u,v)*grad(tau) - F/kappa; *** ASSUMES 1 - A_0 = kappa ***
  see=see+eet*wka
   !Above (1+eta)*wka is added to S_eta, completing its definition.

  do i=1,nt
    do j=1,ng
      sqd(j,i)=sqd(j,i)-(qq(j,i)+cof(j))*wka(j,i)
       !Above q*wka is subtracted from PV source, S_qd
       !***Note: qq contains the PV anomaly q-f throughout the code.
       !wka is no longer needed.

       !Compute flux terms in the delta forcing:
      wka(j,i)=csq*hhp(j,i)*clati(j)*dton(j,i)-aa0i*wkd(j,i)*uu(j,i)
      wkb(j,i)=clat(j)*(csq*hhp(j,i)*dtat(j,i)-aa0i*wkd(j,i)*vv(j,i))
       !(wka,wkb/r) = c^2*h*grad(tau) - F(u,v)/A_0 where r = cos(lat);
       !    *** ASSUMES 1 - A_0 = kappa ***
    enddo
  enddo
   !Compute div(wka,wkb/r):
  call latder(wkb,wkc)                !wkc = dwkb/dlat
  call forfft(ng,nt,wka,trig,factors) 
  call deriv(ng,nt,rk,wka,wke)
  call revfft(ng,nt,wke,trig,factors) !wke = dwka/dlon
   !Add this after weighting by 1/r in definition of div to S_d:
  do i=1,nt
    do j=1,ng
      sdd(j,i)=sdd(j,i)+clati(j)*(wkc(j,i)+wke(j,i))
    enddo
  enddo
   !This completes the definition of the delta source, S_d.

   !Redefine wkd as F*(theta_0/bar_theta_0) = F*exp(tau) = F_hat:
  wkd=wkd*exp(ttp)
   !Prepare to calculate curl(F_hat*(u,v)):
  do i=1,nt
    do j=1,ng
      wka(j,i)=clat(j)*wkd(j,i)*uu(j,i)  !wka = r*F_hat*u
      wkb(j,i)=wkd(j,i)*vv(j,i)          !wkb =   F_hat*v
    enddo
  enddo
  call latder(wka,wkc)                !wkc = d(r*F_hat*u)/dlat
  call forfft(ng,nt,wkb,trig,factors) 
  call deriv(ng,nt,rk,wkb,wke)
  call revfft(ng,nt,wke,trig,factors) !wke = d(F_hat*v)/dlon
   !Complete definition of the curl after weighting by 1/r and subtract
   !after division by A_0*(1+eta) from the PV sources, S_qd:
  do i=1,nt
    do j=1,ng
      sqd(j,i)=sqd(j,i)+aa0i*clati(j)*(wke(j,i)-wkc(j,i))/eet(j,i)
    enddo
  enddo

   !Finally add the Jacobian term c^2*J(h,tau)/(1+eta) to the PV source:
  call latder(hhp,wkc)                !wkc = dh/dtau
  wka=hhp !(Need to copy hhp to avoid overwritting in forfft)
  call forfft(ng,nt,wka,trig,factors) 
  call deriv(ng,nt,rk,wka,wke)
  call revfft(ng,nt,wke,trig,factors) !wke = dh/dlon
   !Use J(h,tau) = (dh/dlon*dtau/dlat - dh/dlat*dtau/dlon)/r:
  do i=1,nt
    do j=1,ng
      sqd(j,i)=sqd(j,i)+csq*clati(j)* &
            (wke(j,i)*dtat(j,i)-wkc(j,i)*dton(j,i))/eet(j,i)
    enddo
  enddo
   !This completes the definition of the PV source, S_d.
endif

 !De-aliase tendencies and convert to semi-spectral space:
call dealiase(sdd)
call dealiase(see)
call dealiase(sqs)
call dealiase(sqd)
if (thermal) call dealiase(stt)

return
end subroutine

!========================================================================

subroutine sistep
! Updates eta & delta given R_eta & R_delta (in the arrays wkc & wkd) 
! as well as the original values of eta & delta (in the arrays eeo
! & ddo).  Here, R_eta = (see^n+see^{n+1})/2+(2/dt)*ee^n where n
! refers to the current time level and n+1 refers to the next one.

implicit none

double precision:: avdd
integer:: j

 !Re-define R_delta:
wkd=csqi*(wkd-omv*wkc)
 !omv = 2/dt + longitudinal hyperviscous operator (see adapt).

 !Invert Helmholtz operator on wkd to get (eta^n+eta^{n+1})/2:
call helminv(wkd,wka,wkb)

 !Finish calculation to get new eta and delta fields:
ee=two*wka-eeo
dd=two*(wkc-omv*wka)-ddo

 !Remove global mean value of delta:
avdd=f1112*(dd(1,1)*clat(1)+dd(ng,1)*clat(ng))
do j=2,ngm1
  avdd=avdd+dd(j,1)*clat(j)
enddo
avdd=avdd*rsumi

do j=1,ng
  dd(j,1)=dd(j,1)-avdd
enddo

 !Update eta, h and delta in physical space:
eep=ee
call revfft(ng,nt,eep,trig,factors) 
hhp=(one+eep)**kappa
ddp=dd
call revfft(ng,nt,ddp,trig,factors)

return
end subroutine

!=======================================================================

subroutine adapt
! Adapts the time step to ensure 
!    dt < cflmax*dl/|u|_max, 
!    dt < dtzz/|z|_max, and
!    dt < dtmax,
! where |z|_max is the maximum relative vorticity, dtzz = pi/10,
! cflmax = 0.7, dl = longitudinal grid spacing, and dtmax is 
! provided in params.dat.  Also adjusts the damping rate of
! residual PV and delta.

implicit none

double precision,parameter:: cflmax=0.7d0, dtzz=pi/10.d0
double precision,parameter:: fdtmin=0.9d0, fdtmax=one/fdtmin
! *** The time step is adapted only if dtnew < fdtmin*dt or 
!                                      dtnew > fdtmax*dt
double precision:: ddrms,zzrms,enstro,umax,usq
double precision:: hhmin,hhmax,zzmax,frmax
double precision:: dtcfl,dtacc,dtnew,trem,cfl,sm,pm
double precision:: a0,b0,c0,d0,deti,cphe,dphe
integer:: i,j,k,m
logical:: change

!---------------------------------------------------------------------
 !Get rms delta and zeta, as well as the enstrophy:
call getrms(ddp,ddrms)
call getrms(zz,zzrms)
enstro=f12*zzrms**2

 !Compute cfl time step and other diagnostics:
umax=small
hhmin=one
hhmax=one
zzmax=small
frmax=zero
do i=1,nt
  do j=1,ng
    usq=uu(j,i)**2+vv(j,i)**2
    umax=max(umax,usq)
    hhmin=min(hhmin,hhp(j,i))
    hhmax=max(hhmax,hhp(j,i))
    zzmax=max(zzmax,abs(zz(j,i)))
    frmax=max(frmax,usq/hhp(j,i))
  enddo
enddo
umax=sqrt(umax)
frmax=sqrt(frmax)/cgw

 !CFL time step: appears to work despite shrinking grid lengths
dtcfl=dl*cflmax/umax

 !Accurate advection time step:
dtacc=dtzz/zzmax

!---------------------------------------------------------------------
 !Choose the smallest time step:
dtnew=min(dtmax,dtcfl,dtacc)

if ((dtnew .lt. fdtmin*dt) .or. (dtnew .gt. fdtmax*dt)) then
  dt=dtnew
  change=.true.
else
  change=.false.
endif

 !See how close we are to the end of the simulation:
trem=tmax-t

if (trem .lt. dt .and. trem .gt. 1.d-6*dtmax) then
  dt=trem
  change=.true.
endif

!---------------------------------------------------------------------
 !Increment the integral of max|zz|:
twist=twist+dt*zzmax

 !Record cfl parameter, enstrophy, |Ro|_max, etc in cfl.dat:
cfl=umax*dt/dl
write(18,'(f10.4,1x,f6.4,7(1x,f8.5))') & 
    t,cfl,enstro,zzmax/fpole,hhmin,hhmax,frmax,ddrms/zzrms,dqdthi

if (.not. change) return
!----------------------------------------------------------------
 !Change the time step & define arrays needed for semi-implicit 
 !time-stepping and hyperviscous damping of eta & delta:
dt2=f12*dt
hdti=one/dt2
 !om = hdti below

 !Initialise Helmholtz block-tridiagonal inversion:
do m=1,nt
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  do j=1,ng
    omv(j,m)=hdti+dislon(j,m)
    simp(j)=d2lon(j,m)+(fsq(j)+omv(j,m)**2)*csqi
  enddo

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*simp(1)
  d0=-pm*tlat(1)

  deti=one/(a0*d0-b0*c0)
  ahe0(1,m)=a0*deti
  bhe0(1,m)=b0*deti
  che0(1,m)=c0*deti
  dhe0(1,m)=d0*deti

  cphe=pc12-0.1d0*simp(2)
  dphe=-0.1d0*tlat(2)
  ehe1(1,m)= cphe*bhe0(1,m)- apla*dhe0(1,m)
  ehe2(1,m)= dphe*bhe0(1,m)-0.2d0*dhe0(1,m)
  fhe1(1,m)= apla*che0(1,m)- cphe*ahe0(1,m)
  fhe2(1,m)=0.2d0*che0(1,m)- dphe*ahe0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-simp(j)
    d0=-tlat(j)
    cmhe(j,m)=pc12-0.1d0*simp(j-1)
    dmhe(j)=-0.1d0*tlat(j-1)

    ahe0(j,m)=   ehe1(j-1,m)*amla      +fhe1(j-1,m)*0.2d0
    bhe0(j,m)=b0+ehe2(j-1,m)*amla      +fhe2(j-1,m)*0.2d0
    che0(j,m)=c0+ehe1(j-1,m)*cmhe(j,m)+fhe1(j-1,m)*dmhe(j)
    dhe0(j,m)=d0+ehe2(j-1,m)*cmhe(j,m)+fhe2(j-1,m)*dmhe(j)
    deti=one/(ahe0(j,m)*dhe0(j,m)-bhe0(j,m)*che0(j,m))
    ahe0(j,m)=ahe0(j,m)*deti
    bhe0(j,m)=bhe0(j,m)*deti
    che0(j,m)=che0(j,m)*deti
    dhe0(j,m)=dhe0(j,m)*deti

    cphe=pc12-0.1d0*simp(j+1)
    dphe=-0.1d0*tlat(j+1)
    ehe1(j,m)= cphe*bhe0(j,m)- apla*dhe0(j,m)
    ehe2(j,m)= dphe*bhe0(j,m)-0.2d0*dhe0(j,m)
    fhe1(j,m)= apla*che0(j,m)- cphe*ahe0(j,m)
    fhe2(j,m)=0.2d0*che0(j,m)- dphe*ahe0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*simp(ng)
  d0=-pm*tlat(ng)

  cmhe(ng,m)=pc12-0.1d0*simp(ngm1)
  dmhe(ng)=-0.1d0*tlat(ngm1)

  ahe0(ng,m)=a0+ehe1(ngm1,m)*amla       +fhe1(ngm1,m)*0.2d0
  bhe0(ng,m)=b0+ehe2(ngm1,m)*amla       +fhe2(ngm1,m)*0.2d0
  che0(ng,m)=c0+ehe1(ngm1,m)*cmhe(ng,m)+fhe1(ngm1,m)*dmhe(ng)
  dhe0(ng,m)=d0+ehe2(ngm1,m)*cmhe(ng,m)+fhe2(ngm1,m)*dmhe(ng)

  deti=one/(ahe0(ng,m)*dhe0(ng,m)-bhe0(ng,m)*che0(ng,m))
  ahe0(ng,m)=ahe0(ng,m)*deti
  bhe0(ng,m)=bhe0(ng,m)*deti
  che0(ng,m)=che0(ng,m)*deti
  dhe0(ng,m)=dhe0(ng,m)*deti
enddo

 !Latitudinal hyperviscous damping used on eta & delta:
do k=1,nt
  hyplat(k)=one/(one+dt*dislat(k))
enddo

return
end subroutine

!========================================================================

subroutine genpv
! Computes the PV field by combining the contour PV with that 
! in the fields qs & qd.

! On return, qq contains the full gridded PV anomaly BUT
! with an incorrect mean value.  The mean value of qq is 
! found in getzet by the requirement that the global mean
! absolute vorticity remain zero.

implicit none

integer:: i,j

!--------------------------------

 !Grid the PV contours as qq (the PV anomaly, qc - f):
call con2grid(qc)

 !Subtract Coriolis frequency to define anomaly (qq):
do i=1,nt
  do j=1,ng
    qq(j,i)=qc(j,i)-cof(j)
  enddo
enddo

 !FFT qq in longitude:
call forfft(ng,nt,qq,trig,factors) 
 !qq contains the PV anomaly arising from the contours in 
 !semi-spectral space (called "qc" in the comments below)

 !Store qs-qc in wka in preparation for filtering:
wka=qs-qq

 !Apply low-pass filter F to qs-qc:
call lofilter(wka)
 !wka = F(qs-qc) (now semi-spectral)

 !Define PV anomaly from qq = F(qs-qc)+qd+qc:
qq=wka+qd+qq
 !Return qq to physical space:
call revfft(ng,nt,qq,trig,factors)

return
end subroutine

!========================================================================

subroutine reset
! Resets qs & qd at the beginning of each time step, or
! just before recontouring at the end of a run.

! On return, qq contains the full gridded PV anomaly BUT
! with an incorrect mean value.  The mean value of qq is 
! found in getzet by the requirement that the global mean
! absolute vorticity remains zero.

implicit none

integer:: i,j

!------------------------------------------------------------------------
 !Grid the PV contours as qq (the PV anomaly, qc - f):
call con2grid(qc)

 !Subtract Coriolis frequency to define anomaly (qq):
do i=1,nt
  do j=1,ng
    qq(j,i)=qc(j,i)-cof(j)
  enddo
enddo

 !FFT qq in longitude:
call forfft(ng,nt,qq,trig,factors) 
 !qq contains the PV anomaly arising from the contours in 
 !semi-spectral space (called "qc" in the comments below)

 !Store qs-qc in wka in preparation for filtering:
wka=qs-qq

 !Apply low-pass filter F to qs-qc:
call lofilter(wka)
 !wka = F(qs-qc) (now semi-spectral)

 !Add on qd and obtain reset value of qs <-- F(qs-qc)+qd+qc:
qd=wka+qd
qs=qd+qq
qq=qs
 !Note: now qd = *new* qs-qc

 !Return qq (the full PV anomaly) to physical space:
call revfft(ng,nt,qq,trig,factors)

 !Reset qd = (1-F)(qs-qc):
call hifilter(qd)

return
end subroutine

!========================================================================

subroutine dump
! Writes various fields to main job directory and also computes 
! and writes the total energy and angular momentum to ene-ang.asc.

implicit none

double precision:: powhh(ngp1),powdd(ngp1),powqs(ngp1),powqd(ngp1)
double precision:: r1,r2,hht,eke,epe,ang
real:: qqr4(ng,nt),tr4
integer:: i,j,m

!------------------------------------------------------------
 !Increment counter for dumping data:
idump=idump+1

 !Store single precision PV field for output below:
do i=1,nt
  do j=1,ng
    qqr4(j,i)=real(qq(j,i)+cof(j))
  enddo
enddo
tr4=real(t)

 !Write fields for later imaging or diagnostics:
write(40,rec=idump) tr4,real(hhp)
if (thermal) write(41,rec=idump) tr4,real(ttp)
write(42,rec=idump) tr4,real(uu(1:ng,:))
write(43,rec=idump) tr4,real(vv(1:ng,:))
write(44,rec=idump) tr4,real(zz)
write(45,rec=idump) tr4,qqr4
write(46,rec=idump) tr4,real(ddp)

!------------------------------------------------------------------
 !Compute and write energy and angular momentum per unit area:
do i=1,nt
  do j=1,ng
    r1=clat(j)
    r2=r1*r1
    wka(j,i)=r1*hhp(j,i)*(uu(j,i)**2+vv(j,i)**2)
    wkb(j,i)=r1*csq*(hhp(j,i)-one)**2
    wkc(j,i)=r2*(hhp(j,i)*uu(j,i)+omega*r1*(hhp(j,i)-one))
  enddo
enddo

eke=zero
epe=zero
ang=zero
do i=1,nt
  eke=eke+f1112*(wka(1,i)+wka(ng,i))
  epe=epe+f1112*(wkb(1,i)+wkb(ng,i))
  ang=ang+f1112*(wkc(1,i)+wkc(ng,i))
  do j=2,ngm1
    eke=eke+wka(j,i)
    epe=epe+wkb(j,i)
    ang=ang+wkc(j,i)
  enddo
enddo
eke=f12*eke*dsumi
epe=f12*epe*dsumi
ang=    ang*dsumi

write(17,'(f12.5,4(1x,f15.11))') t,eke,epe,eke+epe,ang

!----------------------------------------------------------------------
 !Compute and write longitudinal power spectra for h, delta, q_s & q_d:
call forfft(ng,nt,wka,trig,factors) !wka now contains h in spectral space
call longspec(wka,powhh)
call longspec(dd,powdd)
call longspec(qs,powqs)
call longspec(qd,powqd)

write(50,'(f12.5)') t
do m=1,ngp1
  write(50,'(5(1x,f12.8))') log10(dble(m)),log10(powhh(m)+1.d-32), &
                    log10(powdd(m)+1.d-32),log10(powqs(m)+1.d-32), &
                    log10(powqd(m)+1.d-32)
enddo

return
end subroutine

!========================================================================

subroutine writedata
! This routine writes the current contours and the residual PV
! to the cont subdirectory and eta, u, v etc to grid subdirectory

implicit none

real:: qqr4(ng,nt),tr4
integer:: i,j,loop,irec
character(len=3):: ofile

!----------------------------------------------------------------
 !Re-initialise qs & qd before recontouring:
 !    Reset qs = q = F*qs + (1-F)*qc + qd
 !      and qd = (1-F)*(q-qc)
call reset
 !qq contains the full gridded PV anomaly after this call.
 !qd contains the *spectral* residual PV.

 !Correct qq so that relative vorticity has zero mean:
call getzet

 !Define residual needed for recontouring:
do i=1,nt
  do j=1,ng
    qqr4(j,i)=real(qq(j,i)-qc(j,i)+cof(j))
  enddo
enddo
 !Note: qc is the total PV coming from the contours.

!------------------------------------------------------------------------

 !This is the end of a normal period; just write data:
loop=nint(t/tsim)
irec=loop+1
ofile='000'
write(ofile(1:3),'(i3.3)') loop

tr4=real(t)
 !Write contours to the cont subdirectory:
write(80,'(1x,f12.5,1x,i9,1x,i10)') t,n,npt

open(81,file='cont/index'//ofile,form='unformatted',status='replace')
write(81) np(1:n),i1(1:n),ind(1:n)
close(81)

open(82,file='cont/nodes'//ofile,form='unformatted',status='replace')
write(82) x(1:npt),y(1:npt),z(1:npt)
close(82)

 !Write residual PV, F(qs-qc)+qd contained in qqr4:
write(83,rec=irec) tr4,qqr4

 !Store single precision array of full PV field:
do i=1,nt
  do j=1,ng
    qqr4(j,i)=real(qq(j,i)+cof(j))
  enddo
enddo

 !Write fields for later imaging or diagnostics:
write(30,rec=irec) tr4,real(hhp)
if (thermal) write(31,rec=irec) tr4,real(ttp)
write(32,rec=irec) tr4,real(uu(1:ng,:))
write(33,rec=irec) tr4,real(vv(1:ng,:))
write(34,rec=irec) tr4,real(zz)
write(35,rec=irec) tr4,qqr4
write(36,rec=irec) tr4,real(ddp)

 !Write final data for imaging with mgrid, llgrid etc:
if (loop .eq. nperiod) call dump

return
end subroutine

!=======================================================================

 !Main end module
end module
