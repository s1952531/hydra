module evolution

! Module contains subroutines to evolve fields as detailed in spe.f90.

use common

implicit none

 !PV anomaly & relative vorticity:
double precision:: qq(ng,nt),zz(ng,nt)

 !PV associated with the contours:
double precision:: qc(ng,nt)

 !Target PV anomaly at the end of a ramp period (if ramping is used):
double precision:: qr(ng,nt)

 !Velocity fields and contour velocities:
double precision:: uu(0:ng+1,nt),vv(0:ng+1,nt)
double precision:: u(npm),v(npm),w(npm)

 !Semi-implicit hyperviscous damping operator (semi-spectral):
double precision:: omv(ng,nt)

 !Global scalar variables needed only in this module:
double precision:: twist,hdti,dt2

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
if (.not. ramping) then 
  call init
else
   !Define qr & initialise qq to zero:
  do i=1,nt
    do j=1,ng
      qr(j,i)=qs(j,i)-cof(j)
      qq(j,i)=zero
    enddo
  enddo
endif  

 !Initialise counter for regularisations:
nreg=0

 !Used for regularising contours:
twist=zero

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (nint(tmax/tsim) .le. nper) 
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

implicit double precision(a-h,o-z)
implicit integer(i-n)

!------------------------------------------------------------------------
 !Update gridded version of contour PV:
call con2grid(qc)

 !Take away Coriolis frequency from qs & qc and prepare 
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
 !qs contains the entire PV anomaly in semi-spectral space

 !Calculate the residual PV from qd = (1-F)*(q-qc):
do m=1,nt
  do j=1,ng
    qd(j,m)=qs(j,m)-qq(j,m)
    qq(j,m)=qs(j,m)
  enddo
enddo

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
! Prepares PV fields for recontouring

implicit double precision(a-h,o-z)
implicit integer(i-n)

!-------------------------------------------
 !Re-initialise qs & qd before recontouring:
 !    Reset qs = q = F*qs + (1-F)*qc + qd
 !      and qd = (1-F)*(q-qc)
call reset
 !qq contains the full gridded PV anomaly after this call.
 !qd contains the *spectral* residual PV.

 !Correct qq so that relative vorticity has zero mean:
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

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define local parameters and arrays:
integer,parameter:: niter=2
!   niter+1: number of interations of the implicit scheme

 !Source terms:
double precision:: shh(ng,nt),sdd(ng,nt),sqs(ng,nt),sqd(ng,nt)

 !Intermediate fields:
double precision:: qso(ng,nt),qdo(ng,nt),hho(ng,nt),ddo(ng,nt)

 !Source terms in semi-implicit time stepping of h & delta:
double precision:: rhh(ng,nt),rhho(ng,nt),rdd(ng,nt),rddo(ng,nt)

 !Intermediate contour positions:
double precision:: xo(npm),yo(npm),zo(npm)

!------------------------------------------------------------------------
 !Re-initialise qs & qd at the beginning of the time step:
 !    Reset qs = q = F*qs + (1-F)*qc + qd
 !      and qd = (1-F)*(q-qc)
if (.not. ramping) call reset
 !qq contains the full spectral PV anomaly after this call.

 !Find gridded velocity and interpolate it at each contour node:
call invert

 !Save various fields approximately every tdsave time units:
itime=int(t/tdsave+small)
tcrit=tdsave*dble(itime)
if ((t-dt-tcrit)*(t+small*tdsave-tcrit) .lt. zero) then
  if (itime .gt. iref) then
    call dump
    iref=itime
  endif
endif

 !Possibly adapt the time step and adjust damping rates:
call adapt

 !Find tendencies:
call tends(shh,sdd,sqs,sqd)

!--------------------------------------------------------------------
!$omp parallel

 !Store quantities at the beginning of a time step needed in the
 !iteration below and take first forward time step:
if (.not. ramping) then

   !Contour nodes:
!$omp do private(hdx,hdy,hdz,xnew,ynew,znew,fac)
  do i=1,npt
    hdx=dt2*u(i)
    hdy=dt2*v(i)
    hdz=dt2*w(i)
    xo(i)=x(i)+hdx
    yo(i)=y(i)+hdy
    zo(i)=z(i)+hdz
    xnew=xo(i)+hdx
    ynew=yo(i)+hdy
    znew=zo(i)+hdz
    fac=one/sqrt(xnew**2+ynew**2+znew**2)
    x(i)=fac*xnew
    y(i)=fac*ynew
    z(i)=fac*znew
  enddo
!$omp enddo nowait

   !Semi-spectral PV fields:
!$omp do private(hdqs,hdqd)
  do m=1,nt
    do j=1,ng
      hdqs=dt2*sqs(j,m)
      qso(j,m)=qs(j,m)+hdqs
      qs(j,m)=qso(j,m)+hdqs

      hdqd=dt2*sqd(j,m)
      qdo(j,m)=qd(j,m)+hdqd
      qd(j,m)=qdo(j,m)+hdqd
    enddo
  enddo
!$omp enddo nowait

endif

 !Height and divergence:
!$omp do private(hshh,hsdd)
do m=1,nt
  do j=1,ng
    hho(j,m)=hh(j,m)
    hshh=f12*shh(j,m)
    rhho(j,m)=hshh+hdti*hh(j,m)
    rhh(j,m)=hshh+rhho(j,m)

    ddo(j,m)=dd(j,m)
    hsdd=f12*sdd(j,m)
    rddo(j,m)=hsdd+hdti*dd(j,m)
    rdd(j,m)=hsdd+rddo(j,m)
  enddo
enddo
!$omp enddo
!$omp end parallel

 !Get approximate height & divergence at t + dt:
call sistep(hho,ddo,rhh,rdd)

!--------------------------------------------------------
 !Iterate now to solve the implicit problem in time:
do iter=1,niter

   !Generate the PV anomaly at t+dt from the contours & PV fields:
  call genpv

   !Update velocity at this time and interpolate it to the contours:
  call invert

   !Update tendencies at this time:
  call tends(shh,sdd,sqs,sqd)

!$omp parallel

   !Update contours and PV fields if not ramping:
  if (.not. ramping) then
!$omp do private(xnew,ynew,znew,fac)
    do i=1,npt
      xnew=xo(i)+dt2*u(i)
      ynew=yo(i)+dt2*v(i)
      znew=zo(i)+dt2*w(i)
      fac=one/sqrt(xnew**2+ynew**2+znew**2)
      x(i)=fac*xnew
      y(i)=fac*ynew
      z(i)=fac*znew
    enddo
!$omp enddo
!$omp do
    do m=1,nt
      do j=1,ng
        qs(j,m)=qso(j,m)+dt2*sqs(j,m)
        qd(j,m)=qdo(j,m)+dt2*sqd(j,m)
      enddo
    enddo
!$omp enddo
  endif

   !Define R_h & R_d:
!$omp do
  do m=1,nt
    do j=1,ng
      rhh(j,m)=f12*shh(j,m)+rhho(j,m)
      rdd(j,m)=f12*sdd(j,m)+rddo(j,m)
    enddo
  enddo
!$omp enddo
!$omp end parallel
 
   !Update height & divergence at t + dt:
  call sistep(hho,ddo,rhh,rdd)
enddo

 !Apply hyperviscous damping in latitude only to h & delta:
call latdamp(hh)
call latdamp(dd)

 !Update height & divergence in physical space:
do m=1,nt
  do j=1,ng
    hhp(j,m)=hh(j,m)
    ddp(j,m)=dd(j,m)
  enddo
enddo
call revfft(ng,nt,hhp,trig,factors) 
call revfft(ng,nt,ddp,trig,factors) 

 !Increment time:
t=t+dt

!-----------------------------------------------------------------------
if (ramping) then
   !If there is PV ramping, see if we are at the end of it:
  if (t/tramp .gt. oms) then
     !We've reached the end of the ramp period, t = tramp;
     !Prepare for restart:
    t=zero
    tdsave=tsavori
    tmax =tmaxori
    iref=-1
    nper=nperiod
     !Ensure ramping is off:
    ramping=.false.
    do i=1,nt
      do j=1,ng
        qs(j,i)=qr(j,i)+cof(j)
      enddo
    enddo
    call init
  endif

else

   !Add any stochastic forcing to qd here:
  if (stoch) then
     !rms value of added vorticity:
    dzrms=sqrt(two*esr*dt)
     !generate random field qso with this rms value:
    call ranspec(qso,dzrms)

     !Convert to PV forcing:
    do i=1,nt
      do j=1,ng
        qso(j,i)=qso(j,i)/(one+hhp(j,i))
      enddo
    enddo

     !FFT qso and add to qd in semi-spectral space:
    call forfft(ng,nt,qso,trig,factors) 
    do m=1,nt
      do j=1,ng
        qd(j,m)=qd(j,m)+qso(j,m)
      enddo
    enddo
  endif

endif

return
end subroutine

!==========================================================================

subroutine invert
! Calculates the velocity field (uu,vv) given the height (hh), 
! divergence (dd) and PV anomaly field (qq).  The velocity field 
! is then interpolated at all nodes as (u,v).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Work arrays:
double precision:: pp(ng,nt),wka(ng,nt),wkb(ng,nt),wku(ng,nt),wkv(ng,nt)

!---------------------------------------------------------
 !First compute divergent part of the velocity field:

 !Invert Laplace's operator on the divergence:
call laplinv(dd,pp,wkv)
 !Here the divergence potential is pp while d(pp)/dlat = wkv,
 !the divergent meridional velocity.

 !Compute and store divergent zonal velocity (times cos(lat)) in wku:
call deriv(ng,nt,rk,pp,wku) 
 !wku = d(pp)/dlon

 !Compute relative vorticity zz = (1+hh)*(qq+f) - f:
call getzet

 !Convert the vorticity to semi-spectral space as wka:
do i=1,nt
  do j=1,ng
    wka(j,i)=zz(j,i)
  enddo
enddo
call forfft(ng,nt,wka,trig,factors)  

 !Invert Laplace's operator:
call laplinv(wka,pp,wkb)
 !Here the streamfunction is pp while d(pp)/dlat = wkb.

 !Compute d(pp)/dlon = wka:
call deriv(ng,nt,rk,pp,wka)

 !Complete calculation of zonal and meridional velocity, wku & wkv:
do m=1,nt
  do j=1,ng
    wku(j,m)=clati(j)*wku(j,m)-wkb(j,m)
    wkv(j,m)=clati(j)*wka(j,m)+wkv(j,m)
  enddo
enddo

 !Get physical space velocity:
call revfft(ng,nt,wku,trig,factors)
call revfft(ng,nt,wkv,trig,factors)  

 !Copy into larger arrays for use in interpolation of contour velocities:
do i=1,nt
  do j=1,ng
    uu(j,i)=wku(j,i)
    vv(j,i)=wkv(j,i)
  enddo
enddo

!----------------------------------------------------
 !Interpolate the velocity at the contour nodes:
if (.not. ramping) call velint(uu,vv,u,v,w)
 !(This is not done during initialisation using PV anomaly ramping)

return
end subroutine

!==========================================================================

subroutine getzet
! Gets the relative vorticity (zz) given the PV anomaly (qq) and 
! the height fields (hh).  Ensures that <zz> = 0.  This requires
! qq to be adjusted by a constant.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Work arrays:
double precision:: wka(ng,nt),wkb(ng,nt)

!----------------------------------------------------------------
 !Store work arrays:
do i=1,nt
  do j=1,ng
    wka(j,i)=one+hhp(j,i)
    wkb(j,i)=clat(j)*wka(j,i)*(qq(j,i)+cof(j))
  enddo
enddo

 !Carry out 4th-order average:
zsum=zero
do i=1,nt
  zsum=zsum+f1112*(wkb(1,i)+wkb(ng,i))
  do j=2,ngm1
    zsum=zsum+wkb(j,i)
  enddo
enddo
zsum=zsum*dsumi

 !Adjust PV anomaly (qq) and define vorticity (zz):
do i=1,nt
  do j=1,ng
    qq(j,i)=qq(j,i)-zsum
    zz(j,i)=wka(j,i)*(qq(j,i)+cof(j))-cof(j)
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine tends(shh,sdd,sqs,sqd)
! Computes the tendencies for qs & qd (sqs & sqd) as well as
! for the height & divergence (shh & sdd, omitting the linear
! part used in semi-implicit time stepping).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Source terms:
double precision:: shh(ng,nt),sdd(ng,nt),sqs(ng,nt),sqd(ng,nt)

 !Work arrays:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt),wkd(ng,nt)
double precision:: wku(ng,nt),wkv(ng,nt)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! First compute divergence and height tendencies:

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
 !Do not re-use wku

 !Get semi-spectral coefficients of wka & wkb:
call forfft(ng,nt,wka,trig,factors)  
call forfft(ng,nt,wkb,trig,factors)  

 !Compute longitudinal derivatives:
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
 !*** Do not re-use wka or wkb ***

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 !Calculate latitudinal and longitudinal derivatives of divergence:
call deriv(ng,nt,rk,dd,wkv)
call revfft(ng,nt,wkv,trig,factors)  
 !wkv = d(dd)/dlon
call latder(ddp,wkc)
 !wkc = d(dd)/dlat

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 !Compute the velocity divergence tendency (sdd), omitting
 !(c^2*Lap - f^2)h:
do i=1,nt
  do j=1,ng
    sdd(j,i)=cof(j)*(zz(j,i)-cof(j)*hhp(j,i)) - bet(j)*uu(j,i) &
          &  -uu(j,i)**2-vv(j,i)**2 &
          &  +two*(wka(j,i)*( zz(j,i)-wka(j,i)) &
          &       -wkb(j,i)*(wkb(j,i)+ddp(j,i))) &
          &  -ddp(j,i)**2-wku(j,i)*wkv(j,i)-vv(j,i)*wkc(j,i)
!    sdd=      f*(zeta-f*h) - beta*u
!                    - u^2 - v^2
!                    + 2*[A*(zeta-A) - B*(delta+B)]
!                    - div((u,v)*delta)
  enddo
enddo
 !wku = u/r on rhs above

if (topogr) then
 !Add on the topographic contribution:
  topoamp=atopo+btopo*sin(ftopo*t)
  do i=1,nt
    do j=1,ng
      sdd(j,i)=sdd(j,i)-topoamp*ttdd(j,i)
    enddo
  enddo
endif

 !De-aliase tendency and convert to semi-spectral space:
call dealiase(sdd)

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 !Compute the height tendency (shh), omitting delta:

 !First compute div((u,v)*h):
call latder(hhp,wkc)
 !wkc = dh/dlat
call deriv(ng,nt,rk,hh,wkv)
call revfft(ng,nt,wkv,trig,factors) 
 !wkv = dh/dlon

if (thermal) then
   !Define possibly time dependent thermal relaxation coefficient:
  rth=rtherm*(atherm+rtherm*t)/(one+rtherm*t)

   !Set shh = rth*(h_eq - h) - div((u,v)*h):
  do i=1,nt
    do j=1,ng
      shh(j,i)=rth*(hhe(j,i)-hhp(j,i)) &
             & -ddp(j,i)*hhp(j,i) &
             & -wku(j,i)*wkv(j,i)-vv(j,i)*wkc(j,i)
    enddo
  enddo

else
   !Set shh = -div((u,v)*h):
  do i=1,nt
    do j=1,ng
      shh(j,i)=-ddp(j,i)*hhp(j,i) &
             & -wku(j,i)*wkv(j,i)-vv(j,i)*wkc(j,i)
    enddo
  enddo
endif

 !De-aliase tendency and convert to semi-spectral space:
call dealiase(shh)

if (ramping) return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !Get source term (sqs) for spectral PV field qs_t:
do m=1,nt
  do j=1,ng
    wkc(j,m)=qs(j,m)
  enddo
enddo

 !Calculate longitudinal derivative spectrally:
call deriv(ng,nt,rk,wkc,wka)
 !Return to physical space:
call revfft(ng,nt,wka,trig,factors) 
 !Calculate latitudinal derivative in physical space:
call revfft(ng,nt,wkc,trig,factors) 
call latder(wkc,wkb)

 !Get NL tendency for qs:
do i=1,nt
  do j=1,ng
    sqs(j,i)=-uu(j,i)*clati(j)*wka(j,i)-vv(j,i)*(wkb(j,i)+bet(j))
  enddo
enddo

 !De-aliase tendency and convert to semi-spectral space:
call dealiase(sqs)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 !Get source term (sqd) for spectral residual vorticity field qd_t:
do m=1,nt
  do j=1,ng
    wkc(j,m)=qd(j,m)
  enddo
enddo

 !Calculate longitudinal derivative spectrally:
call deriv(ng,nt,rk,wkc,wka)
 !Return to physical space:
call revfft(ng,nt,wka,trig,factors) 

 !Calculate latitudinal derivative in physical space:
call revfft(ng,nt,wkc,trig,factors) 
call latder(wkc,wkb)

 !Get NL tendency for qd:
if (thermal) then
   !Add thermal relaxation:
  do i=1,nt
    do j=1,ng
      sqd(j,i)=-uu(j,i)*clati(j)*wka(j,i)-vv(j,i)*wkb(j,i) &
             & +rth*(hhp(j,i)-hhe(j,i))*(zz(j,i)+cof(j))/ &
             &      (one+hhp(j,i))**2
    enddo
  enddo
else if (friction) then
   !Add Ekman friction:
  do i=1,nt
    do j=1,ng
      sqd(j,i)=-uu(j,i)*clati(j)*wka(j,i)-vv(j,i)*wkb(j,i) &
             & -rekman*zz(j,i)/(one+hhp(j,i))
    enddo
  enddo
else
   !Conservative case:
  do i=1,nt
    do j=1,ng
      sqd(j,i)=-uu(j,i)*clati(j)*wka(j,i)-vv(j,i)*wkb(j,i)
    enddo
  enddo
endif

 !De-aliase tendency and convert to semi-spectral space:
call dealiase(sqd)

return
end subroutine

!========================================================================

subroutine sistep(hho,ddo,rhh,rdd)
! Updates height & divergence given R_h & R_delta (in the arrays rhh & rdd) 
! as well as the original values of height & divergence (in the arrays hho
! & ddo).  Here, R_h = (shh^n+ssh^{n+1})/2+(2/dt)*hh^n where n
! refers to the current time level and n+1 refers to the next one.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: hho(ng,nt),ddo(ng,nt),rhh(ng,nt),rdd(ng,nt)

 !Work arrays:
double precision:: wka(ng,nt),wkb(ng,nt)

!-------------------------------------------------------------------
 !Re-define R_delta:
do m=1,nt
  do j=1,ng
    rdd(j,m)=csqi*(rdd(j,m)-omv(j,m)*rhh(j,m))
  enddo
enddo
 !omv = 2/dt + longitudinal hyperviscous operator (see adapt).

 !Invert Helmholtz operator on rdd to get (h^n+h^{n+1})/2:
call helminv(rdd,wka,wkb)

 !Finish calculation to get new height and divergence fields:
do m=1,nt
  do j=1,ng
    hh(j,m)=two*wka(j,m)-hho(j,m)
    dd(j,m)=two*(rhh(j,m)-omv(j,m)*wka(j,m))-ddo(j,m)
  enddo
enddo

 !Remove global mean values:
avdd=f1112*(dd(1,1)*clat(1)+dd(ng,1)*clat(ng))
avhh=f1112*(hh(1,1)*clat(1)+hh(ng,1)*clat(ng))
do j=2,ngm1
  avdd=avdd+dd(j,1)*clat(j)
  avhh=avhh+hh(j,1)*clat(j)
enddo
avdd=avdd*rsumi
avhh=avhh*rsumi

do j=1,ng
  dd(j,1)=dd(j,1)-avdd
  hh(j,1)=hh(j,1)-avhh
enddo

 !Update height and divergence in physical space:
do m=1,nt
  do j=1,ng
    hhp(j,m)=hh(j,m)
    ddp(j,m)=dd(j,m)
  enddo
enddo
call revfft(ng,nt,hhp,trig,factors) 
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
! residual PV and divergence.

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision,parameter:: cflmax=0.7d0, dtzz=pi/10.d0
double precision,parameter:: fdtmin=0.9d0, fdtmax=one/fdtmin
! *** The time step is adapted only if dtnew < fdtmin*dt or 
!                                      dtnew > fdtmax*dt
logical change

!---------------------------------------------------------------------
 !Get rms divergence, vorticity & enstrophy:
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
    hht=one+hhp(j,i)
    hhmin=min(hhmin,hht)
    hhmax=max(hhmax,hht)
    zzmax=max(zzmax,abs(zz(j,i)))
    frmax=max(frmax,usq/hht)
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
if (.not. ramping) twist=twist+dt*zzmax

 !Record cfl parameter, enstrophy, |Ro|_max, etc in cfl.dat:
cfl=umax*dt/dl
write(18,'(f10.4,1x,f6.4,6(1x,f8.5))') & 
    t,cfl,enstro,zzmax/fpole,hhmin,hhmax,frmax,ddrms/zzrms

if (.not. change) return
!----------------------------------------------------------------
 !Change the time step & define arrays needed for semi-implicit 
 !time-stepping:
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
    alpha(j)=d2lon(j,m)+(fsq(j)+omv(j,m)**2)*csqi
  enddo

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(1)
  d0=-pm*tlat(1)

  deti=one/(a0*d0-b0*c0)
  ahe0(1,m)=a0*deti
  bhe0(1,m)=b0*deti
  che0(1,m)=c0*deti
  dhe0(1,m)=d0*deti

  cphe=pc12-0.1d0*alpha(2)
  dphe=-0.1d0*tlat(2)
  ehe1(1,m)= cphe*bhe0(1,m)- apla*dhe0(1,m)
  ehe2(1,m)= dphe*bhe0(1,m)-0.2d0*dhe0(1,m)
  fhe1(1,m)= apla*che0(1,m)- cphe*ahe0(1,m)
  fhe2(1,m)=0.2d0*che0(1,m)- dphe*ahe0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-alpha(j)
    d0=-tlat(j)
    cmhe(j,m)=pc12-0.1d0*alpha(j-1)
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

    cphe=pc12-0.1d0*alpha(j+1)
    dphe=-0.1d0*tlat(j+1)
    ehe1(j,m)= cphe*bhe0(j,m)- apla*dhe0(j,m)
    ehe2(j,m)= dphe*bhe0(j,m)-0.2d0*dhe0(j,m)
    fhe1(j,m)= apla*che0(j,m)- cphe*ahe0(j,m)
    fhe2(j,m)=0.2d0*che0(j,m)- dphe*ahe0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(ng)
  d0=-pm*tlat(ng)

  cmhe(ng,m)=pc12-0.1d0*alpha(ngm1)
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

 !For latitudinal hyperviscous damping of h & delta:
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
! vorticity and height anomalies remain zero.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Work array:
double precision:: wka(ng,nt)

!-------------------------------------------------------------
if (ramping) then
   !Apply ramp function to fixed field qr to define PV anomaly qq:
  ramp=f12*(one-cos(pi*t/tramp))
  do i=1,nt
    do j=1,ng
      qq(j,i)=ramp*qr(j,i)
    enddo
  enddo

else
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
  do m=1,nt
    do j=1,ng
      wka(j,m)=qs(j,m)-qq(j,m)
    enddo
  enddo

   !Apply low-pass filter F to qs-qc:
  call lofilter(wka)
   !wka = F(qs-qc) (now semi-spectral)

   !Define PV anomaly from qq = F(qs-qc)+qd+qc:
  do m=1,nt
    do j=1,ng
      qq(j,m)=wka(j,m)+qd(j,m)+qq(j,m)
    enddo
  enddo
   !Return qq to physical space:
  call revfft(ng,nt,qq,trig,factors)
endif

return
end subroutine

!========================================================================

subroutine reset
! Resets qs & qd at the beginning of each time step, or
! just before recontouring at the end of a run.

! On return, qq contains the full gridded PV anomaly BUT
! with an incorrect mean value.  The mean value of qq is 
! found in getzet by the requirement that the global mean
! vorticity and height anomalies remain zero.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Work array:
double precision:: wka(ng,nt)

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
do m=1,nt
  do j=1,ng
    wka(j,m)=qs(j,m)-qq(j,m)
  enddo
enddo

 !Apply low-pass filter F to qs-qc:
call lofilter(wka)
 !wka = F(qs-qc) (now semi-spectral)

 !Add on qd and obtain reset value of qs <-- F(qs-qc)+qd+qc:
do m=1,nt
  do j=1,ng
    qd(j,m)=wka(j,m)+qd(j,m)
    qs(j,m)= qd(j,m)+qq(j,m)
    qq(j,m)=         qs(j,m)
  enddo
enddo
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

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Work arrays:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt)

double precision:: powhh(ngp1),powdd(ngp1),powqs(ngp1),powqd(ngp1)
real:: hhr4(ng,nt),uur4(ng,nt),vvr4(ng,nt)
real:: zzr4(ng,nt),qqr4(ng,nt),ddr4(ng,nt),tr4

!------------------------------------------------------------
 !Increment counter for dumping data:
idump=idump+1

 !Define single precision values of various fields for output:
do i=1,nt
  do j=1,ng
    hhr4(j,i)=real(hhp(j,i))
    uur4(j,i)=real(uu(j,i))
    vvr4(j,i)=real(vv(j,i))
    zzr4(j,i)=real(zz(j,i))
    qqr4(j,i)=real(qq(j,i)+cof(j))
    ddr4(j,i)=real(ddp(j,i))
  enddo
enddo
tr4=real(t)

 !Write fields for later imaging or diagnostics:
write(41,rec=idump) tr4,hhr4
write(42,rec=idump) tr4,uur4
write(43,rec=idump) tr4,vvr4
write(44,rec=idump) tr4,zzr4
write(45,rec=idump) tr4,qqr4
write(46,rec=idump) tr4,ddr4

!------------------------------------------------------------------
 !Compute and write energy and angular momentum per unit area:
do i=1,nt
  do j=1,ng
    r1=clat(j)
    r2=r1*r1
    hht=one+hhp(j,i)
    wka(j,i)=r1*hht*(uu(j,i)**2+vv(j,i)**2)
    wkb(j,i)=r1*csq*hhp(j,i)**2
    wkc(j,i)=r2*(hht*uu(j,i)+omega*r1*hhp(j,i))
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

!------------------------------------------------------------------
 !Compute and write longitudinal power spectra for hh, dd, qs & qd:
call longspec(hh,powhh)
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
! This routine writes the current contours and the residual 
! PV to the cont subdirectory and h, u, v etc to grid subdirectory

implicit double precision(a-h,o-z)
implicit integer(i-n)

real:: hhr4(ng,nt),uur4(ng,nt),vvr4(ng,nt)
real:: zzr4(ng,nt),qqr4(ng,nt),ddr4(ng,nt),tr4
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

 !Write h, u, v, zeta, q & delta to grid subdirectory
 !Define single precision values of various fields for output:
do i=1,nt
  do j=1,ng
    hhr4(j,i)=real(hhp(j,i))
    uur4(j,i)=real(uu(j,i))
    vvr4(j,i)=real(vv(j,i))
    zzr4(j,i)=real(zz(j,i))
    qqr4(j,i)=real(qq(j,i)+cof(j))
    ddr4(j,i)=real(ddp(j,i))
  enddo
enddo

 !Write fields for later imaging or diagnostics:
write(31,rec=irec) tr4,hhr4
write(32,rec=irec) tr4,uur4
write(33,rec=irec) tr4,vvr4
write(34,rec=irec) tr4,zzr4
write(35,rec=irec) tr4,qqr4
write(36,rec=irec) tr4,ddr4

 !Write final data for imaging with mgrid, llgrid etc:
if (loop .eq. nperiod) call dump

return
end subroutine

!=======================================================================

 !Main end module
end module
