module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

double precision:: qt(ng,nt),qc(ng,nt),qd(ng,nt) !Various PVs

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================
subroutine advect

! Main subroutine for advecting fields and contours

implicit none

 !Local variables:
double precision,parameter:: twistmax=2.5d0
!      twistmax: the maximum value of the time integral of |zeta|_max
!                between regularisations of the contours.
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, the code stops
!      to rebuild the contours in a separate memory space.
integer:: ireg,itime,jtime
logical:: ggen

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

 !Used for regularising contours:
twist=zero

 !Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Save data periodically:
  itime=nint(t/dt)
  jtime=itime/ngsave
  if (ngsave*jtime .eq. itime) then
    call inversion(1)
    call savegrid(jtime+1)
    ggen=.false.
  else
    ggen=.true.
  endif
   !ggen is used to indicate if calling inversion is needed in advance below
  jtime=itime/ncsave
  if (ncsave*jtime .eq. itime) call savecont(jtime+1)

   !Perform contour surgery or recontouring when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1

     !Don't continue if maximum number of regularisations reached:
    if (ireg .eq. nregmax) then
       !Prepare PV residual qr for recontouring (and preserve qs):
      call prepare
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the PV contours (surgery + node redistribution):
    call surgery
     !Record active contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i9,1x,i10)') t,n,npt

    twist=twist-twistmax
  endif

   !Advect flow from time t to t + dt:
  call advance(ggen)
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Possibly save final data:
itime=nint(t/dt)
jtime=itime/ngsave
if (ngsave*jtime .eq. itime) then
  call inversion(1)
  call savegrid(jtime+1)
endif
jtime=itime/ncsave
if (ncsave*jtime .eq. itime) call savecont(jtime+1)

return
end subroutine advect

!=======================================================================

subroutine init

! Initialises residual PV for normal time integration following 
! contour regeneration

implicit none
integer:: i

!------------------------------------------------------------------
 !Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,n,npt

!------------------------------------------------------------
 !Re-instate qs as the PV anomaly (here qs is on the grid):
do i=1,nt
  qs(:,i)=qs(:,i)-cof
enddo

 !Convert qs to semi-spectral space:
call forfft(ng,nt,qs,trig,factors) 

 !Convert PV contours (x,y,z) to gridded PV anomaly as qc:
call con2grid(qc)

 !Convert qc to semi-spectral space:
call forfft(ng,nt,qc,trig,factors) 

 !Define (semi-spectral) residual PV qd = (1-F)*(qs-qc):
qd=qs-qc
call hipass(qd)

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current PV anomaly field is stored in qs (semi-spectral), and
! the full and residual PV fields needed in congen.f90 is stored in qq
! and qr (physical).

implicit none

! Local variables:
double precision:: qa(ng,nt),aqa
integer:: i

!------------------------------------------------------------
 !Re-initialise qs & qd before recontouring:
 !          Reset  qs = F*(qs-qc) + qc + qd
 !            and  qd = (1-F)*(qs-qc)
 !Here F & (1-F) are low & high pass filter functions (see spectral.f90)
call reset

 !Obtain full PV anomaly (qs) on the grid:
call revfft(ng,nt,qs,trig,factors)

 !Add f to define full PV:
do i=1,nt
  qs(:,i)=qs(:,i)+cof
enddo

 !Correct average PV by enforcing zero average vorticity:
qa=qs*(one+hh)
aqa=average(qa)
qs=qs-aqa

 !Obtain PV anomaly due to contours (qc) on the grid:
call revfft(ng,nt,qc,trig,factors)

 !Define residual PV (qr) needed for recontouring:
do i=1,nt
  qr(:,i)=qs(:,i)-qc(:,i)-cof
enddo

 !Note: We are leaving this module after this; qd will be redefined
 !      upon re-entry in subroutine init.

return
end subroutine prepare

!=======================================================================

subroutine advance(ggen)

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method 
! (for all semi-spectral fields, namely qs, qd, ds & gs).

! Uses an iterative implicit method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}-F^n)/2] + N[(F^{n+1}-F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms, and N refers to the nonlinear source
! terms.  We start with a guess for F^{n+1} in N and iterate 
! niter times (see parameter statement below).

implicit none

 !Passed variable:
logical:: ggen

 !Local variables:
integer,parameter:: niter=2

 !Semi-spectral fields needed in time stepping:
double precision:: qsi(ng,nt),sqs(ng,nt)
double precision:: qdi(ng,nt),qdm(ng,nt),sqd(ng,nt)
double precision:: dsi(ng,nt),sds(ng,nt),nds(ng,nt)
double precision:: gsi(ng,nt),sgs(ng,nt),ngs(ng,nt)
 !Contour positions needed in time stepping:
double precision:: xi(npt),yi(npt),zi(npt)
 !Contour velocities:
double precision:: u(npt),v(npt),w(npt)
 !Other local quantities:
double precision:: xtmp,ytmp,ztmp,fac
integer:: i,iter

!-------------------------------------------------------------------
 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset  qs = F*(qs-qc) + qc + qd
 !            and  qd = (1-F)*(qs-qc)
 !Here F & (1-F) are low & high pass filter functions (see spectral.f90)
call reset

!------------------------------------------------------------------
 !Invert PV and compute velocity at current time level, say t=t^n:
if (ggen) call inversion(0)
 !If ggen is false, inversion was called previously at this time level.
 !A zero argument above means qc does not need to be regenerated.

 !Compute twist parameter and save various diagnostics each time step:
call diagnose

!------------------------------------------------------------------
 !Start with a guess for F^{n+1} for all contours and fields:

 !Contours:
call velint(uu,vv,u,v,w)
 !Here (u,v,w) stands for the Cartesian velocity on the contours
 !at time level n, i.e. u(x^n,t^n)
do i=1,npt
   !Store x^n+0.5*dt*u(x^n,t^n) for efficiency in the iteration below:
  xtmp=x(i)+dt2*u(i)
  ytmp=y(i)+dt2*v(i)
  ztmp=z(i)+dt2*w(i)
   !Ensure point remains on unit sphere:
  fac=one/sqrt(xtmp**2+ytmp**2+ztmp**2)
  xi(i)=fac*xtmp
  yi(i)=fac*ytmp
  zi(i)=fac*ztmp
   !Preliminary guess for x^{n+1}:
  xtmp=x(i)+dt*u(i)
  ytmp=y(i)+dt*v(i)
  ztmp=z(i)+dt*w(i)
   !Ensure point remains on unit sphere:
  fac=one/sqrt(xtmp**2+ytmp**2+ztmp**2)
  x(i)=fac*xtmp
  y(i)=fac*ytmp
  z(i)=fac*ztmp
enddo

 !Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
 !(sds,sgs) for divergence and linearised acceleration divergence (ds,gs):
call source(sqs,sqd,sds,sgs)

 !Update PV fields:
qdi=qd
qdm=two*qd+dt2*sqd
qd=diss*(qdm+dt2*sqd)-qdi
qsi=qs+dt2*sqs
qs=qs+dt*sqs

 !Update divergence and linearised acceleration divergence:
dsi=ds
gsi=gs
nds=sds+dt4i*dsi
ngs=sgs+dt4i*gsi
sds=nds+sds                !2*N_tilde_delta
sgs=ngs+sgs                !2*N_tilde_gamma
call simp(sds,sgs,ds,gs)   !Solves the semi-implicit system
ds=diss*ds-dsi             !ds = delta^{n+1}
gs=gs-gsi                  !gs = gamma_l^{n+1}

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call inversion(1)

   !Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
   !(sds,sgs) for divergence and linearised acceleration divergence (ds,gs):
  call source(sqs,sqd,sds,sgs)

   !Interpolate gridded velocity (uu,vv) at contour nodes as (u,v):
  call velint(uu,vv,u,v,w)
  do i=1,npt
    xtmp=xi(i)+dt2*u(i)
    ytmp=yi(i)+dt2*v(i)
    ztmp=zi(i)+dt2*w(i)
    fac=one/sqrt(xtmp**2+ytmp**2+ztmp**2)
    x(i)=fac*xtmp
    y(i)=fac*ytmp
    z(i)=fac*ztmp
  enddo
   !Now (x,y,z) contain a new guess for x^{n+1}.

   !Update PV fields:
  qd=diss*(qdm+dt2*sqd)-qdi
  qs=qsi+dt2*sqs

   !Update divergence and linearised acceleration divergence:
  sds=nds+sds                !2*N_tilde_delta
  sgs=ngs+sgs                !2*N_tilde_gamma
  call simp(sds,sgs,ds,gs)   !Solves the semi-implicit system
  ds=diss*ds-dsi             !ds = delta^{n+1}
  gs=gs-gsi                  !gs = gamma_l^{n+1}
enddo

 !Apply latitudinal hyperviscous damping to qd, ds & gs:
call latdamp(qd)
call latdamp(ds)
call latdamp(gs)
 !The longitudinal part for qd is incorporated above in diss,
 !that for gs in simp, while no such damping is applied to ds.

 !Advance time:
t=t+dt

return
end subroutine advance

!=======================================================================

subroutine source(sqs,sqd,sds,sgs)

! Gets the source terms (sqs,sqd) for the PV (qs,qd) as well as those
! (sds,sgs) for divergence and linearised acceleration divergence (ds,gs)
! --- all in semi-spectral space.

! Note that  sds = d(ds)/dt - gs  while  sgs = d(gs)/dt - G[ds],
! where G = c^2*Lap - f^2 is the SW gravity-wave operator.

! This follows the semi-implicit treatment used for the SW
! equations - see appendix B of Dritschel & Jalali, JFM (2020).

! The semi-spectral fields ds, gs, qd and qs are all spectrally truncated.
! Note, hh, uu, vv & zz obtained by main_invert before calling this 
! routine are all spectrally truncated.

implicit none

double precision,parameter:: toler=1.d-11
 !toler: maximum error in iteration below to find the vertically-
 !       integrated non-hydrostatic pressure P_n (ppn below).
 !       Note P_n = \bar{P}_n/H in the notes.

 !Passed variables:
double precision:: sqs(ng,nt),sqd(ng,nt),sds(ng,nt),sgs(ng,nt)

 !Local variables:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt),wkd(ng,nt)
double precision:: wke(ng,nt),wkf(ng,nt),wkg(ng,nt),wkh(ng,nt)
double precision:: wkp(ng,nt),wkq(ng,nt),wku(ng,nt),wkv(ng,nt)
double precision:: dd(ng,nt),htot(ng,nt),hinv(ng,nt),hinvsq(ng,nt)
double precision:: dra(ng,nt),drb(ng,nt),dha(ng,nt),dhb(ng,nt)
double precision:: rhs(ng),avgval,errpn
integer:: i

!-----------------------------------------------------------------
 !Compute A = z*U + dV/d(lon) -> wka & B = z*V - dU/d(lon) -> wkb,
 !where U = u/r & V = v/r, while z = sin(lat) and r = cos(lat):
do i=1,nt
  wku(:,i)=clati*uu(:,i)
  wkv(:,i)=clati*vv(:,i)
  wka(:,i)=wkv(:,i)
  wkb(:,i)=wku(:,i)
enddo
 !*** Do not re-use wku ***

 !Convert wka & wkb to semi-spectral space:
call forfft(ng,nt,wka,trig,factors) 
call forfft(ng,nt,wkb,trig,factors) 

 !Compute longitudinal derivatives:
call deriv(ng,nt,rk,wka,wkc)
call deriv(ng,nt,rk,wkb,wkd)

 !Recover physical fields:
call revfft(ng,nt,wkc,trig,factors) 
call revfft(ng,nt,wkd,trig,factors) 

 !Complete definition of A & B, and define wkf = f*zeta & wkd = 2*f*beta*v:
do i=1,nt
  wka(:,i)=slat*wku(:,i)+wkc(:,i) ! slat = z = sin(lat)
  wkb(:,i)=slat*wkv(:,i)-wkd(:,i)
  wkf(:,i)=cof*zz(:,i)            ! cof = f
  wkd(:,i)=dfb*vv(:,i)            ! dfb = 2*f*beta
enddo
 !*** Do not re-use wka, wkb, wkf or wkd ***

!---------------------------------------------------------------
 !Get physical space velocity divergence -> dd:
dd=ds
call revfft(ng,nt,dd,trig,factors) 
 !*** Do not re-use dd ***

 !Compute part of divergence tendency not involving non-hydrostatic pressure:
call deriv(ng,nt,rk,ds,wkv)
call revfft(ng,nt,wkv,trig,factors) ! wkv = d(delta)/d(lon)
call latder(dd,wkc)                 ! wkc = d(delta)/d(lat)
sds=dd**2-wku*wkv-vv*wkc
 !Here sds = delta^2 - (u,v)*grad(delta) in physical space.
 !The rest of the divergence tendency is obtained below.

!---------------------------------------------------------------
 !Form gamma_tilde = gamma_l + 2*(J(u,v) - delta^2) - u^2 - v^2:
wke=uu**2+vv**2
 !*** do not re-use wke ***
wkg=two*(wka*(zz-wka)-wkb*(wkb+dd)-dd**2)-wke
 !*** wka & wkb are now safe to re-use below ***
call dealias(wkg) !wkg is now de-aliased in semi-spectral space
wkg=gs+wkg        !gs contains gamma_l in semi-spectral space
call revfft(ng,nt,wkg,trig,factors)
 !Now wkg contains gamma_tilde on the grid.

 !Multiply next by (1+h) and store in wkg (spectral) the fixed rhs
 !in the non-hydrostatic pressure iteration below:
htot=one+hh
wkg=htot*wkg
call dealias(wkg)
 !Now wkg contains (1+h)*gamma_tilde in semi-spectral space.

!---------------------------------------------------------------
 !Calculate the vertically-integrated non-hydrostatic pressure:

 !Prepare by calculating wkq = (3/H^2)*(1/(1+h)^2 - 1) needed
 !in the pressure iteration below:
hinv=one/htot
call dealias(hinv)
call deriv(ng,nt,rk,hinv,dra)
call revfft(ng,nt,hinv,trig,factors)  ! hinv = (1+h)^{-1} (de-aliased)
call revfft(ng,nt,dra,trig,factors)   ! dra = d(1+h)^{-1}/d(lon)
call latder(hinv,drb)                 ! drb = d(1+h)^{-1}/d(lat)

hinvsq=hinv**2
call dealias(hinvsq)
call revfft(ng,nt,hinvsq,trig,factors)
wkq=hbsq3i*(hinvsq-one)               ! wkq is fully de-aliased

 !Calculate also (1+h)^{-1}*grad(h) and store in dha & dhb:
do i=1,nt
  dha(:,i)=-clatsqi*dra(:,i)*htot(:,i)
  dhb(:,i)=        -drb(:,i)*htot(:,i)
enddo
 !An extra factor of 1/cos(lat) is included in the definition of dha
 !to avoid extra work in the iteration below (when grad(P_n) is needed):

call dealias(dha)
call revfft(ng,nt,dha,trig,factors)
call dealias(dhb)
call revfft(ng,nt,dhb,trig,factors)   ! dha & dhb are now fully de-aliased

 !Now iterate to find P_n (in ppn), starting from the guess
 !ppn = (grad^2 - 3/H^2)^{-1}[(1+h)*gamma_tilde]:
call pope(wkg,ppn,wkb)                ! pope inverts grad^2 - 3/H^2 on wkg
 !Here wkb = d(ppn)/d(lat)
call deriv(ng,nt,rk,ppn,wka)
 !Here wka = d(ppn)/d(lon)
call revfft(ng,nt,ppn,trig,factors)
call revfft(ng,nt,wka,trig,factors)
call revfft(ng,nt,wkb,trig,factors)

errpn=two*toler
do while (errpn .gt. toler)
  wkc=dha*wka+dhb*wkb+wkq*ppn         ! wkc = (1+h)^{-1}*grad(h)*grad(P_n)
  call dealias(wkc)                   !     + (3/H^2)*(1/(1+h)^2 - 1)*P_n
  wkc=wkc+wkg                         ! Now wkc = rhs of elliptic P_n equation
  call pope(wkc,wkp,wkb)              ! pope inverts grad^2 - 3/H^2 on wkc
   !Now wkp is the correction to ppn while wkb = d(ppn)/d(lat)
  call deriv(ng,nt,rk,wkp,wka)
   !Now wka = d(ppn)/d(lon)
  call revfft(ng,nt,wkp,trig,factors)
  call revfft(ng,nt,wka,trig,factors)
  call revfft(ng,nt,wkb,trig,factors)
   !Compute rms error in P_n:
  wkc=(wkp-ppn)**2
  errpn=average(wkc)
  ppn=wkp
enddo
 !*** wkq is now safe to re-use below ***

 !Store non-hydrostatic part of gs tendency, f*J(P_n,1/(1+h)), in sgs:
do i=1,nt
  sgs(:,i)=fri*(wka(:,i)*drb(:,i)-wkb(:,i)*dra(:,i))
enddo
 !fri = f/cos(lat) above.

 !Store (1+h)^{-1}*d(ppn)/d(lon) in wkh for use in sgs below:
wkh=hinv*wka
 !Convert to semi-spectral space and de-alias:
call dealias(wkh)
!*** wka & wkb are now safe to re-use below ***

!---------------------------------------------------------------
 !qd source --- only NL advection term is needed:
call deriv(ng,nt,rk,qd,wkp)
call revfft(ng,nt,wkp,trig,factors)   ! wkp = d(qd)/d(lon)
wka=qd
call revfft(ng,nt,wka,trig,factors) 
call latder(wka,wkq)                  ! wkq = d(qd)/d(lat)
sqd=-wku*wkp-vv*wkq
 !Convert to semi-spectral space and de-alias:
call dealias(sqd)

!---------------------------------------------------------------
 !qs source --- only NL advection term is needed:
call deriv(ng,nt,rk,qs,wkp)
call revfft(ng,nt,wkp,trig,factors)   ! wkp = d(qs)/d(lon)
wka=qs
call revfft(ng,nt,wka,trig,factors) 
call latder(wka,wkq)                  ! wkq = d(qs)/d(lat)
do i=1,nt
  sqs(:,i)=-wku(:,i)*wkp(:,i)-vv(:,i)*(wkq(:,i)+bet) !include beta*v
enddo
 !Convert to semi-spectral space and de-alias:
call dealias(sqs)

!---------------------------------------------------------------
 !Nonlinear part of ds source:

 !Compute 1/(1+h)^3 -> wkp and de-alias:
wkp=hinv*hinvsq
call dealias(wkp)
call revfft(ng,nt,wkp,trig,factors)

 !Form delta^2 - (u,v)*grad(delta) - (3/H^2)*P_n/(1+h)^3:
wkp=sds-hbsq3i*wkp*ppn
 !Here sds = delta^2 - (u,v)*grad(delta) from above.

 !De-alias wkp and subtract P^{-1}[gs] to complete ds source:
call dealias(wkp)

call pope(gs,wkq,wkb)        ! pope inverts grad^2 - 3/H^2 on gs
sds=wkp+hbsq3i*wkq           ! Note: grad^2 - 3/H^2 = -(3/H^2)*P
 !wkq = (grad^2 - 3/H^2)^{-1}[gs] is needed below for sgs.

!---------------------------------------------------------------
 !Nonlinear part of gs source:

 !Compute div((u,v)*Z) (immediately below wkf = Z = f*zeta):
wka=wkf
call forfft(ng,nt,wka,trig,factors)
call deriv(ng,nt,rk,wka,wkv)
call revfft(ng,nt,wkv,trig,factors) ! wkv = dZ/d(lon)
call latder(wkf,wkc)                ! wkc = dZ/d(lat)

 !Store div((u,v)*Z) + 2*f*beta*v - f*J(P_n,1/(1+h)) in wkg:
wkg=dd*wkf+wku*wkv+vv*wkc+wkd-sgs
 !Note  sgs = f*J(P_n,1/(1+h))  &  wkd = 2*f*beta*v  from above.
 !De-alias wkg and convert to semi-spectral space:
call dealias(wkg)
 !*** wkf can now be re-used ***

 !Define B = c^2*h + (u^2+v^2)/2:
wkb=csq*hh+f12*wke
 !*** wke can now be re-used ***
 !De-alias B and convert to semi-spectral space:
call dealias(wkb)
 !Subtract (grad^2-3/H^2)^{-1}[gamma_l] for the SI treatment:
wkb=wkb-wkq
 !Calculate dB/d(lon) -> wkf (keep this semi-spectral):
call deriv(ng,nt,rk,wkb,wkf)

 !Compute div((u,v)*h) -> wke:
wka=hh
call forfft(ng,nt,wka,trig,factors) 
call deriv(ng,nt,rk,wka,wkv)
call revfft(ng,nt,wkv,trig,factors) ! wkv = dh/d(lon)
call latder(hh,wkc)                 ! wkc = dh/d(lat)
wke=dd*hh+wku*wkv+vv*wkc
 !Note wku = u/r on rhs above

 !Compute Laplacian of div((u,v)*h) in wke after de-aliasing:
call dealias(wke)     ! Now wke is in semi-spectral space
call laplace(wke,wka) ! wka is returned in semi-spectral space

 !Complete calculation of nonlinear part of gamma_l tendency:
sgs=csq*wka+fpole*(wkf+wkh)-wkg
 !Note: each part has been separately de-aliased.

!-----------------------------------------------------
 !Remove global mean values of sds & sgs:
rhs=sds(:,1)*clat
avgval=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
sds(:,1)=sds(:,1)-avgval

rhs=sgs(:,1)*clat
avgval=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
sgs(:,1)=sgs(:,1)-avgval

return
end subroutine source

!=======================================================================

subroutine inversion(iopt)

! Finds the gridded dimensionless height anomaly (hh) and velocity
! field (uu,vv) from the PV contours and the (semi-spectral)
! divergence (ds) and linearised acceleration divergence (gs).

! if iopt = 1, qc (semi-spectral) is found from the PV contours.
  
implicit none
integer:: iopt

!------------------------------------------------------------
if (iopt .eq. 1) then
   !Convert PV contours (x,y,z) to gridded PV anomaly as qc:
  call con2grid(qc)

   !Convert qc to semi-spectral space:
  call forfft(ng,nt,qc,trig,factors) 
endif

 !Combine fields to update qt with full (semi-spectral) field,
 !qt = F*(qs-qc) + qc + qd, where F is a low pass filter:
qt=qs-qc
call lopass(qt)
qt=qt+qc+qd

 !Invert PV, divergence and linearised acceleration divergence to obtain
 !the dimensionless height anomaly and velocity field, as well as the
 !gridded PV anomaly and relative vorticity (see spectral.f90):
call main_invert(qt,ds,gs,hh,uu,vv,qq,zz)
 !Note: qt, ds & gs are in semi-spectral space while 
 !      hh, uu, vv, qq and zz are in physical space.

return
end subroutine inversion

!=======================================================================

subroutine reset
! Resets qs & qd at the beginning of each time step, or
! just before recontouring at the end of a run.

implicit none

 !Axillary PV array:
double precision:: qa(ng,nt)

!-------------------------------------------------------------------
 !Reset  qs = F*(qs-qc) + qc + qd  and  qd = (1-F)*(qs-qc)
 !where F & (1-F) are low & high pass filter functions
 !(see spectral.f90)

 !Convert PV contours (x,y,z) to gridded PV anomaly as qc:
call con2grid(qc)

 !Convert qc to semi-spectral space:
call forfft(ng,nt,qc,trig,factors) 
qa=qs-qc
call lopass(qa)
qs=qa+qc+qd
qd=qs-qc
call hipass(qd)

return
end subroutine

!=======================================================================

subroutine diagnose

! Computes the twist parameter, the time integral of |zeta|_max, and
! various quantities at every time step to monitor the flow evolution.

implicit none

 !Local variables:
double precision:: wkp(ng,nt)
double precision:: zmax,frmax,hmax,hmin
double precision:: zrms,drms

!----------------------------------------------------------------------
 !Increment the integral of |zeta|_max:
zmax=maxval(abs(zz))
twist=twist+dt*zmax

 !Compute diagnostics:
frmax=sqrt(csqi*maxval((uu**2+vv**2)/(one+hh)))
call getrms(zz,zrms)
wkp=ds
call revfft(ng,nt,wkp,trig,factors)
call getrms(wkp,drms)
hmax=maxval(hh)
hmin=minval(hh)

 !Record |zeta|_max/f_pole, Fr_max, h_min, h_max, zeta_rms, delta_rms
 !to monitor.asc:
write(17,'(1x,f12.5,6(1x,f14.8))') t,zmax/fpole,frmax,hmin,hmax,zrms,drms

return
end subroutine diagnose

!=======================================================================

subroutine savegrid(igrids)

! Saves various fields, energy components and angular momentum,
! as well as various longitudinal spectra at the desired save time

implicit none

 !Passed variable:
integer:: igrids

 !Local variables:
double precision:: htot(ng,nt),dd(ng,nt),wkp(ng,nt)
double precision:: zspec(ngp1),dspec(ngp1),gspec(ngp1)
double precision:: ekin,epot,ediv,etot,angm
integer:: i,m

!-----------------------------------------------------------------
 !Compute energy components and total as well as angular momentum:
htot=one+hh
wkp=htot*(uu**2+vv**2)
ekin=twopi*average(wkp)
wkp=hh**2
epot=twopi*csq*average(wkp)
 !Get physical space velocity divergence -> dd:
dd=ds
call revfft(ng,nt,dd,trig,factors) 
wkp=htot*(htot*dd)**2
ediv=twopi*hbsq3*average(wkp)
etot=ekin+epot+ediv
do i=1,nt
  wkp(:,i)=clat*(htot(:,i)*uu(:,i)+omega*clat*hh(:,i))
enddo
angm=fourpi*average(wkp)

 !Write energies & angular momentum to ecomp.asc:
write(15,'(f13.6,6(1x,f16.9))') t,ediv,ekin,ekin+ediv,epot,etot,angm
write(*,'(a,f13.6,a,f14.7,a,f14.7)') ' t = ',t,'     E_tot = ',etot, &
                                     '     A = ',angm

!-----------------------------------------------------------------
 !Compute and write longitudinal power spectra for vorticity,
 !divergence & linearised acceleration divergence:
wkp=zz
call forfft(ng,nt,wkp,trig,factors)
call longspec(wkp,zspec)
call longspec(ds, dspec)
call longspec(gs, gspec)

write(51,'(f13.6,1x,i5)') t,ngp1
do m=1,ngp1
  write(51,'(4(1x,f12.8))') log10(dble(m)),log10(zspec(m)+1.d-32), &
                    log10(dspec(m)+1.d-32),log10(gspec(m)+1.d-32)
enddo

!-----------------------------------------------------------------
 !Write various gridded fields to direct access files:
write(31,rec=igrids) real(t),real(qq)
wkp=ds
call revfft(ng,nt,wkp,trig,factors)
write(32,rec=igrids) real(t),real(wkp)
wkp=gs
call revfft(ng,nt,wkp,trig,factors)
write(33,rec=igrids) real(t),real(wkp)
write(34,rec=igrids) real(t),real(hh)
write(35,rec=igrids) real(t),real(zz)
write(36,rec=igrids) real(t),real(ppn)

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont(irec)

! Saves PV contours for post-processing and imaging

implicit none

 !Passed variable:
integer:: irec

 !Local variables:
double precision:: qa(ng,nt)
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') irec

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') n,npt,t,qjump

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
qa=qt-qc
call revfft(ng,nt,qa,trig,factors)
write(83,rec=irec) real(t),real(qa)

 !Save PV contours:
open(81,file='cont/qqindex'//pind,form='unformatted',status='replace')
write(81) np(1:n),i1(1:n),ind(1:n)
close(81)

open(82,file='cont/qqnodes'//pind,form='unformatted',status='replace')
write(82) x(1:npt),y(1:npt),z(1:npt)
close(82)

return
end subroutine savecont

!=======================================================================

 !Main end module
end module evolution
