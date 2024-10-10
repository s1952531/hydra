!#########################################################################
!                   The Doubly-Periodic Single-Layer 
!                 Green-Naghdi Pseudo-Spectral Method
!#########################################################################

!       Code adapted from ~/hydra/ps/plane/sw/swps codes in April 2019
!       by D G Dritschel @ St Andrews.

!       This code simulates the unforced Green-Naghdi Equations (GNE) 
!       in variables (ql,delta,gamma), where ql is the linearised 
!       potential vorticity anomaly, delta is the velocity divergence, 
!       and gamma is the acceleration divergence (called ageostrophic 
!       vorticity).  Note: ql = zeta - f*h, where zeta is the relative
!       vorticity, f is the Coriolis frequency and h is the dimensionless
!       depth anomaly.

!       The full algorithm consists of the following modules:

!       gnps.f90      : This source - main program to evolve fields;
!       parameters.f90: User defined parameters for a simulation;
!       constants.f90 : Fixed constants used throughout the other modules;
!       spectral.f90  : Fourier transform common storage and routines;

!----------------------------------------------------------------------------
program gnps

 !Import contants, parameters and common arrays:
use constants
use spectral

implicit none

 !Define common space:

 !Velocity field, dimensionless height anomaly & vorticity (physical):
double precision:: uu(ng,ng),vv(ng,ng),hh(ng,ng),zz(ng,ng)

 !Prognostic fields (spectral):
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng)

 !Time:
double precision:: t

 !Number of time steps between field saves and related indices:
integer:: ngsave,itime,jtime

 !Logical for use in calling inversion routine:
logical:: ggen

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Save data periodically:
  itime=nint(t/dt)
  jtime=itime/ngsave
  if (ngsave*jtime .eq. itime) then
     !Invert PV, divergence and acceleration divergence to obtain the
     !dimensionless height anomaly and velocity field, as well as the
     !relative vorticity (see spectral.f90):
    call main_invert(qs,ds,gs,hh,uu,vv,zz)
     !Note: qt, ds & gs are in spectral space while 
     !      hh, uu, vv and zz are in physical space.
    call savegrid(jtime+1)
    ggen=.false.
  else
    ggen=.true.
  endif
   !ggen is used to indicate if calling inversion is needed in advance below

   !Advect flow from time t to t + dt:
  call advance(ggen)
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Possibly save final data:
itime=nint(t/dt)
jtime=itime/ngsave
if (ngsave*jtime .eq. itime) then
  call main_invert(qs,ds,gs,hh,uu,vv,zz)
  call savegrid(jtime+1)
endif

 !Close all files:
call finalise


 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!----------------------------------------------------------------------
 !Read linearised PV anomaly and convert to spectral space as qs:
open(11,file='ql_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,qs,xfactors,yfactors,xtrig,ytrig)
 !Ensure domain average qs is zero (this does not matter):
qs(1,1)=zero

!----------------------------------------------------------------------
 !Read divergence and convert to spectral space as ds:
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,ds,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
ds(1,1)=zero

!----------------------------------------------------------------------
 !Read acceleration divergence and convert to spectral space as gs:
open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)

call ptospc(ng,ng,zz,gs,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
gs(1,1)=zero

!----------------------------------------------------------------------
 !Spectrally-truncate all fields for use in de-aliasing:
qs=filt*qs
ds=filt*ds
gs=filt*gs

 !Obtain initial dimensionless height anomaly (hh), velocity (uu,vv),
 !and relative vorticity (zz):
call main_invert(qs,ds,gs,hh,uu,vv,zz)
 !Note: qs, ds & gs are in spectral space while 
 !      hh, uu, vv and zz are in physical space.

!--------------------------------------
 !Open all plain text diagnostic files:
open(15,file='ecomp.asc',status='replace')
open(16,file='ubar.asc',status='replace')
open(17,file='monitor.asc',status='replace')

 !Open file for 1d vorticity & divergence spectra:
open(51,file='spectra.asc',status='replace')

 !Open files for coarse grid saves:
open(31,file='qq.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(32,file='dd.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(33,file='gg.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(34,file='hh.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(35,file='zz.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)

 !Define number of time steps between grid and contour saves:
ngsave=nint(tgsave/dt)
 !*** WARNING: tgsave should be an integer multiple of dt

return
end subroutine

!=======================================================================

subroutine advance(ggen)

! Advances fields from time t to t+dt using an iterative implicit 
! method of the form
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

 !Spectral fields needed in time stepping:
double precision:: qsi(ng,ng),qsm(ng,ng),sqs(ng,ng)
double precision:: dsi(ng,ng),sds(ng,ng),nds(ng,ng)
double precision:: gsi(ng,ng),sgs(ng,ng),ngs(ng,ng)

 !Other local quantities:
integer:: iter

!------------------------------------------------------------------
 !Invert PV and compute velocity at current time level, say t=t^n:
if (ggen) call main_invert(qs,ds,gs,hh,uu,vv,zz)
 !If ggen is false, main_invert was called previously at this time level.

 !Save various diagnostics each time step:
call diagnose

!------------------------------------------------------------------
 !Start with a guess for F^{n+1} for all fields:

 !Calculate the source terms (sqs,sds,sgs) for linearised PV (qs), 
 !divergence (ds) and acceleration divergence (gs):
call source(sqs,sds,sgs)

 !Update PV field:
qsi=qs
qsm=qs+dt4*sqs
qs=diss*(qsm+dt4*sqs)-qsi

 !Update divergence and acceleration divergence:
dsi=ds
gsi=gs
nds=sds+dt4i*dsi
ngs=sgs+pope*gsi
sds=nds+sds          !2*N_tilde_delta
sgs=ngs+sgs          !2*N_tilde_gamma
ds=simp*(pdis*sds+sgs)-dsi       !simp = 1/(P*R^2-G);  pdis = P*R
gs=simp*(rdis*sgs+opak*sds)-gsi  !rdis = R;  opak = G

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call main_invert(qs,ds,gs,hh,uu,vv,zz)

   !Calculate the source terms (sqs,sds,sgs) for linearised PV (qs), 
   !divergence (ds) and acceleration divergence (gs):
  call source(sqs,sds,sgs)

   !Update PV field:
  qs=diss*(qsm+dt4*sqs)-qsi

   !Update divergence and acceleration divergence:
  sds=nds+sds          !2*N_tilde_delta
  sgs=ngs+sgs          !2*N_tilde_gamma
  ds=simp*(pdis*sds+sgs)-dsi       !simp = 1/(P*R^2-G);  pdis = P*R
  gs=simp*(rdis*sgs+opak*sds)-gsi  !rdis = R;  opak = G
enddo

 !Advance time:
t=t+dt

return
end subroutine

!=======================================================================

subroutine source(sqs,sds,sgs)

! Gets the nonlinear source terms for linearised PV, divergence and 
! acceleration divergence  --- all in spectral space.  These are 
! returned in sqs, sds and sgs respectively.

! Note that (sds,sgs) only include the nonlinear terms for a 
! semi-implicit treatment, closely analogous to that described in 
! the appendix of Mohebalhojeh & Dritschel (2004).

! The spectral fields qs, ds and gs are all spectrally truncated.
! Note, hh, uu, vv & zz obtained by main_invert before calling this 
! routine are all spectrally truncated.

implicit none

 !Passed variables:
double precision:: sqs(ng,ng),sds(ng,ng),sgs(ng,ng)

 !Local variables (physical):
double precision:: htot(ng,ng),hx(ng,ng),hy(ng,ng),ht(ng,ng)
double precision:: dd(ng,ng),ddt(ng,ng),aa(ng,ng),bb(ng,ng)
double precision:: ut(ng,ng),vt(ng,ng),ztn(ng,ng)
double precision:: wkg(ng,ng),wkh(ng,ng),wkm(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng),wkz(ng,ng)
double precision:: errgt

 !Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng)

double precision,parameter:: tolgt=1.d-10*cof**3
 !tolgt: maximum error in iteration below to find gamma_t; the f^3 factor
 !       comes from the units of gamma (1/T^2) and d/dt (1/T).
 !       Note: gamma_t is returned in sgs (spectral).

!---------------------------------------------------------------
 !qs source --- compute div(ql*u,ql*v) (wka in spectral space):
wka=qs
call spctop(ng,ng,wka,wkq,xfactors,yfactors,xtrig,ytrig)
 !wkq contains the linearised PV in physical space
wkp=wkq*uu
wkq=wkq*vv
 !Compute spectral divergence from physical fields:
call divs(wkp,wkq,wke)
 !Do not overwrite wke until sqs is determined below.

 !Obtain a copy of the divergence in physical space (dd):
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)

 !Define h_tot = 1 + hh:
htot=one+hh

 !Compute J(u,v):
call jacob(uu,vv,wkp)

 !Form gamma_tilde = gamma + 2J(u,v) - 2*delta^2:
wkq=two*(wkp-dd**2)
 !Convert to spectral space, add on gs to get gamma_tilde (spectral):
call ptospc(ng,ng,wkq,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka+gs
 !Return to physical space:
call spctop(ng,ng,wka,wkq,xfactors,yfactors,xtrig,ytrig)
 !Compute (H^2/3)*(1+h)*J(h,gamma_tilde):
call jacob(hh,wkq,wkz)
 !Filter J before multiplying it by (1+h):
call ptospc(ng,ng,wkz,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka
call spctop(ng,ng,wka,wkz,xfactors,yfactors,xtrig,ytrig)
 !Finish q_l source:
wkz=hbsq3*htot*wkz
call ptospc(ng,ng,wkz,wka,xfactors,yfactors,xtrig,ytrig)
sqs=filt*(wka-wke)

!---------------------------------------------------------------
 !Nonlinear part of ds source --- 2J(u,v) - div(delta*(u,v)):
call jacob(uu,vv,wkp)
 !Convert J(u,v) (in wkp) to spectral space as wkc:
call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)

 !Compute div(delta*u,delta*v) (to be put into wka, in spectral space):
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)
 !dd contains the divergence in physical space
hx=dd*uu
hy=dd*vv
 !Compute spectral divergence from physical fields:
call divs(hx,hy,wka)

 !Add everything up to define delta_t - gamma (spectral, filtered)
sds=filt*(two*wkc-wka)
 !Add gamma and convert to physical space as ddt for use below:
wkd=gs+sds
call spctop(ng,ng,wkd,ddt,xfactors,yfactors,xtrig,ytrig)

!---------------------------------------------------------------
 !gs source --- find iteratively from implicit equation:

 !First compute all fixed quantities in the iteration:
 !Define aa = gamma + 2J(u,v) - 2*delta^2 in physical space:
wka=gs+two*wkc
call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)
aa=aa-two*dd**2
 !Spectrally truncate aa:
call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka
call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)

 !Obtain derivatives of h:
wkp=hh
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call xderiv(ng,ng,hrkx,wka,wkd)
call spctop(ng,ng,wkd,hx,xfactors,yfactors,xtrig,ytrig)
 !hx is h_x in physical space (partial derivative wrt x)
call yderiv(ng,ng,hrky,wka,wkd)
call spctop(ng,ng,wkd,hy,xfactors,yfactors,xtrig,ytrig)
 !hy is h_y in physical space (partial derivative wrt y)

 !Obtain h_t = N_h - delta:
wkp=hh*uu
wkq=hh*vv
 !Compute div(h*u,h*v) = -N_h spectrally and filter for use below:
call divs(wkp,wkq,wkd)
wkd=filt*wkd
 !For use below, calculate -c^2*Lap(N_h):
wke=c2g2*wkd
 !Above, c2g2 = c^2*Lap in spectral space;
 !Convert spectral -N_h in wkd to physical space:
call spctop(ng,ng,wkd,ht,xfactors,yfactors,xtrig,ytrig)
 !Flip sign and subtract delta to finish calculation of h_t (filtered):
ht=-ht-dd

 !Obtain zeta_t -f*delta (nonlinear part, N_zeta, only):
wkz=zz+cof
wkp=zz*uu
wkq=zz*vv
 !Compute div(zeta*u,zeta*v) spectrally:
call divs(wkp,wkq,wkd)
 !Spectrally truncate:
wkd=filt*wkd
 !Convert to physical space:
call spctop(ng,ng,wkd,ztn,xfactors,yfactors,xtrig,ytrig)
 !Compute J(h,aa):
call jacob(hh,aa,wkp)
 !Add everything up to define N_zeta (** unfiltered **):
ztn=hbsq3*wkp*htot-ztn

 !Obtain u_t & v_t (local rate of change):
wkp=f12*(uu**2+vv**2)
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
wka=filt*wka
call xderiv(ng,ng,hrkx,wka,wkd)
call spctop(ng,ng,wkd,wkp,xfactors,yfactors,xtrig,ytrig)
 !wkp is K_x in physical space where K = (u^2+v^2)/2
call yderiv(ng,ng,hrky,wka,wkd)
call spctop(ng,ng,wkd,wkq,xfactors,yfactors,xtrig,ytrig)
 !wkq is K_y in physical space where K = (u^2+v^2)/2
 !Compute (1+h)^3*aa = (1+h)^2 * (1+h)*a (ensure de-aliasing):
wkh=htot**2
call ptospc(ng,ng,wkh,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka
call spctop(ng,ng,wka,wkh,xfactors,yfactors,xtrig,ytrig)
wkm=htot*aa
call ptospc(ng,ng,wkm,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka
call spctop(ng,ng,wka,wkm,xfactors,yfactors,xtrig,ytrig)
wkh=wkh*wkm
call ptospc(ng,ng,wkh,wka,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
wka=filt*wka
 !Compute derivatives of (1+h)^3*aa (in the array wka):
call xderiv(ng,ng,hrkx,wka,wkd)
call spctop(ng,ng,wkd,wkg,xfactors,yfactors,xtrig,ytrig)
 !wkg is ((1+h)^3*aa)_x in physical space
call yderiv(ng,ng,hrky,wka,wkd)
call spctop(ng,ng,wkd,wkh,xfactors,yfactors,xtrig,ytrig)
 !wkh is ((1+h)^3*aa)_y in physical space
 !Add all terms up in physical space to obtain u_t & v_t:
wkm=hbsq3/htot
 !Spectrally truncate wkm:
call ptospc(ng,ng,wkm,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka
call spctop(ng,ng,wka,wkm,xfactors,yfactors,xtrig,ytrig)
ut=wkm*wkg-wkp+wkz*vv-csq*hx
vt=wkm*wkh-wkq-wkz*uu-csq*hy

 !Obtain bb = 2(1+h)*(J(u_t,v)+J(u,v_t)-2*delta*delta_t):
call jacob(ut,vv,wkp)
call jacob(uu,vt,wkq)
bb=two*htot*(wkp+wkq-two*dd*ddt)
 !Spectrally truncate bb:
call ptospc(ng,ng,bb,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call spctop(ng,ng,wkb,bb,xfactors,yfactors,xtrig,ytrig)

 !Obtain f*N_zeta in spectral space:
call ptospc(ng,ng,ztn,wkb,xfactors,yfactors,xtrig,ytrig)

 !Obtain f*N_zeta - c^2*Lap{N_h} + (H^2/3)*Lap{(1+h)*(2*h_t*aa+bb)}:
wkp=htot*(two*ht*aa+bb)
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
 !Apply spectral Laplacian, -(k^2+l^2), to wka and add to 
 !wke = -c^2*Lap(N_h) computed above:
wka=cof*wkb+wke-hbsq3*rksq*wka

 !Obtain (H^2/3)*div(px,py) where px = aa*((1+h)*h_t)_x+bb*h_x
 !                            and py = aa*((1+h)*h_t)_y+bb*h_y
wkp=htot*ht
call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
wkc=filt*wkc
call xderiv(ng,ng,hrkx,wkc,wkd)
call spctop(ng,ng,wkd,wkg,xfactors,yfactors,xtrig,ytrig)
 !wkg is ((1+h)*h_t)_x in physical space
call yderiv(ng,ng,hrky,wkc,wkd)
call spctop(ng,ng,wkd,wkh,xfactors,yfactors,xtrig,ytrig)
 !wkh is ((1+h)*h_t)_y in physical space
wkp=aa*wkg+bb*hx
 !This is px in physical space
wkq=aa*wkh+bb*hy
 !This is py in physical space
call divs(wkp,wkq,wkb)
 !Multiply spectral div(px,py) by H^2/3 and add to 
 !f*N_zeta - c^2*Lap{N_h} + (H^2/3)*Lap{(1+h)*(2*h_t*aa+bb)}
 !in wka above to complete *nonlinear part* of sbar_gamma:
wkd=filt*(wka+hbsq3*wkb)
 !Add linear part to define Sbar_gamma for use in iteration below:
wka=wkd+opak*ds
 !Note: opak*ds = (c^2*Lap - f^2)[delta]

 !Filter and store other fixed arrays for use in iteration below:
wkz=htot**2-one
call ptospc(ng,ng,wkz,wke,xfactors,yfactors,xtrig,ytrig)
wke=filt*wke
call spctop(ng,ng,wke,wkz,xfactors,yfactors,xtrig,ytrig)
aa=htot*hx
call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
wke=filt*wke
call spctop(ng,ng,wke,aa,xfactors,yfactors,xtrig,ytrig)
bb=htot*hy
call ptospc(ng,ng,bb,wke,xfactors,yfactors,xtrig,ytrig)
wke=filt*wke
call spctop(ng,ng,wke,bb,xfactors,yfactors,xtrig,ytrig)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Begin iteration to find gamma_t:

sgs=wka !Do not use wka below in iteration!
 !Convert first guess for gamma_t (in sgs) to physical space:
call spctop(ng,ng,sgs,wkg,xfactors,yfactors,xtrig,ytrig)
 !Note: sgs is overwritten here; only wkg (physical gamma_t) is used below:

errgt=two*tolgt
do while (errgt .gt. tolgt)
  wkp=wkz*wkg !wkz = (1+h)^2-1 = h^2+2*h from above
  call ptospc(ng,ng,wkp,wke,xfactors,yfactors,xtrig,ytrig)
   !wke contains spectral (1+h)^2*gamma_t
  wkp=aa*wkg !aa = (1+h)*h_x from above
  wkq=bb*wkg !bb = (1+h)*h_y from above
   !Compute div((1+h)*h_x*gamma_t,(1+h)*h_y*gamma_t) spectrally:
  call divs(wkp,wkq,wkb)
   !Add up and apply operator adop to get corrected gamma_t in spectral space:
  wkc=hbsq3*(wkb-rksq*wke)
  sgs=adop*(wka+wkc)
   !Convert to physical space to calculate error:
  wkb=sgs
  call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
   !Compute rms error in gamma_t:
  errgt=sqrt(sum((wkp-wkg)**2))
   !Store physical gamma_t in wkg for use in next iteration:
  wkg=wkp
enddo

 !End of iteration; gamma_t (spectral) is now in sgs.
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Redefine sgs as the (spectrally-filtered) nonlinear part needed in
 !the semi-implicit time stepping:
sgs=filt*(wkd+wkc)

return
end subroutine

!=======================================================================

subroutine diagnose

! Computes various quantities every time step to monitor the flow evolution.

implicit none

 !Local variables:
double precision:: umax,zrms,zmax
double precision:: uio,vio

!----------------------------------------------------------------------
 !Compute diagnostics:
umax=sqrt(maxval(uu**2+vv**2))
zmax=maxval(abs(zz))
zrms=sqrt(dsumi*sum(zz**2))

 !Record various diagnostics to monitor.asc:
write(17,'(1x,f12.5,4(1x,f12.6))') t,f12*zrms**2,zrms,zmax,umax

 !Compute and write mean flow:
uio=sum(uu)*dsumi
vio=sum(vv)*dsumi
write(16,'(1x,f12.5,2(1x,f14.10))') t,uio,vio

return
end subroutine

!=======================================================================

subroutine savegrid(igrids)

! Saves PV, energy and various spectra at the desired save time

implicit none

 !Passed variable:
integer:: igrids

 !Local variables:
double precision:: dd(ng,ng),htot(ng,ng),wkp(ng,ng) !Physical
double precision:: wks(ng,ng)                       !Spectral
double precision:: zspec(0:ng),dspec(0:ng),gspec(0:ng)
double precision:: ekin,epot,ediv,etot
integer:: k

!---------------------------------------------------------------
 !Compute energy components and total:
wks=ds
call spctop(ng,ng,wks,dd,xfactors,yfactors,xtrig,ytrig)
htot=one+hh
ekin=f12*garea*sum(htot*(uu**2+vv**2))
epot=f12*garea*csq*sum(hh**2)
dd=htot*dd
ediv=f12*garea*hbsq3*sum(htot*dd**2)
etot=ekin+epot+ediv

 !Write energies to ecomp.asc:
write(15,'(f9.2,5(1x,f16.9))') t,ediv,ekin,ekin+ediv,epot,etot
write(*,'(a,f9.2,a,f13.6)') ' t = ',t,'  E_tot = ',etot

!---------------------------------------------------------------
 !Compute 1d vorticity, divergence & acceleration divergence spectra:
wkp=zz
call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
call spec1d(wks,zspec)
call spec1d( ds,dspec)
call spec1d( gs,gspec)
 !Normalise to take into account uneven sampling of wavenumbers 
 !in each shell [k-1/2,k+1/2]:
zspec=spmf*zspec
dspec=spmf*dspec
gspec=spmf*gspec

write(51,'(f13.6,1x,i5)') t,kmaxred
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(51,'(4(1x,f12.8))') alk(k),log10(zspec(k)),log10(dspec(k)+1.d-32), &
                                                   log10(gspec(k)+1.d-32)
enddo
 !Note: alk(k) = log_10(k)

!---------------------------------------------------------------
 !Write various gridded fields to direct access files:
 !PV field:
wkp=(zz+cof)/(one+hh)-cof
write(31,rec=igrids) real(t),real(wkp)
 !Divergence:
wks=ds
call spctop(ng,ng,wks,wkp,xfactors,yfactors,xtrig,ytrig)
write(32,rec=igrids) real(t),real(wkp)
 !Acceleration divergence:
wks=gs
call spctop(ng,ng,wks,wkp,xfactors,yfactors,xtrig,ytrig)
write(33,rec=igrids) real(t),real(wkp)
 !Dimensionless depth anomaly:
write(34,rec=igrids) real(t),real(hh)
 !Relative vorticity:
write(35,rec=igrids) real(t),real(zz)

return
end subroutine

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(15)
close(16)
close(17)
close(31)
close(32)
close(33)
close(34)
close(35)
close(51)

return
end subroutine

 !End main program
end program gnps

!=======================================================================
