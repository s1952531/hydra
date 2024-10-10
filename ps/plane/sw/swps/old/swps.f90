!#########################################################################
!                   The Doubly-Periodic Single-Layer 
!                 Shallow-Water Pseudo-Spectral Method
!#########################################################################

!       Code adapted from ~/hydra/ca/plane/sw/caps codes in March 2019
!       by D G Dritschel @ St Andrews.

!       This code simulates the unforced Shallow-Water Equations (SWE) 
!       in variables (ql,delta,gamma), where ql is the linearised 
!       potential vorticity anomaly, delta is the velocity divergence, 
!       and gamma is the acceleration divergence (called ageostrophic 
!       vorticity).  Note: ql = zeta - f*h, where zeta is the relative
!       vorticity, f is the Coriolis frequency and h is the dimensionless
!       depth anomaly.

!       The full algorithm consists of the following modules:

!       swps.f90      : This source - main program to evolve fields;
!       parameters.f90: User defined parameters for a simulation;
!       constants.f90 : Fixed constants used throughout the other modules;
!       spectral.f90  : Fourier transform common storage and routines;

!----------------------------------------------------------------------------
program swps

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
double precision:: gsi(ng,ng),sgs(ng,ng),ndg(ng,ng)

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

 !PV fields:
 !Store qs^n & qs^n+0.25*dt*sqs(t^n) for efficiency in the iteration below:
qsi=qs
qsm=qs+dt4*sqs
 !Preliminary guess for qs^{n+1}:
qs=two*dissi*(qsm+dt4*sqs)-qsi

 !Divergence and acceleration divergence:
dsi=ds
gsi=gs
 !Intermediate arrays used in semi-implicit time stepping below:
nds=f12*sds+hdti*dsi
ndg=f12*sgs+hdti*(gsi+nds)
 !Preliminary guess for ds^{n+1} & gs^{n+1}:
sgs=simp*(ndg+f12*(sgs+hdti*sds))
ds=two*sgs-dsi
gs=two*(hdti*sgs-nds-f12*sds)-gsi

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call main_invert(qs,ds,gs,hh,uu,vv,zz)

   !Calculate the source terms (sqs,sds,sgs) for linearised PV (qs), 
   !divergence (ds) and acceleration divergence (gs):
  call source(sqs,sds,sgs)

   !Obtain new guess for the linearised PV (in spectral space):
  qs=two*dissi*(qsm+dt4*sqs)-qsi

   !Obtain new guesses for the divergence and acceleration divergence
   !(in spectral space) at time t^{n+1} using a semi-implicit scheme:
  sgs=simp*(ndg+f12*(sgs+hdti*sds))
  ds=two*sgs-dsi
  gs=two*(hdti*sgs-nds-f12*sds)-gsi
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
double precision:: dd(ng,ng),wkp(ng,ng),wkq(ng,ng)

 !Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng)

!---------------------------------------------------------------
 !qs source --- compute div(ql*u,ql*v) (wka in spectral space):
wka=qs
call spctop(ng,ng,wka,wkq,xfactors,yfactors,xtrig,ytrig)
 !wkq contains the linearised PV in physical space
wkp=wkq*uu
wkq=wkq*vv
 !Compute spectral divergence from physical fields:
call divs(wkp,wkq,wka)
sqs=-filt*wka

!---------------------------------------------------------------
 !Nonlinear part of ds source --- 2J(u,v) - div(delta*(u,v)):
call jacob(uu,vv,wkp)
 !Convert J(u,v) (in wkp) to spectral space as wkb:
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)

 !Compute div(delta*u,delta*v) (to be put into wka, in spectral space):
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)
 !dd contains the divergence in physical space
wkp=dd*uu
wkq=dd*vv
 !Compute spectral divergence from physical fields:
call divs(wkp,wkq,wka)

 !Add everything up to define delta_t - gamma (spectral, filtered)
sds=filt*(two*wkb-wka)

!---------------------------------------------------------------
 !Nonlinear part of gs source --- fN_zeta - c^2*Lap(N_h) where
 !N_zeta and N_h are the nonlinear parts of the vorticity and
 !depth sources:

 !Obtain N_h:
wkp=hh*uu
wkq=hh*vv
 !Compute div(h*u,h*v) = -N_h spectrally and filter for use below:
call divs(wkp,wkq,wka)
 !For use below, calculate -c^2*Lap(N_h):
wka=c2g2*wka

 !Obtain N_zeta:
wkp=zz*uu
wkq=zz*vv
 !Compute div(zeta*u,zeta*v) spectrally:
call divs(wkp,wkq,wkb)

 !Add everything up to define S_gamma:
sgs=wka-cof*filt*wkb

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
double precision:: wkp(ng,ng) !Physical
double precision:: wks(ng,ng) !Spectral
double precision:: zspec(0:ng),dspec(0:ng),gspec(0:ng)
double precision:: ekin,epot,etot
integer:: k

!---------------------------------------------------------------
 !Compute energy components and total:
ekin=f12*garea*sum((one+hh)*(uu**2+vv**2))
epot=f12*garea*csq*sum(hh**2)
etot=ekin+epot

 !Write energies to ecomp.asc:
write(15,'(f13.6,5(1x,f16.9))') t,zero,ekin,ekin,epot,etot
write(*,'(a,f13.6,a,f13.6)') ' t = ',t,'  E_tot = ',etot

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
end program swps

!=======================================================================
