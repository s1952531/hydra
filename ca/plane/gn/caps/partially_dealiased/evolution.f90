module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

double precision:: qt(ng,ng),qc(ng,ng),qd(ng,ng) !Various PVs

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
    call inversion
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
    write(14,'(1x,f12.5,1x,i9,1x,i10)') t,nq,nptq

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

 !Local variables:
double precision:: qa(ng,ng)

!------------------------------------------------------------------
 !Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

!-------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qa:
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(ng,ng,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Define (spectral) residual PV qd = (1-F)[qs-qc]:
qd=bfhi*(qs-qc)
 !Here bfhi = 1-F is a high-pass spectral filter

return
end subroutine

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current (spectral) PV anomaly field is stored in qs, and the 
! residual PV needed in congen.f90 is stored in qr.

implicit none

 !Local variables:
double precision:: qa(ng,ng)

!-----------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qa:
call con2grid(qa)

 !Convert qa to spectral space as qc (note, qa is overwritten):
call ptospc(ng,ng,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Define spectral PV anomaly (qs) and PV residual (qd):
qs=bflo*qs+bfhi*qc+qd
qd=qs-qc

 !Convert qd to physical space as qr (used in recontouring):
call spctop(ng,ng,qd,qr,xfactors,yfactors,xtrig,ytrig)
 !Note: qd is overwritten, but we are leaving this module next
 !      and qd will be redefined upon re-entry in subroutine init.

return
end subroutine

!=======================================================================

subroutine advance(ggen)

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method 
! (for all spectral fields, namely qs, qd, ds & gs).

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

 !Spectral fields needed in time stepping:
double precision:: qsi(ng,ng),sqs(ng,ng)
double precision:: qdi(ng,ng),qdm(ng,ng),sqd(ng,ng)
double precision:: dsi(ng,ng),sds(ng,ng),nds(ng,ng)
double precision:: gsi(ng,ng),sgs(ng,ng),ngs(ng,ng)
 !Physical fields:
double precision:: qx(ng,ng),qy(ng,ng)
 !Contour positions needed in time stepping:
double precision:: xqi(nptq),yqi(nptq)
 !Contour velocities:
double precision:: uq(nptq),vq(nptq)
 !Other local quantities:
double precision:: xx,yy
integer:: i,iter

!-------------------------------------------------------------------
 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*(qs-qc) + qc + qd
 !            and qd = (1-F)*(qs-qc)
 !Here F is a low pass filter (see spectral.f90)
qs=bflo*qs+bfhi*qc+qd
qd=bfhi*(qs-qc)

!------------------------------------------------------------------
 !Invert PV and compute velocity at current time level, say t=t^n:
if (ggen) call inversion
 !If ggen is false, inversion was called previously at this time level.

 !Compute twist parameter and save various diagnostics each time step:
call diagnose

!------------------------------------------------------------------
 !Start with a guess for F^{n+1} for all contours and fields:

 !Contours:
call velint(uu,vv,uq,vq)
 !Here (uq,vq) stands for the velocity at time level n, i.e. u(x^n,t^n)
do i=1,nptq
   !Store x^n+0.5*dt*u(x^n,t^n) for efficiency in the iteration below:
  xx=xq(i)+dt2*uq(i)
  yy=yq(i)+dt2*vq(i)
  xqi(i)=oms*(xx-twopi*dble(int(xx*pinv)))
  yqi(i)=oms*(yy-twopi*dble(int(yy*pinv)))
   !Preliminary guess for x^{n+1}:
  xx=xq(i)+dt*uq(i)
  yy=yq(i)+dt*vq(i)
  xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
  yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
enddo

 !Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
 !(sds,sgs) for divergence and acceleration divergence (ds,gs):
call source(sqs,sqd,sds,sgs)

 !Update PV fields:
qsi=qs+dt2*sqs
qs=qs+dt*sqs
qdi=qd
qdm=qd+dt4*sqd
qd=diss*(qdm+dt4*sqd)-qdi

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
  call inversion

   !Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
   !(sds,sgs) for divergence and acceleration divergence (ds,gs):
  call source(sqs,sqd,sds,sgs)

   !Interpolate gridded velocity (uu,vv) at contour nodes as (uq,vq):
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqi(i)+dt2*uq(i)
    yy=yqi(i)+dt2*vq(i)
    xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
    yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
  enddo
   !Now (xq,yq) contain a new guess for x^{n+1}.

   !Update PV fields:
  qs=qsi+dt2*sqs
  qd=diss*(qdm+dt4*sqd)-qdi

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

subroutine source(sqs,sqd,sds,sgs)

! Gets the source terms (sqs,sqd) for the PV (qs,qd) as well as those
! (sds,sgs) for divergence and acceleration divergence (ds,gs) ---
! --- all in spectral space.  Note that (sds,sgs) only include the
! nonlinear terms for a semi-implicit treatment, closely analogous
! to that described in the appendix of Mohebalhojeh & Dritschel (2004).

! The spectral fields ds, gs, qd and qs are all spectrally truncated.
! Note, hh, uu, vv & zz obtained by main_invert before calling this 
! routine are all spectrally truncated.

implicit none

 !Passed variables:
double precision:: sqs(ng,ng),sqd(ng,ng),sds(ng,ng),sgs(ng,ng)

 !Local variables (physical):
double precision:: htot(ng,ng),hx(ng,ng),hy(ng,ng),ht(ng,ng)
double precision:: dd(ng,ng),ddt(ng,ng),aa(ng,ng),bb(ng,ng)
double precision:: ut(ng,ng),vt(ng,ng),ztn(ng,ng)
double precision:: wkg(ng,ng),wkh(ng,ng),wkm(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng),wkz(ng,ng)
double precision:: errgt

 !Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng),wkf(ng,ng)

double precision,parameter:: tolgt=1.d-10*cof**3
 !tolgt: maximum error in iteration below to find gamma_t; the f^3 factor
 !       comes from the units of gamma (1/T^2) and d/dt (1/T).
 !       Note: gamma_t is returned in sgs (spectral).

!---------------------------------------------------------------
 !qd source --- only NL advection term is needed:
call gradient(qd,hx,hy)
wkp=-uu*hx-vv*hy
 !Convert to spectral space:
call ptospc(ng,ng,wkp,sqd,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sqd=filt*sqd

!---------------------------------------------------------------
 !qs source --- only NL advection term is needed:
call gradient(qs,hx,hy)
wkp=-uu*hx-vv*hy
 !Convert to spectral space:
call ptospc(ng,ng,wkp,sqs,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sqs=filt*sqs

 !Define h_tot = 1 + hh:
htot=one+hh

!---------------------------------------------------------------
 !Nonlinear part of ds source --- 2J(u,v) - div(delta*(u,v)):
call jacob(uu,vv,wkp)
 !Convert J(u,v) (in wkp) to spectral space as wkc:
call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)

 !Get divergence, delta, in physical space as dd:
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)

 !Compute div(delta*u,delta*v) (to be put into wka, in spectral space):
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
 !Define gamma_tilde = gamma + 2*(J(u,v) - delta^2) in physical space as aa:
wka=gs+two*wkc
call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)
aa=aa-two*dd**2
 !Spectrally truncate gamma_tilde:
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
wkp=zz*uu
wkq=zz*vv
 !Compute div(zeta*u,zeta*v) spectrally:
call divs(wkp,wkq,wkd)
 !Spectrally truncate:
wkd=filt*wkd
 !Convert to physical space as ztn:
call spctop(ng,ng,wkd,ztn,xfactors,yfactors,xtrig,ytrig)
 !Compute J(gamma_tilde,h):
call jacob(aa,hh,wkp)
 !Add (H^2/3)*(1+h)*J(gamma_tilde,h) to -div(zeta*u,zeta*v):
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
 !Compute (1+h)^3*gamma_tilde = (1+h)^2 * (1+h)*gamma_tilde (for de-aliasing):
wkh=htot**2
call ptospc(ng,ng,wkh,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka
call spctop(ng,ng,wka,wkh,xfactors,yfactors,xtrig,ytrig)
 !wkh=(1+h)^2 (spectrally truncated)
wkf=htot*aa
call ptospc(ng,ng,wkf,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka
call spctop(ng,ng,wka,wkf,xfactors,yfactors,xtrig,ytrig)
 !wkf=(1+h)*gamma_tilde (spectrally truncated) *** needed below ****
wkh=wkh*wkf
call ptospc(ng,ng,wkh,wka,xfactors,yfactors,xtrig,ytrig)
 !Spectrally truncate:
wka=filt*wka
 !Compute derivatives of (1+h)^3*gamma_tilde (in the array wka):
call xderiv(ng,ng,hrkx,wka,wkd)
call spctop(ng,ng,wkd,wkg,xfactors,yfactors,xtrig,ytrig)
 !wkg is ((1+h)^3*gamma_tilde)_x in physical space
call yderiv(ng,ng,hrky,wka,wkd)
call spctop(ng,ng,wkd,wkh,xfactors,yfactors,xtrig,ytrig)
 !wkh is ((1+h)^3*gamma_tilde)_y in physical space
 !Add all terms up in physical space to obtain u_t & v_t:
wkm=hbsq3/htot
 !Spectrally truncate wkm:
call ptospc(ng,ng,wkm,wka,xfactors,yfactors,xtrig,ytrig)
wka=filt*wka
call spctop(ng,ng,wka,wkm,xfactors,yfactors,xtrig,ytrig)
wkz=zz+cof
ut=wkm*wkg-wkp+wkz*vv-csq*hx
vt=wkm*wkh-wkq-wkz*uu-csq*hy

 !Obtain B = 2(1+h)*(J(u_t,v)+J(u,v_t)-2*delta*delta_t) (in the array bb)
call jacob(ut,vv,wkp)
call jacob(uu,vt,wkq)
bb=two*htot*(wkp+wkq-two*dd*ddt)
 !Spectrally truncate B:
call ptospc(ng,ng,bb,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call spctop(ng,ng,wkb,bb,xfactors,yfactors,xtrig,ytrig)

 !Obtain f*N_zeta in spectral space:
call ptospc(ng,ng,ztn,wkb,xfactors,yfactors,xtrig,ytrig)

 !Obtain f*N_zeta - c^2*Lap{N_h} + (H^2/3)*Lap{(1+h)*(2*h_t*gamma_tilde+B)}:
wkp=two*ht*wkf+htot*bb ! wkf = (1+h)*gamma_tilde here
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
 !Apply spectral Laplacian, -(k^2+l^2), to wka and add to 
 !wke = -c^2*Lap(N_h) computed above:
wka=cof*wkb+wke-hbsq3*rksq*wka

 !Obtain (H^2/3)*div(px,py) where px = gamma_tilde*((1+h)*h_t)_x+B*h_x
 !                            and py = gamma_tilde*((1+h)*h_t)_y+B*h_y
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
 !f*N_zeta - c^2*Lap{N_h} + (H^2/3)*Lap{(1+h)*(2*h_t*gamma_tilde+B)}
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
 !wkz = (1+h)^2-1 = h^2+2*h here and below
aa=htot*hx
call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
wke=filt*wke
call spctop(ng,ng,wke,aa,xfactors,yfactors,xtrig,ytrig)
 !aa = (1+h)*h_x here and below
bb=htot*hy
call ptospc(ng,ng,bb,wke,xfactors,yfactors,xtrig,ytrig)
wke=filt*wke
call spctop(ng,ng,wke,bb,xfactors,yfactors,xtrig,ytrig)
 !bb = (1+h)*h_y here and below

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

subroutine inversion

! Finds the gridded dimensionless height anomaly (hh) and velocity
! field (uu,vv) from the PV contours and the (spectral) divergence (ds)
! and acceleration divergence (gs).

implicit none

 !Local variables:
double precision:: qa(ng,ng)

!------------------------------------------------------------
 !Call con2grid to get updated contour PV (qc):
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is overwritten):
call ptospc(ng,ng,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Combine fields to update qt with full (spectral) field,
 !qt = F[qs-qc]+qc+qd, where F is a low pass filter:
qt=bflo*qs+bfhi*qc+qd

 !Invert PV, divergence and acceleration divergence to obtain the
 !dimensionless height anomaly and velocity field, as well as the
 !gridded PV anomaly and relative vorticity (see spectral.f90):
call main_invert(qt,ds,gs,hh,uu,vv,qq,zz)
 !Note: qt, ds & gs are in spectral space while 
 !      hh, uu, vv, qq and zz are in physical space.

return
end subroutine

!=======================================================================

subroutine diagnose

! Computes the twist parameter, the time integral of |zeta|_max, and
! various quantities at every time step to monitor the flow evolution.

implicit none

 !Local variables:
double precision:: umax,zrms,zmax
double precision:: uio,vio

!----------------------------------------------------------------------
 !Compute diagnostics:
umax=sqrt(maxval(uu**2+vv**2))
zmax=maxval(abs(zz))
zrms=sqrt(dsumi*sum(zz**2))

 !Increment the integral of |zeta|_max:
twist=twist+dt*zmax

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

write(51,'(f9.2,1x,i5)') t,kmaxred
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(51,'(4(1x,f12.8))') alk(k),log10(zspec(k)),log10(dspec(k)+1.d-32), &
                                                   log10(gspec(k)+1.d-32)
enddo
 !Note: alk(k) = log_10(k)

!---------------------------------------------------------------
!Write various gridded fields to direct access files:
write(31,rec=igrids) real(t),real(qq)
wks=ds
call spctop(ng,ng,wks,wkp,xfactors,yfactors,xtrig,ytrig)
write(32,rec=igrids) real(t),real(wkp)
wks=gs
call spctop(ng,ng,wks,wkp,xfactors,yfactors,xtrig,ytrig)
write(33,rec=igrids) real(t),real(wkp)
write(34,rec=igrids) real(t),real(hh)
write(35,rec=igrids) real(t),real(zz)

return
end subroutine

!=======================================================================
      
subroutine savecont(irec)

! Saves PV contours for post-processing and imaging

implicit none

 !Passed variable:
integer:: irec

 !Local variables:
double precision:: ss(ng,ng)
double precision:: qa(ng,ng)
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') irec

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
ss=qt-qc
call spctop(ng,ng,ss,qa,xfactors,yfactors,xtrig,ytrig)
write(83,rec=irec) real(t),real(qa)

 !Save PV contours:
open(81,file='cont/qqindex'//pind,form='unformatted',status='replace')
write(81) npq(1:nq),i1q(1:nq),indq(1:nq)
close(81)

open(82,file='cont/qqnodes'//pind,form='unformatted',status='replace')
write(82) xq(1:nptq),yq(1:nptq)
close(82)

return
end subroutine

!=======================================================================

 !Main end module
end module
