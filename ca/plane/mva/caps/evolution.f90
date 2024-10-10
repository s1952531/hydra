module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

 !Spectral (linearised) PV fields (total, contour part, residual):
double precision:: qt(ng,ng),qc(ng,ng),qd(ng,ng)

 !Spectral operators used in time stepping (adapted):
double precision:: pdis(ng,ng),simp(ng,ng),disq(ng,ng),disb(ng,ng)

 !Physical space velocity divergence & horizontal magnetic field:
double precision:: dd(ng,ng),bx(ng,ng),by(ng,ng)

 !Logicals for saving gridded & contour data
logical:: gsave,csave

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
double precision,parameter:: qratmax=0.2d0
!      qratmax:  the maximum ratio r of the mean-square residual PV qd
!                to the mean-square PV qc in the contours
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, the code stops
!      to rebuild the contours in a separate memory space.
double precision:: wka(ng,ng),qa(ng,ng),qrat
integer:: ireg

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

   !Perform contour surgery or recontouring when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1

     !Compute ratio of mean-square residual and contour PV:
    call con2grid(qa)
    wka=qd
    call spctop(ng,ng,wka,qr,xfactors,yfactors,xtrig,ytrig)
    qrat=sum(qr**2)/sum(qa**2)
    
     !Don't continue if maximum number of regularisations reached:
    if (ireg .eq. nregmax .or. qrat .gt. qratmax) then
       !Prepare PV residual qr for recontouring (and preserve qs):
      call prepare
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the PV contours (surgery + node redistribution):
    call surgery
     !Record active contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

    twist=twist-twistmax
  endif

   !Advect flow from time t to t + dt:
  call advance
  
enddo
 !End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Save final data if not done already:
if (int(t/tgsave) .eq. igrids) call savegrid
if (int(t/tcsave) .eq. iconts) call savecont

return
end subroutine advect

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
 !Logicals used for saving gridded fields and contours:
gsave=.false.
csave=.false.

 !Convert PV contours (xq,yq) to gridded values as qa:
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(ng,ng,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Define (spectral) residual PV qd = (1-F)[qs-qc]:
qd=bfhi*(qs-qc)
 !Here bfhi = 1-F is a high-pass spectral filter

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current (spectral) PV anomaly field is stored in qs, and the 
! residual PV needed in congen.f90 is stored in qr.

implicit none

 !Local variables:
double precision:: qa(ng,ng)
double precision:: wka(ng,ng),qmin,ql1,ql2
integer:: ix,iy

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

!-------------------------------------------------------------------
 !Update PV contour interval for re-contouring; choose contour
 !interval based on <q^2>/<|q|> for |q| > q_rms:
wka=qs
call spctop(ng,ng,wka,qa,xfactors,yfactors,xtrig,ytrig)
qmin=sqrt(dsumi*sum(qa**2))
ql1=zero
ql2=zero
do ix=1,ng
  do iy=1,ng
    if (abs(qa(iy,ix)) .gt. qmin) then
      ql1=ql1+abs(qa(iy,ix))
      ql2=ql2+qa(iy,ix)**2
    endif
  enddo
enddo
qjump=ql2/(ql1*dble(ncont))

return
end subroutine prepare

!=======================================================================

subroutine advance

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method 
! (for all spectral fields, namely qs, qd, ds, gs, bxs & bys).

! Uses an iterative implicit method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}-F^n)/2] + N[(F^{n+1}-F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms, and N refers to the nonlinear source
! terms.  We start with a guess for F^{n+1} in N and iterate 
! niter times (see parameter statement below).

implicit none

 !Local variables:
integer,parameter:: niter=2

 !Spectral fields needed in time stepping:
double precision:: qsi(ng,ng),sqs(ng,ng)
double precision:: qdi(ng,ng),qdm(ng,ng),sqd(ng,ng)
double precision:: bxsi(ng,ng),bxsm(ng,ng),sbxs(ng,ng)
double precision:: bysi(ng,ng),bysm(ng,ng),sbys(ng,ng)
double precision:: dsi(ng,ng),sds(ng,ng),nds(ng,ng)
double precision:: gsi(ng,ng),sgs(ng,ng),ngs(ng,ng)
 !Contour positions needed in time stepping:
double precision:: xqi(nptq),yqi(nptq)
 !Contour velocities:
double precision:: uq(nptq),vq(nptq)
 !Other local quantities:
double precision:: xx,yy
integer:: i,iter

!-------------------------------------------------------------------
 !Invert PV and compute velocity at current time level, say t=t^n:
call inversion
 !This also returns qc, needed immediately below.

 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*(qs-qc) + qc + qd
 !            and qd = (1-F)*(qs-qc)
 !Here F is a low pass filter (see spectral.f90)
qs=bflo*qs+bfhi*qc+qd
qd=bfhi*(qs-qc)

 !Adapt the time step and save various diagnostics each time step:
call adapt

 !Possibly save data (gsave & csave set by adapt):
if (gsave) call savegrid
if (csave) call savecont

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
 !(sds,sgs) for divergence and acceleration divergence (ds,gs), and for
 !the horizontal magnetic field (sbxs,sbys):
call source(sqs,sqd,sds,sgs,sbxs,sbys)

 !Update PV fields:
qsi=qs+dt2*sqs
qs=qs+dt*sqs
qdi=qd
qdm=qd+dt4*sqd
qd=disq*(qdm+dt4*sqd)-qdi

 !Update magnetic field:
bxsi=bxs
bxsm=bxs+dt4*sbxs
bxs=disb*(bxsm+dt4*sbxs)-bxsi
bysi=bys
bysm=bys+dt4*sbys
bys=disb*(bysm+dt4*sbys)-bysi

 !Update divergence and acceleration divergence:
dsi=ds
gsi=gs
nds=sds+dt4i*dsi
ngs=sgs+dt4i*gsi
sds=nds+sds          !2*N_tilde_delta
sgs=ngs+sgs          !2*N_tilde_gamma
ds=simp*(pdis*sds+sgs)-dsi       !simp = 1/(P*R^2-G);  pdis = P*R
gs=simp*(pdis*sgs+pgop*sds)-gsi  !pgop = P*G

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call inversion

   !Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
   !(sds,sgs) for divergence and acceleration divergence (ds,gs), and for
   !the horizontal magnetic field (sbxs,sbys):
  call source(sqs,sqd,sds,sgs,sbxs,sbys)

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
  qd=disq*(qdm+dt4*sqd)-qdi

   !Update magnetic field:
  bxs=disb*(bxsm+dt4*sbxs)-bxsi
  bys=disb*(bysm+dt4*sbys)-bysi

   !Update divergence and acceleration divergence:
  sds=nds+sds          !2*N_tilde_delta
  sgs=ngs+sgs          !2*N_tilde_gamma
  ds=simp*(pdis*sds+sgs)-dsi       !simp = 1/(P*R^2-G);  pdis = P*R
  gs=simp*(pdis*sgs+pgop*sds)-gsi  !pgop = P*G
enddo

 !Advance time:
t=t+dt

return
end subroutine advance

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

 !Invert linearised PV, divergence and linearised acceleration divergence
 !to obtain the dimensionless height anomaly and velocity field, as well
 !as the relative vorticity (see spectral.f90):
call main_invert(qt,ds,gs,hh,uu,vv,zz)
 !Note: qt, ds & gs are in spectral space while 
 !      hh, uu, vv and zz are in physical space.

return
end subroutine inversion

!=======================================================================

subroutine adapt

! Adapts the time step, computes the twist parameter (the time integral
! of |zeta|_max) and various quantities every time step to monitor the
! flow evolution, and updates spectral operators needed in the time
! stepping of delta and gamma_l.

implicit none

 !Local variables:
double precision,parameter:: dtfac=pi/10.d0, cflmax=0.75d0
! The time step dt is chosen to be less than or equal to the minimum
! of dtfac/max(|zeta|_max) and cflmax*dx/max(|u|_max,|B|_max,c_gw)
! where dx is the grid spacing, u is the horizontal velocity field,
! B is the horizontal magnetic field, and c_gw is the short-scale
! inertia-gravity wave speed.

double precision:: rdis(ng,ng) !R operator (only needed temporarily)
double precision:: wka(ng,ng)
double precision:: bmax,umax,zmax,zrms,ro,fr,hmin,hmax

!----------------------------------------------------------------------
 !Get horizontal magnetic field in physical space:
wka=bxs
call spctop(ng,ng,wka,bx,xfactors,yfactors,xtrig,ytrig)
wka=bys
call spctop(ng,ng,wka,by,xfactors,yfactors,xtrig,ytrig)

 !Maximum horizontal magnetic field:
bmax=sqrt(maxval(bx**2+by**2))

 !Maximum horizontal velocity:
wka=uu**2+vv**2
umax=sqrt(maxval(wka))

 !Froude number:
fr=sqrt(maxval(wka/(one+hh)))/cgw

 !Maximum vorticity:
zmax=maxval(abs(zz))

 !R.m.s. vorticity:
zrms=sqrt(dsumi*sum(zz**2))

 !Rossby number:
ro=zmax/cof

 !Min/max dimensionless height anomaly:
hmin=minval(hh)
hmax=maxval(hh)

 !Write data:
write(16,'(1x,f12.5,4(1x,f12.8))') t,ro,fr,hmin,hmax
write(17,'(1x,f12.5,4(1x,f12.8))') t,bmax,umax,zmax,zrms

!---------------------------------------------------------------------
 !Compute new time step:
dt=min(gl*cflmax/max(umax,bmax,cgw),dtfac/(zmax+small),dtmax)
 !Note, dtmax is specified in constants.f90

 !Time step parameters:
dt2=dt*f12
dt4=dt*f14
dt2i=one/dt2
dt4i=one/dt4

 !Update operators needed in time stepping:
rdis=dt2i+(cof+zrms)*qdis
pdis=rdis/pope
simp=filt*pope/(rdis**2-opak*pope)
disq=two/(dt2*rdis)
disb=two/(one+dt2*bdis)

 !Increment the integral of |zeta|_max:
twist=twist+dt*zmax

!---------------------------------------------------------------------
 !Set flags to save data:
gsave=(int(t/tgsave) .eq. igrids)
 !Gridded data will be saved at time t if gsave is true.
csave=(int(t/tcsave) .eq. iconts)
 !Contour data will be saved at time t if csave is true.

return
end subroutine adapt

!=======================================================================

subroutine source(sqs,sqd,sds,sgs,sbxs,sbys)

! Gets the source terms (sqs,sqd) for the PV (qs,qd) as well as those
! (sds,sgs,sbxs,sbys) for divergence ds, acceleration divergence gs
! and the magnetic field (bxs,bys) --- all in spectral space.
! Note that sds & sgs only include the nonlinear terms for a semi-implicit
! treatment, closely analogous to that described in the appendix of
! Mohebalhojeh & Dritschel (2004); see also Dritschel & Jalali (2020).

! The fields ds, gs, qd, qs, bxs & bys are all spectrally truncated.
! Note, hh, uu, vv & zz obtained by inversion before calling this 
! routine are all spectrally truncated (but are in physical space).

! The horizontal magnetic field is available in physical space as
! (bx,by) due to a prior call to subroutine adapt.

implicit none

double precision,parameter:: toler=1.d-11
 !toler: maximum error in iteration below to find the vertically-
 !       integrated non-hydrostatic pressure \bar{P}_n (ppn below).
double precision,parameter:: w=0.75d0, wc=one-w
 !w: weight used in the pressure iteration to accelerate convergence

 !Passed variables:
double precision:: sqs(ng,ng),sqd(ng,ng),sds(ng,ng),sgs(ng,ng)
double precision:: sbxs(ng,ng),sbys(ng,ng)

 !Local variables (physical):
double precision:: ux(ng,ng),uy(ng,ng),vx(ng,ng),vy(ng,ng)
double precision:: rx(ng,ng),ry(ng,ng),jz(ng,ng)
double precision:: hinv(ng,ng),hinvsq(ng,ng)
double precision:: pnx(ng,ng),pny(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng),wkr(ng,ng)
double precision:: errpn

 !Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng),wkg(ng,ng)
double precision:: sql(ng,ng)

!---------------------------------------------------------------
 !Prepare for source calculations:

 !Get divergence in physical space as dd:
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)

 !Get (B_x,B_y) in physical space as (bx,by):
wka=bxs
call spctop(ng,ng,wka,bx,xfactors,yfactors,xtrig,ytrig)
wka=bys
call spctop(ng,ng,wka,by,xfactors,yfactors,xtrig,ytrig)

 !Get j_z = dB_y/dx - dB_x/dy in physical space as jz:
call xderiv(ng,ng,hrkx,bys,wka)
call yderiv(ng,ng,hrky,bxs,wkb)
wka=wka-wkb
call spctop(ng,ng,wka,jz,xfactors,yfactors,xtrig,ytrig)

 !Get q_l in physical space as wkq:
sql=qt
call spctop(ng,ng,sql,wkq,xfactors,yfactors,xtrig,ytrig)

 !Compute div(j_z(B_x,B_y) - q_l*(u,v)) in spectral space as sql:
rx=jz*bx-wkq*uu
ry=jz*by-wkq*vv
call divs(rx,ry,sql)

 !Get tau = dB_x/dx + dB_y/dy in physical space as wkp, compute
 !2*delta^2 + (B_x,B_y)*grad(tau) and store in sds for use below:
call xderiv(ng,ng,hrkx,bxs,wka)
call yderiv(ng,ng,hrky,bys,wkb)
wka=wka+wkb
 !Here wka = tau (spectral); compute grad(tau) -> (rx,ry):
call gradient(wka,rx,ry)
 !Store tau in wkp (physical):
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
 !Compute 2*delta^2 + (B_x,B_y)*grad(tau) -> sds (spectral):
wkr=bx*rx+by*ry !(B_x,B_y)*grad(tau) -> wkr (physical)
wkq=two*dd**2+wkr
call ptospc(ng,ng,wkq,sds,xfactors,yfactors,xtrig,ytrig)
 !Append -div(delta(u,v)) to sds (spectral):
rx=dd*uu
ry=dd*vv
call divs(rx,ry,wkb)
sds=sds-wkb
 !*** Do not re-use sds until it is used in the divergence tendency

 !Compute (1/2)*tau^2 and de-alias:
wkp=f12*wkp**2
 !Convert to spectral space:
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
wka=filt*wka
 !Compute grad{(1/2)*tau^2} as (rx,ry) in physical space:
call gradient(wka,rx,ry)
 !Return (1/2)*tau^2 to physical space as wkp:
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)

 !Compute (H^2/3)*(1+h)^2 and de-alias:
wkq=hbsq3*(one+hh)**2
call dealias(wkq)

 !Compute gamma_tilde_b in physical space as hinv (temporary), where
 !gamma_tilde_b = div{F_b} - (B_x,B_y)*grad(tau) and where
 !F_b = j_z(-B_y,B_x) - (H^2/6)*div{(1+h)^2*grad{tau^2}:
rx=-jz*by-wkq*rx
ry= jz*bx-wkq*ry
call divs(rx,ry,wka) !wka = div{F_b} in spectral space
call spctop(ng,ng,wka,hinv,xfactors,yfactors,xtrig,ytrig)
hinv=hinv-wkr !Note: wkr = (B_x,B_y)*grad(tau) on rhs here

 !Compute (H^2/6)*J(tau^2,(1+h)^2) in physical space as jz (overwrite):
call jacob(wkp,wkq,jz)
 !*** do not use jz again until it is used to define sql below

!---------------------------------------------------------------
 !Calculate the vertically-integrated non-hydrostatic pressure.

 !First compute gamma_tilde in physical space as wkp, where
 !gamma_tilde = gamma_l + 2*(J(u,v) - delta^2) + gamma_tilde_b:
wkp=uu
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,ux,uy) !ux = du/dx & uy = du/dy (physical)
vx=zz+uy !vx = dv/dx using definition of vorticity
vy=dd-ux !vy = dv/dy using definition of diverence
 !Form 2*(J(u,v) - delta^2) + gamma_tilde_b:
wkp=two*(ux*vy-uy*vx-dd**2)+hinv !hinv = gamma_tilde_b from above
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
wkg=gs+filt*wkg
call spctop(ng,ng,wkg,wkp,xfactors,yfactors,xtrig,ytrig)
 !wkp now contains gamma_tilde in physical space (de-aliased)

 !Multiply gamma_tilde next by (1+h) and re-store in wkg (spectral)
 !as the fixed rhs in the pressure iteration below:
wkq=one+hh
wkp=wkq*wkp
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
 !wkg is not de-aliased but the prop operator below takes care of this.

 !Next calculate wkq = 3/H^2*(1/(1+h)^2 - 1) needed for the pressure
 !iteration below:
hinv=one/wkq !wkq = 1+h from above
call dealias(hinv)
hinvsq=hinv**2
call dealias(hinvsq)
wkq=hbsq3i*(hinvsq-one)

 !Calculate also (1+h)^{-1}*grad(h) and store in rx & ry:
wkp=hh
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,rx,ry)
rx=hinv*rx
call dealias(rx)
ry=hinv*ry
call dealias(ry)

 !Now iterate to find \bar{P}_n (in ppn) starting from the guess
 !ppn = (grad^2 - 3/H^2)^{-1}((1+h)*gamma_tilde):
wka=prop*wkg !Here, prop = (grad^2 - 3/H^2)^{-1}
call spctop(ng,ng,wka,ppn,xfactors,yfactors,xtrig,ytrig)
errpn=two*toler
do while (errpn .gt. toler)
  wkp=ppn
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call gradient(wka,pnx,pny)
  wkp=rx*pnx+ry*pny+wkq*ppn
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  wka=prop*(wkg+wkb)
  call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
  errpn=sqrt(dsumi*sum((wkp-ppn)**2))
  ppn=w*wkp+wc*ppn
enddo

 !Store rhs of pressure iteration, minus gs, in wkb:
wkb=wkg+wkb-gs
 !*** do not re-use wkb (needed below for the divergence tendency)

!---------------------------------------------------------------
 !qd source --- S_q + (u,v)*grad{q_l - q_d}:

 !Complete S_q; first compute J(ppn,1/(1+h)) and store in wkp:
call jacob(ppn,hinv,wkp)

 !Add (H^2/6)*J((1+h)^2,tau^2) in the array jz from above:
wkp=wkp+jz

 !Convert to spectral space as wka:
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)

 !Complete definition of S_q (in sql) after de-aliasing:
sql=filt*(sql+wka)

 !Form grad{q_l - q_d} as (rx,ry) in physical space:
wka=qt-qd
call gradient(wka,rx,ry)

 !Take inner product with velocity field:
wkp=uu*rx+vv*ry
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)

 !Filter and add S_q to complete qd source:
sqd=filt*wka+sql

!---------------------------------------------------------------
 !qs source --- only NL advection term is needed:

call gradient(qs,rx,ry)
wkp=-uu*rx-vv*ry
 !Convert to spectral space:
call ptospc(ng,ng,wkp,sqs,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sqs=filt*sqs

!---------------------------------------------------------------
 !Nonlinear part of ds source (part of this is in sds already):

 !Compute 1/(1+h)^3 -> wkp and de-alias:
wkp=hinv*hinvsq
call dealias(wkp)
 !Form \bar{P}_n*(1/(1+h)^3 - 1):
wkp=ppn*(wkp-one)
 !Transform to spectral space as wkg:
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)

 !Add everything up to define delta_t - pope*gamma_l (spectral, filtered)
sds=filt*(sds-hbsq3i*wkg+pope*wkb)
 !Here sds = 2*delta^2 + (B_x,B_y)*grad(tau) - div(delta(u,v)),
 !pope = (1 - (H^2/3)*grad^2)^{-1}, and wkb is the rhs of the pressure
 !iteration, minus gamma_l (all in spectral space).

!---------------------------------------------------------------
 !Nonlinear part of gs source --- G{div{h*(u,v)}} + f*S_q, where
 !G = c^2*grad^2 - f^2 is the spectral gravity-wave operator and
 !S_q = dq_l/dt (partial time derivative).

 !Compute div(h*(u,v)) -> wka (spectral):
rx=hh*uu
ry=hh*vv
call divs(rx,ry,wka)

 !Apply G and add f*S_q (S_q is in sql from above) to define S_gamma:
sgs=gwop*wka+cof*sql
 !Here gwop = G (spectrally filtered)

!---------------------------------------------------------------
 !bxs source (except diffusion), (B_x,B_y)*grad(u)-(u,v)*grad(B_x):

call gradient(bxs,wkp,wkq)
 !(wkp,wkq) = grad{B_x} in physical space
wkp=bx*ux+by*uy-wkp*uu-wkq*vv
call ptospc(ng,ng,wkp,sbxs,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sbxs=filt*sbxs

!---------------------------------------------------------------
 !bys source (except diffusion), (B_x,B_y)*grad(v)-(u,v)*grad(B_y):

call gradient(bys,wkp,wkq)
 !(wkp,wkq) = grad{B_y} in physical space
wkp=bx*vx+by*vy-wkp*uu-wkq*vv
call ptospc(ng,ng,wkp,sbys,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
sbys=filt*sbys

return
end subroutine source

!=======================================================================

subroutine nhpsolve

! Finds the NH pressure (used only diagnostically in data writes).
! The physical space diverence is assumed to be available in dd.

! The horizontal magnetic field is available in physical space as
! (bx,by) due to a prior call to subroutine adapt.
  
implicit none

double precision,parameter:: toler=1.d-11
 !toler: maximum error in iteration below to find the vertically-
 !       integrated non-hydrostatic pressure \bar{P}_n (ppn below).
double precision,parameter:: w=0.75d0, wc=one-w
 !w: weight used in the pressure iteration to accelerate convergence
double precision:: rx(ng,ng),ry(ng,ng),jz(ng,ng)
double precision:: hinv(ng,ng),hinvsq(ng,ng)
double precision:: pnx(ng,ng),pny(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng),wkr(ng,ng)
double precision:: errpn

 !Local variables (spectral):
double precision:: wka(ng,ng),wkb(ng,ng),wkg(ng,ng)

!---------------------------------------------------------------
 !Get j_z = dB_y/dx - dB_x/dy in physical space as jz:
call xderiv(ng,ng,hrkx,bys,wka)
call yderiv(ng,ng,hrky,bxs,wkb)
wka=wka-wkb
call spctop(ng,ng,wka,jz,xfactors,yfactors,xtrig,ytrig)

 !Get tau = dB_x/dx + dB_y/dy in physical space as wkp, compute
 !2*delta^2 + (B_x,B_y)*grad(tau) and store in sds for use below:
call xderiv(ng,ng,hrkx,bxs,wka)
call yderiv(ng,ng,hrky,bys,wkb)
wka=wka+wkb
 !Here wka = tau (spectral); compute grad(tau) -> (rx,ry):
call gradient(wka,rx,ry)
 !Store tau in wkp (physical):
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
 !Compute (B_x,B_y)*grad(tau) -> wkr (physical):
wkr=bx*rx+by*ry

 !Compute (1/2)*tau^2 and de-alias:
wkp=f12*wkp**2
 !Convert to spectral space:
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
wka=filt*wka
 !Compute grad{(1/2)*tau^2} as (rx,ry) in physical space:
call gradient(wka,rx,ry)
 !Return (1/2)*tau^2 to physical space as wkp:
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)

 !Compute (H^2/3)*(1+h)^2 and de-alias:
wkq=hbsq3*(one+hh)**2
call dealias(wkq)

 !Compute gamma_tilde_b in physical space as hinv (temporary), where
 !gamma_tilde_b = div{F_b} - (B_x,B_y)*grad(tau) and where
 !F_b = j_z(-B_y,B_x) - (H^2/6)*div{(1+h)^2*grad{tau^2}:
rx=-jz*by-wkq*rx
ry= jz*bx-wkq*ry
call divs(rx,ry,wka) !wka = div{F_b} in spectral space
call spctop(ng,ng,wka,hinv,xfactors,yfactors,xtrig,ytrig)
hinv=hinv-wkr !Note: wkr = (B_x,B_y)*grad(tau) on rhs here

!---------------------------------------------------------------
 !Calculate the vertically-integrated non-hydrostatic pressure.

 !First compute gamma_tilde in physical space as wkp, where
 !gamma_tilde = gamma_l + 2*(J(u,v) - delta^2) + gamma_tilde_b:
call jacob(uu,vv,wkp)
 !wkp contains J(u,v) in physical space temporarily
wkp=two*(wkp-dd**2)+hinv !hinv = gamma_tilde_b from above
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
wkg=gs+filt*wkg
call spctop(ng,ng,wkg,wkp,xfactors,yfactors,xtrig,ytrig)
 !wkp now contains gamma_tilde in physical space (de-aliased)

 !Multiply next by (1+h) and re-store in wkg (spectral) as the 
 !fixed rhs in the pressure iteration below:
wkq=one+hh
wkp=wkq*wkp
call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
 !wkg is not de-aliased but the prop operator below takes care of this.

 !Next calculate wkq = 3/H^2*(1/(1+h)^2 - 1) needed for the pressure
 !iteration below:
hinv=one/wkq
call dealias(hinv)
hinvsq=hinv**2
call dealias(hinvsq)
wkq=hbsq3i*(hinvsq-one)

 !Calculate also (1+h)^{-1}*grad(h) and store in rx & ry:
wkp=hh
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,rx,ry)
rx=hinv*rx
call dealias(rx)
ry=hinv*ry
call dealias(ry)

 !Now iterate to find \bar{P}_n (in ppn) starting from the guess
 !ppn = (grad^2 - 3/H^2)^{-1}((1+h)*gamma_tilde):
wka=prop*wkg
call spctop(ng,ng,wka,ppn,xfactors,yfactors,xtrig,ytrig)
errpn=two*toler
do while (errpn .gt. toler)
  wkp=ppn
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call gradient(wka,pnx,pny)
  wkp=rx*pnx+ry*pny+wkq*ppn
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  wka=prop*(wkg+wkb)
  call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
  errpn=sqrt(dsumi*sum((wkp-ppn)**2))
  ppn=w*wkp+wc*ppn
enddo

return
end subroutine nhpsolve

!=======================================================================

subroutine savegrid

! Saves PV, energy and various spectra at the desired save time

implicit none

 !Local variables:
double precision:: wkp(ng,ng),wkq(ng,ng),tt(ng,ng)  !Physical
double precision:: wka(ng,ng),wkb(ng,ng)            !Spectral
double precision:: spec(0:ng)                       !For 1D spectra
double precision:: ekin,epot,emag,etot              !Energy components
integer:: k

!---------------------------------------------------------------
 !Increment counter for direct file access:
igrids=igrids+1

 !Update fields:
call inversion

!Get delta = div(u,v) in physical space as dd:
wka=ds
call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)

 !Get tau = div(B_x,B_y) in physical space as tt:
call xderiv(ng,ng,hrkx,bxs,wka)
call yderiv(ng,ng,hrky,bys,wkb)
wka=wka+wkb
call spctop(ng,ng,wka,tt,xfactors,yfactors,xtrig,ytrig)

 !Compute energy components and total:
wkq=one+hh
ekin=f12*garea*sum(wkq*(uu**2+vv**2+hbsq3*(wkq*dd)**2))
emag=f12*garea*sum(wkq*(bx**2+by**2+hbsq3*(wkq*tt)**2))
epot=f12*garea*csq*sum(hh**2)
etot=ekin+epot+emag

 !Write energies to ecomp.asc:
write(15,'(f9.2,4(1x,f16.9))') t,ekin,epot,emag,etot
write(*,'(a,f9.2,a,f13.6)') ' t = ',t,'  E_tot = ',etot

!---------------------------------------------------------------
 !Write various gridded fields to direct access files:
wka=qt
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
write(31,rec=igrids) t,wkp !linearised PV
write(32,rec=igrids) t,dd  !divergence
wka=gs
call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
write(33,rec=igrids) t,wkp !linearised acceleration divergence
write(34,rec=igrids) t,hh  !dimensionless height anomaly
write(35,rec=igrids) t,zz  !relative vorticity
call nhpsolve
write(36,rec=igrids) t,ppn !non-hydrostatic pressure

!---------------------------------------------------------------
 !Compute 1d spectra for various fields:
wkp=zz
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call spec1d(wka,spec)
spec=log10(spmf*spec+1.d-32)
write(51,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(51,'(2(1x,f12.8))') alk(k),spec(k)
enddo

call spec1d(ds,spec)
spec=log10(spmf*spec+1.d-32)
write(52,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(52,'(2(1x,f12.8))') alk(k),spec(k)
enddo

call spec1d(gs,spec)
spec=log10(spmf*spec+1.d-32)
write(53,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(53,'(2(1x,f12.8))') alk(k),spec(k)
enddo

wkp=hh
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call spec1d(wka,spec)
spec=log10(spmf*spec+1.d-32)
write(54,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(54,'(2(1x,f12.8))') alk(k),spec(k)
enddo

wkp=ppn
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call spec1d(wka,spec)
spec=log10(spmf*spec+1.d-32)
write(55,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(55,'(2(1x,f12.8))') alk(k),spec(k)
enddo

 !j_z = dB_y/dx - dB_x/dy:
call xderiv(ng,ng,hrkx,bys,wka)
call yderiv(ng,ng,hrky,bxs,wkb)
wkb=wka-wkb
call spec1d(wkb,spec)
spec=log10(spmf*spec+1.d-32)
write(55,'(f12.5,1x,i5)') t,kmaxred
do k=1,kmaxred
  write(55,'(2(1x,f12.8))') alk(k),spec(k)
enddo

 !Convert j_z to physical space and write:
call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
write(37,rec=igrids) t,wkp
write(38,rec=igrids) t,tt  !tau = div(B_x,B_y)

 !spmf takes into account uneven sampling of wavenumbers in each
 !shell [k-1/2,k+1/2].

 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated

 !alk(k) = log_10(k)

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont

! Saves PV contours for post-processing and imaging

implicit none

 !Local variables:
double precision:: ss(ng,ng)
double precision:: qa(ng,ng)
character(len=3):: pind

!---------------------------------------------------------------
 !Increment counter for direct file access:
iconts=iconts+1

write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') iconts-1

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
ss=qt-qc
call spctop(ng,ng,ss,qa,xfactors,yfactors,xtrig,ytrig)
write(83,rec=iconts) t,qa

 !Save PV contours:
open(81,file='contours/qqindex'//pind,form='unformatted',status='replace')
write(81) npq(1:nq),i1q(1:nq),indq(1:nq)
close(81)

open(82,file='contours/qqnodes'//pind,form='unformatted',status='replace')
write(82) xq(1:nptq),yq(1:nptq)
close(82)

return
end subroutine savecont

!=======================================================================

 !Main end module
end module evolution
