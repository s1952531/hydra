module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

 !PV field after conversion from contours and residual as in CLAM:
double precision:: qc(0:ny,0:nxm1),qd(0:ny,0:nxm1)
 !h_tot = 1 + hh and physical-space copy of cc:
double precision:: hht(0:ny,0:nxm1),ccp(0:ny,0:nxm1)
 !Gradient of the dimensionless height anomaly:
double precision:: hx(0:ny,0:nxm1),hy(0:ny,0:nxm1)
 !Relative vorticity and velocity divergence:
double precision:: zz(0:ny,0:nxm1),dd(0:ny,0:nxm1)

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
end subroutine advect

!=======================================================================

subroutine init

! Initialises residual PV for normal time integration following 
! contour regeneration

implicit none

! Work array:
double precision:: wk(0:ny,0:nxm1)

!------------------------------------------------------------------
 !Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

!-------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)
 !Convert qc to semi-spectral space:
call forfft(nyp1,nx,qc,xtrig,xfactors)

 !Define (semi-spectral) residual PV qd = (1-F)[qs-qc]:
qd=qs-qc
wk=qd
call lowpass(wk,1)
 !wk now contains F[qs-qc] in semi-spectral space.
qd=qd-wk

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current (semi-spectral) PV field is stored in qs, and the
! residual PV needed in congen.f90 is stored in qr.

implicit none

! Work array:
double precision:: wk(0:ny,0:nxm1)

!-----------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)

 !Convert qc to semi-spectral space:
call forfft(nyp1,nx,qc,xtrig,xfactors)

 !Define spectral PV (qs) and PV residual (qr):
wk=qs-qc
call lowpass(wk,1)
qr=qd+wk
qs=qc+qr

 !Convert qr to physical space:
call revfft(nyp1,nx,qr,xtrig,xfactors)
 !Note: we are leaving this module next and qd will be redefined
 !upon re-entry in subroutine init.

return
end subroutine prepare

!=======================================================================

subroutine advance(ggen)

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-semi-spectral method 
! (for all semi-spectral fields, namely qs, qd, cc, gg, uum & uup).

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
double precision:: qsi(0:ny,0:nxm1),sqs(0:ny,0:nxm1)
double precision:: qdi(0:ny,0:nxm1),qdm(0:ny,0:nxm1),sqd(0:ny,0:nxm1)
double precision:: cci(0:ny,0:nxm1),scc(0:ny,0:nxm1),ncc(0:ny,0:nxm1)
double precision:: ggi(0:ny,0:nxm1),sgg(0:ny,0:nxm1),ngg(0:ny,0:nxm1)
double precision:: uumi(0:nxm1),suum(0:nxm1),auum(0:nxm1),uumm(0:nxm1)
double precision:: uupi(0:nxm1),suup(0:nxm1),auup(0:nxm1),uupm(0:nxm1)

 !Contour positions and velocities needed in time stepping:
double precision:: xqi(nptq),yqi(nptq),uq(nptq),vq(nptq)

 !Other local quantities:
double precision:: xx,yy
integer:: i,j,iter

!-------------------------------------------------------------------
 !Compute height and velocity at current time level, say t=t^n:
if (ggen) call inversion
 !If ggen is false, inversion was called previously at this time level.

 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*(qs-qc) + qc + qd
 !            and qd = (1-F)*(qs-qc)
 !Here F is a low pass filter (see spectral.f90)
qsi=qs-qc
call lowpass(qsi,1)
qdi=qsi+qd
qs=qdi+qc
qd=qdi
call lowpass(qdi,1)
qd=qd-qdi

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
  xqi(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  yqi(i)=min(ymax,max(ymin,yy))
   !Preliminary guess for x^{n+1}:
  xx=xq(i)+dt*uq(i)
  yy=yq(i)+dt*vq(i)
  xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  yq(i)=min(ymax,max(ymin,yy))
enddo

 !Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
 !(scc,sgg) for momentum divergence and acceleration divergence (cc,gg):
call source(sqs,sqd,scc,sgg,suum,suup)

 !Update PV fields:
qsi=qs+dt2*sqs
qs=qs+dt*sqs
qdi=qd
qdm=qd+dt4*sqd
do j=0,ny
  qd(j,:)=diss*(qdm(j,:)+dt4*sqd(j,:))-qdi(j,:)
enddo

 !Update zonal velocity at y boundaries:
uumi=uum
uumm=uum+dt4*suum
uum=diss*(uumm+dt4*suum)-uumi
auum=dt2i*(uum-uumi)
uupi=uup
uupm=uup+dt4*suup
uup=diss*(uupm+dt4*suup)-uupi
auup=dt2i*(uup-uupi)

 !Update cc = div((1+h)*(u,v)) and gg = div(D(u,v)/Dt):
cci=cc
ggi=gg
ncc=scc+dt4i*cci
ngg=sgg+dt4i*ggi
scc=ncc+scc          !2*N_tilde_delta
sgg=ngg+sgg          !2*N_tilde_gamma
call sistep(scc,sgg,auum,auup)
 !On return, scc = cc^n + cc^{n+1} & sgg = gg^n + gg^{n+1}
cc=scc-cci
gg=sgg-ggi

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call inversion

   !Calculate the source terms (sqs,sqd) for PV (qs,qd) as well as those
   !(scc,sgg) for momentum divergence and acceleration divergence (cc,gg):
  call source(sqs,sqd,scc,sgg,suum,suup)

   !Interpolate gridded velocity (uu,vv) at contour nodes as (uq,vq):
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqi(i)+dt2*uq(i)
    yy=yqi(i)+dt2*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=min(ymax,max(ymin,yy))
  enddo
   !Now (xq,yq) contain a new guess for x^{n+1}.

   !Update PV fields:
  qs=qsi+dt2*sqs
  do j=0,ny
    qd(j,:)=diss*(qdm(j,:)+dt4*sqd(j,:))-qdi(j,:)
  enddo

   !Update zonal velocity at y boundaries:
  uum=diss*(uumm+dt4*suum)-uumi
  auum=dt2i*(uum-uumi)
  uup=diss*(uupm+dt4*suup)-uupi
  auup=dt2i*(uup-uupi)

   !Update cc = div((1+h)*(u,v)) and gg = div(D(u,v)/Dt):
  scc=ncc+scc          !2*N_tilde_delta
  sgg=ngg+sgg          !2*N_tilde_gamma
  call sistep(scc,sgg,auum,auup)
   !On return, scc = cc^n + cc^{n+1} & sgg = gg^n + gg^{n+1}
  cc=scc-cci
  gg=sgg-ggi
enddo

 !Advance time:
t=t+dt

return
end subroutine advance

!=======================================================================

subroutine inversion

! Finds the gridded dimensionless height anomaly (hh), velocity
! field (uu,vv), relative vorticity (zz) and divergence (dd) from
! the PV contours (qc) and auxillary PV fields (qs,qd), the momentum
! divergence (cc), the acceleration divergence (gg), and the zonal 
! velocities (uum & uup) at y = y_min & y_max.

! hh, uu, vv, zz and dd are all returned in physical space.
! A physical space copy of cc (in ccp) is also returned, together
! with 1 + hh (hht) and the gradient of hh (hx,hy).

implicit none

 !Local variables:
double precision,parameter:: utol=1.e-10 !error tolerance for (uu,vv)
double precision:: wkd(0:ny,0:nxm1),hti(0:ny,0:nxm1)
double precision:: wku(0:ny,0:nxm1),wkv(0:ny,0:nxm1)
double precision:: uub(0:nxm1)
double precision:: uerr

!---------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)

 !Convert qc to semi-spectral space:
call forfft(nyp1,nx,qc,xtrig,xfactors)

 !Combine fields to update qq with full (semi-spectral) field,
 !qq = F[qs-qc]+qc+qd, where F is a low pass filter:
qq=qs-qc
call lowpass(qq,1)
qq=qq+qc+qd

 !Convert qq to physical space:
call revfft(nyp1,nx,qq,xtrig,xfactors)

!---------------------------------------------------------------
 !Obtain the dimensionless height anomaly (hh) from the PV, the
 !acceleration divergence (gg), and the zonal velocities
 !(uum & uup) at y = y_min and y_max:
call height(qq,gg,uum,uup,qoff,hh,zz)
 !The (semi-spectral) relative vorticity (zz) is also returned.

 !Calculate 1 + hh -> hht:
hht=one+hh
 !Calculate the gradient of hh -> (hx,hy):
wkd=hh
call forfft(nyp1,nx,wkd,xtrig,xfactors)
call gradient(wkd,hx,hy)

!---------------------------------------------------------------
 !Iterate to find the divergence (dd) and velocity field (uu,vv)
 !from the momentum divergence (cc), relative vorticity (zz),
 !and the boundary zonal velocities (uum & uup):

 !Obtain a physical space copy of cc = div((1+h)(u,v)) -> ccp:
ccp=cc
call revfft(nyp1,nx,ccp,xtrig,xfactors)

 !Calculate and de-alias 1/(1 + hh) -> hti:
hti=one/hht
call dealias(hti,0)

 !Start with the geostrophic velocity as a guess (geo = c^2/f):
uu=-geo*hy
 !Insert boundary values from reverse FFTs of uum & uup:
uub=uum
call revfft(1,nx,uub,xtrig,xfactors)
uu(0 ,:)=uub
uub=uup
call revfft(1,nx,uub,xtrig,xfactors)
uu(ny,:)=uub
 !y velocity component:
vv=geo*hx

 !Begin iteration to find dd, uu & vv:
uerr=one
do while (uerr .gt. utol)
   !Get new estimate for dd from (cc - (u,v)*grad(h))/(1+h):
  wkd=uu*hx+vv*hy
  call dealias(wkd,0)
   !Now wkd = (u,v)*grad(h)
  dd=(ccp-wkd)*hti
   !hti = 1/(1 + h) from above
  call forfft(nyp1,nx,dd,xtrig,xfactors)
   !Now dd, like zz, is in semi-spectral space for use below
  call dealias(dd,1)

   !Get new estimate for the velocity field -> (wku,wkv):
  call velocity(dd,zz,uum,uup,wku,wkv)

   !Compute maximum error:
  wkd=(wku-uu)**2+(wkv-vv)**2
  uerr=sqrt(maxval(wkd))

   !Store updated velocity field:
  uu=wku
  vv=wkv
enddo

 !FFT divergence and vorticity back to physical space:
call revfft(nyp1,nx,dd,xtrig,xfactors)
call revfft(nyp1,nx,zz,xtrig,xfactors)

return
end subroutine inversion

!=======================================================================

subroutine source(sqs,sqd,scc,sgg,suum,suup)

! Gets the source terms (sqs,sqd) for the PV field (qs,qd), as well as
! (scc,sgg) for momentum divergence and acceleration divergence (cc,gg),
! and (suum,suup) for the boundary zonal velocities (uum,uup) ---
! all in semi-spectral space.  Note that (scc,sgg) only include the
! nonlinear terms for a semi-implicit treatment, closely analogous
! to that described in the appendix of Mohebalhojeh & Dritschel (2004).

! The spectral fields cc, gg, qd and qs are all de-aliased.  Note,
! hh, uu, vv, dd & zz, obtained by calling subroutine inversion
! before calling this routine, are all de-aliased, but are in
! physical space.

! A physical space copy of cc (in ccp) is also available, together
! with 1 + hh (hht) and the gradient of hh (hx,hy).
  
implicit none

 !Passed variables:
double precision:: sqs(0:ny,0:nxm1),sqd(0:ny,0:nxm1)
double precision:: scc(0:ny,0:nxm1),sgg(0:ny,0:nxm1)
double precision:: suum(0:nxm1),suup(0:nxm1)

 !x derivative of velocity field (u_x,v_x):
double precision:: ux(0:ny,0:nxm1),vx(0:ny,0:nxm1)
 !y derivative of velocity field (u_y,v_y):
double precision:: uy(0:ny,0:nxm1),vy(0:ny,0:nxm1)

 !Work arrays:
double precision:: wkp(0:ny,0:nxm1),wkq(0:ny,0:nxm1)
double precision:: wka(0:ny,0:nxm1),wkb(0:ny,0:nxm1)

 !Other variables:
integer:: m,l

!---------------------------------------------------------------
 !qs source --- only nonlinear advection term is needed:
call gradient(qs,wkp,wkq)
sqs=-uu*wkp-vv*wkq
 !Convert to spectral space:
call forfft(nyp1,nx,sqs,xtrig,xfactors)
 !Apply de-aliasing filter:
call dealias(sqs,1)

!---------------------------------------------------------------
 !qd source --- only nonlinear advection term is needed:
call gradient(qd,wkp,wkq)
sqd=-uu*wkp-vv*wkq
 !Convert to spectral space:
call forfft(nyp1,nx,sqd,xtrig,xfactors)
 !Apply de-aliasing filter:
call dealias(sqd,1)

!---------------------------------------------------------------
 !Nonlinear part of cc = div((1+h)*(u,v)) source:

 !Compute x derivative of velocity field -> (ux,vx):
wkp=uu
call forfft(nyp1,nx,wkp,xtrig,xfactors)
call xderiv(wkp,ux)
call revfft(nyp1,nx,ux,xtrig,xfactors)
wkp=vv
call forfft(nyp1,nx,wkp,xtrig,xfactors)
call xderiv(wkp,vx)
call revfft(nyp1,nx,vx,xtrig,xfactors)

 !Compute y derivative of velocity field -> (uy,vy):
uy=vx-zz
vy=dd-ux

 !Compute gradient of divergence -> (wka,wkb):
wkp=dd
call forfft(nyp1,nx,wkp,xtrig,xfactors)
call gradient(wkp,wka,wkb)

 !Compute 2*J(u,v) - (u,v)*grad(delta) -> wkp:
wkp=two*(ux*vy-uy*vx)-uu*wka-vv*wkb
 !De-alias:
call dealias(wkp,0)

 !Compute grad(delta_tilde) -> (wka,wkb):
call gradient(cc,wka,wkb)

 !Compute gamma in physical space -> wkq:
wkq=gg
call revfft(nyp1,nx,wkq,xtrig,xfactors)

 !Form nonlinear part of cc = div((1+h)*(u,v)) source, apart from
 !terms coming from grad(h)*[(u_t,v_t) + delta*(u,v)]:
scc=hh*wkq+hht*wkp-uu*wka-vv*wkb-two*dd*ccp
 !All work arrays are now free to use. scc is de-aliased below.

 !Compute u_t + delta*u -> wkp:
wkp=uu*vy-vv*uy
call dealias(wkp,0)
wkp=wkp+cof*vv-csq*hx

 !Compute v_t + delta*v -> wkq:
wkq=vv*ux-uu*vx
call dealias(wkq,0)
wkq=wkq-cof*uu-csq*hy
 !Ensure zero boundary values:
wkq(0 ,:)=zero
wkq(ny,:)=zero

 !Complete definition of nonlinear source term:
scc=scc+hx*wkp+hy*wkq
 !Convert to semi-spectral space and de-alias:
call forfft(nyp1,nx,scc,xtrig,xfactors)
call dealias(scc,1)

!---------------------------------------------------------------
 !Nonlinear part of gg = div(Du/Dt,Dv/Dt) source:
wkq=zz-cof*hh !This is the linearised PV, q_l
wkp=wkq
call forfft(nyp1,nx,wkp,xtrig,xfactors)
call gradient(wkp,wka,wkb)
 !Here grad(q_l) = (wka,wkb)

 !Complete definition of nonlinear source term:
sgg=-cof*(dd*wkq+uu*wka+vv*wkb)
 !Convert to semi-spectral space and de-alias:
call forfft(nyp1,nx,sgg,xtrig,xfactors)
call dealias(sgg,1)

!---------------------------------------------------------------
 !Source terms for the boundary zonal velocities:
suum=-csq*hx(0 ,:)-uu(0 ,:)*ux(0 ,:)
call forfft(1,nx,suum,xtrig,xfactors)
suum=filt*suum

suup=-csq*hx(ny,:)-uu(ny,:)*ux(ny,:)
call forfft(1,nx,suup,xtrig,xfactors)
suup=filt*suup

return
end subroutine source

!=======================================================================

subroutine diagnose

! Computes the twist parameter, the time integral of |zeta|_max, and
! various quantities at every time step to monitor the flow evolution.

implicit none

 !Local variables:
double precision:: wka(0:ny,0:nxm1)
double precision:: ro,fr,hmin,hmax

!----------------------------------------------------------------------
 !Compute Rossby number:
ro=maxval(abs(zz))/cof

 !Compute Froude number:
wka=(uu**2+vv**2)/hht
 !Above, hht = 1 + hh.
fr=sqrt(maxval(wka))/cgw

 !Compute h_min and h_max:
hmin=minval(hh)
hmax=maxval(hh)

 !Write data:
write(16,'(1x,f12.5,4(1x,f12.8))') t,ro,fr,hmin,hmax

 !Increment the integral of |zeta|_max:
twist=twist+dt*maxval(abs(zz))

return
end subroutine diagnose

!=======================================================================

subroutine savegrid(igrids)

! Saves various fields and energy components at the desired save time

implicit none

 !Passed variable:
integer:: igrids

 !Local variables:
double precision:: wka(0:ny,0:nxm1)
double precision:: ekin,epot,etot

!---------------------------------------------------------------
 !Compute energy components and total:
wka=hht*(uu**2+vv**2)
call average(wka,ekin)
ekin=f12*domarea*ekin
wka=hh**2
call average(wka,epot)
epot=f12*domarea*csq*epot
etot=ekin+epot

 !Write energies to ecomp.asc:
write(15,'(f13.6,3(1x,f16.9))') t,ekin,epot,etot
write(*,'(a,f13.6,a,f13.6)') ' t = ',t,'  E_tot = ',etot

!---------------------------------------------------------------
 !Write various gridded fields to direct access files:
write(31,rec=igrids) t,qq
write(32,rec=igrids) t,dd
wka=gg
call revfft(nyp1,nx,wka,xtrig,xfactors)
write(33,rec=igrids) t,wka
write(34,rec=igrids) t,hh
write(35,rec=igrids) t,zz
write(36,rec=igrids) t,uu
write(37,rec=igrids) t,vv

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont(irec)

! Saves PV contours for post-processing and imaging

implicit none

 !Passed variable:
integer:: irec

 !Local variables:
double precision:: wka(0:ny,0:nxm1)
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') irec

 !Write contours to the contours subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,dq

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
call con2grid(qc)
wka=qq-qc
write(83,rec=irec) t,wka

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
