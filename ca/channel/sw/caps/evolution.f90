module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

 !PV field after conversion from contours and residual as in CLAM:
double precision:: qc(0:ny,0:nxm1),qd(0:ny,0:nxm1)

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
integer:: ireg
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
  call savedata(ggen)
   
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
call savedata(ggen)

return
end subroutine advect

!=======================================================================

subroutine init

! Initialises residual PV for normal time integration following 
! contour regeneration

implicit none

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
qc=qd !Re-use qc as a work array:
call lowpass(qc)
 !qc now contains F[qs-qc] in semi-spectral space.
qd=qd-qc

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current (semi-spectral) PV field is stored in qs, and the
! residual PV needed in congen.f90 is stored in qr.

implicit none

!-----------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)

 !Convert qc to semi-spectral space:
call forfft(nyp1,nx,qc,xtrig,xfactors)

 !Define spectral PV (qs) and PV residual (qr):
qr=qs-qc
call lowpass(qr)
qr=qr+qd           !qr = F(qs-qc) + qd
qs=qr+qc           !qs = F(qs-qc) + qd + qc

 !Convert the PV residual to physical space for contouring purposes:
call revfft(nyp1,nx,qr,xtrig,xfactors)
 !Note: we are leaving this module next and qd will be redefined
 !upon re-entry in subroutine init.

return
end subroutine prepare

!=======================================================================

subroutine savedata(ggen)

! Saves data periodically

implicit none

 !Passed variable:
logical:: ggen

 !Local variables:
integer:: itime,jtime

!-----------------------------------------------------
itime=nint(t/dt)
jtime=itime/ngsave

if (ngsave*jtime .eq. itime) then
  call inversion
  call savegrid(jtime+1)
  ggen=.false.
else
  ggen=.true.
endif
 !ggen is used to work out if calling inversion is needed in advance

jtime=itime/ncsave
if (ncsave*jtime .eq. itime) call savecont(jtime+1)

return
end subroutine savedata

!=======================================================================

subroutine advance(ggen)

! Advances PV from time t to t+dt by a combination of contour
! advection (for PV contours) and the pseudo-semi-spectral method 
! for all remaining fields, namely qs, qd, vs, bs, uum & uup.

! Uses an iterative implicit (trapezoidal) method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}+F^n)/2] + N[(F^{n+1}+F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms (if present), and N refers to the
! nonlinear source terms.  We start with a guess for F^{n+1}
! in N and iterate niter times (see parameter statement below).

implicit none

 !Passed variable:
logical:: ggen

 !Local variables:
integer,parameter:: niter=2

 !Semi-spectral fields needed in PV time stepping:
double precision:: qsi(0:ny,0:nxm1),sqs(0:ny,0:nxm1)
double precision:: qdi(0:ny,0:nxm1),qdm(0:ny,0:nxm1),sqd(0:ny,0:nxm1)

 !Semi-spectral fields needed in v' & Dv/Dt (= b) time stepping:
double precision:: vsi(0:ny,0:nxm1),svs(0:ny,0:nxm1),fvs(0:ny,0:nxm1)
double precision:: bsi(0:ny,0:nxm1),sbs(0:ny,0:nxm1),fbs(0:ny,0:nxm1)

 !Spectral fields needed in boundary zonal velocity time stepping:
double precision:: uumi(0:nxm1),suum(0:nxm1),fuum(0:nxm1),duum(0:nxm1)
double precision:: uupi(0:nxm1),suup(0:nxm1),fuup(0:nxm1),duup(0:nxm1)

 !Contour positions and velocities needed in time stepping:
double precision:: xqi(nptq),yqi(nptq),uq(nptq),vq(nptq)

 !Other local quantities:
double precision:: xx,yy
integer:: i,j,iter

!-------------------------------------------------------------------
 !Compute height and velocity at current time level, say t=t^n:
if (ggen) call inversion
 !If ggen is false, inversion was called previously at this time level.
 !Note, inversion returns qc in semi-spectral space, needed just below.

 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*(qs-qc) + qc + qd
 !            and qd = (1-F)*(qs-qc)
 !Here F is a low pass filter (see spectral.f90)
qsi=qs-qc
call lowpass(qsi)
qdi=qsi+qd          !qdi = F*(qs-qc) + qd
qs=qdi+qc           !qs now reset; qdi above is now the new qs-qc
qd=qdi
call lowpass(qdi)   !This is F(qs-qc)
qd=qd-qdi           !qd now reset (high pass filter of qs-qc)

 !Compute twist parameter and save various diagnostics each time step:
call diagnose

!------------------------------------------------------------------
 !Start with a guess for F^{n+1} for all contours and fields:

 !Contours:
call velint(uu,vv,uq,vq)
 !Here (uq,vq) is the contour velocity at time level n, i.e. u(x^n,t^n)
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

 !Calculate the source terms for the PV fields, the meridional velocity,
 !the meridional acceleration and the boundary zonal velocity:
call source(sqs,sqd,svs,sbs,suum,suup)

 !Update PV fields:
qsi=qs+dt2*sqs
qs=qs+dt*sqs
qdi=qd
qdm=qd+dt4*sqd
do j=0,ny
  qd(j,:)=diss*(qdm(j,:)+dt4*sqd(j,:))-qdi(j,:)
enddo
 !diss is a hyperdiffusion operator (see spectral.f90)

 !Update zonal velocity at y boundaries (below dt4i = 4/dt):
uumi=uum
uupi=uup
fuum=suum+dt4i*uumi
fuup=suup+dt4i*uupi
suum=fuum+suum           !2*S^{-}_check
call xderiv1d(suum,duum) !x derivative of the above
suup=fuup+suup           !2*S^{+}_check
call xderiv1d(suup,duup) !x derivative of the above
uum=divu*(rdis*suum+betu*duup-alpu*duum)-uumi
uup=divu*(rdis*suup+alpu*duup-betu*duum)-uupi

 !Update interior v' & Dv/Dt fields (vs,bs):
vsi=vv-vb
call forfft(nyp1,nx,vsi,xtrig,xfactors)
bsi=bs
fvs=svs+dt4i*vsi
fbs=sbs+dt4i*bsi
svs=fvs+svs          !2*\check{N}_v
sbs=fbs+sbs          !2*\check{N}_b
 !Solve semi-implicit system:
call sistep(svs,sbs)
 !On return, svs = vs^n + vs^{n+1} & sbs = bs^n + bs^{n+1}
vv=svs-vsi
bs=sbs-bsi
 !Create v' + vb in physical space -> vv:
call revfft(nyp1,nx,vv,xtrig,xfactors)
vv=vv+vb

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call inversion

   !Calculate the source terms for the PV fields, the meridional velocity,
   !the meridional acceleration and the boundary zonal velocity:
  call source(sqs,sqd,svs,sbs,suum,suup)

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
  suum=fuum+suum           !2*S^{-}_check
  call xderiv1d(suum,duum) !x derivative of the above
  suup=fuup+suup           !2*S^{+}_check
  call xderiv1d(suup,duup) !x derivative of the above
  uum=divu*(rdis*suum+betu*duup-alpu*duum)-uumi
  uup=divu*(rdis*suup+alpu*duup-betu*duum)-uupi

   !Update interior v' & Dv/Dt fields (vs,bs):
  svs=fvs+svs          !2*\check{N}_v
  sbs=fbs+sbs          !2*\check{N}_b
   !Solve semi-implicit system:
  call sistep(svs,sbs)
   !On return, svs = vs^n + vs^{n+1} & sbs = bs^n + bs^{n+1}
  vv=svs-vsi
  bs=sbs-bsi
   !Create v' + vb in physical space -> vv:
  call revfft(nyp1,nx,vv,xtrig,xfactors)
  vv=vv+vb
enddo

 !Advance time:
t=t+dt

return
end subroutine advance

!====================================================================
subroutine sistep(svs,sbs)

! Solves the semi-implicit system for v' & Dv/Dt (= b). Returns the
! sums of these quantities at time levels n and n+1, from which the
! values at level n+1 are easily computed in the calling routine.

! Here svs = 2*\check{N}_v and sbs = 2*\check{N}_b, (twice) the
! the source terms needed to find 2*\bar{v'} = v'^n + v'^{n+1}
! and 2*\bar{b} = b^n + b^{n+1}.

! All quantities are in semi-spectral space.

 !Declarations:
implicit none

 !Passed variables:
double precision:: svs(0:ny,0:nxm1),sbs(0:ny,0:nxm1)

 !Local variables:
double precision:: wks(0:ny,0:nxm1)
integer:: j,m

!---------------------------------------------------------------------
 !Form r.h.s. of the equation for v'^n + v'^{n+1}:
do j=1,nym1
  sbs(j,:)=-csqi*(sbs(j,:)+rdis*svs(j,:))
enddo
 !sbs = 2*\check{S}_v in the notes.

 !Solve the tri-diagonal system for v'^n + v'^{n+1} -> wks:
 !(Below ap = 1/dy^2)
do m=0,nxm1
  wks(1,m)=sbs(1,m)*htdt(1,m)
  do j=2,nym1
    wks(j,m)=(sbs(j,m)-ap*wks(j-1,m))*htdt(j,m)
  enddo
  do j=nym2,1,-1
    wks(j,m)=etdt(j,m)*wks(j+1,m)+wks(j,m)
  enddo
enddo
 !Note, the boundary values for j = 0 and ny are zero (not needed)

 !Obtain b^n + b^{n+1} -> sbs, and store v'^n + v'^{n+1} -> svs:
do j=1,nym1
  sbs(j,:)=rdis*wks(j,:)-svs(j,:)
  svs(j,:)=wks(j,:)
enddo

return
end subroutine sistep

!=======================================================================

subroutine inversion

! Finds the gridded dimensionless height anomaly (hh), horizontal
! velocity component (uu) and relative vorticity (zz) from the PV
! contours (qc) and auxillary PV fields (qs,qd), the meridional
! velocity component (vv), the meridional acceleration (bs) and
! the zonal velocities (uum & uup) at y = y_min & y_max.

! hh, uu and zz are all returned in physical space.

implicit none

!---------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)

 !Convert qc to semi-spectral space:
call forfft(nyp1,nx,qc,xtrig,xfactors)

 !Combine fields to update qq with full (semi-spectral) field,
 !qq = F[qs-qc]+qc+qd, where F is a low pass filter:
qq=qs-qc
call lowpass(qq)
qq=qq+qc+qd

 !Convert qq to physical space:
call revfft(nyp1,nx,qq,xtrig,xfactors)

!---------------------------------------------------------------
 !Obtain the dimensionless height anomaly (hh), horizontal
 !velocity (uu).relative vorticity (zz), linearised PV anomaly
 !zeta - f*h (qt) and balanced y velocity (vb) by inversion
 !(see spectral.f90 for full details):
call invert(qq,vv,bs,uum,uup,qoff,hh,uu,zz,qt,vb)

return
end subroutine inversion

!=======================================================================

subroutine source(sqs,sqd,svs,sbs,suum,suup)

! Gets the source terms (sqs,sqd) for the PV fields (qs,qd), as well
! as (svs,sbs) for v' & Dv/Dt (= b) and (suum,suup) for the boundary
! zonal velocities (uum,uup) --- all in semi-spectral space.
! Note that svs ... suup only include the nonlinear terms for a
! semi-implicit treatment.

! In the channel SW notes (nse.pdf), svs & sbs correspond to 
! \check{N}_v & \check{N}_b.

! The spectral fields qs, qd, bs, uum & uup are all de-aliased.
! Furthermore, the physical fields hh, uu, vv, zz, qt & vb obtained
! before calling this routine, are all de-aliased.
  
implicit none

 !Passed variables:
double precision:: sqs(0:ny,0:nxm1),sqd(0:ny,0:nxm1)
double precision:: svs(0:ny,0:nxm1),sbs(0:ny,0:nxm1)
double precision:: suum(0:nxm1),suup(0:nxm1)

 !x derivative of velocity field (u_x,v_x):
double precision:: ux(0:ny,0:nxm1),vx(0:ny,0:nxm1)
 !y derivative of velocity field (u_y,v_y):
double precision:: uy(0:ny,0:nxm1),vy(0:ny,0:nxm1)

 !x derivative of acceleration field (a_x,b_x):
double precision:: ax(0:ny,0:nxm1),bx(0:ny,0:nxm1)
 !y derivative of acceleration field (a_y,a_y):
double precision:: ay(0:ny,0:nxm1),by(0:ny,0:nxm1)

 !Geopotential anomaly:
double precision:: phi(0:ny,0:nxm1)

 !Work arrays:
double precision:: wka(0:ny,0:nxm1),wkb(0:ny,0:nxm1)
double precision:: bwk(0:nxm1)

 !Loop indices:
integer:: j,m

!---------------------------------------------------------------
 !qs source --- only nonlinear advection term is needed:
call gradient(qs,wka,wkb)
sqs=-uu*wka-vv*wkb
 !Convert to spectral space:
call forfft(nyp1,nx,sqs,xtrig,xfactors)
 !Apply de-aliasing filter:
call dealias(sqs,1)

!---------------------------------------------------------------
 !qd source --- only nonlinear advection term is needed:
call gradient(qd,wka,wkb)
sqd=-uu*wka-vv*wkb
 !Convert to spectral space:
call forfft(nyp1,nx,sqd,xtrig,xfactors)
 !Apply de-aliasing filter:
call dealias(sqd,1)

!---------------------------------------------------------------
 !Nonlinear part of v' source --- c^2 G^{-1}{w_x} - (u,v)*grad(v)
 !where w = div(qt(u,v)) and qt = zeta - f*h:
wka=qt*uu
call forfft(nyp1,nx,wka,xtrig,xfactors)
call xderiv(wka,wkb) !(qt*u)_x
wka=qt*vv
call forfft(nyp1,nx,wka,xtrig,xfactors)
call yderiv(wka,phi) !(qt*v)_y
wkb=wkb+phi          !w = div(qt(u,v))
call xderiv(wkb,phi) !w_x
do m=0,nxm1
  phi(0,m)=zero
  phi(1,m)=phi(1,m)*htdv(1,m)
  do j=2,nym1
    phi(j,m)=(phi(j,m)-ap*phi(j-1,m))*htdv(j,m)
  enddo
  phi(ny,m)=zero
  do j=nym2,1,-1
    phi(j,m)=etdv(j,m)*phi(j+1,m)+phi(j,m)
  enddo
enddo
 !phi now contains c^2 G^{-1}{w_x} in semi-spectral space

 !Calculate grad(v):
svs=vv
call forfft(nyp1,nx,svs,xtrig,xfactors)
call gradient(svs,wka,wkb)
svs=uu*wka+vv*wkb  !(u,v)*grad(v)
 !Convert to spectral space:
call forfft(nyp1,nx,svs,xtrig,xfactors)
 !Add c^2 G^{-1}{w} in phi from above:
svs=phi-svs
 !Apply de-aliasing filter:
call dealias(svs,1)

!---------------------------------------------------------------
 !Nonlinear part of b source --- T_x + (phi*v)_yy - f*zeta*v
 !where T = (f/2)*(u^2+v^2) + (phi*u)_y

 !Compute geopotential anomaly -> phi:
phi=csq*hh
 !Compute T -> wkb:
wka=phi*uu
call yderiv(wka,wkb) !Note: y derivative at boundary is set to zero
wkb=f12*cof*(uu**2+vv**2)+wkb
call forfft(nyp1,nx,wkb,xtrig,xfactors)
 !Compute T_x -> sbs (semi-spectral):
call xderiv(wkb,sbs)

 !Add (phi*v)_yy - f*zeta*v after conversion to semi-spectral space:
wkb=phi*vv
 !Use centred differences for the second derivative (ap = 1/dy^2):
do j=1,nym1
  wka(j,:)=ap*(wkb(j+1,:)-two*wkb(j,:)+wkb(j-1,:))-cof*zz(j,:)*vv(j,:)
enddo
call forfft(nyp1,nx,wka,xtrig,xfactors)
sbs=sbs+wka
 !Apply de-aliasing filter:
call dealias(sbs,1)

!---------------------------------------------------------------
 !Source terms for the boundary zonal velocities:
bwk=phi(0,:)+f12*uu(0,:)**2
call forfft(1,nx,bwk,xtrig,xfactors)
bwk=alpu*uum-betu*uup-bwk
call xderiv1d(bwk,suum)
suum=filt*suum

bwk=phi(ny,:)+f12*uu(ny,:)**2
call forfft(1,nx,bwk,xtrig,xfactors)
bwk=betu*uum-alpu*uup-bwk
call xderiv1d(bwk,suup)
suup=filt*suup

return
end subroutine source

!=======================================================================

subroutine diagnose

! Computes the twist parameter, the time integral of |zeta|_max, and
! various quantities at every time step to monitor the flow evolution.

implicit none

 !Local variables:
double precision:: wk(0:ny,0:nxm1)
double precision:: ro,fr,hmin,hmax

!----------------------------------------------------------------------
 !Compute Rossby number:
ro=maxval(abs(zz))/cof

 !Compute Froude number:
wk=(uu**2+vv**2)/(one+hh)
fr=sqrt(maxval(wk))/cgw

 !Compute h_min and h_max:
hmin=minval(hh)
hmax=maxval(hh)

 !Write data to ro-fr-hm.asc in the evolution subdirectory:
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
double precision:: wk(0:ny,0:nxm1)
double precision:: ekin,epot,etot

!---------------------------------------------------------------
 !Compute kinetic energy:
wk=(one+hh)*(uu**2+vv**2)
call average(wk,ekin)
ekin=f12*domarea*ekin

 !Compute potential energy:
wk=hh**2
call average(wk,epot)
epot=f12*domarea*csq*epot

 !Compute total energy:
etot=ekin+epot

 !Write data to ecomp.asc in evolution subdirectory:
write(15,'(f13.6,3(1x,f16.9))') t,ekin,epot,etot
write(*,'(a,f13.6,a,f13.6)') ' t = ',t,'  E_tot = ',etot

!---------------------------------------------------------------
 !Write various gridded fields to direct access files:
write(32,rec=igrids) t,qq
wk=bs
call revfft(nyp1,nx,wk,xtrig,xfactors)
write(33,rec=igrids) t,wk
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
double precision:: wk(0:ny,0:nxm1)
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') irec

 !Write contours to the contours subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,dq

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
call con2grid(qc)
wk=qq-qc
write(83,rec=irec) t,wk

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
