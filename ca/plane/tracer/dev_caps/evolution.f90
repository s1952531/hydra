module evolution

! Module contains subroutines to evolve PV contours and all fields 
! according to the algorithm detailed in caps.f90.

use common

implicit none

 !Velocity field:
double precision:: uu(ny,nx),vv(ny,nx)

 !Spectral fields (note array order):
double precision:: pp(nx,ny),qq(nx,ny),qc(nx,ny),qd(nx,ny)
double precision:: qspre(nx,ny)
double precision:: emq(nx,ny),epq(nx,ny)

 !Energies:
double precision:: eneupre

 !Logicals to indicate presence of contours or data saves:
logical:: contex,gsave,csave

!Internal subroutine definitions (inherit global variables):
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

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

 !Used for regularising contours:
twist=zero

 !Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tfin)

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
    write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

    twist=twist-twistmax
  endif

   !Advect flow from time t to t + dt:
  call advance
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

return
end subroutine

!=======================================================================

subroutine init

! Initialises quantities needed for normal time integration following 
! contour regeneration

implicit none

 !Local variables:
double precision:: qa(ny,nx)
integer:: kx,ky,nc,nptc,i,j

!------------------------------------------------------------------
 !Logical to indicate presence of contours:
contex=(nptq .gt. 0)

 !Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

!-------------------------------------------------------------
 !Logicals used for saving gridded fields and contours:
gsave=.false.
csave=.false.

 !Convert PV contours to gridded values (qa):
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Define residual PV qd = (1-F)[qs-qc], and copy current gridded 
 !field (qs) for use in time interpolation:
do ky=1,ny
  do kx=1,nx
    qd(kx,ky)=fhi(kx,ky)*(qs(kx,ky)-qc(kx,ky))
    qspre(kx,ky)=qs(kx,ky) 
  enddo
enddo
 !Here fhi = 1-F is a high-pass spectral filter

return
end subroutine

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current PV anomaly field is stored in qs, the residual PV needed 
! in congen.f90 is stored in qr, and (if present) tracer contours are
! stored in tracer.bin

implicit none

 !Local variables:
double precision:: qa(ny,nx)
integer:: kx,ky,i,j

!-----------------------------------------------------------------
 !Convert PV contours to gridded values (qa):
call con2grid(qa)

 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Put current PV anomaly (spectral in qs) and define residual (qd):
do ky=1,ny
  do kx=1,nx
    qs(kx,ky)=flo(kx,ky)*(qs(kx,ky)-qc(kx,ky))+qc(kx,ky)+qd(kx,ky)
    qd(kx,ky)=qs(kx,ky)-qc(kx,ky)
  enddo
enddo

 !Convert qd to physical space as qr (used in recontouring):
call spctop(nx,ny,qd,qr,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!=======================================================================

subroutine advance

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method (for all
! gridded fields, i.e. qs & qd).

! *** Uses a 4th-order Runge-Kutta method ***

implicit none

 !Local variables:

 !Spectral fields needed in Runge-Kutta time stepping (note array order):
double precision:: qsi(nx,ny),qsf(nx,ny),sqs(nx,ny)
double precision:: qdi(nx,ny),qdf(nx,ny),sqd(nx,ny)
 !Contour positions needed in Runge-Kutta time stepping:
double precision:: xqi(nptq),yqi(nptq),xqf(nptq),yqf(nptq)
 !Contour velocities:
double precision:: uq(nptq),vq(nptq)
 !Other local quantities:
double precision:: xx,yy
integer:: kx,ky,i

!-------------------------------------------------------------------
 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*(qs-qc) + qc + qd
 !            and qd = (1-F)*(qs-qc)
 !after qs is reset; here F is a low pass filter (see spectral.f90)
do ky=1,ny
  do kx=1,nx
    qs(kx,ky)=flo(kx,ky)*(qs(kx,ky)-qc(kx,ky))+qc(kx,ky)+qd(kx,ky)
    qd(kx,ky)=fhi(kx,ky)*(qs(kx,ky)-qc(kx,ky))
  enddo
enddo

!------------------------------------------------------------------
 !RK4 predictor step to time t0 + dt/2:

 !Invert PV and compute velocity:
call inversion

 !Possibly save data (gsave & csave set by adapt in the previous time step):
if (gsave) call savegrid
if (csave) call savecont

 !Adjust timestep (dt) on maximum vorticity magnitude or CFL:
call adapt

 !Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,0)

if (contex) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xqi(i)=xq(i)
    yqi(i)=yq(i)
    xx=xqi(i)+dt2*uq(i)
    yy=yqi(i)+dt2*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
    xx=xqi(i)+dt6*uq(i)
    yy=yqi(i)+dt6*vq(i)
    xqf(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yqf(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

do ky=1,ny
  do kx=1,nx
    qsi(kx,ky)=qs(kx,ky)
    qs(kx,ky)=qsi(kx,ky)+dt2*sqs(kx,ky)
    qsf(kx,ky)=qsi(kx,ky)+dt6*sqs(kx,ky)
    qdi(kx,ky)=qd(kx,ky)
    qd(kx,ky)=emq(kx,ky)*(qdi(kx,ky)+dt2*sqd(kx,ky))
    qdf(kx,ky)=qdi(kx,ky)+dt6*sqd(kx,ky)
  enddo
enddo

!------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt/2:
t=t+dt2

 !Invert PV and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,1)

if (contex) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqi(i)+dt2*uq(i)
    yy=yqi(i)+dt2*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
    xx=xqf(i)+dt3*uq(i)
    yy=yqf(i)+dt3*vq(i)
    xqf(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yqf(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

do ky=1,ny
  do kx=1,nx
    qs(kx,ky)=qsi(kx,ky)+dt2*sqs(kx,ky)
    qsf(kx,ky)=qsf(kx,ky)+dt3*sqs(kx,ky)
    qd(kx,ky)=emq(kx,ky)*(qdi(kx,ky)+dt2*sqd(kx,ky))
    qdf(kx,ky)=qdf(kx,ky)+dt3*sqd(kx,ky)
  enddo
enddo

!------------------------------------------------------------------
 !RK4 predictor step at time t0 + dt:

 !Invert PV and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,1)

if (contex) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqi(i)+dt*uq(i)
    yy=yqi(i)+dt*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
    xx=xqf(i)+dt3*uq(i)
    yy=yqf(i)+dt3*vq(i)
    xqf(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yqf(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

do ky=1,ny
  do kx=1,nx
    qs(kx,ky)=qsi(kx,ky)+dt*sqs(kx,ky)
    qsf(kx,ky)=qsf(kx,ky)+dt3*sqs(kx,ky)
    emq(kx,ky)=emq(kx,ky)**2
    qd(kx,ky)=emq(kx,ky)*(qdi(kx,ky)+dt*sqd(kx,ky))
    qdf(kx,ky)=qdf(kx,ky)+dt3*sqd(kx,ky)
  enddo
enddo

!------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt:
t=t+dt2

 !Invert PV and compute velocity:
call inversion

 !Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,2)

if (contex) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqf(i)+dt6*uq(i)
    yy=yqf(i)+dt6*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

do ky=1,ny
  do kx=1,nx
    qs(kx,ky)=qsf(kx,ky)+dt6*sqs(kx,ky)
    qd(kx,ky)=emq(kx,ky)*(qdf(kx,ky)+dt6*sqd(kx,ky))
  enddo
enddo

 !Add any stochastic forcing (if present) to qd:
if (stoch) call pvforce

return
end subroutine

!=======================================================================

subroutine pvforce

! Adds stochastic forcing to PV in the array qd at the end of a time step

implicit none

 !Local variables:
double precision:: wkp(ny,nx) !Physical work array
double precision:: wks(nx,ny) !Spectral work array
double precision:: dnx,dny,svor,xx,pxc,px,yy,pyc,py
double precision:: theta,cth,sth
integer:: nvor,k,ix,iy,ix0,ix1,iy0,iy1,kx,ky

!-----------------------------------------------------------------
 !Initialise array to uptake added PV:
do ix=1,nx
  do iy=1,ny
    wkp(iy,ix)=zero
  enddo
enddo

 !Add point vortices here to force the flow (they are converted to 
 !gridded PV values and added to qd):
if (ivor .eq. 1) then
   !Add nvor vortices (as +/- arbitrarily separated pairs):
  nvor=2*nint(f12*(dnvor*t-totnvor))

   !The vortices are place randomly in the domain:
  dnx=dble(nx)
  dny=dble(ny)
  svor=-vorvor
  do k=1,nvor
    svor=-svor

     !Divide each vortex's circulation among the corners
     !of the local grid box (inverse bi-linear interpolation);
    xx=dnx*rand(0)
    ix0=1+int(xx)
    pxc=dble(ix0)-xx
    px=one-pxc
    ix1=ixp(ix0)

    yy=dny*rand(0)
    iy0=1+int(yy)
    pyc=dble(iy0)-yy
    py=one-pyc
    iy1=iyp(iy0)

    wkp(iy0,ix0)=wkp(iy0,ix0)+svor*pyc*pxc
    wkp(iy0,ix1)=wkp(iy0,ix1)+svor*pyc*px
    wkp(iy1,ix0)=wkp(iy1,ix0)+svor*py*pxc
    wkp(iy1,ix1)=wkp(iy1,ix1)+svor*py*px
  enddo

else
   !Add nvor dipoles:
  nvor=nint(dnvor*t-totnvor)

   !The dipoles are place randomly in the domain with random
   !orientation:
  dnx=dble(nx)
  dny=dble(ny)
  do k=1,nvor

    theta=twopi*rand(0)
    cth=cos(theta)
    sth=sin(theta)

    xx=dnx*rand(0)
    ix0=1+int(xx)
    pxc=dble(ix0)-xx
    px=one-pxc
    ix1=ixp(ix0)

    yy=dny*rand(0)
    iy0=1+int(yy)
    pyc=dble(iy0)-yy
    py=one-pyc
    iy1=iyp(iy0)

    wkp(iy0,ix0)=wkp(iy0,ix0)-vorvor*(pyc*cth+pxc*sth)
    wkp(iy0,ix1)=wkp(iy0,ix1)+vorvor*(pyc*cth -px*sth)
    wkp(iy1,ix0)=wkp(iy1,ix0)+vorvor*(pxc*sth -py*cth)
    wkp(iy1,ix1)=wkp(iy1,ix1)+vorvor*( px*sth +py*cth)
  enddo

endif

 !Total number of vortices/dipoles added so far:
totnvor=totnvor+dble(nvor)

 !Convert wkp to spectral space as wks:
call ptospc(nx,ny,wkp,wks,xfactors,yfactors,xtrig,ytrig)

 !Add wks to qd in spectral space:
do ky=1,ny
  do kx=1,nx
    qd(kx,ky)=qd(kx,ky)+wks(kx,ky)
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine source(sqs,sqd,lev)

! Gets the source terms (sqs,sqd) for the PV (qs,qd)
! evolution equation (all in spectral space):

implicit none

 !Passed variables:
double precision:: sqs(nx,ny),sqd(nx,ny)
integer:: lev

 !Local variables:
double precision:: qqx(ny,nx),qqy(ny,nx)
double precision:: wkp(ny,nx),x,skx
integer:: ix,iy,kx,ky

!---------------------------------------------------------------
 !qd source:
call gradient(qd,qqx,qqy)
do ix=1,nx
  x=xmin+glx*dble(ix-1)
  skx=exp(-srcalp*t)*sin(kf*x)
  do iy=1,ny
    wkp(iy,ix)=skx-uu(iy,ix)*qqx(iy,ix)-vv(iy,ix)*qqy(iy,ix)
  enddo
enddo
 !Convert to spectral space:
call ptospc(nx,ny,wkp,sqd,xfactors,yfactors,xtrig,ytrig)

!---------------------------------------------------------------
 !qs source - only NL term is needed:
call gradient(qs,qqx,qqy)
do ix=1,nx
  do iy=1,ny
    wkp(iy,ix)=-uu(iy,ix)*qqx(iy,ix)-vv(iy,ix)*qqy(iy,ix)
  enddo
enddo
 !Convert to spectral space:
call ptospc(nx,ny,wkp,sqs,xfactors,yfactors,xtrig,ytrig)

!---------------------------------------------------------------
 !Implement Ekman and/or thermal damping (add to sqd if present):
if (heating) then
  if (friction) then 
   !Use thermal and Ekman damping:
    do ky=1,ny
      do kx=1,nx
        sqd(kx,ky)=sqd(kx,ky)+therm*(pp(kx,ky)-ppeq(kx,ky)) &
                          & -rekman*(qq(kx,ky)+kdsq*pp(kx,ky))
      enddo
    enddo
  else
    do ky=1,ny
      do kx=1,nx
        sqd(kx,ky)=sqd(kx,ky)+therm*(pp(kx,ky)-ppeq(kx,ky))
      enddo
    enddo
  endif
else
   !Only use Ekman damping
  do ky=1,ny
    do kx=1,nx
      sqd(kx,ky)=sqd(kx,ky)-rekman*(qq(kx,ky)+kdsq*pp(kx,ky))
    enddo
  enddo
endif

!----------------------------------------------------------------
if (lev .eq. 0) return

 !Apply exponential integrating factors:
if (lev .eq. 1) then
  do ky=1,ny
    do kx=1,nx
      sqd(kx,ky)= epq(kx,ky)*sqd(kx,ky)
    enddo
  enddo
else
  do ky=1,ny
    do kx=1,nx
      sqd(kx,ky)= epq(kx,ky)**2*sqd(kx,ky)
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine inversion

! Inverts Laplace's operator on PV anomaly (PV - beta*y) to obtain 
! the streamfunction (pp) ***in spectral space*** and the velocity 
! (uu,vv) = (-dpp/dy,dpp/dx) ***in physical space***

implicit none

 !Local variables:
double precision:: qa(ny,nx)
integer:: kx,ky

!------------------------------------------------------------
 !Call con2grid to get updated contour PV (qc):
call con2grid(qa)
 !Convert qa to spectral space as qc (note, qa is modified):
call ptospc(nx,ny,qa,qc,xfactors,yfactors,xtrig,ytrig)

 !Combine fields to update qq with full field,
 !qq = F[qs-qc]+qc+qd, where F is a low pass filter:
do ky=1,ny
  do kx=1,nx
    qq(kx,ky)=flo(kx,ky)*(qs(kx,ky)-qc(kx,ky))+qc(kx,ky)+qd(kx,ky)
  enddo
enddo

 !Invert PV to obtain velocity field: 
!call main_invert(qq,uu,vv,pp)

if (nper .ne. int(time/tper)) then
  theta1=rand(0)*twopi
  theta2=rand(0)*twopi
  nper=nper+1
endif

if (nper .eq. nint(time/tper)) then 
  do ix=1,nx
    do iy=1,ny
      y=ymin+gly*dble(iy-1)
      uu(iy,ix)=uscale*cos(ku*y+theta1)
      vv(iy,ix)=zero
    enddo
  enddo
else 
  do ix=1,nx
    do iy=1,ny
      x=xmin+glx*dble(ix-1)
      uu(iy,ix)=zero
      vv(iy,ix)=uscale*cos(ku*x+theta2)
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine adapt

! Adapts the time step dt to ensure that it is less than or equal to
! the minimum of dtfac/|zeta|_max and C*dx/|u|_max
! where dx is the grid spacing, u is the vector velocity field.
! C = cfl_max is specified below.

implicit none

 !Local variables:
double precision,parameter:: dtfac=pi/10.d0, cflmax=0.7d0
double precision:: zz(ny,nx) !Physical
double precision:: ss(nx,ny) !Spectral
double precision:: umax,zzrms,zzmax,dtacc,tcont,epot,cfl
integer:: kx,ky,ix,iy,itime

!----------------------------------------------------------------------
 !Compute the gridded relative vorticity (zz):
if (btropic) then
   !Here kd = 0:
  do ky=1,ny
    do kx=1,nx
      ss(kx,ky)=qq(kx,ky)
    enddo
  enddo
else
   !Here kd > 0:
  do ky=1,ny
    do kx=1,nx
      ss(kx,ky)=qq(kx,ky)+kdsq*pp(kx,ky)
    enddo
  enddo
endif
 !Above, qq = PV anomaly (q - beta*y) in spectral space.
call spctop(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)

 !Compute accurate advection time step:
umax=small
zzrms=zero
zzmax=small
do ix=1,nx
  do iy=1,ny
    umax=max(umax,uu(iy,ix)**2+vv(iy,ix)**2)
    zzrms=zzrms+zz(iy,ix)**2
    zzmax=max(zzmax,abs(zz(iy,ix)))
  enddo
enddo
umax=sqrt(umax)
zzrms=sqrt(zzrms*dsumi)
dtacc=min(glx*cflmax/umax,dtfac/max(zzmax,srwfm))
 !The restriction on the maximum Rossby wave frequency (srwfm)
 !ensures that the fastest Rossby wave frequency is resolved.

!---------------------------------------------------------------------
 !Choose a new time step, dt:
if (dt .gt. zero) then
  dt=min(dtacc,dtmax)
  if (dt .gt. dtacc) write(*,'(a,f9.5)') 'Warning! dt/dt_acc= ',dt/dtacc
else
   !Limit max timestep to a data save time
  dtmax=min(tgsave,tcsave)
  dt=min(dtacc,dtmax)
endif
 !Fractional time steps used in 4th-order Runge-Kutta time stepping:
dt2=dt*f12
dt3=dt*f13
dt6=dt*f16

 !Increment the integral of max|zz|:
twist=twist+dt*zzmax

!---------------------------------------------------------------------
 !Record various diagnostics to monitor.asc:
cfl=umax*dt/glx
write(17,'(1x,f12.5,1x,f6.4,4(1x,f12.6),1x,f5.3)') & 
     & t,cfl,f12*zzrms**2,zzrms,zzmax,umax,twist

!---------------------------------------------------------------------
 !Set flag to save gridded data every tgsave time units:
itime=int((t+dt)/tgsave)
tgrid=tgsave*dble(itime)+small
 !small = 1.d-12 is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  gsave=.true.

   !Copy current gridded fields for use in time interpolation:
  do ky=1,ny
    do kx=1,nx
      qspre(kx,ky)=qs(kx,ky) 
    enddo
  enddo

   !Compute energy:
  call energy(eneupre)
endif

 !Set flag to save contour data every tcsave time units:
itime=int((t+dt)/tcsave)
tcont=tcsave*dble(itime)
if (t .le. tcont .and. t+dt .gt. tcont) then
   !The save time is between t & t+dt; set flag to save data:
  csave=(sign(1,2*itime-1) .gt. 0)
   !This construction avoids saving the contours at t = 0.
endif

!---------------------------------------------------------------------
 !Define spectral integrating factors used in Runge-Kutta integration:
do ky=1,ny
  do kx=1,nx
    epq(kx,ky)=exp(dt2*qdiss(kx,ky))
    emq(kx,ky)=one/epq(kx,ky)
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine energy(eu)

! This routine computes the total energy from uu, vv & pp (if kd > 0).

implicit none

 !Passed variable:
double precision:: eu

 !Local variables:
double precision:: ss(nx,ny) !Spectral
double precision:: zz(nx,ny) !Physical
double precision:: epot
integer:: ix,iy,kx,ky

!----------------------------------------------------------------------
if (btropic) then
   !kd = 0 and hence there is no potential energy:
  epot=zero
else
   !Get streamfunction pp in physical space as zz:
  do ky=1,ny
    do kx=1,nx
      ss(kx,ky)=pp(kx,ky)
    enddo
  enddo
  call spctop(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)
   !Compute L_2 norm:
  epot=zero
  do ix=1,nx
    do iy=1,ny
      epot=epot+zz(iy,ix)**2
    enddo
  enddo
  epot=f12*garea*kdsq*epot
endif

 !Compute kinetic energy, and add to PE from above:
eu=zero
do ix=1,nx
  do iy=1,ny
    eu=eu+uu(iy,ix)**2+vv(iy,ix)**2
  enddo
enddo
eu=f12*garea*eu+epot

return
end subroutine

!=======================================================================

subroutine savegrid

! Saves PV, energy and various spectra at the desired save time

implicit none

 !Local variables:
double precision:: wka(ny,nx),wkb(ny,nx) !Physical
double precision:: qqs(nx,ny) !Spectral
double precision:: qspec(0:max(nx,ny))
double precision:: qql2,pt,ptc,eneu
double precision:: eneupost,sumqspec
real:: qqr4(ny,nx),tr4
integer:: ix,iy,kx,ky,k

!---------------------------------------------------------------
 !Increment counter for direct file access:
igrids=igrids+1

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Interpolate PV anomaly at save time:
do ky=1,ny
  do kx=1,nx
    qqs(kx,ky)=pt*qspre(kx,ky)+ptc*qq(kx,ky)
  enddo
enddo

!---------------------------------------------------------------
 !Compute energy:
call energy(eneupost)

 !Compute time interpolated energy:
eneu=pt*eneupre+ptc*eneupost

 !Write energy to ene.asc:
write(15,'(f7.2,1x,f16.9)') tgrid,eneu

!---------------------------------------------------------------
 !Compute 1d PV spectrum:
call spec1d(qqs,qspec,0)
sumqspec=zero
do k=1,kmax
  sumqspec=sumqspec+qspec(k)
   !Normalise to take into account uneven sampling of wavenumbers 
   !in each shell [k-1/2,k+1/2]:
  qspec(k)=spmf(k)*qspec(k)
enddo
sumqspec=8.d0*sumqspec*dsumi

!---------------------------------------------------------------
!Write full PV to qq.r4:
call spctop(nx,ny,qqs,wka,xfactors,yfactors,xtrig,ytrig)

tr4=real(tgrid)
do ix=1,nx
  do iy=1,ny
    qqr4(iy,ix)=real(wka(iy,ix)+bety(iy))
  enddo
enddo
write(31,rec=igrids) tr4,qqr4

 !Compute domain integral of (q-beta*y)^2 and write to norms.asc:
call l2norm(wka,qql2)
write(16,'(f7.2,1x,f16.9)') tgrid,qql2

write(*,'(a,f7.2,2(a,f13.6))') ' t = ',tgrid,' <q^2> = ',qql2,' E_tot = ',eneu

 !Write out spectrum to file:
write(51,'(f7.2,2(1x,f16.9),1x,i5)') tgrid,sumqspec,qql2,kmaxred
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(51,'(2(1x,f12.8))') alk(k),log10(qspec(k))
enddo
 !Note: alk(k) = log_10(k)

 !Unset flag for saving data:
gsave=.false.

return
end subroutine

!=======================================================================
      
subroutine savecont

! Saves PV contours for post-processing and imaging

implicit none

 !Local variables:
double precision:: ss(nx,ny)
double precision:: qa(ny,nx)
real:: qdr4(ny,nx),tr4
integer:: irec,kx,ky,ix,iy
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
irec=nint(t/tcsave)
write(pind(1:3),'(i3.3)') irec

tr4=real(t)
 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
do ky=1,ny
  do kx=1,nx
    ss(kx,ky)=qq(kx,ky)-qc(kx,ky)
  enddo
enddo
call spctop(nx,ny,ss,qa,xfactors,yfactors,xtrig,ytrig)
do ix=1,nx
  do iy=1,ny
    qdr4(iy,ix)=real(qa(iy,ix))
  enddo
enddo
write(83,rec=irec) tr4,qdr4

 !Save PV contours if any exist:
if (contex) then
  open(81,file='cont/qqindex'//pind,form='unformatted', &
      & access='direct',status='replace',recl=12*nq)
  write(81,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq)
  close(81)

  open(82,file='cont/qqnodes'//pind,form='unformatted', &
      & access='direct',status='replace',recl=16*nptq)
  write(82,rec=1) xq(1:nptq),yq(1:nptq)
  close(82)
endif

 !Unset flag for saving data:
csave=.false.

return
end subroutine

!=======================================================================

 !Main end module
end module
