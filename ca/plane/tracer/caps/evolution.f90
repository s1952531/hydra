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

subroutine getvel(time)

! Sets the current velocity based on the formula:
! u=   {U*cos(ku*y+theta1(n)     n*T < t < (n+1/2)*T
!      {U*cos(ku*x+theta2(n)     (n+1/2)*T < t < (n+1)*T
! See paper Ahmadi etal (in prep.) for more details.

implicit double precision(a-h,o-z)
implicit integer(i-n)

!------------------------------------------------------------

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
call getvel(t+small)

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
call getvel(t)

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
call getvel(t)

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
call getvel(t)

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
double precision:: wkp(ny,nx),x
integer:: ix,iy,kx,ky

!---------------------------------------------------------------
 !qd source:
call gradient(qd,qqx,qqy)
do ix=1,nx
  do iy=1,ny
    wkp(iy,ix)=-uu(iy,ix)*qqx(iy,ix)-vv(iy,ix)*qqy(iy,ix)
  enddo
enddo
 !Convert to spectral space:
call ptospc(nx,ny,wkp,sqd,xfactors,yfactors,xtrig,ytrig)

!---------------------------------------------------------------
 !qs source - only NL term is needed:
call gradient(qs,qqx,qqy)
do ix=1,nx
  x=xmin+glx*dble(ix-1)
  do iy=1,ny
    wkp(iy,ix)=exp(-srcalp*t)*sin(kf*x)-uu(iy,ix)*qqx(iy,ix)-vv(iy,ix)*qqy(iy,ix)
  enddo
enddo
 !Convert to spectral space:
call ptospc(nx,ny,wkp,sqs,xfactors,yfactors,xtrig,ytrig)

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

subroutine adapt

! Adapts the time step dt to ensure that it is less than or equal to
! the minimum of dtfac/|zeta|_max and C*dx/|u|_max
! where dx is the grid spacing, u is the vector velocity field.
! C = cfl_max is specified below.

implicit none

 !Local variables:
double precision:: zzmax,tcont
integer:: kx,ky,ix,iy,itime

!----------------------------------------------------------------------
 !Increment the integral of max|zz|:
zzmax=ku*uscale
twist=twist+dt*zzmax

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
endif

 !Set flag to save contour data every tcsave time units:
itime=int((t+dt)/tcsave)
tcont=tcsave*dble(itime)
if (t .le. tcont .and. t+dt .gt. tcont) then
   !The save time is between t & t+dt; set flag to save data:
  csave=(sign(1,2*itime-1) .gt. 0)
   !This construction avoids saving the contours at t = 0.
endif

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
    qqr4(iy,ix)=real(wka(iy,ix))
  enddo
enddo
write(31,rec=igrids) tr4,qqr4

 !Compute domain integral of (q-beta*y)^2 and write to norms.asc:
call l2norm(wka,qql2)

write(*,'(a,f7.2,1(a,f13.6))') ' t = ',tgrid,' <q^2> = ',qql2

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
