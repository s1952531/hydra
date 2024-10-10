module evolution

! Module contains subroutines to evolve PV according
! to the algorithm detailed in casl.f90.

use common

implicit none

 !Energies:
double precision:: enepre,enepost

 !Physical fields:
double precision:: qq(ny,nx),qc(ny,nx)
double precision:: qqpre(ny,nx),qspre(ny,nx),qdpre(ny,nx)
double precision:: uu(ny,nx),vv(ny,nx),pp(ny,nx)
 !Semi-Lagrangian advection arrays:
double precision:: x0(ny,nx), y0(ny,nx)
double precision:: ula(ny,nx),vla(ny,nx)
 !Allocatable source arrays:
double precision,allocatable,dimension(:,:):: sq,sqpre,sqstoch

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

double precision:: qqmin,qqmax
integer:: ireg,igsave,icsave,ix,iy

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

   !Perform surgery & field reset 
  if (twist .gt. twistmax) then
    ireg=ireg+1
     !Don't continue if maximum number of regularisations reached;
     !it is time to recontour (pass control back to main program):
    if (ireg .eq. nregmax) then
       !Prepare PV residual for recontouring:
      do ix=1,nx
        do iy=1,ny
          qd(iy,ix)=qq(iy,ix)-qc(iy,ix)
          qs(iy,ix)=qq(iy,ix)
        enddo
      enddo
      if (vardq) then
         !Re-compute PV contour interval (not done if beta .ne. 0):
        qqmax=zero
        qqmin=zero
        do ix=1,nx
          do iy=1,ny
            qqmax=max(qqmax,qq(iy,ix))
            qqmin=min(qqmin,qq(iy,ix))
          enddo
        enddo
        qjump=(qqmax-qqmin)/dble(ncontq)
         !ncontq is set in parameters.f90. 
      endif
       !De-allocate memory for the PV source arrays if they are non-zero:
      if (damping)  deallocate(sq,sqpre)
      if (stoch)    deallocate(sqstoch)
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the PV contours (surgery + node redistribution):
    call surgery
     !Record contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

     !Convert PV contours to gridded values (qc):
    call con2grid(qc)
     !Reset qs and qd:
    call reset(qc,qs,qd)

     !Copy gridded PV fields to old time level:
    do ix=1,nx
      do iy=1,ny
        qqpre(iy,ix)=qs(iy,ix) 
        qdpre(iy,ix)=qd(iy,ix) 
        qspre(iy,ix)=qs(iy,ix) 
      enddo
    enddo
    twist=twist-twistmax
  endif

   !Adjust timestep (dt) on maximum vorticity magnitude:
  call adapt(igsave,icsave)
   !Advect PV from time t to t + dt:
  call advance

   !Update the time:
  t=t+dt

   !Possibly save PV & energy at chosen save time (tgrid):
  if (igsave .eq. 1) call savegrid

   !Possibly save contours for post processing:
  if (icsave .eq. 1) call savecont
  
   !Copy new fields into previous time:
  do ix=1,nx
    do iy=1,ny
      qqpre(iy,ix)=qq(iy,ix)
      qdpre(iy,ix)=qd(iy,ix)
      qspre(iy,ix)=qs(iy,ix)
    enddo
  enddo

  if (damping) then  
    do ix=1,nx
      do iy=1,ny
        sqpre(iy,ix)=sq(iy,ix)
      enddo
    enddo
  endif
enddo

!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
return
end subroutine

!=======================================================================

subroutine init

! Initialises quantities needed for normal time integration following 
! contour regeneration

implicit double precision(a-h,o-z)
implicit integer(i-n)

!---------------------------------------------
 !Record contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

!---------------------------------------------------
 !Allocate memory for the PV source arrays if they are non-zero:
if (damping)  allocate(sq(ny,nx),sqpre(ny,nx))
if (stoch)    allocate(sqstoch(ny,nx))

 !Convert PV contours to gridded values (qc):
call con2grid(qc)

!-----------------------------------------------------
 !Define residual PV qd = qs-qc-F[qs-qc]
do ix=1,nx
  do iy=1,ny
    ula(iy,ix)=qs(iy,ix)-qc(iy,ix)
    vla(iy,ix)=ula(iy,ix)
  enddo
enddo

call filter(ula,2)

do ix=1,nx
  do iy=1,ny
    qd(iy,ix)=vla(iy,ix)-ula(iy,ix)
  enddo
enddo

 !Copy gridded PV fields to old time level:
do ix=1,nx
  do iy=1,ny
    qqpre(iy,ix)=qs(iy,ix) 
    qdpre(iy,ix)=qd(iy,ix) 
    qspre(iy,ix)=qs(iy,ix) 
  enddo
enddo

!------------------------------------------------------
 !Calculate the source term for PV (sqpre) if non-zero:
if (damping) call getqqsrc(sqpre)
 !Note: qqpre and sqpre are needed by subroutine advance.

 !Get the initial velocity field (uu,vv):
if (qjump .gt. zero) then 
  call inversion
else
   !Here PV is zero identically - no need to do inversion:
  do ix=1,nx
    do iy=1,ny
      uu(iy,ix)=zero
      vv(iy,ix)=zero 
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine getqqsrc(qqsrc)

! Gets the source term for the PV evolution equation.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed array:
double precision:: qqsrc(ny,nx)

!-----------------------------------------------------------------
 !Implement Ekman and/or thermal damping:
if (heating) then
  if (friction) then 
   !Use thermal and Ekman damping:
    do ix=1,nx
      do iy=1,ny
        qqsrc(iy,ix)=therm*(pp(iy,ix)-ppeq(iy,ix))-&
               & rekman*(qq(iy,ix)+kdsq*pp(iy,ix))
      enddo
    enddo
  else
    do ix=1,nx
      do iy=1,ny
        qqsrc(iy,ix)=therm*(pp(iy,ix)-ppeq(iy,ix))
      enddo
    enddo
  endif
else
   !Only use Ekman damping
  do ix=1,nx
    do iy=1,ny
      qqsrc(iy,ix)=-rekman*(qq(iy,ix)+kdsq*pp(iy,ix))
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine inversion

! Inverts Laplace's operator on PV anomaly (PV - beta*y) to obtain 
! the streamfunction (pp) and the velocity (uu,vv) = (-dpp/dy,dpp/dx).

implicit double precision(a-h,o-z)
implicit integer(i-n)

!------------------------------------------------------------
 !Call con2grid to get updated contour PV (qc):
call con2grid(qc)

 !Combine fields to update qq with full field:
call combine(qq,qc,qs,qd)

 !Invert PV to obtain velocity field: 
call main_invert(qq,uu,vv,pp)

return
end subroutine

!=======================================================================
      
subroutine advance

! Computes qq(t+dt) by contour advection by trajectory 
! integration using a standard semi-Lagrangian (SL) scheme.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define local parameters and arrays:
integer,parameter:: niter=2
double precision:: uq(nptq),vq(nptq),xqm(nptq),yqm(nptq)

!------------------------------------------------------------------------
 !Copy current velocity field into (ula,vla) for use in subroutines below:
do ix=1,nx
  do iy=1,ny
    ula(iy,ix)=uu(iy,ix)
    vla(iy,ix)=vv(iy,ix)
  enddo
enddo

 !Increments in grid units needed below for trajectory integration:
gcx=dt*glxi
gcy=dt*glyi
hgcx=f12*gcx
hgcy=f12*gcy

 !Compute Euler backward predictor for departure grid point (x0,y0):
do ix=1,nx
  do iy=1,ny
    x0(iy,ix)=mod(xigmax+xig(ix)-gcx*uu(iy,ix),xigmax)
    y0(iy,ix)=mod(yigmax+yig(iy)-gcy*vv(iy,ix),yigmax)
  enddo
enddo
 !Note, (uu,vv) is used since we have no other velocity field available

 !Prepare contour evolution; get velocity on PV contour nodes:
if (nptq .gt. 0) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xq(i)+hfdt*uq(i)
    yy=yq(i)+hfdt*vq(i)
    xqm(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yqm(i)=oms*(yy-elly*dble(int(yy*hlyi)))
    xx=xq(i)+dt*uq(i)
    yy=yq(i)+dt*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

 !Iterate to converge on implicit trapezoidal integration:
do iter=1,niter
   !Obtain qs & qd at time t + dt:
  call sl_step

   !Obtain qq & hence uu & vv at time t + dt:
  call inversion

   !Correct departure grid point (x0,y0):
  do ix=1,nx
    do iy=1,ny
       !Obtain old velocity (time t) at the departure grid point using
       !bi-linear interpolation of (ula,vla):
      ix0=1+int(x0(iy,ix))
      pxc=dble(ix0)-x0(iy,ix)
      px=one-pxc
      ix1=ixp(ix0)

      iy0=1+int(y0(iy,ix))
      pyc=dble(iy0)-y0(iy,ix)
      py=one-pyc
      iy1=iyp(iy0)

      uod=pyc*(pxc*ula(iy0,ix0)+px*ula(iy0,ix1)) &
      &   +py*(pxc*ula(iy1,ix0)+px*ula(iy1,ix1))

      vod=pyc*(pxc*vla(iy0,ix0)+px*vla(iy0,ix1)) &
      &   +py*(pxc*vla(iy1,ix0)+px*vla(iy1,ix1))

      x0(iy,ix)=mod(xigmax+xig(ix)-hgcx*(uod+uu(iy,ix)),xigmax)
      y0(iy,ix)=mod(yigmax+yig(iy)-hgcy*(vod+vv(iy,ix)),yigmax)
       !(uu,vv) is the new velocity (time t+dt) at the arrival grid point
    enddo
  enddo

   !Update the PV contour points:
  if (nptq .gt. 0) then
    call velint(uu,vv,uq,vq)
    do i=1,nptq
      xx=xqm(i)+hfdt*uq(i)
      yy=yqm(i)+hfdt*vq(i)
      xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
      yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
    enddo
  endif

enddo

 !Obtain final corrected qs & qd at time t + dt:
call sl_step

 !Obtain final corrected uu & vv at time t + dt from qq:
call inversion

 !Update the PV contour points:
if (nptq .gt. 0) then
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqm(i)+hfdt*uq(i)
    yy=yqm(i)+hfdt*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo
endif

 !Call con2grid to get updated contour PV (qc):
call con2grid(qc)

 !Combine fields to update qq with full field:
call combine(qq,qc,qs,qd)

return
end subroutine

!=======================================================================
      
subroutine sl_step

! Carries out semi-Lagrangian integration of qs & qd (including any 
! sources in qd) from t to t+dt

implicit double precision(a-h,o-z)
implicit integer(i-n)

!-------------------------------------------------------------------
 !Integrate qs using bi-cubic Lagrange interpolation of qspre at x0,y0:
call interpol(qspre,qs,x0,y0)
 !This approximates Dqs/Dt = 0 when beta = 0

if (beffect) then
   !Add displacement in y to correct for beta:
  do ix=1,nx
    do iy=1,ny
      yy=gly*(y0(iy,ix)-yig(iy))
      qs(iy,ix)=qs(iy,ix)+beta*(yy-elly*dble(int(yy*hlyi)))
    enddo
  enddo
endif

 !Integrate also qd:
call interpol(qdpre,qd,x0,y0)
 !Here, we obtain only the adiabatic evolution, Dqd/Dt = 0; add source next:

!-------------------------------------------------------------------
if (damping) then
  
   !Add diabatic source, Dqd/Dt = sq:
  call getqqsrc(sq)

  do ix=1,nx
    do iy=1,ny
       !Obtain old source (at time t) at the departure grid point using
       !bi-linear interpolation of sqpre:
      ix0=1+int(x0(iy,ix))
      pxc=dble(ix0)-x0(iy,ix)
      px=one-pxc
      ix1=ixp(ix0)

      iy0=1+int(y0(iy,ix))
      pyc=dble(iy0)-y0(iy,ix)
      py=one-pyc
      iy1=iyp(iy0)

      sqod=pyc*(pxc*sqpre(iy0,ix0)+px*sqpre(iy0,ix1)) &
        &  +py*(pxc*sqpre(iy1,ix0)+px*sqpre(iy1,ix1))
  
       !Integrate in time using trapezoidal rule (sq is the new source
       !at time t+dt) and add to qd:
      qd(iy,ix)=qd(iy,ix)+hfdt*(sqod+sq(iy,ix))
    enddo
  enddo
endif

!---------------------------------------------------------------------
 !Add any point vortices here to force the flow (they are converted
 !to gridded PV values and added to qd in spectral space):
if (stoch) then
   !Initialise array to uptake the vorticity of the vortices:
  do ix=1,nx
    do iy=1,ny
      sqstoch(iy,ix)=zero
    enddo
  enddo

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

      sqstoch(iy0,ix0)=sqstoch(iy0,ix0)+svor*pyc*pxc
      sqstoch(iy0,ix1)=sqstoch(iy0,ix1)+svor*pyc*px
      sqstoch(iy1,ix0)=sqstoch(iy1,ix0)+svor*py*pxc
      sqstoch(iy1,ix1)=sqstoch(iy1,ix1)+svor*py*px
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

      sqstoch(iy0,ix0)=sqstoch(iy0,ix0)-vorvor*(pyc*cth+pxc*sth)
      sqstoch(iy0,ix1)=sqstoch(iy0,ix1)+vorvor*(pyc*cth -px*sth)
      sqstoch(iy1,ix0)=sqstoch(iy1,ix0)+vorvor*(pxc*sth -py*cth)
      sqstoch(iy1,ix1)=sqstoch(iy1,ix1)+vorvor*( px*sth +py*cth)
    enddo

  endif

   !Total number of vortices/dipoles added so far:
  totnvor=totnvor+dble(nvor)

   !Finally add forcing to qd:
  do ix=1,nx
    do iy=1,ny
      qd(iy,ix)=qd(iy,ix)+sqstoch(iy,ix)
    enddo
  enddo
endif

return
end subroutine

!=======================================================================

subroutine adapt(igsave,icsave)

! Adapts the time step to ensure dt < dtfac/max(|zeta|_max)

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision,parameter:: dtfac=pi/40.d0

!----------------------------------------------------------------------

 !Compute accurate advection time step:
zzrms=zero
zzmax=small
do ix=1,nx
  do iy=1,ny
    zz=qq(iy,ix)+kdsq*pp(iy,ix)
    zzrms=zzrms+zz**2
    zzmax=max(zzmax,abs(zz))
     !Here, qq = gridded total PV anomaly (q - beta*y).
  enddo
enddo
zzrms=sqrt(zzrms*dsumi)
dtacc=dtfac/max(zzmax,srwfm)

 !The restriction on the maximum Rossby wave frequency (srwfm)
 !ensures that the fastest Rossby wave frequency is resolved.

!---------------------------------------------------------------------
 !Choose a new time step: 
if (dt .gt. zero) then 
  dt=min(dtacc,dtmax)
  if (dt .gt. dtacc) write(*,'(a,f9.5)') 'Warning! dt/dt_acc= ',dt/dtacc
else
   !Limit max timestep to a data save time
  dtmax=min(tgsave,tcsave)
  dt=min(dtacc,dtmax)
endif 
hfdt=dt/two

 !Increment the integral of max|zz|:
twist=twist+dt*zzmax

!---------------------------------------------------------------------
 !Record various diagnostics to monitor.asc:
write(16,'(1x,f16.8,4(1x,f12.6))') t,f12*zzrms**2,zzrms,zzmax,twist

!---------------------------------------------------------------------

 !Set flag to save gridded data every tgsave time units:
itime=int((t+dt)/tgsave)
tgrid=tgsave*dble(itime)+small
 !small = 1.d-12 is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  igsave=1
   !Compute energy:
  call binorm(pp,qq,enepre)
  enepre=-f12*enepre
else
   !Do not save data:
  igsave=0
endif

 !Set flag to save contour data every tcsave time units:
itime=int((t+dt)/tcsave)
tcont=tcsave*dble(itime)
if (t .le. tcont .and. t+dt .gt. tcont) then
   !The save time is between t & t+dt; set flag to save data:
  icsave=sign(1,2*itime-1)
   !This construction avoids saving the contours at t = 0.
else
   !Do not save data:
  icsave=0
endif

return
end subroutine

!=======================================================================
      
subroutine savegrid

! Saves PV, energy and enstrophy spectrum at the desired save time

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: qspec(0:max(nx,ny))
real:: qqr4(ny,nx),tr4

 !Increment counter for direct file access:
igrids=igrids+1

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Compute energy
call binorm(pp,qq,enepost)
enepost=-f12*enepost

 !Compute time interpolated energies:
ene=pt*enepre+ptc*enepost

 !Interpolate PV anomaly at save time:
do ix=1,nx
  do iy=1,ny
    ula(iy,ix)=pt*qqpre(iy,ix)+ptc*qq(iy,ix)
  enddo
enddo

!Write full PV to qq.r4:
do ix=1,nx
  do iy=1,ny
    qqr4(iy,ix)=real(ula(iy,ix)+bety(iy))
  enddo
enddo
tr4=real(tgrid)

write(31,rec=igrids) tr4,qqr4

 !Compute domain integral of (q-beta*y)^2:
call l2norm(ula,qql2)

 !Write energy to ene.asc:
write(15,'(f7.2,1x,f16.9)') tgrid,ene

write(*,'(a,f7.2,2(a,f13.6))') &
    & ' t = ',tgrid,'  <q^2> = ',qql2,'  E = ',ene

 !Compute 1d PV spectrum:
call spec1d(ula,qspec)
sumqspec=zero
do k=1,kmax
  sumqspec=sumqspec+qspec(k)
   !Normalise to take into account uneven sampling of wavenumbers 
   !in each shell [k-1/2,k+1/2]:
  qspec(k)=spmf(k)*qspec(k)
enddo
sumqspec=8.d0*sumqspec*dsumi
 !Write out spectrum to file:
write(50,'(f7.2,2(1x,f16.9),1x,i5)') tgrid,sumqspec,qql2,kmaxred
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(50,'(2(1x,f12.8))') alk(k),log10(qspec(k))
enddo
 !Note: alk(k) = log_10(k)

return
end subroutine

!=======================================================================
      
subroutine savecont

! Saves PV contours and residual PV for post-processing and imaging

implicit double precision(a-h,o-z)
implicit integer(i-n)

real:: qdr4(ny,nx),tr4
character(len=3):: pind

write(*,'(a,f12.5)') ' Saving contours at t = ',t
irec=nint(t/tcsave)
write(pind(1:3),'(i3.3)') irec

tr4=real(t)
 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,qjump

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
do ix=1,nx
  do iy=1,ny
    qdr4(iy,ix)=real(qq(iy,ix)-qc(iy,ix))
  enddo
enddo
write(83,rec=irec) tr4,qdr4

 !Save PV contours if any exist:
if (nq .gt. 0) then
  open(81,file='cont/qqindex'//pind,form='unformatted', &
      & access='direct',status='replace',recl=12*nq)
  write(81,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq)
  close(81)

  open(82,file='cont/qqnodes'//pind,form='unformatted', &
      & access='direct',status='replace',recl=16*nptq)
  write(82,rec=1) xq(1:nptq),yq(1:nptq)
  close(82)
endif

return
end subroutine

!=======================================================================

 !Main end module
end module
