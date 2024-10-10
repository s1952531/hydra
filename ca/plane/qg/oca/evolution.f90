module evolution

! Module contains subroutines to evolve PV according
! to the algorithm detailed in casl.f90.

use common

implicit none

 !Energies:
double precision:: enepre,enepost

 !Physical fields:
double precision:: uu(ny,nx),vv(ny,nx),pp(ny,nx),qqpre(ny,nx)

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

integer:: igsave,icsave,ix,iy

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

 !Used for regularising contours:
twist=zero

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tfin)

   !Perform surgery & field reset 
  if (twist .gt. twistmax) then
     !Regularise the PV contours (surgery + node redistribution):
    call surgery
     !Record contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

     !Convert PV contours to gridded values (qq):
    call con2grid(qq)

    twist=twist-twistmax
  endif

   !Adjust timestep (dt) on maximum PV magnitude:
  call adapt(igsave,icsave)
   !Advect PV from time t to t + dt:
  call advance

   !Update the time:
  t=t+dt

   !Possibly save PV & energy at chosen save time (tgrid):
  if (igsave .eq. 1) call savegrid

   !Possibly save contours for post processing:
  if (icsave .eq. 1) call savecont
  
   !Copy qq into previous time:
  do ix=1,nx
    do iy=1,ny
      qqpre(iy,ix)=qq(iy,ix)
    enddo
  enddo
enddo

!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
return
end subroutine

!=======================================================================

subroutine init

! Initialises quantities needed for normal time integration:

implicit double precision(a-h,o-z)
implicit integer(i-n)

!---------------------------------------------
 !Record contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

!---------------------------------------------------
 !Convert PV contours to gridded values (qq):
call con2grid(qq)

 !Copy qq into previous time:
do ix=1,nx
  do iy=1,ny
    qqpre(iy,ix)=qq(iy,ix)
  enddo
enddo

 !Get the initial velocity field (uu,vv):
call inversion(.true.)

return
end subroutine

!=======================================================================

subroutine inversion(ppflag)

! Inverts Laplace's operator on PV (qq) to obtain the 
! streamfunction (pp) and the velocity (uu,vv) = (-dpp/dy,dpp/dx).
! pp is returned only if ppflag is true.

implicit double precision(a-h,o-z)
implicit integer(i-n)

logical:: ppflag

!------------------------------------------------------------
 !Call con2grid to get updated PV (qq):
call con2grid(qq)

 !Invert PV to obtain velocity field: 
call main_invert(qq,uu,vv,pp,ppflag)      

return
end subroutine

!=======================================================================
      
subroutine advance

! Computes contour positions at t+dt by an implicit trapezoidal scheme

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define local parameters and arrays:
integer,parameter:: niter=2
double precision:: uq(nptq),vq(nptq),xqm(nptq),yqm(nptq)

!------------------------------------------------------------------------
 !Prepare contour evolution; get velocity on PV contour nodes:
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

 !Iterate to converge on implicit trapezoidal integration:
do iter=1,niter

   !Obtain qq & hence uu & vv at time t + dt:
  call inversion(.false.)

   !Update the PV contour points:
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqm(i)+hfdt*uq(i)
    yy=yqm(i)+hfdt*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
  enddo

enddo

 !Obtain final corrected uu & vv at time t + dt from qq:
call inversion(.true.)

 !Update the PV contour points:
call velint(uu,vv,uq,vq)
do i=1,nptq
  xx=xqm(i)+hfdt*uq(i)
  yy=yqm(i)+hfdt*vq(i)
  xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  yq(i)=oms*(yy-elly*dble(int(yy*hlyi)))
enddo

 !Call con2grid to get updated PV (qq):
call con2grid(qq)

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
write(16,'(1x,f12.5,4(1x,f10.5))') t,f12*zzrms**2,zzrms,zzmax,twist

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

 !---------------------------------------------
 !Compute rotation rate of flow and steadiness:
call dynomega(uu,vv,pp,omega,epsb)
 !Record in steady.asc:
write(17,'(1x,f12.5,1x,f14.10,1x,f16.12)') t,omega,epsb

return
end subroutine

!=======================================================================
      
subroutine savegrid

! Saves qq at the desired save time to files:

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: wka(ny,nx),qspec(0:max(nx,ny))
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

 !Store PV at save time:
do ix=1,nx
  do iy=1,ny
    wka(iy,ix)=pt*qqpre(iy,ix)+ptc*qq(iy,ix)
  enddo
enddo

!Write gridded fields to file:             
do ix=1,nx
  do iy=1,ny
    qqr4(iy,ix)=real(wka(iy,ix))
  enddo
enddo
tr4=real(tgrid)

write(31,rec=igrids) tr4,qqr4

 !Compute domain integral of q^2:
call l2norm(wka,qql2)

 !Write energy to ene.asc:
write(15,'(f7.2,1x,f14.9)') tgrid,ene

write(*,'(a,f7.2,2(a,f10.6))') &
    & ' t = ',tgrid,'  <q^2> = ',qql2,'  E = ',ene

 !Compute 1d PV spectrum:
call spec1d(wka,qspec)
sumqspec=zero
do k=1,kmax
  sumqspec=sumqspec+qspec(k)
   !Normalise to take into account uneven sampling of wavenumbers 
   !in each shell [k-1/2,k+1/2]:
  qspec(k)=spmf(k)*qspec(k)
enddo
sumqspec=8.d0*sumqspec*dsumi
 !Write out spectrum to file:
write(50,'(f7.2,2(1x,f14.9),1x,i5)') tgrid,sumqspec,qql2,kmaxred
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

! Saves qq contours and residual qq for post-processing via
! congen.f90

implicit double precision(a-h,o-z)
implicit integer(i-n)

character(len=3):: pind

write(*,'(a,f12.5)') ' Saving contours at t = ',t
irec=nint(t/tcsave)
write(pind(1:3),'(i3.3)') irec

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5)') nq,nptq,t

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
