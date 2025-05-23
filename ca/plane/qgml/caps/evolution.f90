module evolution

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! This module contains subroutines to evolve the PV contours and fields
! according to the algorithm detailed in caps.f90.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 !Import common areas:
use common

implicit none

! Streamfunction and velocity field:
double precision:: pp(ng,ng,nz),uu(ng,ng,nz),vv(ng,ng,nz)

! Relative vertical vorticity field:
double precision:: zz(ng,ng,nz)

! Spectral integrating factors used for residual PV evolution:
double precision:: emq(ng,ng),epq(ng,ng)

! Time step, time step fractions and "twist" parameter:
double precision:: dt,dt2,dt3,dt6,twist

! Logicals for saving gridded & contour data
logical:: gsave,csave

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================

subroutine advect

! Main subroutine for evolving PV contours and fields.

implicit none

! Local variables:
double precision,parameter:: twistmax=2.5d0
!      twistmax: the maximum value of the time integral of |zeta|_max
!                between regularisations of the contours.
double precision,parameter:: qratmax=0.2d0
!      qratmax:  the maximum ratio r of the mean-square residual PV qd
!                to the mean-square PV qc in the contours
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, or when r > qratmax, 
!      the code rebuild the PV contours in a separate memory space.
integer:: ireg,iz

double precision:: qc(ng,ng,nz)
double precision:: wkp(ng,ng),wks(ng,ng),qrat

!------------------------------------------------------------------
! Initialise after contour regeneration:
call init

! Used for regularising contours:
twist=zero

! Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Start the time loop:
do while (t <= tsim)

   ! Regularise contours when the twist parameter is large enough:
   if (twist > twistmax) then
      ireg=ireg+1

      ! Compute ratio of mean-square residual and contour PV:
      call con2grid(qc)
      qrat=zero
      do iz=1,nz
         ! Only process layers having contours:
         if (jl2q(iz) > 0) then
            wks=qd(:,:,iz)
            call spctop(ng,ng,wks,wkp,xfactors,yfactors,xtrig,ytrig)
            qrat=max(qrat,sum(wkp**2)/sum(qc(:,:,iz)**2))
         endif
      enddo

      ! Don't continue if maximum number of regularisations reached:
      if (ireg == nregmax .or. qrat > qratmax) then
         ! Prepare PV residual for recontouring:
         call prepare
         ! Exit module and go to recontouring:
         return
      endif

      ! Regularise the PV contours (surgery + node redistribution):
      call surgery
      ! Record active contour complexity to complexity.asc:
      write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

      twist=twist-twistmax
   endif

   ! Advect flow from time t to t + dt:
   call advance

enddo
! End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! Save final data if not done already:
if (int(t/tgsave) == igrids) call savegrid
if (int(t/tcsave) == iconts) call savecont

return
end subroutine advect

!=======================================================================

subroutine init

! Initialises residual PV for normal time integration following 
! contour regeneration

implicit none

! Local variables:
double precision:: qc(ng,ng,nz)
double precision:: wkp(ng,ng),wks(ng,ng)
integer:: iz

!-------------------------------------------------------------
! Logicals used for saving gridded fields and contours:
gsave=.false.
csave=.false.

! Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)

! Define (spectral) residual PV qd = (1-F)[qs-qc]:
do iz=1,nz
   wkp=qc(:,:,iz)
   call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
   qd(:,:,iz)=bfhi*(qs(:,:,iz)-wks)
enddo
! Here bfhi = 1-F is a high-pass spectral filter

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current (spectral) PV field is stored in qs, and the residual PV
! needed in congen.f90 is stored in qq.

implicit none

! Local variables:
double precision:: qc(ng,ng,nz)
double precision:: wkp(ng,ng),wks(ng,ng)
integer:: iz

!-------------------------------------------------------------
! Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)

! Define (spectral) PV qs and residual PV qd:
do iz=1,nz
   wkp=qc(:,:,iz)
   call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
   qs(:,:,iz)=bflo*qs(:,:,iz)+bfhi*wks+qd(:,:,iz)
   ! Here bflo = F and bfhi = 1-F are low- and high-pass filters
   qd(:,:,iz)=qs(:,:,iz)-wks
   ! Convert qd to physical space as qq (used in recontouring):
   call spctop(ng,ng,qd(:,:,iz),qq(:,:,iz),xfactors,yfactors,xtrig,ytrig)
   ! Note: qd is overwritten, but we are leaving this module next and
   !       qd will be redefined upon re-entry in subroutine init.
enddo

return
end subroutine prepare

!=======================================================================

subroutine inversion(level)

! Inverts the PV definition to obtain the streamfunction (pp) and
! the velocity (uu,vv) = (-dpp/dy,dpp/dx).

! The parameter level = 0 at the beginning of a time step in order
! to reset the PV fields qs & qd.  Otherwise, we merely combine qc
! (the contour PV field), qs & qd to define the total PV for inversion.

implicit none

! Passed variable:
integer:: level

! Local variables:
double precision:: qc(ng,ng,nz)
double precision:: wkp(ng,ng),wks(ng,ng)
integer:: iz

!-------------------------------------------------------------
! Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)

if (level == 0) then
   ! Reset qs & qd, and combine qs, qd & qc to update
   ! the full PV field qq:
   do iz=1,nz
      wkp=qc(:,:,iz)
      call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
      ! Here wks is qc in spectral space.
      qs(:,:,iz)=bflo*qs(:,:,iz)+bfhi*wks+qd(:,:,iz)
      qd(:,:,iz)=bfhi*(qs(:,:,iz)-wks)
      wks=qs(:,:,iz)
      call spctop(ng,ng,wks,qq(:,:,iz),xfactors,yfactors,xtrig,ytrig)
   enddo
   ! bflo & bfhi = 1 - bflo are low- & high-pass filters (see spectral.f90)

   ! Reset also the mean value of the PV in the bottom layer (the only
   ! layer in which the mean PV can change due to Ekman friction):
   qavg(nz)=danorm*sum(qq(:,:,nz))
   ! Here, danorm = 1/ng^2.
else
   ! Combine qs, qd & qc to update the full PV field qq:
   do iz=1,nz
      wkp=qc(:,:,iz)
      call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
      ! Here wks is qc in spectral space.
      wks=bflo*qs(:,:,iz)+bfhi*wks+qd(:,:,iz)
      call spctop(ng,ng,wks,qq(:,:,iz),xfactors,yfactors,xtrig,ytrig)
   enddo
endif

! Invert PV field:
call main_invert(qq,pp,uu,vv)

return
end subroutine inversion

!=======================================================================
      
subroutine advance

! Advanced PV contours and fields one time step using a 4th-order
! Runge-Kutta time-integration method.

implicit none

! Local variables:
double precision:: qsi(ng,ng,nz),qsf(ng,ng,nz),sqs(ng,ng,nz)
double precision:: qdi(ng,ng,nz),qdf(ng,ng,nz),sqd(ng,ng,nz)
double precision:: xqi(nptq),yqi(nptq),xqf(nptq),yqf(nptq)
double precision:: uq(nptq),vq(nptq)
double precision:: xx,yy
integer:: i,iz

!------------------------------------------------------------------
! RK4 predictor step to time t0 + dt/2:

! Invert PV and compute velocity:
call inversion(0)

! Adapt the time step and save various diagnostics each time step:
call adapt

! Possibly save data (gsave & csave set by adapt):
if (gsave) call savegrid
if (csave) call savecont

! Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,0)

call velint(uu,vv,uq,vq)
do i=1,nptq
   xqi(i)=xq(i)
   yqi(i)=yq(i)
   xx=xqi(i)+dt2*uq(i)
   xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
   yy=yqi(i)+dt2*vq(i)
   yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
   xx=xqi(i)+dt6*uq(i)
   yy=yqi(i)+dt6*vq(i)
   xqf(i)=oms*(xx-twopi*dble(int(xx*pinv)))
   yqf(i)=oms*(yy-twopi*dble(int(yy*pinv)))
enddo

qsi=qs
qs=qsi+dt2*sqs
qsf=qsi+dt6*sqs
qdi=qd
do iz=1,nz
   qd(:,:,iz)=emq*(qdi(:,:,iz)+dt2*sqd(:,:,iz))
enddo
qdf=qdi+dt6*sqd

!------------------------------------------------------------------
! RK4 corrector step at time t0 + dt/2:
t=t+dt2

! Invert PV and compute velocity:
call inversion(1)

! Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,1)

call velint(uu,vv,uq,vq)
do i=1,nptq
   xx=xqi(i)+dt2*uq(i)
   xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
   yy=yqi(i)+dt2*vq(i)
   yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
   xx=xqf(i)+dt3*uq(i)
   yy=yqf(i)+dt3*vq(i)
   xqf(i)=oms*(xx-twopi*dble(int(xx*pinv)))
   yqf(i)=oms*(yy-twopi*dble(int(yy*pinv)))
enddo

qs=qsi+dt2*sqs
qsf=qsf+dt3*sqs
do iz=1,nz
   qd(:,:,iz)=emq*(qdi(:,:,iz)+dt2*sqd(:,:,iz))
enddo
qdf=qdf+dt3*sqd

!------------------------------------------------------------------
! RK4 predictor step at time t0 + dt:

! Invert PV and compute velocity:
call inversion(1)

! Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,1)

call velint(uu,vv,uq,vq)
do i=1,nptq
   xx=xqi(i)+dt*uq(i)
   xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
   yy=yqi(i)+dt*vq(i)
   yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
   xx=xqf(i)+dt3*uq(i)
   yy=yqf(i)+dt3*vq(i)
   xqf(i)=oms*(xx-twopi*dble(int(xx*pinv)))
   yqf(i)=oms*(yy-twopi*dble(int(yy*pinv)))
enddo

qs=qsi+dt*sqs
qsf=qsf+dt3*sqs
emq=emq**2
do iz=1,nz
   qd(:,:,iz)=emq*(qdi(:,:,iz)+dt*sqd(:,:,iz))
enddo
qdf=qdf+dt3*sqd

!------------------------------------------------------------------
! RK4 corrector step at time t0 + dt:
t=t+dt2

! Invert PV and compute velocity:
call inversion(1)

! Calculate the source terms (sqs,sqd) for PV (qs,qd):
call source(sqs,sqd,2)

call velint(uu,vv,uq,vq)
do i=1,nptq
   xx=xqf(i)+dt6*uq(i)
   xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
   yy=yqf(i)+dt6*vq(i)
   yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
enddo

qs=qsf+dt6*sqs
do iz=1,nz
   qd(:,:,iz)=emq*(qdf(:,:,iz)+dt6*sqd(:,:,iz))
enddo

return
end subroutine advance

!=======================================================================

subroutine source(sqs,sqd,lev)

! Gets the source terms (sqs,sqd) for the PV (qs,qd), including
! wind-stress forcing in the uppermost layer (1) and Ekman drag
! in the lowest layer (nz). Here, lev is the level of the Runge
! Kutta integration, used for defining the integrating factors.

! --- All sources are in spectral space.

implicit none

! Passed variables:
double precision:: sqs(ng,ng,nz),sqd(ng,ng,nz)
integer:: lev

! Local variables:
double precision:: qx(ng,ng,nz),qy(ng,ng,nz)
double precision:: wkp(ng,ng),wks(ng,ng)
integer:: iz

!---------------------------------------------------------------
! qs source --- only NL advection term is needed:
call gradient(qs,qx,qy)
qx=-uu*qx-vv*qy

do iz=1,nz
   call ptospc(ng,ng,qx(:,:,iz),sqs(:,:,iz),xfactors,yfactors,xtrig,ytrig)
   ! Mean source must remain zero since Dq_s/Dt = 0 & div(u,v) = 0
   sqs(1,1,iz)=zero
enddo

!---------------------------------------------------------------
! qd source --- include wind-stress forcing and Ekman damping:
call gradient(qd,qx,qy)
if (beffect) then
   !Include beta effect if beta > 0 (evolve PV anomaly):
  qx=-uu*qx-vv*(qy+beta)
else
  qx=-uu*qx-vv*qy
endif
do iz=1,nz
   call ptospc(ng,ng,qx(:,:,iz),sqd(:,:,iz),xfactors,yfactors,xtrig,ytrig)
   ! Mean advective source must remain zero
   sqd(1,1,iz)=zero
enddo

if (wind) then
   ! Add wind-stress forcing to the upper layer (iz = 1), i.e.
   ! Dq_1/Dt = fwind*sin(y):
   sqd(:,:,1)=sqd(:,:,1)+sfwind
   ! This is done in spectral space; sfwind is proportional to fwind
   ! and is defined in subroutine initialise in caps.f90.
   ! *** Note, this has zero average so won't change the mean value
   !     of qd; be careful when using other wind-stress forcing
   !     functions - it may be best to avoid any mean tendency.
endif

if (friction) then
   ! Add contribution of Ekman damping to the lowest layer (iz = nz):
   wkp=qq(:,:,nz)-bety+kkm(nz)*(pp(:,:,nz)-pp(:,:,nz-1))
   ! wkp is the lowest layer relative vorticity in physical space;
   ! kkm(iz) = f^2/(H_iz b'_{iz-1}), see init_spectral in spectral.f90.
   if (bath) wkp=wkp-qb
   ! Need to remove bathymetry if present.
   call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)   
   ! Add -r_ekman * vorticity to lowest layer PV tendency:
   sqd(:,:,nz)=sqd(:,:,nz)-rekman*wks
   ! *** Note, the mean tendency in sqd(1,1,nz) may change as a
   !     result of Ekman friction.  This is accounted for when
   !     resetting qs & qd at the beginning of each time step:
   !     see the call to inversion(0).
endif

if (lev == 0) then
   ! Spectrally truncate sources:
   do iz=1,nz
      sqs(:,:,iz)=filt*sqs(:,:,iz)
      sqd(:,:,iz)=filt*sqd(:,:,iz)
   enddo
   ! Apply exponential integrating factors and spectrally truncate sources:
else if (lev == 1) then
   do iz=1,nz
      sqs(:,:,iz)=filt*sqs(:,:,iz)
      sqd(:,:,iz)= epq*sqd(:,:,iz)
   enddo
else
   ! Apply exponential integrating factors and spectrally truncate sources:
   epq=epq**2
   do iz=1,nz
      sqs(:,:,iz)=filt*sqs(:,:,iz)
      sqd(:,:,iz)= epq*sqd(:,:,iz)
   enddo
endif

return
end subroutine source

!=======================================================================

subroutine adapt
! Adapts the time step, computes the twist parameter (the time integral
! of the maximum absolute vorticity), and various quantities every
! time step to monitor the flow evolution. Also updates the spectral
! integrating factors used in the evolution of the PV residual.

implicit none

! Local parameter used for setting time step:
double precision,parameter:: dtfac=pi/40.d0

! Other local variables:
double precision:: wka(ng,ng)
double precision:: zzmax,zzrms,dtacc,dfac
integer:: iz

!-------------------------------------------------------------------
! Compute vertical (relative) vorticity:
call vorticity(qq,pp,zz)
! Note: pp is available from a prior call to inversion.

! Compute maximum and rms vorticity:
zzmax=zero
zzrms=zero
do iz=1,nz
   zzmax=max(zzmax,maxval(abs(zz(:,:,iz))))
   zzrms=zzrms+hhat(iz)*sum(zz(:,:,iz)**2)
enddo
zzrms=danorm*zzrms
! Here, danorm = 1/ng^2.

! Save data to monitor.asc:
write(16,'(f12.5,2(1x,f15.7))') t,zzmax,zzrms

! Time step needed for accuracy:
dtacc=dtfac/max(zzmax,srwfm)
! The restriction on the maximum Rossby wave frequency (srwfm)
! ensures that the fastest Rossby wave frequency is resolved.

! Choose a new time step, limiting it to next gridded data save:
if (t > zero) then
   dt=min(dtacc,tgsave*(dble(igrids)+small)-t)
   ! small is added to ensure data are saved after this time step.
else
   dt=min(dtacc,tgsave*(one+small))
endif

! Time step fractions for 4th-order Runge-Kutta method:
dt2=f12*dt
dt3=f13*dt
dt6=f16*dt

! Increment the integral of max|zz|:
twist=twist+dt*zzmax

!---------------------------------------------------------------------
! Define spectral integrating factors used in Runge-Kutta integration:
dfac=dt2*zzrms
emq=exp(-dfac*diss)
epq=filt/emq
! diss & filt are defined in spectral.f90

!---------------------------------------------------------------------
! Set flags to save data:
gsave=(int(t/tgsave) == igrids)
! Gridded data will be saved at time t if gsave is true.
csave=(int(t/tcsave) == iconts)
! Contour data will be saved at time t if csave is true.

return
end subroutine adapt

!=======================================================================
      
subroutine savegrid

! Saves PV field at the desired save time, and computes and saves
! the kinetic and available potential energy per unit mass.

implicit none

! Local variables:
double precision:: tke,ape,tot
integer:: iz

!------------------------------------------------------------------
! Increment counter for direct file access:
igrids=igrids+1

! Update PV, streamfunction and velocity fields:
call inversion(1)

! Write PV field (single precision):
write(31,rec=igrids) real(t),real(qq)

! Write streamfunction (single precision):
write(32,rec=igrids) real(t),real(pp)

! Compute vertical vorticity for post-processing:
call vorticity(qq,pp,zz)

! Write vorticity (single precision):
write(33,rec=igrids) real(t),real(zz)

! Compute kinetic energy per unit mass:
tke=zero
do iz=1,nz
   tke=tke+hhat(iz)*sum((uu(:,:,iz)**2+vv(:,:,iz)**2))
enddo
tke=f12*danorm*tke

! Compute available potential energy per unit mass:
ape=zero
do iz=1,nz-1
   ape=ape+kdsq(iz)*sum((pp(:,:,iz)-pp(:,:,iz+1))**2)
enddo
ape=f12*danorm*ape
! Above, danorm = 1/ng^2.

! Total energy:
tot=tke+ape

! Write energy components and total:
write(15,'(f10.3,3(1x,f14.9))') t,tke,ape,tot

write(*,'(a,f10.3,3(a,f10.7))') &
    & ' t = ',t,'  K = ',tke,'  P = ',ape,'  K + P = ',tot

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont

! Saves PV contours and residual PV for fine-scale imaging:

implicit none

! Local variables
double precision:: qc(ng,ng,nz)
integer:: j
character(len=3):: pind

!---------------------------------------------------------------
! Increment counter for direct file access:
iconts=iconts+1

write(pind(1:3),'(i3.3)') iconts-1

! Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f9.2)') nq,nptq,t

! Write PV jumps for potentially restoring PV fields from contours
write(80,*) qjump

! Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)

! Save residual needed to build ultra-fine-grid vorticity with congen:
qc=qq-qc
write(83,rec=iconts) real(t),real(qc)

! Save PV contours:
open(81,file='cont/index'//pind,form='unformatted', &
      access='direct',status='replace',recl=16*max(1,nq))
write(81,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq),layq(1:nq)
close(81)

open(82,file='cont/nodes'//pind,form='unformatted', &
      access='direct',status='replace',recl=16*max(1,nptq))
write(82,rec=1) xq(1:nptq),yq(1:nptq)
close(82)

return
end subroutine savecont

!=======================================================================

! Main end module
end module evolution
