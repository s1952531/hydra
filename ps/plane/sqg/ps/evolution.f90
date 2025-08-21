module evolution

! Module contains subroutines to evolve the buoyancy field 
! according to the algorithm detailed in ps.f90.

use common

implicit none

 !Velocity field:
double precision:: uu(ny,nx),vv(ny,nx)

 !Spectral fields (note array order):
double precision:: emq(nx,ny),epq(nx,ny)

 !Energy & enstrophy for time interpolation:
double precision:: enepre,enspre

 !Logical to indicate data saves:
logical:: gsave

!Internal subroutine definitions (inherit global variables):
contains

!=============================================================
subroutine advect

! Main subroutine for advection

implicit none

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tfin)

   !Advect flow from time t to t + dt:
  call advance
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!Write final data; first invert buoyancy to compute velocity:
call inversion
call savegrid

return
end subroutine advect

!=======================================================================

subroutine init

! Initialises various quantities for the time integration

implicit none

!-------------------------------------------------------------
 !Logical used for saving gridded fields:
gsave=.false.

return
end subroutine init

!=======================================================================

subroutine advance

! Advances buoyancy from time t to t+dt by the pseudo-spectral method

! *** Uses a 4th-order Runge-Kutta method ***

implicit none

 !Local variables:

 !Spectral fields needed in Runge-Kutta time stepping (note array order):
double precision:: qsi(nx,ny),qsf(nx,ny),sqs(nx,ny)
double precision,allocatable,dimension(:,:):: csi,csf,scs

!------------------------------------------------------------------
 !RK4 predictor step to time t0 + dt/2:

 !Invert buoyancy and compute velocity:
call inversion

 !Possibly save data (gsave is set by adapt in the previous time step):
if (gsave) call savegrid

 !Adjust timestep (dt) on maximum vorticity magnitude or CFL:
call adapt

 !Calculate the source term (sqs) for buoyancy (qs):
call source(sqs,0)

qsi=qs
qs=emq*(qsi+dt2*sqs)
qsf=qsi+dt6*sqs

!------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt/2:
t=t+dt2

 !Invert buoyancy and compute velocity:
call inversion

 !Calculate the source term (sqs) for buoyancy (qs):
call source(sqs,1)

qs=emq*(qsi+dt2*sqs)
qsf=qsf+dt3*sqs

!------------------------------------------------------------------
 !RK4 predictor step at time t0 + dt:

 !Invert buoyancy and compute velocity:
call inversion

 !Calculate the source term (sqs) for buoyancy (qs):
call source(sqs,1)

emq=emq**2
qs=emq*(qsi+dt*sqs)
qsf=qsf+dt3*sqs

!------------------------------------------------------------------
 !RK4 corrector step at time t0 + dt:
t=t+dt2

 !Invert buoyancy and compute velocity:
call inversion

 !Calculate the source term (sqs) for buoyancy (qs):
call source(sqs,2)

qs=emq*(qsf+dt6*sqs)

return
end subroutine advance

!=======================================================================

subroutine source(sqs,lev)

! Computes the source term (sqs) for the buoyancy (qs)
! evolution equation (in spectral space).

implicit none

 !Passed variables:
double precision:: sqs(nx,ny)
integer:: lev

 !Local variables:
double precision:: qqx(ny,nx),qqy(ny,nx)
double precision:: wkp(ny,nx)

!---------------------------------------------------------------
 !qs source - only NL term is needed:
call gradient(qs,qqx,qqy)
wkp=-uu*qqx-vv*qqy
 !Convert to spectral space:
call ptospc(nx,ny,wkp,sqs,xfactors,yfactors,xtrig,ytrig)

!----------------------------------------------------------------
if (lev .eq. 0) then
 !Spectrally truncate source:
  sqs=filt*sqs
else if (lev .eq. 1) then
 !Apply exponential integrating factor:
  sqs=epq*sqs
else
  sqs=epq**2*sqs
endif

return
end subroutine source

!=======================================================================

subroutine inversion

! Inverts SQG operator on buoyancy to obtain the streamfunction (pp)
! ***in spectral space*** and the velocity (uu,vv) = (-dpp/dy,dpp/dx) 
! ***in physical space***

implicit none

 !Local variable (spectral streamfunction):
double precision:: pp(nx,ny)

!------------------------------------------------------------
 !Invert buoyancy to obtain the velocity field (uu,vv): 
call main_invert(qs,uu,vv,pp)

return
end subroutine inversion

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
double precision:: umax,zzrms,zzmax,dtacc,cfl,dfac,eif
integer:: itime

!----------------------------------------------------------------------
 !Compute the gridded relative vorticity (zz):
ss=qs*vorop
 !vorop = |k|, and is defined in spectral.f90
 !Above, qq = buoyancy in spectral space.
call spctop(nx,ny,ss,zz,xfactors,yfactors,xtrig,ytrig)

 !Compute accurate advection time step:
umax=sqrt(maxval(uu**2+vv**2))
zzmax=maxval(abs(zz))
zzrms=sqrt(dsumi*sum(zz**2))
dtacc=min(glx*cflmax/umax,dtfac/zzmax)

!---------------------------------------------------------------------
 !Choose a new time step, dt:
if (dt .gt. zero) then
  dt=min(dtacc,dtmax)
  if (dt .gt. dtacc) write(*,'(a,f9.5)') 'Warning! dt/dt_acc= ',dt/dtacc
else
   !Limit max timestep to a data save time
  dtmax=tgsave
  dt=min(dtacc,dtmax)
endif
 !Fractional time steps used in 4th-order Runge-Kutta time stepping:
dt2=dt*f12
dt3=dt*f13
dt6=dt*f16

!---------------------------------------------------------------------
 !Record various diagnostics to monitor.asc:
cfl=umax*dt/glx
write(17,'(1x,f9.5,1x,f12.7,1x,f14.7,1x,f12.7)') t,zzrms,zzmax,umax

!---------------------------------------------------------------------
 !Set flag to save gridded data every tgsave time units:
itime=int((t+dt)/tgsave)
tgrid=tgsave*dble(itime)+small
 !small = 1.d-12 is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  gsave=.true.

   !Copy current gridded fields for use in time interpolation:
  qspre=qs

   !Compute energy & enstrophy:
  call energy(enepre,enspre)
endif

!---------------------------------------------------------------------
 !Define spectral integrating factors used in Runge-Kutta integration:
if (fixed_viscosity) then
  ss=exp(-dt2*qdiss) !See notes at the end of constants.f90
else
  dfac=dt2*zzrms     !Viscosity varies in proportion to zeta_rms
  ss=exp(dfac*qdiss)
endif
epq=ss*filt
emq=one/ss

return
end subroutine adapt

!=======================================================================

subroutine energy(ene,ens)

! This routine computes the total energy & enstrophy

implicit none

 !Passed variables:
double precision:: ene,ens

 !Local variables:
double precision:: wka(ny,nx) !Physical
double precision:: qqs(nx,ny) !Spectral

!----------------------------------------------------------------------
 !Compute energy & enstrophy:
ene=f12*garea*sum(uu**2+vv**2)

qqs=qs
call spctop(nx,ny,qqs,wka,xfactors,yfactors,xtrig,ytrig)
ens=f12*garea*sum(wka**2)

return
end subroutine energy

!=======================================================================

subroutine savegrid

! Saves buoyancy, energy and various spectra at the desired save time

implicit none

 !Local variables:
double precision:: wka(ny,nx),wkb(ny,nx) !Physical
double precision:: qqs(nx,ny),zzs(nx,ny) !Spectral
double precision:: spec(0:max(nx,ny))
double precision:: pt,ptc,ene,ens
double precision:: enepost,enspost
integer:: k

!---------------------------------------------------------------
 !Increment counter for direct file access:
igrids=igrids+1

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Interpolate buoyancy at save time:
qqs=pt*qspre+ptc*qs

 !Compute vertical vorticity:
zzs=vorop*qqs

!---------------------------------------------------------------
 !Compute energy:
call energy(enepost,enspost)

 !Compute time interpolated energy & enstrophy:
ene=pt*enepre+ptc*enepost
ens=pt*enspre+ptc*enspost

 !Write energy & enstrophy to ene-ens.asc:
write(15,'(f9.5,2(1x,f16.9))') tgrid,ene,ens

!---------------------------------------------------------------
 !Compute 1d spectra for various fields:
call spec1d(qqs,spec)
spec=log10(spmf*spec+1.d-32)
write(51,'(f9.5,1x,i5)') tgrid,kmaxred
do k=1,kmaxred
  write(51,'(2(1x,f12.8))') alk(k),spec(k)
enddo

call spec1d(zzs,spec)
spec=log10(spmf*spec+1.d-32)
write(52,'(f9.5,1x,i5)') tgrid,kmaxred
do k=1,kmaxred
  write(52,'(2(1x,f12.8))') alk(k),spec(k)
enddo

!---------------------------------------------------------------
 !Write scaled buoyancy (b_0/N) to bb.r4:
call spctop(nx,ny,qqs,wka,xfactors,yfactors,xtrig,ytrig)
write(31,rec=igrids) real(tgrid),real(wka)

 !Write vertical vorticity to zz.r4:
call spctop(nx,ny,zzs,wka,xfactors,yfactors,xtrig,ytrig)
write(32,rec=igrids) real(tgrid),real(wka)

write(*,'(a,f9.5,2(a,f13.6))') ' t = ',tgrid,' enstrophy = ',ens, &
                                             ' energy = ',ene

 !Unset flag for saving data:
gsave=.false.

return
end subroutine savegrid

!=======================================================================

 !Main end module
end module evolution
