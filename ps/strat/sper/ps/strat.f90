!#########################################################################
!                The Singly-Periodic Density-Stratified
!                 2D Boussinesq Pseudo-Spectral Method
!#########################################################################

!       Code adapted from ~/hydra/ps/plane/sw/caps codes in January 2020
!       by D G Dritschel @ St Andrews.

!       This code solves: 
!            zeta_t + (u,v)*grad(zeta) = db/dx + D[zeta]
!               b_t + (u,v)*grad(b) = D[b]
!               u_x + v_y = 0
!       in the domain xmin < x < xmax ; ymin < y < ymax
!       (free slip boundary conditions in y, periodic in x).
!       Here D is a hyperviscous dissipation operator.

!       Incompressibility is handled via the introduction of a stream-
!       function such that:
!             Lap(psi) = zeta

!       The full algorithm consists of the following modules:

!       strat.f90      : This source - main program to evolve fields;
!       parameters.f90 : User defined parameters for a simulation;
!       constants.f90  : Fixed constants used throughout the other modules;
!       spectral.f90   : Fourier transform common storage and routines.

!----------------------------------------------------------------------------
program strat

 !Import contants, parameters and common arrays:
use constants
use spectral

implicit none

 !Define common space:

 !Velocity field, vorticity & buoyancy (physical):
double precision:: uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
double precision:: zz(0:ny,0:nxm1),bb(0:ny,0:nxm1)

 !Prognostic fields and dissipation operator (spectral):
double precision:: zs(0:nxm1,0:ny),bs(0:nxm1,0:ny)
double precision:: diss(0:nxm1,0:ny)

 !Time, time step, etc:
double precision:: t,dt,dt4

 !Used for saving data approximately every tgsave units of time:
integer:: igrids

 !Logical for use in saving data:
logical:: gsave

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t < tsim)

   !Advect flow from time t to t + dt:
  call advance
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Save final data if not done already:
if (int(t/tgsave) .eq. igrids) then
  call main_invert(zs,uu,vv,zz)
  call adapt
  call savegrid
endif

 !Close all files:
call finalise


 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none

 !Local variables:
double precision:: qq(0:ny,0:nxm1)
double precision:: bbdif

!----------------------------------------------------------------------
 !Read buoyancy and find min/max value for use in determining viscosity:
open(11,file='bb_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qq
close(11)

bbdif=maxval(qq)-minval(qq)

 !Initialise inversion constants and arrays:
call init_spectral(bbdif)

 !Convert buoyancy (in qq) to spectral space as bs:
call ptospc_fc(nx,ny,qq,bs,xfactors,yfactors,xtrig,ytrig)
 !qq is overwritten here.

 !Read vorticity and convert to spectral space as zs:
open(11,file='zz_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qq
close(11)

call ptospc_fc(nx,ny,qq,zs,xfactors,yfactors,xtrig,ytrig)
 !qq is overwritten here.

!----------------------------------------------------------------------
 !Spectrally-truncate all fields for use in de-aliasing:
bs=filt*bs
zs=filt*zs

 !Obtain initial velocity field (uu,vv):
call main_invert(zs,uu,vv,zz)
 !Note: zs is in spectral space while uu, vv & zz are in physical space.

!----------------------------------------------------------------------
 !Open all plain text diagnostic files:
open(21,file='evolution/ecomp.asc',status='replace')
open(22,file='evolution/monitor.asc',status='replace')
open(23,file='evolution/vorticity.asc',status='replace')

 !Open files to save vorticity and buoyancy fields periodically:
open(31,file='evolution/zz.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(32,file='evolution/bb.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)

 !Open file for 1d vorticity & buoyancy spectra:
open(51,file='spectra/zspec.asc',status='replace') ! zeta
open(52,file='spectra/bspec.asc',status='replace') ! delta

 !Initialise counter for saving gridded data:
igrids=0

 !Logical used for saving gridded fields:
gsave=.false.

return
end subroutine initialise

!=======================================================================

subroutine advance

! Advances fields from time t to t+dt using an iterative implicit 
! trapezoidal method of the form
!
!     (F^{n+1}-F^n)/dt = (L^{n+1}+L^n)/2 + (S^{n+1}+S^n)/2
!
! for a field F, where n refers to the time level, L[F] refers to
! the linear dissipation terms (hyperdiffusion), and S[F] refers to
! the remaining source terms.

! We start with the guess S^{n+1} = S^n and iterate  niter  times
! (see parameter statement below).

implicit none

 !Number of iterations of above scheme:
integer,parameter:: niter=2

 !Spectral fields needed in time stepping:
double precision:: bsi(0:nxm1,0:ny),bsm(0:nxm1,0:ny),sbs(0:nxm1,0:ny)
double precision:: zsi(0:nxm1,0:ny),zsm(0:nxm1,0:ny),szs(0:nxm1,0:ny)

 !Other local quantities:
integer:: iter

!-------------------------------------------------------------------
 !Invert vorticity for velocity at current time level, say t=t^n:
call main_invert(zs,uu,vv,zz) !zs is spectral here

 !Adapt the time step and save various diagnostics each time step:
call adapt

 !Possibly save field data (gsave is set by adapt):
if (gsave) call savegrid

!------------------------------------------------------------------
 !Start with a guess for F^{n+1} for all fields:

 !Calculate the source terms (sbs,szs) for buoyancy (bs) and
 !vorticity (zs) in spectral space:
call source(sbs,szs)

 !Initialise iteration (dt = dt/4 below):
bsi=bs
bsm=bs+dt4*sbs
bs=diss*(bsm+dt4*sbs)-bsi
zsi=zs
zsm=zs+dt4*szs
zs=diss*(zsm+dt4*szs)-zsi
 !diss is related to the hyperdiffusive operator (see end of adapt)

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call main_invert(zs,uu,vv,zz) !zs is spectral here

   !Calculate the source terms (sbs,szs):
  call source(sbs,szs)

   !Update fields:
  bs=diss*(bsm+dt4*sbs)-bsi
  zs=diss*(zsm+dt4*szs)-zsi
enddo

 !Advance time:
t=t+dt

return
end subroutine advance

!=======================================================================

subroutine source(sbs,szs)

! Gets the source terms for vorticity and buoyancy in spectral space.

! The spectral fields bs and zs are all spectrally truncated.
! Note, uu and vv obtained by main_invert before calling this 
! routine are spectrally truncated as well.

implicit none

 !Passed variables (spectral):
double precision:: sbs(0:nxm1,0:ny),szs(0:nxm1,0:ny)

 !Local variables (physical):
double precision:: px(0:ny,0:nxm1),py(ny,0:nxm1)

 !Local variables (spectral):
double precision:: sx(0:nxm1,0:ny),sy(0:nxm1,ny)

!--------------------------------------------------------------
 !Buoyancy source bb_t = -(u,v)*grad(bb):

 !Obtain x & y derivatives of buoyancy -> px, py (physical):
call xderiv_fc(nx,ny,hrkx,bs,sx)
 !Store spectral db/dx in szs for use in vorticity source below:
szs=sx
call spctop_fc(nx,ny,sx,px,xfactors,yfactors,xtrig,ytrig)
call yderiv_fc(nx,ny,rky,bs,sy)
call spctop_fs(nx,ny,sy,py,xfactors,yfactors,xtrig,ytrig)

 !Compute (u,v)*grad(bb) -> px in physical space:
px(0,:)=uu(0,:)*px(0,:)
px(1:nym1,:)=uu(1:nym1,:)*px(1:nym1,:)+vv(1:nym1,:)*py(1:nym1,:)
px(ny,:)=uu(ny,:)*px(ny,:)

 !Convert to spectral space as sbs and apply de-aliasing filter:
call ptospc_fc(nx,ny,px,sbs,xfactors,yfactors,xtrig,ytrig)
sbs=-filt*sbs

!--------------------------------------------------------------
 !Vorticity source zz_t = bb_x - (u,v)*grad(zz):

 !Obtain x & y derivatives of vorticity:
call xderiv_fc(nx,ny,hrkx,zs,sx)
call spctop_fc(nx,ny,sx,px,xfactors,yfactors,xtrig,ytrig)
call yderiv_fc(nx,ny,rky,zs,sy)
call spctop_fs(nx,ny,sy,py,xfactors,yfactors,xtrig,ytrig)

 !Compute (u,v)*grad(zz) -> px in physical space:
px(0,:)=uu(0,:)*px(0,:)
px(1:nym1,:)=uu(1:nym1,:)*px(1:nym1,:)+vv(1:nym1,:)*py(1:nym1,:)
px(ny,:)=uu(ny,:)*px(ny,:)

 !Convert to spectral space as sx and apply de-aliasing filter:
call ptospc_fc(nx,ny,px,sx,xfactors,yfactors,xtrig,ytrig)
szs=szs-filt*sx

return
end subroutine source

!=======================================================================

subroutine adapt

! Adapts the time step and computes various diagnostics

implicit none

 !For defining the max strain & buoyancy frequency based time step:
double precision,parameter:: alpha=0.1d0
 !Note: EPIC-2D paper recommends alpha = 0.2 for ls-rk4 method

 !For controlling numerical stability (CFL_max <= 0.8 recommended):
double precision,parameter:: cflmax=0.8d0
double precision,parameter:: cflpf=cflmax*glmin

 !Local variables (physical):
double precision:: px(0:ny,0:nxm1),py(ny,0:nxm1)
double precision:: fp(0:ny,0:nxm1)

 !Local variables (spectral):
double precision:: sx(0:nxm1,0:ny),sy(0:nxm1,ny)
double precision:: fs(0:nxm1,0:ny)

 !Diagnostic quantities:
double precision:: bfmax,zzmax,zzrms,ggmax,uumax,dfac
double precision:: zztmp,zzl1,zzl2,zzch
integer:: ix,iy

!----------------------------------------------------------
 !Obtain x & y derivatives of buoyancy -> px, py (physical):
call xderiv_fc(nx,ny,hrkx,bs,sx)
call spctop_fc(nx,ny,sx,px,xfactors,yfactors,xtrig,ytrig)
call yderiv_fc(nx,ny,rky,bs,sy)
call spctop_fs(nx,ny,sy,py,xfactors,yfactors,xtrig,ytrig)

 !Compute px^2 + py^2 -> px in physical space:
px(0,:)=px(0,:)**2
px(1:nym1,:)=px(1:nym1,:)**2+py(1:nym1,:)**2
px(ny,:)=px(ny,:)**2

 !Maximum buoyancy frequency:
bfmax=sqrt(sqrt(maxval(px)))

 !Maximum vorticity:
px=zz**2
zzmax=sqrt(maxval(px))

 !R.m.s. vorticity:
zzrms=sqrt(dsumi*(f12*sum(px(0,:)+px(ny,:))+sum(px(1:nym1,:))))

 !Characteristic vorticity, <zz^2>/<|zz|> for |zz| > zz_rms:
zzl1=small
zzl2=zero
do ix=0,nxm1
  do iy=1,ny
    zztmp=f12*(zz(iy-1,ix)+zz(iy,ix))
    if (abs(zztmp) .gt. zzrms) then
      zzl1=zzl1+abs(zztmp)
      zzl2=zzl2+zztmp**2
    endif
  enddo
enddo
zzch=zzl2/zzl1

 !Compute x derivative of velocity components:
fp=uu
call ptospc_fc(nx,ny,fp,fs,xfactors,yfactors,xtrig,ytrig)
call xderiv_fc(nx,ny,hrkx,fs,sx)
call spctop_fc(nx,ny,sx,fp,xfactors,yfactors,xtrig,ytrig)
 !fp = u_x
px=vv
call ptospc_fc(nx,ny,px,fs,xfactors,yfactors,xtrig,ytrig)
call xderiv_fc(nx,ny,hrkx,fs,sx)
call spctop_fc(nx,ny,sx,px,xfactors,yfactors,xtrig,ytrig)
 !px = v_x

 !Strain rate squared, u_x^2 + (v_x - zz/2)^2:
fp=fp**2+(px-f12*zz)**2

 !Maximum strain rate:
ggmax=sqrt(maxval(fp))

 !Maximum speed:
uumax=sqrt(maxval(uu**2+vv**2))

 !Choose new time step:
dt=min(alpha/(ggmax+small),alpha/(bfmax+small),cflpf/(uumax+small),tsim-t)

 !Update value of dt/4:
dt4=dt/four

!---------------------------------------------------------------------
if (nnu .eq. 1) then
  !Update diffusion operator used in time stepping:
  dfac=dt/two
  diss=two/(one+dfac*hdis)
   !hdis = nu*(k_x^2+k_y^2) where nu is the viscosity coefficient
   !(see spectral.f90 and parameters.f90).
else
   !Update hyperdiffusion operator used in time stepping:
  dfac=zzch*dt/two
   !zzch is the characteristic vorticity defined above.
  diss=two/(one+dfac*hdis)
   !hdis = C*(K/K_max)^{2p} where K^2 = k_x^2+k_y^2, p is the order,
   !K_max is the maximum x or y wavenumber and C is a dimensionless
   !prefactor (see spectral.f90 and parameters.f90 where C = prediss).
endif

!---------------------------------------------------------------------
 !Save |u|_max, N_max and gamma_max to monitor.asc:
write(22,'(1x,f13.6,3(1x,1p,e14.7))') t,uumax,bfmax,ggmax

 !Save vorticity diagnostics to vorticity.asc:
write(23,'(1x,f13.6,3(1x,1p,e14.7))') t,zzmax,zzrms,zzch

!---------------------------------------------------------------------
 !Set flag to save data:
gsave=(int(t/tgsave) .eq. igrids)
 !Gridded data will be saved at time t if gsave is true.

return
end subroutine adapt

!=======================================================================

subroutine savegrid

! Saves fields, their spectra and the energy components

implicit none

 !Local variables (spectral):
double precision:: fs(0:nxm1,0:ny)
double precision:: spec(0:max(nx,ny))

!Other variables:
double precision:: ekin,epot,etot,yg
integer:: iy,k

!---------------------------------------------------------------
 !Increment counter for direct file access:
igrids=igrids+1

!---------------------------------------------------------------
 !Kinetic energy:
ekin=f12*garea*(f12*sum(uu(0,:)**2+uu(ny,:)**2)+ &
                sum(uu(1:nym1,:)**2+vv(1:nym1,:)**2))

 !Gridded buoyancy:
fs=bs
call spctop_fc(nx,ny,fs,bb,xfactors,yfactors,xtrig,ytrig)

 !Potential energy:
epot=f12*(ymin*sum(bb(0,:))+ymax*sum(bb(ny,:)))
do iy=0,ny
  yg=ymin+gly*dble(iy)
  epot=epot+yg*sum(bb(iy,:))
enddo
epot=-garea*epot

 !Total energy:
etot=ekin+epot

 !Write energies to ecomp.asc:
write(21,'(1x,f13.6,3(1x,1p,e14.7))') t,ekin,epot,etot
write(*,'(a,f13.6,a,1p,e14.7)') ' t = ',t,'  E_tot = ',etot

!---------------------------------------------------------------
 !Compute 1d vorticity spectrum and write:
call spec1d_fc(zs,spec)
spec(0:kmax)=log10(spmf(0:kmax)*spec(0:kmax)+1.d-32)
write(51,'(f13.6,1x,i5)') t,kmax
do k=1,kmax
  write(51,'(3(1x,f12.8))') alk(k),spec(k)
enddo

 !Compute 1d buoyancy spectrum and write:
call spec1d_fc(bs,spec)
spec(0:kmax)=log10(spmf(0:kmax)*spec(0:kmax)+1.d-32)
write(52,'(f13.6,1x,i5)') t,kmax
do k=1,kmax
  write(52,'(3(1x,f12.8))') alk(k),spec(k)
enddo

 !spmf takes into account uneven sampling of wavenumbers in each
 !shell [K-1/2,K+1/2], where K is proportional to k (see spectral.f90)

 !alk(k) = log_10(K)

!---------------------------------------------------------------
 !Write vorticity and buoyancy fields:
write(31,rec=igrids) real(t),real(zz)
write(32,rec=igrids) real(t),real(bb)
 !Note, zz is available from a previous call to adapt.

return
end subroutine

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(21)
close(22)
close(23)
close(31)
close(32)
close(51)
close(52)

return
end subroutine

 !End main program
end program strat

!=======================================================================
