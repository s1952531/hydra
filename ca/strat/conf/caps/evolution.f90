module evolution

! Module contains subroutines to evolve bb and zz fields according to the 
! algorithm detailed in casl.f90.

use common

implicit none

 !Energies:
double precision:: ekepre,ekepost,epepre,epepost

 !Physical fields:
double precision:: zz(0:ny,0:nx),zc(0:ny,0:nx),szs(0:ny,0:nx),szd(0:ny,0:nx)
double precision:: uu(0:ny,0:nx),vv(0:ny,0:nx)
double precision:: zzpre(0:ny,0:nx),zspre(0:ny,0:nx),zdpre(0:ny,0:nx)
double precision:: bbpre(0:ny,0:nx),szspre(0:ny,0:nx),szdpre(0:ny,0:nx)

 !Spectral arrays:
double precision:: zopi(0:nx,0:ny)

!Internal subroutine definitions (inherit global variables):

contains 

!=============================================================
subroutine advect

! Main subroutine for advecting fields and contours

implicit none

 !Local variables:
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, the code stops 
!      to rebuild the contours in a separate memory space.
integer:: istep,ireg,igsave,icsave,ix,iy
double precision:: zzl1,zzl2

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

 !Counter used for resetting fields and performing contour surgery:
istep=0

 !Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tfin)

   !Perform surgery & field reset every nstep time steps:
  istep=mod(istep,nstep)+1

  if (istep .eq. nstep) then
    ireg=ireg+1

     !Don't continue if maximum number of regularisations reached;
     !it is time to recontour (pass control back to main program):
    if (ireg .eq. nregmax) then
       !Update average vorticity:
      call average(zz,zavg)
       !Prepare vorticity residual for recontouring:
      zd=zz-zc
      zs=zz
       !Compute contour interval for vorticity:
      call l1norm(zz,zzl1)
      call l2norm(zz,zzl2)
      zjump=(zzl2/zzl1)/dble(ncontz)
       !ncontz is set in parameters.f90. 
       !*** Exit module and go to recontouring ***
      return
    endif

     !Regularise the buoyancy contours (surgery + node redistribution):
    call surgery(xb,yb,nextb,indb,npb,i1b,i2b,nb,nptb)
    call surgery(xz,yz,nextz,indz,npz,i1z,i2z,nz,nptz)

     !Record contour complexity to complexity.dat:
    write(14,'(1x,f12.5,2(1x,i8,1x,i9))') t,nb,nptb,nz,nptz

     !Convert vorticity contours to gridded values (zc):
    call con2grid(zc,xz,yz,zjump,nextz,nptz)

     !Reset zs and zd:
    call reset(zc,zs,zd,zavg)

     !Copy gridded vorticity fields to old time level:
    zdpre=zd 
    zspre=zs 
  endif

   !Adapt timestep (dt) on maximum vorticity magnitude & buoyancy gradient:
  call adapt(igsave,icsave)

   !Advect buoyancy & vorticity from time t to t + dt:
  call advance

   !Update the time:
  t=t+dt

   !Possibly save buoyancy, vorticity & energy at chosen save time (tgrid):
  if (igsave .eq. 1) call savegrid

   !Possibly save contours and residual vorticity (zd) for post processing:
  if (icsave .eq. 1) call savecont
  
   !Copy new fields into previous time:
  szspre=szs
  szdpre=szd
  zdpre=zd
  zspre=zs
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

 !Work arrays:
double precision:: wka(0:ny,0:nx),wkb(0:ny,0:nx)

!-----------------------------------------------------------------------
 !Record contour complexity to complexity.dat:
write(14,'(1x,f12.5,2(1x,i8,1x,i9))') t,nb,nptb,nz,nptz

!---------------------------------------------------
 !Convert vorticity contours to gridded values (zc):
call con2grid(zc,xz,yz,zjump,nextz,nptz)

 !Define residual vorticity zd = zs-zc-F[zs-zc]
wka=zs-zc
wkb=wka

call filter(wka,0,2)
zd=wkb-wka

!-----------------------------------------------------
 !Copy gridded vorticity to old time level:
zzpre=zs
zdpre=zd
zspre=zs

!------------------------------------------------------
 !Get the initial velocity field (uu,vv):
if (zjump .gt. zero) then 
   !Invert Laplace operator on zs (which here contains
   !the full vorticity) and obtain velocity field: 
  call main_invert(zs,uu,vv)      
   !Note: zz is not needed yet (it is defined in combine).
else
   !Here vorticity is zero identically - no need to do inversion:
  uu=zero
  vv=zero 
endif

 !Calculate the source terms for zd and zs at time t:
call zztend(szd,szs,t)

return
end subroutine

!=======================================================================

subroutine inversion

! Inverts Laplace's operator on vorticity (zz) to obtain the 
! streamfunction (pp) and the velocity (uu,vv) = (-dpp/dy,dpp/dx).

implicit none

!-----------------------------------------------------------------------
 !Call con2grid to get updated contour vorticity (zc):
call con2grid(zc,xz,yz,zjump,nextz,nptz)

 !Combine fields to update zz with full field:
call combine(zz,zc,zs,zd,zavg)

 !Invert Laplace operator on zz spectrally and obtain velocity field: 
call main_invert(zz,uu,vv)      

return
end subroutine

!=======================================================================
      
subroutine advance

! Computes bb(t+dt) by contour advection and zz(t+dt) by a combination of 
! contour advection and a standard pseudo-spectral scheme (PS)

implicit none

 !Define local parameters and arrays:
integer,parameter:: niter=2
double precision:: ub(nptb),vb(nptb),xbm(nptb),ybm(nptb)
double precision:: uz(nptz),vz(nptz),xzm(nptz),yzm(nptz)
integer:: i,iter

!------------------------------------------------------------------------
 !Calculate the source term for zeta at time t:
call zztend(szd,szs,t)

 !Copy gridded vorticity & form vorticity tendency term at old time level (time t):
zdpre=zd
zspre=zs
szdpre=zd+qudt*szd
szspre=zs+qudt*szs
 !At the end of the subroutine zz will contain vorticity at t + dt

 !Prepare buoyancy contour evolution; get velocity on contour nodes:
call velint(uu,vv,xb,yb,ub,vb,nptb)
do i=1,nptb
  xbm(i)=min(xmax,max(xmin,xb(i)+hfdt*ub(i)))
  ybm(i)=min(ymax,max(ymin,yb(i)+hfdt*vb(i)))
  xb(i) =min(xmax,max(xmin,xb(i)+  dt*ub(i)))
  yb(i) =min(ymax,max(ymin,yb(i)+  dt*vb(i)))
enddo

 !Prepare vorticity contour evolution; get velocity on contour nodes:
if (nptz .gt. 0) then
  call velint(uu,vv,xz,yz,uz,vz,nptz)
  do i=1,nptz
    xzm(i)=min(xmax,max(xmin,xz(i)+hfdt*uz(i)))
    yzm(i)=min(ymax,max(ymin,yz(i)+hfdt*vz(i)))
    xz(i) =min(xmax,max(xmin,xz(i)+  dt*uz(i)))
    yz(i) =min(ymax,max(ymin,yz(i)+  dt*vz(i)))
  enddo
endif

 !Iterate to converge on implicit trapezoidal integration:
do iter=1,niter
   !Obtain zs & zd at time t + dt:
  call ps_step

   !Obtain zz & hence uu & vv at time t + dt:
  call inversion

   !Update the buoyancy contour nodes:
  call velint(uu,vv,xb,yb,ub,vb,nptb)
  do i=1,nptb
    xb(i)=min(xmax,max(xmin,xbm(i)+hfdt*ub(i)))
    yb(i)=min(ymax,max(ymin,ybm(i)+hfdt*vb(i)))
  enddo

   !Update the vorticity contour nodes:
  if (nptz .gt. 0) then
    call velint(uu,vv,xz,yz,uz,vz,nptz)
    do i=1,nptz
      xz(i)=min(xmax,max(xmin,xzm(i)+hfdt*uz(i)))
      yz(i)=min(ymax,max(ymin,yzm(i)+hfdt*vz(i)))
    enddo
  endif

   !Update vorticity tendency:
  call zztend(szd,szs,t+dt)
enddo

 !Obtain final corrected zs & zd at time t + dt:
call ps_step

 !Obtain final corrected uu & vv at time t + dt from zz:
call inversion

 !Update the buoyancy contour nodes:
call velint(uu,vv,xb,yb,ub,vb,nptb)
do i=1,nptb
  xb(i)=min(xmax,max(xmin,xbm(i)+hfdt*ub(i)))
  yb(i)=min(ymax,max(ymin,ybm(i)+hfdt*vb(i)))
enddo

 !Update the vorticity contour nodes:
if (nptz .gt. 0) then
  call velint(uu,vv,xz,yz,uz,vz,nptz)
  do i=1,nptz
    xz(i)=min(xmax,max(xmin,xzm(i)+hfdt*uz(i)))
    yz(i)=min(ymax,max(ymin,yzm(i)+hfdt*vz(i)))
  enddo
endif

 !Call con2grid to get updated contour vorticity (zc):
call con2grid(zc,xz,yz,zjump,nextz,nptz)

 !Combine fields to update zz with full field:
call combine(zz,zc,zs,zd,zavg)

return
end subroutine

!=======================================================================
      
subroutine ps_step

! Evolve zd & zs from t to t + dt using pseudo-spectral method.

implicit none

double precision:: zdbs(0:nx,0:ny)

!-----------------------------------------------------------------------
 !Include source S_zz in zd:
zd=szdpre+qudt*szd

 !Transform to spectral space:
call ptospc_cc(nx,ny,zd,zdbs,xfactors,yfactors,xtrig,ytrig)

 !Apply hyperviscous damping:
zdbs=zopi*zdbs

 !Return to physical space:
call spctop_cc(nx,ny,zdbs,zd,xfactors,yfactors,xtrig,ytrig)

 !Finish:
zd=two*zd-zdpre

 !zs evolution (needs only advective source term and no damping):
zs=two*(szspre+qudt*szs)-zspre

return
end subroutine

!=======================================================================
      
subroutine zztend(vard,vars,time)

! Get vorticity local rate-of-change \pa{zs}/\pa{t} & \pa{zd}/\pa{t}:

implicit none

 !Passed arrays and variables:
double precision:: vard(0:ny,0:nx),vars(0:ny,0:nx)
double precision:: time

 !Local arrays and variables:
double precision:: ztmp(0:ny,0:nx)
double precision:: wks(0:nx,0:ny)
double precision:: dzdxs(nx,0:ny),dzdxp(0:ny,nx)
double precision:: dzdys(0:nx,ny),dzdyp(ny,0:nx)
integer:: kx,ky,ix,iy

!-----------------------------------------------------------------------
 !Get linear source term (Dzz/Dt):
call getzzsrc(vard,time,yox,yoy)

 !copy zd into ztmp to avoid overwriting zd:
ztmp=zd

 !Obtain x & y derivatives of vorticity:
call ptospc_cc(nx,ny,ztmp,wks,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
do ky=0,ny
  do kx=0,nx
    wks(kx,ky)=wks(kx,ky)*dafx(kx)*dafy(ky)
  enddo
enddo
call xderiv_cc(nx,ny,rkx,wks,dzdxs)
call spctop_sc(nx,ny,dzdxs,dzdxp,xfactors,yfactors,xtrig,ytrig)
call yderiv_cc(nx,ny,rky,wks,dzdys)
call spctop_cs(nx,ny,dzdys,dzdyp,xfactors,yfactors,xtrig,ytrig)

 !Subtract u.grad(zd) from the linear source term to get \pa{zd}/\pa{t}:
do ix=1,nxm1
  do iy=1,nym1
    vard(iy,ix)=vard(iy,ix)-uu(iy,ix)*dzdxp(iy,ix)-vv(iy,ix)*dzdyp(iy,ix)
  enddo
enddo

 !Deal with edge values separately:
do ix=1,nxm1
  vard(0, ix)=vard(0 ,ix)-uu(0 ,ix)*dzdxp(0 ,ix)
  vard(ny,ix)=vard(ny,ix)-uu(ny,ix)*dzdxp(ny,ix)
enddo
do iy=1,nym1
  vard(iy, 0)=vard(iy, 0)-vv(iy, 0)*dzdyp(iy, 0)
  vard(iy,nx)=vard(iy,nx)-vv(iy,nx)*dzdyp(iy,nx)
enddo
 !Note corner values are unchanged since uu=vv=0 there

!------------------------------------------------------
 !copy zs into ztmp to avoid overwriting zs:
ztmp=zs

 !Obtain x & y derivatives of vorticity:
call ptospc_cc(nx,ny,ztmp,wks,xfactors,yfactors,xtrig,ytrig)
 !Apply de-aliasing filter:
do ky=0,ny
  do kx=0,nx
    wks(kx,ky)=wks(kx,ky)*dafx(kx)*dafy(ky)
  enddo
enddo
call xderiv_cc(nx,ny,rkx,wks,dzdxs)
call spctop_sc(nx,ny,dzdxs,dzdxp,xfactors,yfactors,xtrig,ytrig)
call yderiv_cc(nx,ny,rky,wks,dzdys)
call spctop_cs(nx,ny,dzdys,dzdyp,xfactors,yfactors,xtrig,ytrig)

 !Subtract u.grad(zs) from the linear source term to get \pa{zs}/\pa{t}:
do ix=1,nxm1
  do iy=1,nym1
    vars(iy,ix)=-uu(iy,ix)*dzdxp(iy,ix)-vv(iy,ix)*dzdyp(iy,ix)
  enddo
enddo

 !Deal with edge values separately:
do ix=1,nxm1
  vars(0, ix)=-uu(0 ,ix)*dzdxp(0 ,ix)
  vars(ny,ix)=-uu(ny,ix)*dzdxp(ny,ix)
enddo
do iy=1,nym1
  vars(iy, 0)=-vv(iy, 0)*dzdyp(iy, 0)
  vars(iy,nx)=-vv(iy,nx)*dzdyp(iy,nx)
enddo

 !Set corner values to zero since uu=vv=0 there:
vars(0,  0)=zero
vars(ny, 0)=zero
vars(0, nx)=zero
vars(ny,nx)=zero

return
end subroutine

!=======================================================================

subroutine adapt(igsave,icsave)

! Adapts the time step to ensure dt < dtfac/max(|zeta|_max,sqrt|S_zz|_max)

implicit none

 !Passed variables:
integer:: igsave,icsave

 !Local variables:
double precision:: wka(0:ny,0:nx),wkb(0:ny,0:nx)
double precision:: zzmax,uumax,dtacc,dtcfl
double precision:: dbmax,dissfac,cfl,tcont
integer:: itime

 !Maximum CFL number:
double precision,parameter:: cflmax=0.7d0,glmin=min(glx,gly)

!-----------------------------------------------------------------------
 !Compute accurate advection time step (dtacc) and a stable time step (dtcfl):
zzmax=maxval(abs(zz))
uumax=maxval(uu**2+vv**2)
dtacc=dtfac/(zzmax+small)
uumax=sqrt(uumax)
dtcfl=cflmax*glmin/(uumax+small)

!---------------------------------------------------------------------
 !Choose a new time step: 
dt=min(dtacc,dtcfl,dtmax)
if (dt .gt. dtacc) write(*,'(a,f9.5)') 'Warning! dt/dt_acc= ',dt/dtacc
 !Define substeps:
hfdt=dt/two
qudt=dt/four

!---------------------------------------------------------------------
 !Define spectral diffusion operator:
dissfac=hfdt*zzmax
zopi=one/(one+dissfac*diss)

!---------------------------------------------------------------------
 !Compute CFL number:
cfl=uumax*dt/glmin

 !Record cfl, max |z|, and max and l1,l2 norms of b to monitor.dat:
write(12,'(1x,f12.5,1x,f6.4,1x,1p,e14.7)') t,cfl,zzmax

!---------------------------------------------------------------------
 !Set flag to save gridded data every tgsave time units:
itime=int((t+dt)/tgsave)
tgrid=tgsave*dble(itime)+small
 !small = 1.d-12 is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  igsave=1 
   !Copy zz into zzpre for output:
  zzpre=zz
   !Compute current gridded buoyancy and store in bbpre:
  call con2grid(bbpre,xb,yb,bjump,nextb,nptb)
   !Restore correct average:
  call restore(bbpre,bavg)
   !Compute kinetic energy:
  call kinetic(uu,vv,ekepre)
   !Compute potential energy:
  call potential(bbpre,epepre)
else
   !Do not save data:
  igsave=0
endif

 !Set flag to save contour data every tgsave time units:
itime=int((t+dt)/tcsave)
tcont=tcsave*dble(itime)-small
if (t .lt. tcont .and. t+dt .ge. tcont) then
   !The save time is between t & t+dt; set flag to save data:
  icsave=1
else
   !Do not save data:
  icsave=0
endif

return
end subroutine

!=======================================================================
      
subroutine savegrid

! Saves zz & bb at the desired save time to files (process with image):

implicit none

 !Local variables:
double precision:: wka(0:ny,0:nx),wkb(0:ny,0:nx)
double precision:: pt,ptc,epot,ekin
double precision:: zzl1,zzl2,bbl1,bbl2

!-----------------------------------------------------------------------
 !Increment counter for direct file access:
igrids=igrids+1

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Update gridded buoyancy field:
call con2grid(bb,xb,yb,bjump,nextb,nptb)

 !Restore correct average:
call restore(bb,bavg)

 !Compute kinetic energy
call kinetic(uu,vv,ekepost)

 !Compute potential energy:
call potential(bb,epepost)

 !Compute time interpolated energies:
epot=pt*epepre+ptc*epepost
ekin=pt*ekepre+ptc*ekepost

 !Store vorticity and buoyancy at save time:
wka=pt*zzpre+ptc*zz
wkb=pt*bbpre+ptc*bb

 !Write gridded fields to file:
write(31,rec=igrids) real(tgrid),real(wka)
write(32,rec=igrids) real(tgrid),real(wkb)

 !Compute domain integrals of zeta, zeta^2, b and b^2:
call average(wka,zzl1)
zzl1=zzl1*domarea
call l2norm(wka,zzl2)
call average(wkb,bbl1)
bbl1=bbl1*domarea
call l2norm(wkb,bbl2)

 !Write diagnostics to the files norms.dat & ene.dat:
write(13,'(f7.2,4(1x,f14.9))') tgrid,bbl1,bbl2,zzl1,zzl2
write(15,'(f7.2,3(1x,f14.9))') tgrid,ekin+epot,ekin,epot

write(*,'(a,f7.2,4(a,f10.7))') &
    & ' t = ',tgrid,'  b_1 = ',bbl1,'  z_1 = ',zzl1,'  K = ',ekin,'  P = ',epot

return
end subroutine

!=======================================================================
      
subroutine savecont

! Saves zz & bb contours and residual zz for post-processing via
! congen.f90

implicit none

 !Local variables
double precision:: wka(0:ny,0:nx)
integer:: iop(max(nz,nb))
integer:: irec,j
character(len=3):: pind

!-----------------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
irec=nint(t/tcsave)
write(pind(1:3),'(i3.3)') irec

 !Write contours to the cont subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,2(1x,f16.12))') nz,nptz,t,zjump,zavg
write(90,'(i8,1x,i9,1x,f12.5,2(1x,f16.12))') nb,nptb,t,bjump,bavg

 !Save residual needed to build ultra-fine-grid vorticity with congen:
wka=zz-zc
 !Ensure zero average:
call restore(wka,zero)
write(83,rec=irec) real(t),real(wka)

 !Save vorticity contours if any exist:
if (nz .gt. 0) then
   !First form iop; open/closed indicator:
  do j=1,nz
    iop(j)=nextz(i2z(j))/i1z(j)
     !iop = 0 for an open contour, and 1 for a closed one
  enddo
  open(81,file='cont/zzindex'//pind,form='unformatted', &
      & access='direct',status='replace',recl=16*nz)
  write(81,rec=1) npz(1:nz),i1z(1:nz),indz(1:nz),iop(1:nz)
  close(81)

  open(82,file='cont/zznodes'//pind,form='unformatted', &
      & access='direct',status='replace',recl=16*nptz)
  write(82,rec=1) xz(1:nptz),yz(1:nptz)
  close(82)
endif

 !Save buoyancy contours:
 !First form iop; open/closed indicator:
do j=1,nb
  iop(j)=nextb(i2b(j))/i1b(j)
   !iop = 0 for an open contour, and 1 for a closed one
enddo
open(91,file='cont/bbindex'//pind,form='unformatted', &
    & access='direct',status='replace',recl=16*nb)
write(91,rec=1) npb(1:nb),i1b(1:nb),indb(1:nb),iop(1:nb)
close(91)

open(92,file='cont/bbnodes'//pind,form='unformatted', &
    & access='direct',status='replace',recl=16*nptb)
write(92,rec=1) xb(1:nptb),yb(1:nptb)
close(92)

return
end subroutine

!=======================================================================

 !Main end module
end module
