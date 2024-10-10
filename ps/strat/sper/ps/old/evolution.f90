module evolution

! Module contains subroutines to evolve bb and zz fields according to the 
! algorithm detailed in ps.f90.

use common

implicit none

 !Energies:
double precision:: ekepre,ekepost,apepre,apepost,dkdtpre,dkdtpost

 !Physical fields:
double precision:: sbb(0:ny,0:nxm1),bbpre(0:ny,0:nxm1),sbbpre(0:ny,0:nxm1)
double precision:: szz(0:ny,0:nxm1),zzpre(0:ny,0:nxm1),szzpre(0:ny,0:nxm1)
double precision:: uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)

 !Spectral arrays:
double precision:: zopi(0:nxm1,0:ny)

!Internal subroutine definitions (inherit global variables):

contains 

!=============================================================
subroutine advect

! Main subroutine for advecting fields

implicit none

 !Local variable:
integer:: igsave

!-----------------------------------------------------------------------
 !Obtain the initial velocity field:
call main_invert(zz,uu,vv)      

 !Calculate the initial source terms for bb & zz:
call tend(sbb,szz)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Adjust timestep (dt) on maximum vorticity magnitude:
  call adapt(igsave)

   !Advect buoyancy & vorticity from time t to t + dt:
  call advance

   !Update the time:
  t=t+dt

   !Possibly save buoyancy, vorticity & energy at chosen save time (tgrid):
  if (igsave .eq. 1) call savegrid

   !Copy new fields into previous time (for saving data):
  bbpre=bb
  sbbpre=sbb
  zzpre=zz
  szzpre=szz

enddo
 !End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 
return
end subroutine advect

!=======================================================================
      
subroutine advance

! Computes bb(t+dt) and zz(t+dt) by a standard pseudo-spectral scheme

implicit none

 !Number of iterations of trapezoidal rule
integer,parameter:: niter=2

integer:: iter

!------------------------------------------------------------------------
 !Calculate the source terms at time t:
call tend(sbb,szz)

 !Copy gridded fields & form tendencies at old time level (time t):
bbpre=bb
zzpre=zz
sbbpre=bb+qudt*sbb
szzpre=zz+qudt*szz
 !At the end of the subroutine bb & zz will contain the fields at t + dt

 !Iterate to converge on implicit trapezoidal integration:
do iter=1,niter
   !Obtain bb & zz at time t + dt:
  call ps_step

   !Obtain uu & vv at time t + dt:
  call main_invert(zz,uu,vv)

   !Update vorticity tendency:
  call tend(sbb,szz)
enddo

 !Obtain final corrected zs & zd at time t + dt:
call ps_step

 !Obtain final corrected uu & vv at time t + dt from zz:
call main_invert(zz,uu,vv)

return
end subroutine advance

!=======================================================================
      
subroutine ps_step

! Evolve bb & zz from t to t + dt using pseudo-spectral method.

implicit none

double precision:: wks(0:nxm1,0:ny)

bb=sbbpre+qudt*sbb
call ptospc_fc(nx,ny,bb,wks,xfactors,yfactors,xtrig,ytrig)
wks=zopi*wks
call spctop_fc(nx,ny,wks,bb,xfactors,yfactors,xtrig,ytrig)
bb=two*(sbbpre+qudt*sbb)-bbpre

zz=szzpre+qudt*szz
call ptospc_fc(nx,ny,zz,wks,xfactors,yfactors,xtrig,ytrig)
wks=zopi*wks
call spctop_fc(nx,ny,wks,zz,xfactors,yfactors,xtrig,ytrig)
zz=two*(szzpre+qudt*szz)-zzpre

return
end subroutine ps_step

!=======================================================================
      
subroutine tend

! Computes bb & zz tendencies (apart from hyperdiffusion)

implicit none

double precision:: wkx(0:ny,0:nxm1),wky(ny,0:nxm1)
double precision:: wks(0:nxm1,0:ny),wkt(0:nxm1,0:ny),wku(0:nxm1,ny)

integer:: kx,ky

!--------------------------------------------------------------
 !Buoyancy tendency:

 !copy bb to avoid overwriting:
wkx=bb

 !Obtain x & y derivatives of buoyancy:
call ptospc_fc(nx,ny,wkx,wks,xfactors,yfactors,xtrig,ytrig)
call xderiv_fc(nx,ny,hrkx,wks,wkt)
call spctop_fc(nx,ny,wkt,wkx,xfactors,yfactors,xtrig,ytrig)
call yderiv_fc(nx,ny,rky,wks,wku)
call spctop_fs(nx,ny,wku,wky,xfactors,yfactors,xtrig,ytrig)

 !Compute bb_t = -(u,v)*grad(bb):
sbb(0,:)=-uu(0,:)*wkx(0,:)
sbb(1:nym1,:)=-uu(1:nym1,:)*wkx(1:nym1,:)-vv(1:nym1,:)*wky(1:nym1,:)
sbb(ny,:)=-uu(ny,:)*wkx(ny,:)

 !Apply de-aliasing filter:
call ptospc_fc(nx,ny,sbb,wks,xfactors,yfactors,xtrig,ytrig)
do ky=0,ny
  do kx=0,nxm1
    wks(kx,ky)=wks(kx,ky)*dafx(kx)*dafy(ky)
  enddo
enddo
call spctop_fc(nx,ny,wks,sbb,xfactors,yfactors,xtrig,ytrig)

!--------------------------------------------------------------
 !Vorticity tendency:

 !Store linear source term (db/dx):
szz=wkx

 !copy zz to avoid overwriting:
wkx=zz

 !Obtain x & y derivatives of vorticity:
call ptospc_fc(nx,ny,wkx,wks,xfactors,yfactors,xtrig,ytrig)
call xderiv_fc(nx,ny,hrkx,wks,wkt)
call spctop_fc(nx,ny,wkt,wkx,xfactors,yfactors,xtrig,ytrig)
call yderiv_fc(nx,ny,rky,wks,wku)
call spctop_fs(nx,ny,wku,wky,xfactors,yfactors,xtrig,ytrig)

 !Compute zz_t = bb_x -(u,v)*grad(zz):
szz(0,:)=szz(0,:)-uu(0,:)*wkx(0,:)
szz(1:nym1,:)=szz(1:nym1,:)-uu(1:nym1,:)*wkx(1:nym1,:) &
                           -vv(1:nym1,:)*wky(1:nym1,:)
szz(ny,:)=szz(ny,:)-uu(ny,:)*wkx(ny,:)

 !Apply de-aliasing filter:
call ptospc_fc(nx,ny,szz,wks,xfactors,yfactors,xtrig,ytrig)
do ky=0,ny
  do kx=0,nxm1
    wks(kx,ky)=wks(kx,ky)*dafx(kx)*dafy(ky)
  enddo
enddo
call spctop_fc(nx,ny,wks,szz,xfactors,yfactors,xtrig,ytrig)

return
end subroutine tend

!=======================================================================

subroutine adapt(igsave)

! Adapts the time step

implicit none

 !For defining the max vorticity & buoyancy frequency based time step:
double precision,parameter:: dtfac=0.1 !alpha in EPIC paper
 !Keep CFL prefactor < 1 (0.8 recommended):
double precision,parameter:: glmin=min(glx,gly),cflpf=0.8d0*glmin

 !Work arrays:
double precision:: wkx(0:ny,0:nxm1),wky(ny,0:nxm1)
double precision:: wks(0:nxm1,0:ny),wkt(0:nxm1,0:ny),wku(0:nxm1,ny)
double precision:: wka(0:ny,0:nx),wkb(0:ny,0:nx)

double precision:: bfmax,zzmax,uumax,dissfac,cfl

!----------------------------------------------------------
 !Copy bb to avoid overwriting:
wkx=bb

 !Obtain x & y derivatives of buoyancy:
call ptospc_fc(nx,ny,wkx,wks,xfactors,yfactors,xtrig,ytrig)
call xderiv_fc(nx,ny,hrkx,wks,wkt)
call spctop_fc(nx,ny,wkt,wkx,xfactors,yfactors,xtrig,ytrig)
call yderiv_fc(nx,ny,rky,wks,wku)
call spctop_fs(nx,ny,wku,wky,xfactors,yfactors,xtrig,ytrig)

 !Maximum buoyancy frequency:
bfmax=sqrt(sqrt(max(maxval(wkx**2),maxval(wky**2))))

 !Maximum vorticity:
zzmax=maxval(abs(zz))+small

 !Maximum velocity:
uumax=sqrt(maxval(uu**2+vv**2))+small

 !Choose new time step:
dt=min(dtfac/zzmax,dtfac/bfmax,cflpf/uumax)

 !Update time step fractions:
hfdt=dt/two
qudt=dt/four

!---------------------------------------------------------------------
 !Hyper-viscous dissipation (proportional to |zz|_max)
dissfac=hfdt*zzmax   
do ky=0,ny
  do kx=0,nxm1
    zopi(kx,ky)=one/(one+dissfac*diss(kx,ky))
  enddo
enddo

!---------------------------------------------------------------------
 !Compute CFL number:
cfl=uumax*dt/glmin

 !Record cfl, max |z|, and max and l1,l2 norms of b to monitor.dat:
write(12,'(1x,f12.5,1x,f6.4,1x,1p,e14.7)') t,cfl,zzmax

!---------------------------------------------------------------------

 !Set flag to save gridded data every tgsave time units:
tgrid=tgsave*dble(int((t+dt)/tgsave))+small
 !small = 1.d-12 is added so that t = 0 is saved correctly.
if (t .lt. tgrid .and. t+dt .ge. tgrid) then
   !The save time is between t & t+dt; set flag to save data:
  igsave=1
  apepre=ape
  iene=1
   !Compute kinetic energy:
  call l2norm(uu,uul2)
  call l2norm(vv,vvl2)
  ekepre=f12*(uul2+vvl2)
   !Compute rate of change of kinetic energy:
  call binorm(vv,bb,dkdtpre)
else
   !Do not save data:
  igsave=0
endif

return
end subroutine adapt

!=======================================================================
      
subroutine savegrid

! Saves zz & bb at the desired save time to files

implicit none

double precision:: zspec(0:max(nx,ny)),bspec(0:max(nx,ny))
real:: qqr4(0:ny,0:nxm1),tr4

double precision:: pt,ptc,apot,ekin,dkdt
double precision:: zzl1,bbl1,zzl2,bbl2
double precision:: sumzspec,sumbspec
integer:: k

!---------------------------------------------------------------
 !Increment counter for direct file access:
igrids=igrids+1

 !Weights for time interpolation:
pt=(t-tgrid)/dt
ptc=one-pt

 !Compute available potential energy:
call getape(apepost)
 !Compute kinetic energy
call l2norm(uu,uul2)
call l2norm(vv,vvl2)
ekepost=f12*(uul2+vvl2)

 !Compute time interpolated energies:
apot=pt*apepre+ptc*apepost
ekin=pt*ekepre+ptc*ekepost

 !Compute rate of change of kinetic energy:
call binorm(vv,bb,dkdtpost)
dkdt=pt*dkdtpre+ptc*dkdtpost

 !Compute vorticity and buoyancy at save time:
tr4=real(tgrid)
qqr4=real(pt*zzpre+ptc*zz)
write(31,rec=igrids) tr4,qqr4
qqr4=real(pt*bbpre+ptc*bb)
write(32,rec=igrids) tr4,qqr4

 !Compute domain integrals of zeta, zeta^2, b and b^2:
call average(wka,zzl1)
zzl1=zzl1*domarea
call l2norm(wka,zzl2)
call average(wkb,bbl1)
bbl1=bbl1*domarea
call l2norm(wkb,bbl2)

 !Write diagnostics to the files norms.dat & ene.dat:
write(13,'(f7.2,4(1x,f14.9))') tgrid,bbl1,bbl2,zzl1,zzl2
write(15,'(f7.2,4(1x,f14.9))') tgrid,ekin+apot,ekin,apot,dkdt

write(*,'(a,f7.2,5(a,f10.7))') &
    & ' t = ',tgrid,'  b_1 = ',bbl1,'  z_1 = ',zzl1,'  K = ',ekin,'  P = ',apot,'  E = ',ekin+apot

 !Compute 1d vorticity & buoyancy spectrum:
call spec1d_fc(wka,zspec)
call spec1d_fc(wkb,bspec)
sumzspec=zero
sumbspec=zero
do k=0,kmax
  sumzspec=sumzspec+zspec(k)
  sumbspec=sumbspec+bspec(k)
   !Normalise to take into account uneven sampling of wavenumbers 
   !in each shell [k-1/2,k+1/2]:
  zspec(k)=spmf(k)*zspec(k)
  bspec(k)=spmf(k)*bspec(k)
enddo
sumzspec=8.d0*sumzspec*dsumi
sumbspec=8.d0*sumbspec*dsumi
 !Write out spectrum to file:
write(50,'(f7.2,3(1x,f14.9),1x,i5)') tgrid,sumzspec,sumbspec,zzl2,kmaxred
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(50,'(3(1x,f12.8))') alk(k),log10(zspec(k)),log10(bspec(k))
enddo
 !Note: alk(k) = log_10(k)

return
end subroutine savegrid

!=======================================================================

 !Main end module
end module evolution
