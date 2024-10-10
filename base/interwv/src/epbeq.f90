program epbeq

!  ==========================================================
!               Boussinesq Equilibria Finder
!  ==========================================================
!
!  Finds wavelike equilibria in the 2D Boussinesq equations in a
!  channel geometry (-L_x/2 <= x < L_x/2, periodic; 0 <= y <= L_y).
!  Here, "y" stands for the vertical coordinate, and below "v" 
!  stands for the vertical velocity component.

use constants
use wavesetup
use sta2dfft

implicit double precision (a-h,o-z)

!---------------------------------------------------------------
 !Parameters and array declarations:
integer,parameter:: nbytes=8+nx*(ny+1)*8
double precision,parameter:: widout=outsmoond*gly
 !Grid arrays:
double precision:: zz(ny,nx),zzold(ny,nx),psiold(ny,nx)
double precision:: zzpre(ny,nx),psipre(ny,nx),bfsloc(ny,nx)
double precision:: bb(ny,nx),uh(0:ny,nx),vh(ny,nx)
double precision:: bbar(0:ny),bfs(0:ny),yg(0:ny),ygf(0:nyf)
double precision:: output(0:ny,nx)
 !Spectral arays:
double precision:: swka(nx,ny),swkb(nx,0:ny),psispc(nx,ny)
double precision:: green(nx,0:ny)
 !Interpolation arrays:
double precision:: balp(0:nyf),bbet(0:nyf),bgam(0:nyf)
double precision:: bfsalp(0:nyf),bfsbet(0:nyf),bfsgam(0:nyf)
 !FFT arrays:
integer:: xfactors(5),yfactors(5)
double precision:: xtrig(2*nx),ytrig(2*ny)
double precision:: hrkx(nx),rkx(nx),hrky(ny)
 !Basic variables:
double precision:: da,amp,biga
!--------------------------------------------------------------------

 !Call into setup module to get weakly nonlinear initial wave
 !to use as an initial state in following iterations:
call wave_init

 !Set up bbar and bfs on the coarse grid:
do iy=0,ny
  iyf=mgf*iy
  bbar(iy)=bbarf(iyf)
  bfs( iy)=bfsf( iyf)
enddo

 !Set up the constants used in cubic interpolation of bfs & bbar:
do iy=0,nyf-1
  delb  =bbarf(iy+1)-bbarf(iy)
  delbfs= bfsf(iy+1)- bfsf(iy)
  balp(iy)  =glyf*bfsf(iy)
  bfsalp(iy)=glyf*dbfsfdy(iy)
  bgam(iy)  =glyf*(bfsf(iy)   +bfsf(iy+1))   -two*delb
  bfsgam(iy)=glyf*(dbfsfdy(iy)+dbfsfdy(iy+1))-two*delbfs
  bbet(iy)  =delb  -balp(iy)  -bgam(iy)
  bfsbet(iy)=delbfs-bfsalp(iy)-bfsgam(iy)
enddo

!----------------------------------------------
 !Initialise FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define x wavenumbers:
nxp2=nx+2
rkx(1)=zero
do kx=2,nwx
  wkx=hrkx(2*(kx-1))
  rkx(     kx)=wkx
  rkx(nxp2-kx)=wkx
enddo
rkx(nwx+1)=hrkx(nx)

 !Set up spectral Laplace inversion operator:
do ky=1,nym1
  do kx=1,nx
    green(kx,ky)=-one/(rkx(kx)**2+hrky(ky)**2)
  enddo
enddo

!----------------------------------------------------------------------

 !Work out r.m.s. psi and scaling to give desired amplitude:
psirms=zero
do ix=1,nx
  do iy=1,nym1
    psirms=psirms+psi(iy,ix)**2
  enddo
enddo
cpi=one/cp
amp=cpi*sqrt(psirms/dble(ny*nx))

write(*,'(a,1x,f12.6)') ' The calculated amplitude is:', amp
write(*,*) ' Enter the amplitude you wish to start with:'
read(*,*) amp

 !Set up increment in amplitude between equilibria:
write(*,*) ' Enter the increment in the r.m.s. displacement'
write(*,*) ' amplitude, da:'
read(*,*) da
      
scalef=amp*cp*sqrt(dble(ny*nx)/psirms)
write(*,*) 'The scale factor for that is: ',scalef
biga=dble(ny*nx)*amp**2

 !Define y on the coarse grid:
do iy=0,ny
  yg(iy)=gly*dble(iy)
enddo

 !Define y on the fine grid:
do iy=0,nyf
  ygf(iy)=glyf*dble(iy)
enddo

 !Scale psi by this:
do ix=1,nx
  do iy=1,nym1
    psi(iy,ix)=scalef*psi(iy,ix)
  enddo
enddo

 !Open files to record equilibrium solutions:
open(10,file='speed.asc',status='unknown')
open(21,file='pp.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(22,file='uu.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(23,file='vv.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(24,file='bb.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(25,file='zz.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

 !Define parameters for bfs smoothing out of domain:
const=sqrt(pi)/two
tautop=const*widout*bfsf(nyf)
taubot=const*widout*bfsf(0)

 !Set relaxation parameter:
relax=orelax
relaxc=one-relax

 !Open the input-params data file and write some values into it:
open(14,file='input-params',status='unknown')
write(14,'(1x,f16.12)') amp
write(14,'(1x,f16.12)') da
close(14)

 !Compute vorticity (zz) and velocity components (uh,vh):
cpisq=cpi**2
do ix=1,nx
  do iy=1,nym1
    zz(iy,ix)=-cpisq*bfs(iy)*psi(iy,ix)
  enddo
enddo

 !Start by writing out the linear guess for cp:
write(*,'(a,1x,1p,e14.7,a)') ' cp, rel. abs. error & rms error = ',cp,'        -               -'

 !Counter for number of states found:
nstate=0

 !Used to stop it niter > nitmax:
niter=0

!--------------------------------------------------------------------
 !Begin main loop to find a nonlinear solution <<<<<<<<<<<<
do while (niter .lt. nitmax) 
  niter=niter+1

   !Save cp in cpold to check convergence below:
  cpold=cp
   !Save psi, zz in psiold, zzold to check convergence below:
  do ix=1,nx
    do iy=1,nym1
      psiold(iy,ix)=psi(iy,ix)
      zzold(iy,ix)=zz(iy,ix)
    enddo
  enddo

   !Calculate new value of cp implied by the stream function:
  psisq=zero
  do ix=1,nx
    do iy=1,nym1
       psisq=psisq+psi(iy,ix)**2
    enddo
  enddo
  cp=sqrt(psisq/biga)
  cpi=one/cp
  cpisq=cpi**2
  
   !Define ybar(iy,ix) function for every point in domain, and then
   !define the vorticity (zz) everywhere by use of local N^2, which
   !is found by cubic interpolation:
  do ix=1,nx
    do iy=1,nym1
      ybar=yg(iy)+cpi*psi(iy,ix)
      if (ybar*(elly-ybar) .ge. zero) then 
         !In an ordinary streaming flow region:
        iyf=int(ybar*glyfi)
        p=(ybar-ygf(iyf))*glyfi
        bfsloc(iy,ix)=bfsf(iyf)+p*(bfsalp(iyf)+p*(bfsbet(iyf)+p*bfsgam(iyf)))
      else if (ybar.gt.elly) then            
         !Inside a closed re-circulation region at the top:
        bfsloc(iy,ix)=bfsf(nyf)*exp(-((ybar-elly)/widout)**2)            
      else           
         !Inside a closed re-circulation region at the bottom:
        bfsloc(iy,ix)=bfsf(0)*exp(-(ybar/widout)**2)                                  
      endif
      zz(iy,ix)=-bfsloc(iy,ix)*psi(iy,ix)*cpisq
    enddo
  enddo       

   !Apply relaxation:
  do ix=1,nx
    do iy=1,nym1
      zz(iy,ix)=relax*zz(iy,ix)+relaxc*zzold(iy,ix)
      psi(iy,ix)=zz(iy,ix)
    enddo
  enddo

   !Solve for psi by wavenumber division in Fourier space:
  call ptospc_fs(nx,ny,psi,psispc,xfactors,yfactors,xtrig,ytrig)
  do ky=1,nym1   
    do kx=1,nx
      swka(kx,ky)  =psispc(kx,ky)*green(kx,ky)
      psispc(kx,ky)=  swka(kx,ky)
    enddo
  enddo

   !Return psi to physical space:
  call spctop_fs(nx,ny,swka,psi,xfactors,yfactors,xtrig,ytrig)

   !Check error in psi:
  rmserr=zero
  abserr=zero
  psimax=zero
  psirms=zero
  do ix=1,nx
    do iy=1,nym1
      dpsi=psi(iy,ix)-psiold(iy,ix)
      psimax=max(psimax,abs(psi(iy,ix)))
      psirms=psirms+psi(iy,ix)**2
      abserr=max(abserr,abs(dpsi))
      rmserr=rmserr+dpsi**2
    enddo
  enddo
  abserr=abserr/psimax
  rmserr=sqrt(rmserr/psirms)
  
   !Set cp error:
  cperr=abs(cp-cpold)/abs(cp)

  write(*,'(a,4(1x,1p,e12.5))') ' cperr,rmserr,abserr,cp_nd=',cperr,rmserr,abserr,cp/cplong

  relax=orelax*(one-dble(niter)/dble(2*nitmax))
  relaxc=one-relax      

  if (niter .gt. nitmax) exit
  if (rmserr .gt. derrmax) exit
  if (rmserr .gt. tol) cycle


   !Past this point, we've converged!  Save data:

   !Define ybar(iy,ix) function for every point in domain, and then
   !define the buoyancy pert. everywhere by use of bbar(ybar), which
   !is found by cubic interpolation:
  do ix=1,nx
    do iy=1,nym1
      ybar=yg(iy)+cpi*psi(iy,ix)
      if (ybar*(elly-ybar) .ge. zero) then 
         !In an ordinary streaming flow region:
        iyf=int(ybar*glyfi)
        p=(ybar-ygf(iyf))*glyfi
        bb(iy,ix)=bbarf(iyf)+p*(balp(iyf)+p*(bbet(iyf)+p*bgam(iyf)))-bbar(iy)
      else if (ybar .gt. elly) then            
         !Inside a closed re-circulation region at the top:
         !Set bb to be consistent with bfs in the main loop above:
        bb(iy,ix)=bbarf(nyf)+tautop*erf((ybar-elly)/widout)-bbar(iy)        
      else
         !Inside a closed re-circulation region at the bottom:
         !Set bb to be consistent with bfs in the main loop above:
        bb(iy,ix)=bbarf(0)+taubot*erf(ybar/widout)-bbar(iy)           
      endif
    enddo
  enddo       

   !Calculate u = -dpsi/dy:
  call yderiv_fs(nx,ny,hrky,psispc,swkb)

   !Flip sign of derivative:
  do ky=0,ny
    do kx=1,nx
       swkb(kx,ky)=-swkb(kx,ky)
    enddo
  enddo

   !Calculate v = +dpsi/dx:
  call xderiv_fs(nx,ny,hrkx,psispc,swka)

   !Return swkb & vh to physical space as uh & vh:
  call spctop_fc(nx,ny,swkb,uh,xfactors,yfactors,xtrig,ytrig)
  call spctop_fs(nx,ny,swka,vh,xfactors,yfactors,xtrig,ytrig)

   !Increment the number of states found: 
  nstate=nstate+1

  write(*,'(a,i3,a,i5,a)') ' *** Found state, ',nstate,' in ',niter,' iterations.'

   !Prepare for next iteration:
  niter=0
  relax=orelax
  relaxc=one-relax

   !Output data to various files:
  write(10,'(2(1x,f14.9))') cp,cp/cplong

  do ix=1,nx
    output(0,ix)=zero
    do iy=1,nym1
      output(iy,ix)=psi(iy,ix)+cp*yg(iy)
    enddo
    output(ny,ix)=cp*elly
  enddo
  write(21,rec=nstate) amp,output
  write(22,rec=nstate) amp,uh
  do ix=1,nx
    output(0,ix)=zero
    do iy=1,nym1
      output(iy,ix)=vh(iy,ix)
    enddo
    output(ny,ix)=zero
  enddo
  write(23,rec=nstate) amp,output
  do ix=1,nx
    output(0,ix)=bbar(0)
    do iy=1,nym1
      output(iy,ix)=bb(iy,ix)+bbar(iy)
    enddo
    output(ny,ix)=bbar(ny)
  enddo
  write(24,rec=nstate) amp,output
  do ix=1,nx
    output(0,ix)=zero
    do iy=1,nym1
      output(iy,ix)=zz(iy,ix)
    enddo
    output(ny,ix)=zero
  enddo
  write(25,rec=nstate) amp,output

  !-----------------------------------------------------------------
   !Decide how much of the domain is filled and whether to continue: 
  vmin=vh(1,1)
  vmax=vh(1,1)
  do ix=1,nx
    do iy=1,nym1
      if (vh(iy,ix) .lt. vmin) then 
        ixmin=ix
        vmin=vh(iy,ix)
      endif
      if (vh(iy,ix) .gt. vmax) then 
        ixmax=ix
        vmax=vh(iy,ix)
      endif
    enddo
  enddo
   !If wave fills more than half of the domain then stop:
  if (abs(ixmax-ixmin) .gt. nx/2) exit

   !Increment the amplitude and look for another solution:
  amp=amp+da
  biga=dble(nx*ny)*amp**2

  if (nstate .eq. 1) then      
     !Scale eta and psi by scalef since previous solution not known:
    scalef=amp/(amp-da)
    cppre=cp
    do ix=1,nx
      do iy=1,nym1
        zzpre(iy,ix)=zz(iy,ix)
        zz(iy,ix)=scalef*zz(iy,ix)
        psipre(iy,ix)=psi(iy,ix)
        psi(iy,ix)=scalef*psi(iy,ix)
      enddo
    enddo
  else
     !Extrapolate cp, eta and zz from previous two states:
    cptmp=cp
    cp=two*cp-cppre
    cpi=one/cp
    cppre=cptmp

    do ix=1,nx
      do iy=1,ny
        zztmp=zz(iy,ix)
        zz(iy,ix)=two*zz(iy,ix)-zzpre(iy,ix)
        zzpre(iy,ix)=zztmp
        psitmp=psi(iy,ix)
        psi(iy,ix)=two*psi(iy,ix)-psipre(iy,ix)
        psipre(iy,ix)=psitmp
      enddo
    enddo
  endif
  
enddo
 !End of main do loop <<<<<<<<<<<<<<<<

 !Here, we can't find any more solutions; close files and stop:
write(*,*) ' Code finished normally'
close(10)
close(21)
close(22)
close(23)
close(24)
close(25)

end program
