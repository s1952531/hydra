!#########################################################################
! Generates r, (u,v), zeta, delta and gamma_tilde spectra, in three ways:

! (1) a theta average, then compute spectra
! (2) compute spectra for each layer, then average over theta (usual way)
! (3) Use differences between field and their theta averages,
!     compute spectra for each layer, then average over theta

! These are written to the files rspec.asc, kspec.asc, zspec.asc,
! dspec.asc and gspec.asc for r, (u,v), zeta, delta and gamma_tilde,
! respectively.

! Use dec_spec_view to image the data (allows one to select from above).

!          Written 24/8/2019 by D G Dritschel @ St Andrews
!#########################################################################

program decomp_spec

 !Import spectral module:
use spectral

implicit none

 !Various arrays needed below:
double precision:: u(ng,ng,0:nz),v(ng,ng,0:nz),r(ng,ng,0:nz)
double precision:: d(ng,ng,0:nz),g(ng,ng,0:nz),zeta(ng,ng,0:nz)
double precision:: qs(ng,ng,0:nz),ds(ng,ng,0:nz),gs(ng,ng,0:nz)

double precision:: r2d(ng,ng),u2d(ng,ng),v2d(ng,ng)
double precision:: d2d(ng,ng),zeta2d(ng,ng),gt2d(ng,ng)
double precision:: wkp(ng,ng),wka(ng,ng)

 !Other local variables:
double precision:: rspec1(0:ng),rspec2(0:ng),rspec3(0:ng)
double precision:: kspec1(0:ng),kspec2(0:ng),kspec3(0:ng)
double precision:: zspec1(0:ng),zspec2(0:ng),zspec3(0:ng)
double precision:: dspec1(0:ng),dspec2(0:ng),dspec3(0:ng)
double precision:: gspec1(0:ng),gspec2(0:ng),gspec3(0:ng)
double precision:: tmpspec(0:ng),dk
real:: tr4,q3dr4(ng,ng,0:nz)
integer:: loop,iread,iz,k

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open input data files:
open(31,file='3d/ql.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
open(32,file= '3d/d.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
open(33,file= '3d/g.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)

 !Open output files:
open(61,file='rspec.asc',status='replace')
open(62,file='kspec.asc',status='replace')
open(63,file='zspec.asc',status='replace')
open(64,file='dspec.asc',status='replace')
open(65,file='gspec.asc',status='replace')

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,q3dr4
  if (iread .ne. 0) exit 
  zeta=dble(q3dr4)
  call ptospc3d(zeta,qs,0,nz)

  read(32,rec=loop,iostat=iread) tr4,q3dr4
  if (iread .ne. 0) exit 
  d=dble(q3dr4)
  call ptospc3d(d,ds,0,nz)

  read(33,rec=loop,iostat=iread) tr4,q3dr4
  if (iread .ne. 0) exit 
  g=dble(q3dr4)
  call ptospc3d(g,gs,0,nz)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

   !Obtain velocity field by inversion:
  call main_invert(qs,ds,gs,r,u,v,zeta)
   !Note: qs, ds & gs are in spectral space while 
   !      r, u, v and zeta are in physical space.

   !Compute gamma_tilde (store back in the array g):
  gt2d=zero
  do iz=0,nz
    gt2d=gt2d+weight(iz)*g(:,:,iz)
    call jacob(u(:,:,iz),v(:,:,iz),wka)
    call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
    wkp=g(:,:,iz)+two*(wkp-d(:,:,iz)**2)
    call deal2d(wkp)
    g(:,:,iz)=wkp
  enddo

   !Compute theta-averaged fields:
  r2d=zero
  u2d=zero
  v2d=zero
  d2d=zero
  zeta2d=zero
!  gt2d=zero
  do iz=0,nz
    r2d=r2d+weight(iz)*r(:,:,iz)
    u2d=u2d+weight(iz)*u(:,:,iz)
    v2d=v2d+weight(iz)*v(:,:,iz)
    d2d=d2d+weight(iz)*d(:,:,iz)
    zeta2d=zeta2d+weight(iz)*zeta(:,:,iz)
!    gt2d=gt2d+weight(iz)*g(:,:,iz)
   enddo

   !Obtain theta-averaged gamma from theta-averaged fields:
  call jacob(u2d,v2d,wka)
  call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
  gt2d=gt2d+two*(wkp-d2d**2)
  call deal2d(gt2d)
   
   !Compute spectra of the theta-averaged fields:
  wkp=r2d
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,rspec1)
  rspec1=spmf*rspec1

  wkp=u2d
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,tmpspec)
  wkp=v2d
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,kspec1)
  kspec1=spmf*(kspec1+tmpspec)

  wkp=d2d
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,dspec1)
  dspec1=spmf*dspec1

  wkp=zeta2d
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,zspec1)
  zspec1=spmf*zspec1

  wkp=gt2d
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,gspec1)
  gspec1=spmf*gspec1

   !Compute vertically averaged spectra:
  rspec2=zero
  kspec2=zero
  dspec2=zero
  zspec2=zero
  gspec2=zero
  do iz=0,nz
    wkp=r(:,:,iz)
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    rspec2=rspec2+weight(iz)*tmpspec

    wkp=u(:,:,iz)
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    kspec2=kspec2+weight(iz)*tmpspec

    wkp=v(:,:,iz)
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    kspec2=kspec2+weight(iz)*tmpspec

    wkp=d(:,:,iz)
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    dspec2=dspec2+weight(iz)*tmpspec

    wkp=zeta(:,:,iz)
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    zspec2=zspec2+weight(iz)*tmpspec

    wkp=g(:,:,iz)
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    gspec2=gspec2+weight(iz)*tmpspec
  enddo
  rspec2=spmf*rspec2
  kspec2=spmf*kspec2
  dspec2=spmf*dspec2
  zspec2=spmf*zspec2
  gspec2=spmf*gspec2

   !Compute vertically averaged spectra of vertically-varying parts:
  rspec3=zero
  kspec3=zero
  dspec3=zero
  zspec3=zero
  gspec3=zero
  do iz=0,nz
    wkp=r(:,:,iz)-r2d
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    rspec3=rspec3+weight(iz)*tmpspec

    wkp=u(:,:,iz)-u2d
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    kspec3=kspec3+weight(iz)*tmpspec

    wkp=v(:,:,iz)-v2d
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    kspec3=kspec3+weight(iz)*tmpspec

    wkp=d(:,:,iz)-d2d
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    dspec3=dspec3+weight(iz)*tmpspec

    wkp=zeta(:,:,iz)-zeta2d
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    zspec3=zspec3+weight(iz)*tmpspec

    wkp=g(:,:,iz)-gt2d
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call spec1d(wka,tmpspec)
    gspec3=gspec3+weight(iz)*tmpspec
  enddo
  rspec3=spmf*rspec3
  kspec3=spmf*kspec3
  dspec3=spmf*dspec3
  zspec3=spmf*zspec3
  gspec3=spmf*gspec3

   !Write data:
  write(61,'(f9.2,1x,i5)') tr4,kmaxred
  write(62,'(f9.2,1x,i5)') tr4,kmaxred
  write(63,'(f9.2,1x,i5)') tr4,kmaxred
  write(64,'(f9.2,1x,i5)') tr4,kmaxred
  write(65,'(f9.2,1x,i5)') tr4,kmaxred
  do k=1,kmaxred
    write(61,'(4(1x,f12.8))') alk(k),log10(rspec1(k)+1.d-32), &
             log10(rspec2(k)+1.d-32),log10(rspec3(k)+1.d-32)
    write(62,'(4(1x,f12.8))') alk(k),log10(kspec1(k)+1.d-32), &
             log10(kspec2(k)+1.d-32),log10(kspec3(k)+1.d-32)
    write(63,'(4(1x,f12.8))') alk(k),log10(zspec1(k)+1.d-32), &
             log10(zspec2(k)+1.d-32),log10(zspec3(k)+1.d-32)
    write(64,'(4(1x,f12.8))') alk(k),log10(dspec1(k)+1.d-32), &
             log10(dspec2(k)+1.d-32),log10(dspec3(k)+1.d-32)
    write(65,'(4(1x,f12.8))') alk(k),log10(gspec1(k)+1.d-32), &
             log10(gspec2(k)+1.d-32),log10(gspec3(k)+1.d-32)
  enddo
enddo

 !Close files:
close(31)
close(32)
close(33)
close(61)
close(62)
close(63)
close(64)
close(65)

write(*,*)
write(*,*) ' Decomposed spectra of r, (u,v), zeta, delta & gamma_tilde'
write(*,*) ' written to rspec.asc, kspec.asc, zspec.asc, dspec.asc and'
write(*,*) ' gspec.asc.'
write(*,*)
write(*,*) ' Image using dec_spec_view'
write(*,*)

 !End main program
end program decomp_spec
!=======================================================================
