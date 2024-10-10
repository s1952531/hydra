!#########################################################################
!  Computes gamma-tilde = gamma = 2J(u,v) - 2*delta^2 and writes gt.r4

!      Also computes the 1D spectrum and writes gt-spectra.asc

!          Revised 22/8/2019 by D G Dritschel @ St Andrews
!#########################################################################

program gammat

 !Import spectral module:
use spectral

implicit none

 !Physical arrays:
double precision:: uu(ng,ng),vv(ng,ng),hh(ng,ng)
double precision:: dd(ng,ng),zz(ng,ng),gg(ng,ng),aa(ng,ng)

 !Spectral arrays:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng)

 !Other local variables:
double precision:: gspec(0:ng),zspec,dspec,ogspec,dk
double precision:: t,uio,vio
real:: tr4,qqr4(ng,ng)
integer:: loop,iread,k

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Open input data files:
open(31,file='hh.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(32,file='dd.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(33,file='zz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(34,file='gg.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(50,file='spectra.asc',status='old')

 !Open output data files:
open(41,file='gt.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(60,file='gt-spectra.asc',status='replace')

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit 
  hh=dble(qqr4)

  read(32,rec=loop,iostat=iread) tr4,qqr4
  dd=dble(qqr4)

  read(33,rec=loop,iostat=iread) tr4,qqr4
  zz=dble(qqr4)

  read(34,rec=loop,iostat=iread) tr4,qqr4
  gg=dble(qqr4)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

  !Obtain velocity field from dd & zz:
  aa=dd
  call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
  wke=rlap*wke
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)

  aa=zz
  call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
  wke=rlap*wke
  call xderiv(ng,ng,hrkx,wke,wkc)
  call yderiv(ng,ng,hrky,wke,wkd)

  wka=wka-wkd
  wkb=wkb+wkc
  call spctop(ng,ng,wka,uu,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkb,vv,xfactors,yfactors,xtrig,ytrig)

  !Add mean flow (uio,vio):
  uio=-sum(hh*uu)*dsumi
  vio=-sum(hh*vv)*dsumi
  uu=uu+uio
  vv=vv+vio

  !Compute J(u,v):
  call jacob(uu,vv,aa)

  !Form gamma_tilde = gamma + 2*(J(u,v) - delta^2) (overwrite gamma = gg):
  gg=gg+two*(aa-dd**2)

  !De-alias:
  call dealias(gg)

  !Write data:
  write(41,rec=loop) tr4,real(gg)

  !Compute 1d spectrum:
  call ptospc(ng,ng,gg,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,gspec)
  gspec=spmf*gspec

  read(50,*)
  write(60,'(f9.2,1x,i5)') tr4,kmaxred
  do k=1,kmaxred
    read(50,*) dk,zspec,dspec,ogspec
    write(60,'(4(1x,f12.8))') alk(k),zspec,dspec,log10(gspec(k)+1.d-32)
  enddo
enddo

 !Close files:
close(31)
close(32)
close(33)
close(34)
close(41)
close(50)
close(60)

write(*,*)
write(*,*) ' *** gamma-tilde written to gt.r4'
write(*,*) ' Spectra of zeta, delta & gamma-tilde written to gt-spectra.asc'

 !End main program
end program gammat
!=======================================================================
