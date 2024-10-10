!#########################################################################
!       Computes the non-hydrostatic pressure and writes pn.r4.

!      Also computes the 1D spectrum and writes pn-spectra.asc

!          Revised 22/8/2019 by D G Dritschel @ St Andrews
!#########################################################################

program pressnh

 !Import spectral module:
use spectral

implicit none

 !Physical arrays:
double precision:: uu(ng,ng),vv(ng,ng),hh(ng,ng),htot(ng,ng)
double precision:: dd(ng,ng),zz(ng,ng),gg(ng,ng),pn(ng,ng)

 !Spectral arrays:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng)

 !Other local variables:
double precision:: pspec(0:ng),zspec,dspec,gspec,dk
double precision:: uio,vio
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

 !Open output data file:
open(41,file='pn.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(60,file='pn-spectra.asc',status='replace')

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
  pn=dd
  call ptospc(ng,ng,pn,wke,xfactors,yfactors,xtrig,ytrig)
  wke=rlap*wke
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)

  pn=zz
  call ptospc(ng,ng,pn,wke,xfactors,yfactors,xtrig,ytrig)
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
  call jacob(uu,vv,pn)

  !Complete definition of gamma-tilde (use pn array):
  pn=gg+two*(pn-dd**2)

  !Filter:
  call ptospc(ng,ng,pn,wke,xfactors,yfactors,xtrig,ytrig)
  wke=filt*wke
  call spctop(ng,ng,wke,pn,xfactors,yfactors,xtrig,ytrig)

  !Define 1+hh:
  htot=one+hh

  !Multiply gamma-tilde (in pn) by htot and filter:
  pn=pn*htot
  call ptospc(ng,ng,pn,wke,xfactors,yfactors,xtrig,ytrig)
  wke=filt*wke
  call spctop(ng,ng,wke,pn,xfactors,yfactors,xtrig,ytrig)

  !Square htot and filter:
  wkb=htot**2
  call ptospc(ng,ng,wkb,wke,xfactors,yfactors,xtrig,ytrig)
  wke=filt*wke
  call spctop(ng,ng,wke,wkb,xfactors,yfactors,xtrig,ytrig)

  !Multiply -H^2/3, htot*gamma-tilde and htot^2 to obtain nh pressure:
  pn=-hbsq3*pn*wkb
  call ptospc(ng,ng,pn,wke,xfactors,yfactors,xtrig,ytrig)
  wke=filt*wke
  call spctop(ng,ng,wke,pn,xfactors,yfactors,xtrig,ytrig)

  !Write data:
  write(41,rec=loop) tr4,real(pn)

  !Compute 1d spectrum:
  call ptospc(ng,ng,pn,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,pspec)
  pspec=spmf*pspec

  read(50,*)
  write(60,'(f9.2,1x,i5)') tr4,kmaxred
  do k=1,kmaxred
    read(50,*) dk,zspec,dspec,gspec
    write(60,'(4(1x,f12.8))') alk(k),zspec,dspec,log10(pspec(k)+1.d-32)
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
write(*,*) ' *** Non-hydrostatic pressure P_n written to pn.r4'
write(*,*) ' Spectra of zeta, delta & P_n written to pn-spectra.asc'

write(*,*)

 !End main program
end program pressnh
!=======================================================================
