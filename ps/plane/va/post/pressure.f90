!#########################################################################
!       Computes the non-hydrostatic pressure and writes pn.r4.

!      Also computes the 1D spectrum and writes pn-spectra.asc

!          Revised 22/8/2019 by D G Dritschel @ St Andrews
!#########################################################################

program pressnh

 !Import spectral module:
use spectral

implicit none

double precision,parameter:: toler=1.d-11
 !toler: maximum error in iteration below to find pressure P_0

 !Physical arrays:
double precision:: dd(ng,ng),zz(ng,ng),gg(ng,ng)
double precision:: uu(ng,ng),vv(ng,ng),hh(ng,ng)
double precision:: htot(ng,ng),hinv(ng,ng),hinvsq(ng,ng)
double precision:: ppn(ng,ng),rx(ng,ng),ry(ng,ng)
double precision:: pnx(ng,ng),pny(ng,ng),wkp(ng,ng),wkq(ng,ng)
double precision:: errpn

 !Spectral arrays:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng),wkf(ng,ng),wkg(ng,ng)

 !Other local variables:
double precision:: pspec(0:ng),zspec,dspec,gspec,dk
double precision:: uio,vio,errp0
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
  wkp=dd
  call ptospc(ng,ng,wkp,wke,xfactors,yfactors,xtrig,ytrig)
  wke=rlap*wke
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)

  wkp=zz
  call ptospc(ng,ng,wkp,wke,xfactors,yfactors,xtrig,ytrig)
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
  call jacob(uu,vv,wkp)

  !Complete definition of gamma-tilde (use gg array):
  gg=gg+two*(wkp-dd**2)

  !Multiply next by (1+h) and re-store in wkg (spectral) as the 
  !fixed rhs in the pressure iteration below:
  htot=one+hh
  wkp=htot*gg
  call ptospc(ng,ng,wkp,wkg,xfactors,yfactors,xtrig,ytrig)
  !wkg is not de-aliased by the prop operator below takes care of this.

  !Next calculate wkq = 3/H^2*(1/(1+h)^2 - 1) needed for the pressure
  !iteration below:
  hinv=one/htot
  call dealias(hinv)
  hinvsq=hinv**2
  call dealias(hinvsq)
  wkq=hbsq3i*(hinvsq-one)

  !Calculate also (1+h)^{-1}*grad(h) and store in rx & ry:
  wkp=hh
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call gradient(wka,rx,ry)
  rx=hinv*rx
  call dealias(rx)
  ry=hinv*ry
  call dealias(ry)

  !Now iterate to find \bar{P}_n (in ppn) starting from the guess
  !ppn = (grad^2 - 3/H^2)^{-1}((1+h)*gamma_tilde):
  wka=prop*wkg
  call spctop(ng,ng,wka,ppn,xfactors,yfactors,xtrig,ytrig)
  errpn=two*toler
  do while (errpn .gt. toler)
    wkp=ppn
    call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
    call gradient(wka,pnx,pny)
    wkp=rx*pnx+ry*pny+wkq*ppn
    call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
    wka=prop*(wkg+wkb)
    call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
    errpn=sqrt(dsumi*sum((wkp-ppn)**2))
    ppn=wkp
  enddo

  !Write data:
  write(41,rec=loop) tr4,real(ppn)

  !Compute 1d spectrum:
  call ptospc(ng,ng,ppn,wka,xfactors,yfactors,xtrig,ytrig)
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
