!#########################################################################
!  Computes Ro = |zeta|_max/f, Fr = (|(u,v)|/c*sqrt(h))_max, h_min and
!  h_max.  
!  Writes ro-fr-hm.asc

!           Written 23/5/2018 by D G Dritschel @ St Andrews
!#########################################################################

program rofrhm

 !Import spectral module:
use spectral

implicit none

 !Physical arrays:
double precision:: uu(ng,ng),vv(ng,ng),hh(ng,ng)
double precision:: dd(ng,ng),zz(ng,ng),aa(ng,ng)

 !Spectral arrays:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng)

 !Other local variables:
double precision:: ro,fr,hmin,hmax,uio,vio
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
open(35,file='zz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)

 !Open output data file:

 !Open output file:
open(41,file='ro-fr-hm.asc',status='replace')

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
  read(35,rec=loop,iostat=iread) tr4,qqr4
  zz=dble(qqr4)

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

  !Compute Rossby number:
  ro=maxval(abs(zz))/cof

  !Compute Froude number:
  fr=sqrt(maxval((uu**2+vv**2)/(one+hh)))/cgw

  !Compute h_min and h_max:
  hmin=minval(hh)
  hmax=maxval(hh)

  !Write data:
  write(41,'(1x,f9.2,4(1x,f12.8))') tr4,ro,fr,hmin,hmax
enddo

 !Close files:
close(31)
close(32)
close(35)
close(41)

write(*,*)
write(*,*) ' Ro, Fr, h_min and h_max written to ro-fr-hm.asc'

 !End main program
end program
!=======================================================================
