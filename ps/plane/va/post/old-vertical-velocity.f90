!#########################################################################
!  Computes the vertical velocity w = Dh/Dt at z = h and writes ww.r4

!           Written 13/5/2019 by D G Dritschel @ St Andrews
!#########################################################################

program vvel

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
double precision:: uio,vio
real:: tr4,qqr4(ng,ng)
integer:: loop,iread

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

 !Open output data file:
open(41,file='ww.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)

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
 
  !Compute divergence of (htot*uu,htot*vv):
  hh=one+hh
  uu=hh*uu
  vv=hh*vv
  call divs(uu,vv,wka)

  !Return to physical space and switch sign to obtain w:
  call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)

  !Write data:
  write(41,rec=loop) tr4,real(-aa)
enddo

 !Close files:
close(31)
close(32)
close(33)
close(34)
close(41)

write(*,*)
write(*,*) ' The vertical velocity w at z = 0 is written to ww.r4'

 !End main program
end program vvel
!=======================================================================
