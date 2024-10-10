!#########################################################################
!  Computes u & v and writes uu.r4 and vv.r4

!  Written 5/8/2019 by D G Dritschel @ New York
!#########################################################################

program velocity

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
open(35,file='zz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)

 !Open output files:
open(41,file='uu.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(42,file='vv.r4',form='unformatted',access='direct', &
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

  !Write data:
  write(41,rec=loop) tr4,real(uu)
  write(42,rec=loop) tr4,real(vv)
enddo

 !Close files:
close(31)
close(32)
close(35)
close(41)
close(42)

write(*,*)
write(*,*) ' u & v written to uu.r4 and vv.r4'

 !End main program
end program velocity
!=======================================================================
