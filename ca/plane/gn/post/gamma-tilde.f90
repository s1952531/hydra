!#########################################################################
!  Computes the gamma-tilde = gamma + 2J(u,v) - 2*delta^2.
!  Writes gt.r4

!           Written 22/5/2018 by D G Dritschel @ St Andrews
!#########################################################################

program gtilde

use common
use evolution

implicit none

 !Physical arrays:
double precision:: dd(ng,ng),gg(ng,ng)
double precision:: aa(ng,ng)

 !Spectral arrays:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng)

 !Other local variables:
real:: tr4,qqr4(ng,ng)
integer:: loop,iread,k

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Open input data files:
open(32,file='dd.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(33,file='gg.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(35,file='zz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)

 !Open output data file:

 !Open files for coarse grid saves:
open(41,file='gt.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(32,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit 
  dd=dble(qqr4)

  read(33,rec=loop,iostat=iread) tr4,qqr4
  gg=dble(qqr4)
  read(35,rec=loop,iostat=iread) tr4,qqr4
  zz=dble(qqr4)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

  !Obtain velocity field from dd & zz:
  aa=dd
  call ptospc(ng,ng,aa,ds,xfactors,yfactors,xtrig,ytrig)
  wke=rlap*ds
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

  !Form gamma-tilde:
  call jacob(uu,vv,aa)
  gg=gg+two*(aa-dd**2)

  !Write data:
  write(41,rec=loop) tr4,real(gg)
enddo

 !Close files:
close(32)
close(33)
close(35)
close(41)

write(*,*)
write(*,*) ' gamma-tilde written to gt.r4'

 !End main program
end program
!=======================================================================
