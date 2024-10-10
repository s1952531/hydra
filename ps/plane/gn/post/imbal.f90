!#########################################################################
!  Computes the imbalanced height, vorticity, divergence and 
!  acceleration divergence given the original and balanced fields.
!  Writes ihh.r4, izz.r4, idd.r4 and igg.r4.

!  *** Must run dgbal first ***

!           Written 8/5/2018 by D G Dritschel @ St Andrews
!#########################################################################

program imbal

 !Import constants and parameters:
use constants

implicit none

 !Physical arrays:
double precision::  hh(ng,ng), zz(ng,ng), dd(ng,ng), gg(ng,ng)
double precision:: bhh(ng,ng),bzz(ng,ng),bdd(ng,ng),bgg(ng,ng)

 !Other local variables:
real:: tr4,qqr4(ng,ng)
integer:: loop,iread

!----------------------------------------------------------------------
 !Open input data files:
open(32,file='dd.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(33,file='gg.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(34,file='hh.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(35,file='zz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(42,file='bdd.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(43,file='bgg.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(44,file='bhh.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(45,file='bzz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)

 !Open output data files:
open(52,file='idd.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(53,file='igg.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(54,file='ihh.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(55,file='izz.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  !Read original fields:
  read(32,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit 
  dd=dble(qqr4)
  read(33,rec=loop,iostat=iread) tr4,qqr4
  gg=dble(qqr4)
  read(34,rec=loop,iostat=iread) tr4,qqr4
  hh=dble(qqr4)
  read(35,rec=loop,iostat=iread) tr4,qqr4
  zz=dble(qqr4)

  !Read balanced fields:
  read(42,rec=loop,iostat=iread) tr4,qqr4
  bdd=dble(qqr4)
  read(43,rec=loop,iostat=iread) tr4,qqr4
  bgg=dble(qqr4)
  read(44,rec=loop,iostat=iread) tr4,qqr4
  bhh=dble(qqr4)
  read(45,rec=loop,iostat=iread) tr4,qqr4
  bzz=dble(qqr4)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

  !Compute differences and save data:
  write(52,rec=loop) tr4,real(dd-bdd)
  write(53,rec=loop) tr4,real(gg-bgg)
  write(54,rec=loop) tr4,real(hh-bhh)
  write(55,rec=loop) tr4,real(zz-bzz)

enddo

 !Close files:
close(32)
close(33)
close(34)
close(35)

close(42)
close(43)
close(44)
close(45)

close(52)
close(53)
close(54)
close(55)

write(*,*)
write(*,*) ' Imbalanced fields ready in ihh.r4, izz.r4, idd.r4 & igg.r4'

 !End main program
end program
!=======================================================================
