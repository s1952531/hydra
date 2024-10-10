!#########################################################################
!  Computes the scaled vertical velocity w/H where w = Dh/Dt at z = h
!  and writes ww.r4.  Note, w = -h*delta where delta is the horizontal
!  divergence.

!           Written 13/5/2019 by D G Dritschel @ St Andrews
!#########################################################################

program vvel

 !Import constants and parameters:
use constants

implicit none

 !Physical arrays:
double precision:: hh(ng,ng),dd(ng,ng),ww(ng,ng)

 !Other local variables:
real:: tr4,qqr4(ng,ng)
integer:: loop,iread

!----------------------------------------------------------------------
 !Open input data files:
open(31,file='hh.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(32,file='dd.r4',form='unformatted',access='direct', &
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

   !Vertical velocity scaled by H at z = h:
  ww=-dd*(one+hh)

   !Write data:
  write(41,rec=loop) tr4,real(ww)
enddo

 !Close files:
close(31)
close(32)
close(41)

write(*,*)
write(*,*) ' The vertical velocity w at z = 0 is written to ww.r4'

 !End main program
end program vvel
!=======================================================================
