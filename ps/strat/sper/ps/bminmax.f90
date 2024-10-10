program bminmax

use constants

! This routine reads evolution/bb.r4 and computes the min & max values of b.
! Writes bminmax.asc in the evolution subdirectory.

implicit none

real:: t,bb(0:ny,0:nxm1)
integer:: loop,iread

!---------------------------------------------------------------------
 !Open input file:
open(32,file='evolution/bb.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)

 !Open output file:
open(55,file='evolution/bminmax.asc',status='replace')

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(32,rec=loop,iostat=iread) t,bb
  if (iread .ne. 0) exit

   !Save diagnostics:
  write(55,'(1x,f13.6,2(1x,1p,e14.7))') t,minval(bb),maxval(bb)
enddo

close(31)
close(32)
close(55)

write(*,*) ' t vs b_min and b_max is available in evolution/bminmax.asc'

end program bminmax
