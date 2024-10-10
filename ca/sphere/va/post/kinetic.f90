program kinetic
!  ------------------------------------------------------------------
!  |   Computes the kinetic energy corresponding to the velocity    |
!  |   field available in uu.r4 and vv.r4.  Writes kk.r4            |
!  ------------------------------------------------------------------

use constants

implicit none

integer:: loop,iread
real:: t,uu(ng,nt),vv(ng,nt),kk(ng,nt)

!---------------------------------------------------------------
! Open files:
open(42,file='uu.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
open(43,file='vv.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
open(44,file='kk.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

!---------------------------------------------------------------
! Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(42,rec=loop,iostat=iread) t,uu
  if (iread .ne. 0) exit 
  read(43,rec=loop) t,vv
  write(*,'(a,f11.4)') ' Processing t = ',t
  kk=f12*(uu**2+vv**2)
  write(44,rec=loop) t,kk
enddo

close(42)
close(43)
close(44)

write(*,*)
write(*,*) ' The kinetic energy is available in kk.r4'
end program
