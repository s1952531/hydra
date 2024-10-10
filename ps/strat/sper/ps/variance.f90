program variance

use constants

! This routine reads evolution/bb.r4 and zz.r4 and computes the variances
! <b^2> and <zeta^2> where < > indicates the domain integral.  Writes var.asc
! in the evolution subdirectory.

implicit none

double precision:: bb(0:ny,0:nxm1),zz(0:ny,0:nxm1)
double precision:: t,bbl2,zzl2
real:: tr4,qqr4(0:ny,0:nxm1)
integer:: loop,iread

!---------------------------------------------------------------------
 !Open input files:
open(31,file='evolution/zz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(32,file='evolution/bb.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)

open(55,file='evolution/variance.asc',status='replace')

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit

  t=dble(tr4)
  zz=dble(qqr4)

  read(32,rec=loop) tr4,qqr4

  bb=dble(qqr4)

   !Compute diagnostics:
  bbl2=garea*(f12*sum(bb(0,:)**2+bb(ny,:)**2)+sum(bb(1:nym1,:)**2))
  zzl2=garea*(f12*sum(zz(0,:)**2+zz(ny,:)**2)+sum(zz(1:nym1,:)**2))

   !Save diagnostics:
  write(55,'(1x,f13.6,2(1x,1p,e14.7))') t,bbl2,zzl2
enddo

close(31)
close(32)
close(55)

write(*,*) ' t vs <b^2> and <zeta^2> is available in evolution/variance.asc'

end program
