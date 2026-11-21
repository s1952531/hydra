program random_top

!=============================================================
!  Sets up a random bottom topography field using a function
!  in module randomise.f90
!=============================================================

! Import necessary modules:
use randomise

! Declarations:
implicit none

double precision:: botTopo(numLatGridPts,numLongGridPts),bmax
integer:: ordbeg,ordend

!--------------------------------------------------------------------------
write(*,*) ' Enter the range in order, n_beg to n_end:'
read(*,*) ordbeg,ordend

write(*,*) ' Enter the maximum absolute amplitude of the topography'
write(*,*) ' scaled on the mean lower-layer depth:'
read(*,*) bmax

!Generate random topography scaled on the mean depth:
call generate_random_field(botTopo, bmax, ordbeg, ordend)

!Write data:
open(11,file = 'bb_init.r8',form = 'unformatted', &
     access = 'direct',status = 'replace',recl = 2*nbytes)
write(11,rec=1) zero,botTopo
close(11)

end program random_top
