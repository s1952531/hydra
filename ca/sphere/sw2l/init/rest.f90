program rest

!===================================================================
!  Sets up PV = planetary vorticity in both layers (state of rest)
!===================================================================

! Import constants are parameters:
use constants

! Declarations:
implicit none

double precision:: qq(numLatGridPts,numLongGridPts), qAnom1(numLatGridPts)
integer:: i, j

!------------------------------------------------------------------
 !Define PV for a rest state:
do j = 1,numLatGridPts
  qAnom1(j) = fpole*sin((dble(j)-f12)*dl-hpi)
enddo

do i = 1,numLongGridPts
   qq(:,i) = qAnom1
enddo

 !Write initial PV field:
open(20,file = 'qq_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)
write(20,rec=1) zero, qq
write(20,rec=2) zero, qq
close(20)

end program rest
