program random_pv

!=========================================================
!  Sets up a random baroclinic PV field using a function
!  in module randomise.f90
!=========================================================

 !Import necessary modules:
use randomise

 !Declarations:
implicit none

double precision:: qq(numLatGridPts,numLongGridPts), qAnom2(numLatGridPts,numLongGridPts), qAnom1(numLatGridPts), qmax
integer:: ordbeg, ordend, i, j

!--------------------------------------------------------------------------
write(*,*) ' Enter the range in order, n_beg to n_end:'
read(*,*) ordbeg,ordend

write(*,*) ' Enter the maximum absolute amplitude of the baroclinic PV:'
read(*,*) qmax

! Generate random baroclinic PV field:
call generate_random_field(qAnom2, qmax, ordbeg, ordend)

! Define resting-state barotropic PV, 2*Omega*sin(latitude):
do j=1,numLatGridPts
   qAnom1(j) = fpole*sin((dble(j)-f12)*dl-hpi)
enddo

! Write full PV (including planetary part qAnom1) in each layer:
open(11,file = 'qq_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)

! Layer 1:
do i=1,numLongGridPts
   qq(:,i) = qAnom1-mrat2*qAnom2(:,i)
enddo
write(11,rec=1) zero,qq

! Layer 2:
do i=1,numLongGridPts
   qq(:,i) = qAnom1+mrat1*qAnom2(:,i)
enddo
write(11,rec=2) zero,qq

! Here, mrat1 = m_1/(m_1+m_2) and mrat2 = m_2/(m_1+m_2) where
! m_j is the mass of layer j --- see constants.f90.
! Note, above qAnom1 is the PV in mode 1 (barotropic) while qAnom2 is
! the PV in mode 2. These are projected onto layers 1 and 2
! in the two loops immediately above.

close(11)

end program random_pv
