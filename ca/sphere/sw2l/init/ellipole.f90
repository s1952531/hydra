program ellipole

!===============================================================
!  Sets up a baroclinic PV field in the shape of an elliptical
!  Gaussian centred on the north pole.

!  Developed to test the initial balancing procedure.
!===============================================================

 !Import constants are parameters:
use constants

 !Declarations:
implicit none

double precision:: qq(numLatGridPts,numLongGridPts),qAnom2(numLatGridPts,numLongGridPts),qAnom1(numLatGridPts)
double precision:: amp,asp,aspi,b,bsq
double precision:: lon,lat,x,y,z
integer:: i,j

!------------------------------------------------------------
 !Define planetary vorticity (qAnom1):
do j=1,numLatGridPts
  qAnom1(j)=fpole*sin((dble(j)-f12)*dl-hpi)
enddo

 !Define the PV anomaly:
write(*,*) ' The baroclinic PV field has the form'
write(*,*) '   A*exp(-(x^2*lambda+y^2/lambda)/(4*b^2*(1-z)).'
write(*,*)
write(*,*) ' Enter A, lambda and b:'
read(*,*) amp,asp,b

aspi=one/asp
bsq=b**2

do i=1,numLongGridPts
   lon=dl*dble(i-1)-pi
   do j=1,numLatGridPts
      lat=(dble(j)-f12)*dl-hpi
      x=cos(lat)*cos(lon)
      y=cos(lat)*sin(lon)
      z=sin(lat)
      qAnom2(j,i)=amp*exp(-(asp*x**2+aspi*y**2)/(bsq*(one+z)))
   enddo
enddo

! Write full PV (including planetary part qAnom1) in each layer:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)

! Layer 1:
do i=1,numLongGridPts
   qq(:,i)=qAnom1-mrat2*qAnom2(:,i)
enddo
write(11,rec=1) zero,qq

! Layer 2:
do i=1,numLongGridPts
   qq(:,i)=qAnom1+mrat1*qAnom2(:,i)
enddo
write(11,rec=2) zero,qq

! Here, mrat1=m_1/(m_1+m_2) and mrat2=m_2/(m_1+m_2) where
! m_j is the mass of layer j --- see constants.f90.
! Note, above qAnom1 is the PV in mode 1 (barotropic) while qAnom2 is
! the PV in mode 2. These are projected onto layers 1 and 2
! in the two loops immediately above.

! Write zero fields of divergence and acceleration divergence:
qq=zero
open(11,file='dd_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
write(11,rec=2) zero,qq
close(11)

open(11,file='gg_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
write(11,rec=2) zero,qq
close(11)

end program ellipole
