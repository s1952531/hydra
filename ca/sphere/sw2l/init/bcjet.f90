program bcjet
!==========================================================================
! Initialises a baroclinic (mode 2) jet with a Gaussian shape PV profile.
! Use balinit to balance the initial flow using zero first derivatives
! of divergence and acceleration divergence.
!==========================================================================

! Import constants & parameters:
use constants

! Declarations:
implicit none

double precision:: qq(numLatGridPts,numLongGridPts), qAnom2(numLatGridPts,numLongGridPts), qAnom1(numLatGridPts), phic(numLongGridPts)
double precision:: robc, wid, phi0, amp2, amp3, fac, lon, lat

integer:: i, j

!--------------------------------------------------------------------------
write(*,*) ' Consider a baroclinic jet with PV anomaly = f*A*exp(-s^2)'
write(*,*) ' where s = (phi-phi_c)/w, f is the Coriolis frequency and'
write(*,*) ' phi_c(lambda) = phi_0+A_2*sin(2*lambda)+A_3*sin(3*lambda).'
write(*,*)

write(*,*) ' Enter the baroclinic Rossby number A and the width w (degrees):'
read(*,*) robc, wid
! Convert to radians:
fac = pi/180.d0
wid = fac*wid

write(*,*) ' Enter phi_0 (degrees):'
read(*,*) phi0
! Convert to radians:
phi0 = fac*phi0

write(*,*) ' Enter A_2 & A_3 (degrees):'
read(*,*) amp2, amp3
! Convert to radians:
amp2 = fac*amp2
amp3 = fac*amp3

! Define phi_c(lambda):
do i = 1,numLongGridPts
   lon = dl*dble(i-1)-pi
   phic(i) = amp2*sin(two*lon)+amp3*sin(three*lon)
enddo

! Define baroclinic PV qAnom2:
fac = one/wid**2
robc = fpole*robc
do i = 1,numLongGridPts
   do j = 1,numLatGridPts
      lat = (dble(j)-f12)*dl-hpi
      qAnom2(j,i) = robc*exp(-fac*(lat-phic(i))**2)
   enddo
enddo

!-------------------------------------------------------------
! Define resting-state barotropic PV, 2*Omega*sin(latitude):
do j = 1,numLatGridPts
   lat = (dble(j)-f12)*dl-hpi
   qAnom1(j) = fpole*sin(lat)
enddo

! Write full PV (including planetary part qAnom1) in each layer:
open(11,file = 'qq_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)

! Layer 1:
do i = 1,numLongGridPts
   qq(:,i) = qAnom1-mrat2*qAnom2(:,i)
enddo
write(11,rec=1) zero,qq

! Layer 2:
do i = 1,numLongGridPts
   qq(:,i) = qAnom1+mrat1*qAnom2(:,i)
enddo
write(11,rec=2) zero,qq

! Here, mrat1 = m_1/(m_1+m_2) and mrat2 = m_2/(m_1+m_2) where
! m_j is the mass of layer j --- see constants.f90.
! Note, above qAnom1 is the PV in mode 1 (barotropic) while qAnom2 is
! the PV in mode 2. These are projected onto layers 1 and 2
! in the two loops immediately above.

! Write zero fields of divergence and acceleration divergence:
qq = zero
open(11,file = 'dd_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)
write(11,rec=1) zero,qq
write(11,rec=2) zero,qq
close(11)

open(11,file = 'gg_init.r8',form = 'unformatted', &
      access = 'direct',status = 'replace',recl = 2*nbytes)
write(11,rec=1) zero,qq
write(11,rec=2) zero,qq
close(11)

end program bcjet
