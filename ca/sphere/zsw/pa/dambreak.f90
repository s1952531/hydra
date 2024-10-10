!#########################################################################
! Initialises pam.f90 with h = 1 - a for z > 0 and h = 1 + a for z < 0
! and v = 0 (no motion).
!#########################################################################

program dambreak

use constants

implicit none

 !Latitudes (phi) sine latitudes (z):
double precision:: phi(0:n),z(0:n)

 !Mass and mean PV times mass in each interval (phi_{i-1},phi_i):
double precision:: mass(n),qm(n)

double precision:: a,wnd,w,b,dp,fac
integer:: j

!---------------------------------------------------------
write(*,*) ' We consider a scaled height field of the form'
write(*,*) '           h = 1 - a*tanh(z/w)'
write(*,*) ' with no initial motion.'
write(*,*)
write(*,*) ' Enter a < 1 and w/L_d:'
read(*,*) a,wnd
w=wnd*rrad
write(*,'(a,f9.7)') ' This corresponds to w = ',w

dp=pi/dble(n)
do j=0,n
  phi(j)=dp*dble(j)-hpi
enddo
z=sin(phi)

 !Rest state distribution of PV times mass (gives u = 0):
do j=1,n
  qm(j)=omega*(z(j)**2-z(j-1)**2)
enddo

 !Mass distribution:
fac=a*w
b=one/w
do j=1,n
  mass(j)=z(j)-z(j-1)-fac*log(cosh(z(j)*b)/cosh(z(j-1)*b))
enddo

open(11,file='init.asc',status='replace')
phi(0)=-hpi !This value at the South Pole is not read in below
do j=1,n
  write(11,*) phi(j),mass(j),qm(j),zero
enddo
close(11)

end program dambreak
