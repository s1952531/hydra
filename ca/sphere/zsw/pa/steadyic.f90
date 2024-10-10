!#########################################################################
! Initialises pam.f90 with the time averaged fields obtained from
! profile.f90 and pam.f90
!#########################################################################

program steadyic

use constants

implicit none

 !Lagrangian labels (a), latitudes (phi) sine latitudes (z):
double precision:: a(0:n),phi(0:n),z(0:n)

 !Mass and mean PV times mass in each interval (phi_{i-1},phi_i):
double precision:: mass(n),qm(n)

integer:: j

open(11,file='init.asc',status='old')
do j=1,n
  read(11,*) phi(j),mass(j),qm(j)
enddo
close(11)

 !Re-write original init file:
open(11,file='ori-init.asc',status='replace')
phi(0)=-hpi !This value at the South Pole is not read in below
do j=1,n
  write(11,*) phi(j),mass(j),qm(j),zero
enddo
close(11)

 !Open and read displacement in phi:
open(11,file='davg.asc',status='old')
do j=0,n
  read(11,*) a(j),phi(j)
enddo
close(11)

phi=a+phi !Displaced phi coordinate

 !Re-write init file:
open(11,file='init.asc',status='replace')
do j=1,n
  write(11,*) phi(j),mass(j),qm(j),zero
enddo
close(11)

end program steadyic
