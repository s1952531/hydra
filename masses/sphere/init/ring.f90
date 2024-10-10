program ring
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Initialises a ring of point masses on a surface

! Output files:
! -------------
!   points.dat:    positions  (x,y,z) of all masses
!   speeds.dat:    velocities (u,v,w) of all masses
!   strengths.dat: "strengths" s = m/(4*pi) of all masses
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 ! Import parameters and constants:
use parameters
use constants

implicit double precision (a-h,o-z)

!--------------------------------------------------------------------
write(*,'(a,i2)') '  Number of masses (n) to be placed on the ring = ',n

write(*,*)
write(*,*) ' Co-latitude (theta) of ring (degrees)?'
read(*,*) thetar
thetar=thetar*pi/180.d0

write(*,*) ' All of the vortices have m/4*pi = (2*pi*sin(theta))^2/(n-1)'
write(*,*) ' so they rotate around the North pole in unit time with an'
write(*,*) ' angular velocity of 2*pi.'
str=(twopi*sin(thetar))**2/dble(n-1)
open(20,file='strengths.dat',status='unknown')
do i=1,n
  write(20,'(1x,f14.10)') str
enddo
close(20)

write(*,*)
write(*,*) ' Perturbation in co-latitude (degrees)?'
read(*,*) dtheta
dtheta=dtheta*pi/180.d0

write(*,*)
write(*,*) ' Random seed (integer)?'
read(*,*) iseed
 ! Initialize random # generator:
do i=1,iseed
  uni=rand(0)
enddo

write(*,*)
write(*,*) ' Angular velocity of ring, divided by 2*pi?'
read(*,*) omega
omega=twopi*omega

 ! Write point mass positions & speeds:
open(11,file='points.dat',status='unknown')
open(22,file='speeds.dat',status='unknown')
write(11,'(f12.5)') zero
write(22,'(f12.5)') zero
dphi=twopi/dble(n)
do i=1,n
  phi=dphi*dble(i-1)-pi
  theta=thetar+two*dtheta*rand(0)-dtheta
  r=sin(theta)
  x=r*cos(phi)
  y=r*sin(phi)
  z=cos(theta)
  write(11,'(3(1x,f12.9))') x,y,z
  u=-omega*y
  v=omega*x
  w=zero
  write(22,'(3(1x,f12.9))') u,v,w
enddo
close(11)


 ! Internal subroutine definitions:
contains 

!=======================================================================

function rand(i)
! returns a random number: i is any integer

implicit none
double precision:: rand,r
integer:: i

call random_number(r)
rand=r

return 
end function


 !End main program
end program
!=======================================================================
