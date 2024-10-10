program ring
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Initialises two point masses a distance 2 apart rotating with a
! unit period of rotation.  A small perturbation to the positions
! of the masses is added to examine stability.

! Output files:
! -------------
!   points.dat:    positions  (x,y) of all masses
!   speeds.dat:    velocities (u,v) of all masses
!   strengths.dat: "strengths" s = m/(4*pi) of all masses
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 ! Import parameters and constants:
use parameters
use constants

implicit double precision (a-h,o-z)

!--------------------------------------------------------------------
if (n .ne. 2) then
  write(*,*) ' There must be just 2 masses!  Stopping!'
  stop
endif

write(*,*)
write(*,*) ' Mass ratio?'
read(*,*) alpha

s1=8.d0*pi**2/(one+alpha)
s2=alpha*s1

x2=two/(one+alpha)
x1=-alpha*x2
v1=-s2/twopi
v2= s1/twopi

open(20,file='strengths.dat',status='unknown')
write(20,'(1x,f14.10)') s1
write(20,'(1x,f14.10)') s2
close(20)

write(*,*) ' The masses are a distance 2 apart in equilibrium.'
write(*,*) ' Maximum random displacement?'
read(*,*) disp

write(*,*)
write(*,*) ' Random seed (integer)?'
read(*,*) iseed
 ! Initialize random # generator:
do i=1,iseed
  uni=rand(0)
enddo

 ! Write point mass positions & speeds:
open(11,file='points.dat',status='unknown')
write(11,'(f12.5)') zero
x=x1+disp*(two*rand(0)-one)
y=   disp*(two*rand(0)-one)
write(11,'(3(1x,f12.9))') x,y
x=x2+disp*(two*rand(0)-one)
y=   disp*(two*rand(0)-one)
write(11,'(3(1x,f12.9))') x,y
close(11)

open(22,file='speeds.dat',status='unknown')
write(22,'(f12.5)') zero
write(22,'(3(1x,f12.9))') zero,v1
write(22,'(3(1x,f12.9))') zero,v2
close(22)

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
