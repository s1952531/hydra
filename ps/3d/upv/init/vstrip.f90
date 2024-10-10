program vstrip
! Initialises a scaled surface pressure field p_0(x,y)/f^2 where f is
! the Coriolis frequency.

! Here, we take p_0 = p_max * e^[-(x - c*sin(y))^2/a^2]

use constants

implicit none

double precision:: p0(ng,ng)
double precision:: pm,x0,c,xfac,x,y
integer:: ix,iy

write(*,*)
write(*,*) ' We take p_0 = p_m * e^{-s^2} where s = (x - c*sin(y))/x_0.'
write(*,*) ' Enter p_m/f^2, x_0 and c:'
read(*,*) pm,x0,c

xfac=one/x0

do ix=1,ng
  x=gl*dble(ix-1)-pi
  do iy=1,ng
    y=gl*dble(iy-1)-pi
    p0(iy,ix)=pm*exp(-((x-c*sin(y))*xfac)**2)
  enddo
enddo

! Write data:
open(11,file='p0_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nhbytes)
write(11,rec=1) zero,p0
close(11)

end program vstrip
