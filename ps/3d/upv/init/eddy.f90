program eddy
! Initialises a scaled surface pressure field p_0(x,y)/f^2 where f is
! the Coriolis frequency.

! Here, we take p_0 = p_max * e^(-x^2/x_0^2-y^2/y_0^2).

use constants

implicit none

double precision:: p0(ng,ng)
double precision:: pm,x0,y0,xfac,yfac,x,y
integer:: ix,iy

write(*,*)
write(*,*) ' We take p_0 = p_m * exp[-(x/x_0)^2 - (y/y_0)^2].'
write(*,*) ' Enter p_m/f^2, x_0 and y_0:'
read(*,*) pm,x0,y0

xfac=one/x0
yfac=one/y0

do ix=1,ng
  x=xfac*(gl*dble(ix-1)-pi)
  do iy=1,ng
    y=yfac*(gl*dble(iy-1)-pi)
    p0(iy,ix)=pm*exp(-x**2-y**2)
  enddo
enddo

! Write data:
open(11,file='p0_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nhbytes)
write(11,rec=1) zero,p0
close(11)

end program eddy
