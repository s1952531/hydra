program eddy
! Initialises an elliptical eddy.

use constants

implicit none

double precision:: qq(ny,nx)
double precision:: pow,x0,y0,xfac,yfac,x,y,ssq,qavg
integer:: ix,iy

write(*,*)
write(*,*) ' We take b_0/N = (1 - s)^p where s = (x/x_0)^2 + (y/y_0)^2'
write(*,*) ' for s < 1, and b_0/N = 0 otherwise. *** Use p = 0 to instead'
write(*,*) ' take b_0/N = e^{-s}.'
write(*,*)
write(*,*) ' Enter p, x_0 and y_0:'
read(*,*) pow,x0,y0

xfac=one/x0
yfac=one/y0

if (pow > zero) then
   do ix=1,nx
      x=xfac*(xmin+glx*dble(ix-1))
      do iy=1,ny
         y=yfac*(ymin+gly*dble(iy-1))
         ssq=x**2+y**2
         if (ssq < one) then
            qq(iy,ix)=(one-ssq)**pow
         else
            qq(iy,ix)=zero
         endif
      enddo
   enddo
else
   do ix=1,nx
      x=xfac*(xmin+glx*dble(ix-1))
      do iy=1,ny
         y=yfac*(ymin+gly*dble(iy-1))
         qq(iy,ix)=exp(-x**2-y**2)
      enddo
   enddo
endif

 !Calculate and remove average qq:
qavg=zero
do ix=1,nx
   do iy=1,ny
      qavg=qavg+qq(iy,ix)
   enddo
enddo
qavg=qavg/dble(nx*ny)

 !Save average for plotting purposes:
open(44,file='average_qq.asc',status='replace')
write(44,*) qavg
close(44)

qq=qq-qavg

open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program eddy
