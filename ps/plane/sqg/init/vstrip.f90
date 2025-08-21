program vstrip
! Initialises a wavy buoyancy strip.

use constants

implicit none

double precision:: qq(ny,nx)
double precision:: pow,x0,c1,c2,xfac,x,y,ssq,qavg
integer:: ix,iy

write(*,*)
write(*,*) ' We take b_0/N = (1 - s^2)^p where s = (x - x_c(y))/x_0'
write(*,*) ' for s < 1, and b_0/N = 0 otherwise.'
write(*,*) ' Here, x_c = c_1*sin(y) + c_2*sin(2y) for |y| <= pi.'
write(*,*) ' *** Use p = 0 for the Gaussian profile, b_0/N = e^(-s^2).'
write(*,*) ' Enter p, x_0, c_1 and c_2:'
read(*,*) pow,x0,c1,c2

xfac=one/x0

if (pow > zero) then
   do ix=1,nx
      x=xmin+glx*dble(ix-1)
      do iy=1,ny
         y=ymin+gly*dble(iy-1)
         ssq=((x-c1*sin(y)-c2*sin(two*y))*xfac)**2
         if (ssq < one) then
            qq(iy,ix)=(one-ssq)**pow
         else
            qq(iy,ix)=zero
         endif
      enddo
   enddo
else
   do ix=1,nx
      x=xmin+glx*dble(ix-1)
      do iy=1,ny
         y=ymin+gly*dble(iy-1)
         qq(iy,ix)=exp(-((x-c1*sin(y)-c2*sin(two*y))*xfac)**2)
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

end program vstrip
