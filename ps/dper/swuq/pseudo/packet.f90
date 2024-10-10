program packet

! Program packet initialises a packet of nonlinear waves for 
! subsequent evolution by the swgw code.

use constants

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: hh(ny,nx),dd(ny,nx)
!------------------------------------------------------------

write(*,*) ' We take the height anomaly to be of the form:'
write(*,*) '  A cos(2*pi*k*x/L_x+2*pi*l*y/L_y)*exp(-(x^2+y^2)/b^2)'
write(*,*) ' Enter A, k, l and b:'
read(*,*) amp,rk,rl,b
k=nint(rk)
l=nint(rl)
dok=two*pi*dble(k)/ellx
dol=two*pi*dble(l)/ellx
bsq=b**2
sig=sqrt(fsq+csq*(dok**2+dol**2))

xcen=zero
ycen=zero
do ix=1,nx 
  xg=glx*dble(ix-1)-pi
  do iy=1,ny
    yg=gly*dble(iy-1)-pi
    gauss=exp(-((xg-xcen)**2+(yg-ycen)**2)/bsq)
    hh(iy,ix)= amp*cos(dok*xg+dol*yg)*gauss
    dd(iy,ix)=-amp*sig*sin(dok*xg+dol*yg)*gauss
  enddo
enddo

open(12,file='hh_init.dat',status='unknown')
open(13,file='dd_init.dat',status='unknown')

do ix=1,nx
  do iy=1,ny
    write(12,'(1x,f20.15)')  hh(iy,ix)
    write(13,'(1x,f20.15)')  dd(iy,ix)
  enddo
enddo
close(12)
close(13)

end program
