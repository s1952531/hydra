program gauss

! Program gauss defines a Gaussian seamount somewhere 
! in the domain.

use constants

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: hb(ny,nx)
!------------------------------------------------------------

write(*,*) ' We take the bottom topography to be of the form:'
write(*,*) '    A exp(-((x-x_cen)^2+(y-y_cen)^2)/b^2)'
write(*,*) ' Enter A, x_cen, y_cen and b:'
read(*,*) amp,xcen,ycen,b
bsq=b**2

do ix=1,nx 
  xg=glx*dble(ix-1)-pi
  do iy=1,ny
    yg=gly*dble(iy-1)-pi
    hb(iy,ix)=amp*exp(-((xg-xcen)**2+(yg-ycen)**2)/bsq)
  enddo
enddo

open(12,file='topo.dat',status='unknown')
do ix=1,nx
  do iy=1,ny
    write(12,'(1x,f20.15)')  hb(iy,ix)
  enddo
enddo
close(12)

end program
