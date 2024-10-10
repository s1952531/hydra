program sine_bath
! |---------------------------------------------------|
! |   This routine sets up sinusoidal bathymetry      |
! |   of the form A sin(k(x-x_min)) sin(k(y-y_min).   |
! |---------------------------------------------------|

 !Import constants & parameters:
use constants

implicit none

double precision:: bb(0:ny,0:nx)
double precision:: amp,k_bath
integer:: ix,iy

write(*,*) ' We consider bathymetry eta_b of the form '
write(*,*) ' f_0*eta_b/H_nz = A*sin(k*(x-x_min))*sin(k*(y-y_min),'
write(*,*) ' where H_nz is the depth of the lowest layer.'
write(*,*)

write(*,*) ' Enter A:'
read(*,*) amp

write(*,*) ' Enter k:'
read(*,*) k_bath

 !Set up bathymetry:
do ix=0,nx
   do iy=0,ny
      bb(iy,ix)=amp*sin(k_bath*(xmin+glx*dble(ix)))* &
                    sin(k_bath*(ymin+gly*dble(iy)))
   enddo
enddo

 !Write bathymetry to a file:
open(11,file='bath.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,bb
close(11)

end program sine_bath
