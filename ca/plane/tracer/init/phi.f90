program phi
! Initialises a concentration field

use constants

implicit double precision(a-h,o-z)

double precision:: qq(ny,nx)

write(*,*) ' Enter the wavenumber of the initial condition k_phi:'
read(*,*) kphi

do ix=1,nx
  x=xmin+glx*dble(ix-1)
  do iy=1,ny
    y=ymin+gly*dble(iy-1)
    qq(iy,ix)=two*sin(kphi*x)
  enddo
enddo

open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program
