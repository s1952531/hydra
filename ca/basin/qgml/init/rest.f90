program rest

use constants

! Sets up initial conditions for a flow at rest in the absence of bathymetry.

implicit none
double precision:: qq(0:ny,0:nx,nz)
double precision:: bety(0:ny)
integer:: ix,iy,iz

! Define beta*y:
do iy=0,ny
   bety(iy)=beta*(ymin+gly*dble(iy))
enddo

! Define PV and make sure beta*y is included:
do iz=1,nz
   do ix=0,nx
      qq(:,ix,iz)=bety
   enddo
enddo

 !Write PV distribution to a file:
open(20,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)
      
end program rest
