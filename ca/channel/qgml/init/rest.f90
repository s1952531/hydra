program rest

use constants

! Sets up initial conditions for a flow at rest.

implicit none
double precision:: qq(0:ny,0:nxm1,nz),qb(0:ny,0:nxm1)
double precision:: bety(0:ny),td
integer:: ix,iy,iz

! Define beta*y:
do iy=0,ny
   bety(iy)=beta*(ymin+gly*dble(iy))
enddo

! Define PV and make sure beta*y is included:
do iz=1,nz
   do ix=0,nxm1
      qq(:,ix,iz)=bety
   enddo
enddo

! Get bathymetry
if (bath) then
   open(11,file='bath.r8',form='unformatted', &
        access='direct',status='old',recl=2*nhbytes)
   read(11,rec=1) td,qb
   close(11)

    do ix=0,nxm1
        do iy=0,ny
            qq(iy,ix,nz)=qq(iy,ix,nz)+qb(iy,ix)
        enddo
    enddo
endif

 !Write PV distribution to a file:
open(20,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program rest
