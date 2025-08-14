program double_grid

use constants

! Reads half-resolution data from qq_final.r8 and linearly interpolates
! these to the current resolution; writes qq_init.r8.

! *** WARNING: nx & ny must be divisible by 2.

implicit none
double precision:: qqlo(0:ny/2,0:nx/2,nz),qq(0:ny,0:nx,nz),t
integer:: ix,iy,iz,nblo

 !Read in half-resolution PV (qqlo) as a double-precision field:
nblo=8*((nx/2+1)*(ny/2+1)*nz+1)
open(11,file='qq_final.r8',form='unformatted', &
      access='direct',status='old',recl=nblo)
read(11,rec=1) t,qqlo
close(11)

 !Interpolate to current resolution as qq:
do iz=1,nz
   do ix=0,nx/2
      do iy=0,ny/2
         qq(2*iy,2*ix,iz)=qqlo(iy,ix,iz)
      enddo
      do iy=0,ny/2-1
         qq(2*iy+1,2*ix,iz)=f12*(qqlo(iy,ix,iz)+qqlo(iy+1,ix,iz))
      enddo
   enddo
   do ix=0,nx/2-1
      do iy=0,ny
         qq(iy,2*ix+1,iz)=f12*(qq(iy,2*ix,iz)+qq(iy,2*ix+2,iz))
      enddo
   enddo
enddo

 !Write PV distribution to a file:
open(20,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program double_grid
