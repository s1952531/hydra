program eddy

use constants

! Sets up initial conditions for a circular eddy with a parabolic PV
! profile in each layer.

implicit none
double precision:: qq(0:ny,0:nxm1,nz)
double precision:: xg(0:nxm1),yg(0:ny)
double precision:: xc,yc,qc,rr,rrsq,rsq
integer:: ix,iy,iz

! Define x & y grid points:
do ix=0,nxm1
   xg(ix)=xmin+glx*dble(ix)
enddo

do iy=0,ny
   yg(iy)=ymin+gly*dble(iy)
enddo

! Define PV and make sure beta*y is included:
do iz=1,nz
   write(*,'(a,i2)') ' Enter the centre (xc,yc) of the PV anomaly in layer ',iz
   read(*,*) xc,yc
   write(*,*) ' Enter the radius and amplitude of the PV anomaly'
   read(*,*) rr,qc
   rrsq=rr**2
   do ix=0,nxm1
      do iy=0,ny
         rsq=(xg(ix)-xc)**2+(yg(iy)-yc)**2
         if (rsq < rrsq) then
            qq(iy,ix,iz)=qc*(one-rsq/rrsq)+beta*yg(iy)
         else
            qq(iy,ix,iz)=beta*yg(iy)
         endif
      enddo
   enddo
enddo

 !Write PV distribution to a file:
open(20,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program eddy
