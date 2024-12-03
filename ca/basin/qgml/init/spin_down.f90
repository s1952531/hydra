program spin_down

use constants

! Sets up initial conditions for a spin-down experiment.

implicit none
double precision:: qb(0:ny,0:nx)
double precision:: qq(0:ny,0:nx,nz)
double precision:: xg(0:nx),yg(0:ny)
double precision:: amp_psi1,k_psi1,l_psi1,k2l2,pv_aux,d1,K1sq,d2,K2,F1,F2,td
integer:: ix,iy,iz

! Define x & y grid points:
do ix=0,nx
   xg(ix)=xmin+glx*dble(ix)
enddo

do iy=0,ny
   yg(iy)=ymin+gly*dble(iy)
enddo

! Define PV and make sure beta*y is included:
write(*,*) ' Enter the Amplitude A of the upper layer streamfunction'
read(*,*) amp_psi1
write(*,*) ' Enter the x-wavenumbers k'
read(*,*) k_psi1
write(*,*) ' Enter the y-wavenumbers l'
read(*,*) l_psi1

! Get stretching term coefficients for the top two layers
open(unit=10,file='vertical.asc',status='old')
read(10,*) d1,K1sq
read(10,*) d2
close(10)

do ix=0,nx
   do iy=0,ny
      pv_aux=amp_psi1*sin(k_psi1*xg(ix))*cos(l_psi1*yg(iy))
      qq(iy,ix,1)=-(k_psi1**2+l_psi1**2+K1sq/d1)*pv_aux+beta*yg(iy)
      qq(iy,ix,2)=K1sq/d2*pv_aux+beta*yg(iy)
   enddo
enddo

do iz=3,nz
   do ix=0,nx
      do iy=0,ny
         qq(iy,ix,iz)=beta*yg(iy)
      enddo
   enddo
enddo

! Get bathymetry
if (bath) then
   open(11,file='bath.r8',form='unformatted', &
        access='direct',status='old',recl=2*nhbytes)
   read(11,rec=1) td,qb
   close(11)

    do ix=0,nx
        do iy=0,ny
            qq(iy,ix,nz)=qq(iy,ix,nz)+qb(iy,ix)
        enddo
    enddo
endif

 !Write PV distribution to a file:
open(20,file='qq_init.r8',form='unformatted',access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program spin_down
