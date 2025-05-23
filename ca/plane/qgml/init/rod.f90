program column
! Initialises a columnar distribution of PV whose vertical
! structure is taken from the first baroclinic mode.

use constants

implicit none

double precision,parameter:: qmax=4.d0*pi
double precision:: qq(ng,ng,nz)
double precision:: vl2m(nz,nz),vm2l(nz,nz)
double precision:: a,b,fac,q0,xc,x,y
integer:: ix,iy,iz,m

write(*,*) ' The PV in each layer has a 1 - (r/a)^2 profile.'
write(*,*) ' Enter a:'
read(*,*) a
write(*,*) ' Enter the BC mode-1 displacement the column:'
read(*,*) b

 !Read vertical eigenmodes:
open(60,file='modes.asc',status='old')
do m=1,nz
   read(60,*)
enddo
do m=1,nz
   do iz=1,nz
      read(60,*) vl2m(iz,m),vm2l(iz,m)
   enddo
enddo
close(60)

 !Define baroclinic PV distribution:
fac=one/a**2
do iz=1,nz
   q0=vm2l(iz,2)*qmax
   xc=vm2l(iz,2)*b
   do ix=1,ng
      x=xc+gl*dble(ix-1)-pi
      do iy=1,ng
         y=gl*dble(iy-1)-pi
         qq(iy,ix,iz)=q0*max(one-fac*(x**2+y**2),zero)
      enddo
   enddo
enddo

 !Write data:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program column
