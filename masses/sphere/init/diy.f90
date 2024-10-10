program diy
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Initialises n masses on a sphere and rotates the configuration
! so that the total angular momentum is zero.

! Output files:
! -------------
! coordinates.asc:    positions  (x,y,z) of all masses
! speeds.asc:         velocities (u,v,w) of all masses
! strengths.asc:      "strengths" s = m/(4*pi) of all masses
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 ! Import parameters and constants:
use parameters
use constants

implicit none

double precision,parameter:: pif=pi/180.d0

double precision:: theta(n),phi(n),utheta(n),uphi(n),ptheta(n),pphi(n)
double precision:: x(n),y(n),z(n),px(n),py(n),pz(n)
double precision:: cth,sth,cph,sph
double precision:: angx,angy,angz,angm
double precision:: clat,slat,clon,slon
double precision:: ex,ey,ez,bx,by,bz,ox,oy,oz
double precision:: ux,uy,uz,tx,ty,tz
integer:: i

!--------------------------------------------------------------------
write(*,'(a,i2)') '  Number of masses (n) to be placed = ',n

write(*,*)
write(*,*) ' For each mass, enter its "strength" (s = mass/4*pi),'
write(*,*) ' co-latitude (theta) and longitude (phi) in degrees,'
write(*,*) ' and the corresponding speeds u_theta and u_phi.'
write(*,*)

do i=1,n
  write(*,'(a,i2)') ' Enter s, theta, phi, u_theta & u_phi for mass ',i
  read(*,*) s(i),theta(i),phi(i),utheta(i),uphi(i)
enddo
close(20)

open(20,file='strengths.asc',status='replace')
do i=1,n
  write(20,'(1x,f14.10)') s(i)
enddo
close(20)

do i=1,n
  theta(i)=pif*theta(i)
  cth=cos(theta(i))
  sth=sin(theta(i))
  phi(i)=pif*phi(i)
  cph=cos(phi(i))
  sph=sin(phi(i))
  x(i)=sth*cph
  y(i)=sth*sph
  z(i)=cth
  ! Use p = s*u where s = m/(4*pi) for scaled momenta:
  ptheta(i)=s(i)*utheta(i)
  pphi(i)=s(i)*uphi(i)
  px(i)=ptheta(i)*cth*cph-pphi(i)*sph
  py(i)=ptheta(i)*cth*sph+pphi(i)*cph
  pz(i)=-ptheta(i)*sth
enddo

angx=zero
angy=zero
angz=zero
do i=1,n
  angx=angx+y(i)*pz(i)-z(i)*py(i)
  angy=angy+z(i)*px(i)-x(i)*pz(i)
  angz=angz+x(i)*py(i)-y(i)*px(i)
enddo
angm=sqrt(angx**2+angy**2+angz**2)

write(*,*) ' Angular momentum = ',angx,angy,angz
write(*,*) ' Magnitude of angular momentum = ',angm

slat=angz/angm
clat=sqrt(one-slat**2)
clon=angx/(angm*clat)
slon=angy/(angm*clat)

ex=-slon
ey= clon
ez= zero

bx=-slat*clon
by=-slat*slon
bz= clat

ox=clat*clon
oy=clat*slon
oz=slat

!  ==== USED FOR CHECKING ====
!write(*,*) ' Rotated M_x = ',angx*ex+angy*ey+angz*ez
!write(*,*) ' Rotated M_y = ',angx*bx+angy*by+angz*bz
!write(*,*) ' Rotated M_z = ',angx*ox+angy*oy+angz*oz

! Rotate all vectors:
do i=1,n
  tx=x(i)*ex+y(i)*ey+z(i)*ez
  ty=x(i)*bx+y(i)*by+z(i)*bz
  tz=x(i)*ox+y(i)*oy+z(i)*oz
  x(i)=tx
  y(i)=ty
  z(i)=tz
  tx=px(i)*ex+py(i)*ey+pz(i)*ez
  ty=px(i)*bx+py(i)*by+pz(i)*bz
  tz=px(i)*ox+py(i)*oy+pz(i)*oz
  px(i)=tx
  py(i)=ty
  pz(i)=tz
enddo

!  ==== USED FOR CHECKING ====
angx=zero
angy=zero
angz=zero
do i=1,n
  angx=angx+y(i)*pz(i)-z(i)*py(i)
  angy=angy+z(i)*px(i)-x(i)*pz(i)
  angz=angz+x(i)*py(i)-y(i)*px(i)
enddo
angm=sqrt(angx**2+angy**2+angz**2)

write(*,*) ' Angular momentum = ',angx,angy,angz
write(*,*) ' Magnitude of angular momentum = ',angm

 ! Write point mass positions & speeds:
open(11,file='coordinates.asc',status='replace')
open(22,file='speeds.asc',status='replace')
write(11,'(f12.5)') zero
write(22,'(f12.5)') zero
do i=1,n
  write(11,'(3(1x,f15.12))') x(i),y(i),z(i)
  write(22,'(3(1x,f15.12))') px(i)/s(i),py(i)/s(i),pz(i)/s(i)
enddo
close(11)
close(22)

 !End main program
end program
!=======================================================================
