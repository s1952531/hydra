program hderiv
!===================================================================
! Tests calculating dh/dphi given only the "mass" in each interval
! (phi_{i-1},phi_i) (the integral of h*cos(phi)*dphi), the boundary
! conditions dh/dphi = 0 at phi = +/-pi/2, and the values of phi_i.
  
!        Written 29/10/2020 by D G Dritschel @ St Andrews
!===================================================================

implicit none

 !Number of grid intervals n:
integer, parameter:: n=20, n1=n-1, n2=n-2

 !Basic constants:
double precision, parameter:: zero=0.d0, one=1.d0, two=2.d0, f12=one/two
double precision, parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0
double precision, parameter:: hpi=pi/two

 !Latitudes plus cos & sin latitudes (r & z):
double precision:: phi(0:n),r(0:n),z(0:n)

 !Mass in each interval (phi_{i-1},phi_i):
double precision:: mass(n)

 !Numerical values of h and dh/dphi:
double precision:: h(0:n),dhdp(0:n)

 !Exact values of h and dh/dphi for testing purposes only:
double precision:: he(0:n),dhdpe(0:n)

!-----------------------------------------------------------------
! Initialise:
call initialise

! Find discrete approximation to h & dh/dphi:
call convert

! Print out results:
call finalise

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine initialise

! Sets up the mass in each latitude interval.

implicit none

! Local variables:
double precision:: s(0:n)
double precision:: eps,fni,fac,x0,a
integer:: j,m

!-----------------------------------------------------------------
write(*,*) ' We choose the discrete latitudes phi_j according to'
write(*,*) '    phi(j) = pi*(j/n + eps*sin(j*m*pi/n)) - pi/2.'
write(*,*)
write(*,*) ' Enter eps and m:'
read(*,*) eps,m

fni=one/dble(n)
fac=pi*dble(m)
do j=0,n
  x0=fni*dble(j)
  phi(j)=pi*(x0+eps*sin(x0*fac))-hpi
enddo
r=cos(phi)
z=sin(phi)

write(*,*) ' To test the method, we take h = 1/(a + z) where z = sin(phi).'
write(*,*) ' Enter a > 1:'
read(*,*) a

! Calculate masses and store exact values of h and dh/dphi:
s=log((a+z)/(a-one))
do j=1,n
  mass(j)=s(j)-s(j-1)
enddo

he=one/(a+z)
dhdpe=-r*he**2

return
end subroutine initialise

!=======================================================================

subroutine convert

! Finds h and dh/dphi from the mass in each phi interval

implicit none

! Local variables:
double precision:: rho(n1),am(n1),a0(n1),ap(n1),etd(n1),htd(n1)
double precision:: dphi(n),dr(n),dz(n)
double precision:: dphic(n),dzc(n)
double precision:: hbar(n),vtd(n),wtd(n),rhs(n1)
double precision:: alp(0:n),bet(0:n1)
integer:: j

!-----------------------------------------------------------------
do j=1,n
  dphi(j)=phi(j)-phi(j-1)
  dr(j)=r(j)-r(j-1)
  dz(j)=z(j)-z(j-1)
enddo
dphic=one/dphi
dzc=one/dz

! Average value of h over interval (phi_{j-1},phi_j):
hbar=mass*dzc
! Integral of xi*dz where xi = (phi-phi_{j-1})/dphi:
vtd=(z(1:n)+dr*dphic)/dz
! Integral of xi^2*dz:
wtd=(z(1:n)+two*dphic*(r(1:n)-dz*dphic))/dz

do j=1,n1
  rho(j)=dphi(j)/dphi(j+1)
  a0(j)=two*vtd(j+1)+(one-wtd(j))*rho(j)-wtd(j+1)
enddo

do j=2,n1
  am(j)=two*(one-vtd(j))+wtd(j)-one
  ap(j-1)=wtd(j)*rho(j)
enddo

htd(1)=1.d0/a0(1)
etd(1)=-ap(1)*htd(1)

do j=2,n2
  htd(j)=1.d0/(a0(j)+am(j)*etd(j-1))
  etd(j)=-ap(j)*htd(j)
enddo
htd(n1)=1.d0/(a0(n1)+am(n1)*etd(n2))

do j=1,n1
  rhs(j)=two*(hbar(j+1)-hbar(j))
enddo

alp(1)=rhs(1)*htd(1)

do j=2,n1
  alp(j)=(rhs(j)-am(j)*alp(j-1))*htd(j)
enddo

do j=n2,1,-1
  alp(j)=etd(j)*alp(j+1)+alp(j)
enddo

!-----------------------------------------------------------------
alp(0)=0.d0
alp(n)=0.d0

do j=0,n2
  bet(j)=f12*(rho(j+1)*alp(j+1)-alp(j))
enddo
bet(n1)=-f12*alp(n1)

dhdp(0)=0.d0
do j=1,n1
  dhdp(j)=alp(j)/dphi(j+1)
enddo
dhdp(n)=0.d0

do j=0,n1
  h(j)=hbar(j+1)-vtd(j+1)*alp(j)-wtd(j+1)*bet(j)
enddo
h(n)=h(n1)+alp(n1)+bet(n1)

return
end subroutine convert

!=======================================================================

subroutine finalise

! Prints out results of the test
  
implicit none

! Local variable:
integer:: j

!-----------------------------------------------------------------
write(*,*)

open(77,file='dhdp.asc',status='replace')
do j=0,n
  write(77,*) phi(j),dhdp(j),dhdpe(j)
enddo
close(77)
write(*,*) ' R.m.s. error in dh/dphi = ',sqrt(sum((dhdp-dhdpe)**2)/dble(n))

open(77,file='h.asc',status='replace')
do j=0,n
  write(77,*) phi(j),h(j),he(j)
enddo
close(77)
write(*,*) ' R.m.s. error in    h    = ',sqrt(sum((h-he)**2)/dble(n))

return
end subroutine finalise

end program hderiv
