program equil
!  -------------------------------------------------------------------------
!  |   Checks how well the time average phi displacement satisfies the     |
!  |   equation for equilibrium.                                           |
!  -------------------------------------------------------------------------

use constants

implicit none

double precision:: a(0:n),phi(0:n),phip(0:n),r(0:n),z(0:n),rbhb(0:n)
double precision:: mass(n),qm(n),uu(n1)
double precision:: dphi(n),dr(n),dz(n)
double precision:: hbar(n),h(0:n),alp(0:n),bet(0:n1)
double precision:: t1(n1),t2(n1),s(0:n)
double precision:: tmp,da,fac
integer:: j

!----------------------------------------------------------------------
 !Read in displaced latitudes phi, masses, etc:
open(11,file='init.asc',status='old')
phi(0)=-hpi !This value at the South Pole is not read in below
a(0)=-hpi
da=pi/dble(n)
do j=1,n
  read(11,*) phi(j),mass(j),qm(j)
  a(j)=da*dble(j)-hpi
enddo
close(11)

 !Ensure that sum(qm) = zero:
tmp=sum(qm)/dble(n)
qm=qm-tmp
 !This ensures the global average vorticity is zero.

 !Accumulate partial sums of -qm for defining angular momentum U:
uu(1)=-qm(1)
do j=2,n1
  uu(j)=uu(j-1)-qm(j)
enddo
 !Note uu(n) = zero; this is not needed.

 !Find background h (h_bar) from masses:
phip=phi
phi=a
call convert
 !Store h_bar(a)*cos(a) in rbhb:
rbhb=h*cos(a)

phi=phip

 !Calculate dphi/da by centred differences:
fac=one/(two*da)
do j=1,n1
  phip(j)=fac*(phi(j+1)-phi(j-1))
enddo

 !Fill boundary values:
phip(0)=two*phip(1)-phip(2)
phip(n)=two*phip(n1)-phip(n2)

write(*,*) ' c = ',cgw,' gamma = ',one/rrad

 !Compute terms on left and right hand side of equilibrium equation:
r=cos(phi)
z=sin(phi)
s(1:n1)=rbhb(1:n1)/(r(1:n1)*phip(1:n1))
 !Fill boundary values:
s(0)=two*s(1)-s(2)
s(n)=two*s(n1)-s(n2)

do j=1,n1
  t1(j)=(fac*csq/phip(j))*(s(j+1)-s(j-1))
  t2(j)=((omega*r(j))**2-(uu(j)/r(j))**2)*z(j)/r(j)
enddo

open(77,file='eq.asc',status='replace')
do j=1,n1
  write(77,*) a(j),t1(j),t2(j),t2(j)-t1(j)
enddo
close(77)

write(*,*)
write(*,*) ' a vs lhs, rhs, and rhs-lhs is in eq.asc'
write(*,*)

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine convert

! Finds the quadratic interpolation of h between phi points

implicit none

! Local variables:
double precision:: rho(n1),am(n1),a0(n1),ap(n1),etd(n1),htd(n1)
double precision:: dphic(n),dzc(n)
double precision:: vtd(n),wtd(n),rhs(n1)
integer:: j

!-----------------------------------------------------------------
! Define cosine and sine of latitude (r & z):
r=cos(phi)
z=sin(phi)

! Compute differences in phi, r and z:
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

do j=0,n1
  h(j)=hbar(j+1)-vtd(j+1)*alp(j)-wtd(j+1)*bet(j)
enddo
h(n)=h(n1)+alp(n1)+bet(n1)

return
end subroutine convert

end program equil

