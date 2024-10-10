program hspectrum
!  ----------------------------------------------------------------------------
!  |   Computes the wavenumber spectrum of the dimensionless height anomaly   |
!  |   h.  Output is to hs.asc, and may be viewed by spec_view.               |
!  ----------------------------------------------------------------------------

use constants
use stafft

implicit none

 !For FFTs:
integer, parameter:: nt=1024
double precision:: trig(2*nt),alk(nt),hs(0:nt)
integer:: factors(5)

 !Other variables:
double precision:: q(n),phi(0:n),r(0:n),z(0:n)
double precision:: mass(n),dphi(n),dr(n),dz(n)
double precision:: hbar(n),h(0:n),alp(0:n),bet(0:n1)
double precision:: t,dum,phi0,dphi0,p
integer:: i,j,k,loop,iread

!----------------------------------------------------------------------
 !Initialise the FFT module:
call initfft(nt,factors,trig)
do k=1,nt
  alk(k)=log10(dble(k))
enddo

 !Grid spacing of regular phi grid:
dphi0=pi/dble(nt)

 !Open input data files:
open(31,file='phi.r8',form='unformatted', access='direct', &
                  status='old',recl=nbytes)
open(32,file='h.r8',form='unformatted', access='direct', &
                  status='old',recl=nbytes)

 !Open output data file:
open(21,file='hs.asc',status='replace')

 !Fix polar values of phi:
phi(0)=-hpi
phi(n)=hpi

 !Read in masses:
open(11,file='init.asc',status='old')
do j=1,n
  read(11,*) dum,mass(j)
enddo
close(11)

 !Read data and process:
loop=0
do
  loop=loop+1
  read(31,rec=loop,iostat=iread) t,q
  if (iread .ne. 0) exit 
  do j=1,n1
    phi(j)=two*q(j)-phi(j-1)
  enddo
  read(32,rec=loop) t,hbar

   !Interpolate h between phi values:
  call convert

   !Interpolate values to a regular grid:
  hs(0)=h(0)  !South pole value
  hs(nt)=h(n) !North pole value
  j=0
  do i=1,nt-1
    phi0=dphi0*dble(i)-hpi
    do while (phi(j) .lt. phi0)
      j=j+1
    enddo
    p=(phi0-phi(j-1))/dphi(j)
    hs(i)=h(j-1)+p*(alp(j-1)+p*bet(j-1))
  enddo

   !Do a cosine transform of hs:
  call dct(1,nt,hs,trig,factors)

   !Save spectrum for k > 0 (it is zero for k = 0):
  write(21,'(f12.5,1x,i5)') t,nt
  do k=1,nt
    write(21,'(2(1x,f12.8))') alk(k),log10(hs(k)**2)
  enddo
enddo

close(21)
close(31)
close(32)

write(*,*)
write(*,*) ' The wavenumber spectra for h are ready in hs.asc.'
write(*,*)
write(*,*) ' View using spec_view'
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

end program hspectrum
