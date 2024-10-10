program relax
!============================================================
! Sets up a rest state consisting of PV varying as 2*Omega*z
! together with a thermal equilibrium depth field.
!============================================================

use constants 

 !Declarations:
implicit none

double precision:: hh(ng,nt),qq(ng,nt)
double precision:: slat(ng),clat(ng)
double precision:: eps,phi0,b,amp,rlat,fm,acmlon,binv,phi
double precision:: rsum,rsumi,dsumi,vsum
integer:: i,j,m

!------------------------------------------------------------------------------
write(*,*) ' We consider a thermal equilibrium height field which increases by'
write(*,*) ' a fraction, epsilon, of the mean height from equator to pole.'
write(*,*) ' The main increase occurs around the latitudes phi = +/-phi_0,'
write(*,*) ' and occurs over a half width of b degrees.'
write(*,*)

write(*,*) ' Enter pole - equator height difference: '
read(*,*) eps

write(*,*) ' Enter phi_0 in degrees: '
read(*,*) phi0
phi0=phi0*pi/180.d0

write(*,*) ' Enter b in degrees: '
read(*,*) b
b=b*pi/180.d0

write(*,*)
write(*,*) ' To induce instability, we displace the latitudes by the function'
write(*,*) ' A*cos(phi)*cos(m*lambda), where lambda is the longitude.'
write(*,*)
write(*,*) ' Enter the amplitude A: ' 
read(*,*) amp
write(*,*) ' Enter the wavenumber m (integer): '
read(*,*) m

!------------------------------------------------------------------------------
 !Define cos and sin(latitude):
do j=1,ng
  rlat=dl*(dble(j)-f12)-hpi
  slat(j)=sin(rlat)
  clat(j)=cos(rlat)
enddo

 !For removing global averages:
rsum=f1112*(clat(1)+clat(ng))
do j=2,ngm1
  rsum=rsum+clat(j)
enddo
 !Note, rsum = 1/sin(dl/2) - sin(dl/2)/6 exactly.
rsumi=one/rsum
dsumi=rsumi/dble(nt)

 !Define hh (used for writing h & d) and PV:
do i=1,nt
  do j=1,ng
    hh(j,i)=zero
    qq(j,i)=fpole*slat(j)
  enddo
enddo

 !Write initial height field:
open(20,file='hh_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

 !Write initial divergence field:
open(20,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

!------------------------------------------------------------------------------
 !Compute and write thermal equilibrium height field:
fm=dble(m)
eps=f12*eps
binv=one/b
do i=1,nt
  acmlon=amp*cos(fm*(dl*dble(i-1)-pi))
  do j=1,ng
    phi=dl*(dble(j)-f12)-hpi-clat(j)*acmlon
    hh(j,i)=eps*(two+tanh(binv*(phi-phi0))-tanh(binv*(phi+phi0)))
  enddo
enddo

 !Remove global mean h:
vsum=zero
do i=1,nt
  vsum=vsum+f1112*(clat(1)*hh(1,i)+clat(ng)*hh(ng,i))
  do j=2,ngm1
    vsum=vsum+clat(j)*hh(j,i)
  enddo
enddo
vsum=vsum*dsumi
do i=1,nt
  do j=1,ng
    hh(j,i)=hh(j,i)-vsum
  enddo
enddo

 !Write data:
open(20,file='hequil.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

end program
