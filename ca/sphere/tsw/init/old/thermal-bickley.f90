program thermbickley

! Sets up the height field corresponding to a perturbed bickley 
! jet just off the equator. This is used as a thermal equilibrium 
! profile. Also sets up a rest state initial condition.

use constants 

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: clat(ng),slat(ng)
double precision:: hbar(ng),dhdl(ng)
double precision:: hh(ng,nt),he(ng,nt),qq(ng,nt)
double precision:: lambda

!--------------------------------------------------------------
write(*,*) ' We consider a planet of radius 1 rotating with'
write(*,*) ' a period of one "day".'
write(*,*)
write(*,*) ' We set up a Bickley jet'
write(*,*) '  u = u_0 sech^2((phi-phi_0)/b),'
write(*,*) ' as the thermal equilibrium profile.'
write(*,*) ' Enter u_0/2*pi, phi_0 and b:'
read(*,*) u0nd,phi0,b
u0=u0nd*twopi
write(*,*)
write(*,*) ' We add a height perturbation of the form'
write(*,*) '  h_0 exp(-(lambda/alpha)^2-((phi-phi_c)/beta)^2).'
write(*,*) ' Enter h_0, phi_c, alpha and beta:'
read(*,*) h0,phic,alpha,beta

!------------------------------------------------------------
hdl=f12*dl
 !dl: the latitude & longitude grid spacing
!------------------------------------------------------------
 !Define cos, sin and tan(latitude):
do j=1,ng
  rlat=dl*(dble(j)-f12)-hpi
  clat(j)=cos(rlat)
  slat(j)=sin(rlat)
enddo

!----------------------------------------------------------
 !Define zonal velocity and absolute vorticity:
z0=two*u0/b
do j=1,ng
  phi=dl*(dble(j)-f12)-hpi
  ubar=u0/cosh((phi-phi0)/b)**2
  tlat=slat(j)/clat(j)
  dhdl(j)=-csqi*ubar*(fpole*slat(j)+ubar*tlat)
enddo

 !Compute undisturbed height anomaly:
hbar(1)=0.d0
do j=2,ng
  hbar(j)=hbar(j-1)+hdl*(dhdl(j-1)+dhdl(j))
enddo

 !Remove global mean:
rsum=f1112*(clat(1)+clat(ng))
hsum=f1112*(clat(1)*hbar(1)+clat(ng)*hbar(ng))
do j=2,ng-1
  rsum=rsum+clat(j)
  hsum=hsum+clat(j)*hbar(j)
enddo
 !Note, rsum = 1/sin(dl/2) - sin(dl/2)/6 exactly.
havg=hsum/rsum
dsumi=one/(rsum*dble(nt))

do j=1,ng
  hbar(j)=hbar(j)-havg
enddo

 !Define perturbation height anomaly and remove global average:
do i=1,nt
  lambda=dl*dble(i-1)-pi
  do j=1,ng
    phi=dl*(dble(j)-f12)-hpi
    he(j,i)=h0*exp(-(lambda/alpha)**2-((phi-phic)/beta)**2)
  enddo
enddo

hsum=zero
do i=1,nt
  hsum=hsum+f1112*(he(1,i)*clat(1)+he(ng,i)*clat(ng))
  do j=2,ngm1
    hsum=hsum+he(j,i)*clat(j)
  enddo
enddo
hsum=hsum*dsumi

do i=1,nt
  do j=1,ng
    he(j,i)=he(j,i)-hsum
  enddo
enddo

 !Define hh (used for writing h & d), h_equil & PV:
do i=1,nt
  do j=1,ng
    hh(j,i)=zero
    he(j,i)=hbar(j)+he(j,i)
    qq(j,i)=fpole*slat(j)
  enddo
enddo

 !Write equilibrium height field:
open(20,file='hequil.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,he
close(20)

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

 !Write initial divergence field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program
