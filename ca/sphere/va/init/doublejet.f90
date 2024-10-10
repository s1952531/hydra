program doublejet

!======================================================================
! Sets up a pair of jets with angular velocity of the form
!   Omega = Omega_0*[sech^2((phi-phi_0)/b) + sech^2((phi+phi_0)/b)]
! where Omega_0 = c*F_r/cos(phi_0), and F_r is the specified Froude 
! number at the latitudes phi = +/-phi_0.  A small non-zonal vorticity
! perturbation is also imposed.

! Written 1 June 2016 by D G Dritschel @ St Andrews
!======================================================================

use constants 

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: clat(ng),slat(ng),cof(ng)
double precision:: hbar(ng),obar(ng),zbar(ng),dhdl(ng)
double precision:: zz(ng,nt),hh(ng,nt),dd(ng,nt),qq(ng,nt)

!-----------------------------------------------------------------------
write(*,*) ' We consider a planet of radius 1 rotating with a period of'
write(*,*) ' one "day", i.e. at a rotation rate of 2*pi.'
write(*,*)
write(*,*) ' We start with a pair of jets with angular velocity of the form'
write(*,*) ' Omega = Omega_0*[sech^2((phi-phi_0)/b) + sech^2((phi+phi_0)/b)]'
write(*,*) ' where Omega_0 = c*F_r/cos(phi_0), and F_r is the Froude number.'
write(*,*) ' Enter phi_0 (degrees), b (degrees) & F_r:'
read(*,*) phi0, b, froude
phi0=phi0*pi/180.d0
b=b*pi/180.d0
write(*,*)
write(*,*) ' We add a vorticity perturbation of the form 2*pi*A*exp(-4*d^2/(4-d^2)),'
write(*,*) ' where d is the chord distance from a point at latitude,longitude (p,0).'
write(*,*) ' *** Note, the global vorticity perturbation mean is removed.'
write(*,*) ' Enter A & p (degrees):'
read(*,*) a,pz
pz=pz*pi/180.d0

!------------------------------------------------------------
hdl=f12*dl
 !dl: the latitude & longitude grid spacing
!------------------------------------------------------------
 !Define cos, sin and tan(latitude):
do j=1,ng
  rlat=dl*(dble(j)-f12)-hpi
  clat(j)=cos(rlat)
  slat(j)=sin(rlat)
  cof(j)=fpole*slat(j)
enddo

 !For computing global means:
rsum=f1112*(clat(1)+clat(ng))
do j=2,ng-1
  rsum=rsum+clat(j)
enddo
 !Note, rsum = 1/sin(dl/2) - sin(dl/2)/6 exactly.

!----------------------------------------------------------
 !Define zonal angular velocity and absolute vorticity
 !from the superposition of the two jets:
omega0=cgw*froude/cos(phi0)
binv=one/b

do j=1,ng
  phi=dl*(dble(j)-f12)-hpi

  argm=(phi-phi0)*binv
  chm=cosh(argm)
  shm=sinh(argm)
  obarm=omega0/chm**2

  argp=(phi+phi0)*binv
  chp=cosh(argp)
  shp=sinh(argp)
  obarp=omega0/chp**2

  obar(j)=obarm+obarp
  zbar(j)=cof(j)+two*obarm*(slat(j)+binv*clat(j)*shm/chm) &
                +two*obarp*(slat(j)+binv*clat(j)*shp/chp) 
enddo

 !Remove global mean angular velocity:
osum=f1112*(clat(1)*obar(1)+clat(ng)*obar(ng))
do j=2,ng-1
  osum=osum+clat(j)*obar(j)
enddo
oavg=osum/rsum
do j=1,ng
  obar(j)=obar(j)-oavg
  zbar(j)=zbar(j)-two*oavg*slat(j)
enddo

 !Define dh/d(latitude):
do j=1,ng
  dhdl(j)=-csqi*obar(j)*clat(j)*(cof(j)+obar(j)*slat(j))
enddo

!-------------------------------------------------------------------------
 !Now compute height anomaly h by integration:
hbar(1)=0.d0
do j=2,ng
  hbar(j)=hbar(j-1)+hdl*(dhdl(j-1)+dhdl(j))
enddo

 !Remove global mean h:
hsum=f1112*(clat(1)*hbar(1)+clat(ng)*hbar(ng))
do j=2,ng-1
  hsum=hsum+clat(j)*hbar(j)
enddo
havg=hsum/rsum
do j=1,ng
  hbar(j)=hbar(j)-havg
enddo

!-------------------------------------------------------------------------
 !Record zonal profiles vs z = sin(latitude) for diagnostic purposes:
open(22,file='huzq_zonal.asc',status='unknown')
do j=1,ng
  qbar=zbar(j)/(one+hbar(j))
  write(22,'(5(1x,f15.9))') slat(j),hbar(j),clat(j)*obar(j),zbar(j)-cof(j),qbar
enddo
close(22)

!-------------------------------------------------------------------------
 !Construct vorticity perturbation:
zamp=twopi*a
cpz=cos(pz)
spz=sin(pz)
do i=1,nt
  rlon=dl*dble(i-1)-pi
  clon=cos(rlon)
  slon=sin(rlon)
  do j=1,ng
    zp=cpz*clat(j)*clon+spz*slat(j)
    arg=four*(one-zp)/(abs(one+zp)+small)
    if (arg .lt. 33.d0) then
      zz(j,i)=zamp*exp(-arg)
    else
       !exp(-arg) < machine double precision:
      zz(j,i)=zero
    endif
  enddo
enddo

 !Remove global mean:
zsum=zero
do i=1,nt
  zsum=zsum+f1112*(zz(1,i)*clat(1)+zz(ng,i)*clat(ng))
  do j=2,ngm1
    zsum=zsum+zz(j,i)*clat(j)
  enddo
enddo
zsum=zsum*dsumi

do i=1,nt
  do j=1,ng
    zz(j,i)=zz(j,i)-zsum
  enddo
enddo

!-------------------------------------------------------------------------
 !Define PV:
do i=1,nt
  do j=1,ng
    hh(j,i)=hbar(j)
    dd(j,i)=zero
    qq(j,i)=(zbar(j)+zz(j,i))/(one+hbar(j))
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
write(20,rec=1) zero,dd
close(20)

 !Write initial divergence field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program
