!###########################################################################
!  Given an initial vorticity field in qq_init.r8, this routine diffuses
!  this vorticity over a specified length ldiff while preserving the
!  angular momentum.

!  Writes dqq_init.r8

!  Written 6 March 2016 by D.G. Dritschel
!###########################################################################

program diffuse

use constants
use stafft
use deriv1d

 !Declarations:
implicit none

 !Number of "time" steps used to diffuse input data over a length ldiff:
integer,parameter:: nsteps=100
 !Fine grid used to compute mu(theta) by numerical quadrature:
integer,parameter:: nthf=8*ng

double precision,parameter:: dthf=pi/dble(nthf), hdthf=dthf/two
double precision,parameter:: p1=0.21132486541d0,p2=0.78867513459d0
 ! p1,p2: points for 2-point Gaussian Quadrature

 !FFT commons:
integer:: factors(5)
double precision:: trig(2*nt),rk(nt)
double precision:: wave(nt)

 !Tridiagonal commons:
double precision:: sigm(2),psig(2)
double precision:: amla,apla,pc12,pc24
double precision:: d2lon(ng,nt)
integer:: ksig(nt)

double precision:: cmla(ng,nt),dmla(ng)
double precision:: ala0(ng,nt),bla0(ng,nt)
double precision:: cla0(ng,nt),dla0(ng,nt)
double precision:: ela1(ng,nt),ela2(ng,nt)
double precision:: fla1(ng,nt),fla2(ng,nt)
double precision:: dl24

 !Other quantities:
double precision:: qq(ng,nt),ss(ng,nt),dd(ng,nt)
double precision:: rho(ng),rhoi,tau(ng),rdt(ng),tsqi(ng)
double precision:: tlat(ng),rtisq(ng),mu(ng),ff(ng)
double precision:: muf(0:nthf)

double precision:: ldiff,hdti,fac,ang
double precision:: rsum,rsumi,dsumi
double precision:: aspsqm1,rlat,slat,wm,rho1,rho2,addmu,t
double precision:: sm,pm,a0,b0,c0,d0,deti
double precision:: cpla,dpla
integer:: i,j,k,m,ith

!------------------------------------------------------------------------
write(*,*) ' Length over which to diffuse the input data?'
read(*,*) ldiff
 !2/dt:
hdti=four*dble(nsteps)/ldiff**2

 !Initialise latitude functions:

 !Find mu by integration on a grid 8 times finer than the inversion grid:
aspsqm1=asp**2-one
muf(0)=zero
do ith=0,nthf-1
  rho1=sin((dble(ith)+p1)*dthf)
  rho2=sin((dble(ith)+p2)*dthf)
  muf(ith+1)=muf(ith)+hdthf*(rho1*sqrt(one+aspsqm1*rho1**2)+ &
                             rho2*sqrt(one+aspsqm1*rho2**2))
enddo
addmu=(muf(nthf)-muf(0))/two
write(*,'(a,f5.3)') ' Surface aspect ratio = ',asp
write(*,'(a,f14.10)') '    Area/(4*pi) = ',addmu

muf=muf-addmu
 !Find mu on half grid used below:
do j=1,ng
  mu(j)=muf(8*j-4)
enddo

 !Remaining latitude functions:
do j=1,ng
  rlat=(dble(j)-f12)*dl-hpi
   !sin(lat):
  slat=sin(rlat)
   !rho = cos(lat):
  rho(j)=cos(rlat)
   !1/rho:
  rhoi=one/rho(j)
   !1/tau^2:
  fac=one+aspsqm1*rho(j)**2
   !tau:
  tau(j)=one/sqrt(fac)
   !rho/tau:
  rdt(j)=rho(j)/tau(j)
   !-(rho*tau)'/(rho*tau) = tau^2*tan(lat): 
  tlat(j)=tau(j)**2*slat/rho(j)
   !1/(rho*tau)^2:
  rtisq(j)=rhoi**2+aspsqm1
   !2/(dt*tau^2):
  tsqi(j)=hdti*fac
   !2/(dt*tau^2) - 2*rho'/(mu*tau):
  ff(j)=tsqi(j)-two*slat/(mu(j)*tau(j))
enddo

!Used in integrals over mu where dmu = rho/tau:
rsum=f1112*(rdt(1)+rdt(ng))
do j=2,ngm1
  rsum=rsum+rdt(j)
enddo
rsumi=one/rsum
dsumi=rsumi/dble(nt)

!---------------------------
 !Initialise the FFT module:
call initfft(nt,factors,trig)
 !Initialise the spectral derivative module:
call init_deriv(nt,twopi,rk)

 !Define longitudinal wavenumbers:
wave(1)=zero
do m=2,ng
  wm=rk(2*(m-1))
  wave(     m)=wm
  wave(ntp2-m)=wm
enddo
wave(ng+1)=rk(nt)

 !m^2/(rho^2*tau^2) + 2/(dt*tau^2) - 2*rho'/(mu*tau):
do m=1,nt
  fac=wave(m)**2
  d2lon(:,m)=fac*rtisq+ff
enddo

!----------------------------------------------------------------
 !Initialise tridiagonal coefficients:

!-------------------------------------------------------------
 !Sign changes for longitudinal wavenumbers in inversion:

 !Odd physical wavenumbers (1, 3, ...):
do m=2,nt,2
  ksig(m)=1
enddo
sigm(1)=-one
psig(1)=0.9d0

 !Even physical wavenumbers (0, 2, ...):
do m=1,ntm1,2
  ksig(m)=2
enddo
sigm(2)= one
psig(2)=1.1d0

 !Constants used immediately below:
pc12=1.2d0/dl**2
pc24=two*pc12
dl24=dl/24.d0
amla=0.6d0/dl
apla=-amla

!-----------------------------------------------------
 !Initialise Laplace block-tridiagonal inversion:
do m=1,nt
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*d2lon(1,m)
  d0=-pm*tlat(1)

  deti=one/(a0*d0-b0*c0)
  ala0(1,m)=a0*deti
  bla0(1,m)=b0*deti
  cla0(1,m)=c0*deti
  dla0(1,m)=d0*deti

  cpla=pc12-0.1d0*d2lon(2,m)
  dpla=-0.1d0*tlat(2)
  ela1(1,m)= cpla*bla0(1,m)- apla*dla0(1,m)
  ela2(1,m)= dpla*bla0(1,m)-0.2d0*dla0(1,m)
  fla1(1,m)= apla*cla0(1,m)- cpla*ala0(1,m)
  fla2(1,m)=0.2d0*cla0(1,m)- dpla*ala0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-d2lon(j,m)
    d0=-tlat(j)
    cmla(j,m)=pc12-0.1d0*d2lon(j-1,m)
    dmla(j)=-0.1d0*tlat(j-1)

    ala0(j,m)=   ela1(j-1,m)*amla     +fla1(j-1,m)*0.2d0
    bla0(j,m)=b0+ela2(j-1,m)*amla     +fla2(j-1,m)*0.2d0
    cla0(j,m)=c0+ela1(j-1,m)*cmla(j,m)+fla1(j-1,m)*dmla(j)
    dla0(j,m)=d0+ela2(j-1,m)*cmla(j,m)+fla2(j-1,m)*dmla(j)
    deti=one/(ala0(j,m)*dla0(j,m)-bla0(j,m)*cla0(j,m))
    ala0(j,m)=ala0(j,m)*deti
    bla0(j,m)=bla0(j,m)*deti
    cla0(j,m)=cla0(j,m)*deti
    dla0(j,m)=dla0(j,m)*deti

    cpla=pc12-0.1d0*d2lon(j+1,m)
    dpla=-0.1d0*tlat(j+1)
    ela1(j,m)= cpla*bla0(j,m)- apla*dla0(j,m)
    ela2(j,m)= dpla*bla0(j,m)-0.2d0*dla0(j,m)
    fla1(j,m)= apla*cla0(j,m)- cpla*ala0(j,m)
    fla2(j,m)=0.2d0*cla0(j,m)- dpla*ala0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*d2lon(ng,m)
  d0=-pm*tlat(ng)

  cmla(ng,m)=pc12-0.1d0*d2lon(ngm1,m)
  dmla(ng)=-0.1d0*tlat(ngm1)

  ala0(ng,m)=a0+ela1(ngm1,m)*amla      +fla1(ngm1,m)*0.2d0
  bla0(ng,m)=b0+ela2(ngm1,m)*amla      +fla2(ngm1,m)*0.2d0
  cla0(ng,m)=c0+ela1(ngm1,m)*cmla(ng,m)+fla1(ngm1,m)*dmla(ng)
  dla0(ng,m)=d0+ela2(ngm1,m)*cmla(ng,m)+fla2(ngm1,m)*dmla(ng)

  deti=one/(ala0(ng,m)*dla0(ng,m)-bla0(ng,m)*cla0(ng,m))
  ala0(ng,m)=ala0(ng,m)*deti
  bla0(ng,m)=bla0(ng,m)*deti
  cla0(ng,m)=cla0(ng,m)*deti
  dla0(ng,m)=dla0(ng,m)*deti
enddo
 !Initialisation complete

!---------------------------------------------------------------------
 !Read input data and apply diffusion:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(20,rec=1) t,qq
close(20)

 !Compute angular momentum of input data:
do i=1,nt
  dd(:,i)=mu(:)*rdt(:)*qq(:,i)
enddo
ang=zero
do i=1,nt
  ang=ang+f1112*(dd(1,i)+dd(ng,i))
  do j=2,ngm1
    ang=ang+dd(j,i)
  enddo
enddo
ang=ang*dsumi
write(*,*)
write(*,'(a,f12.9)') ' Angular momentum of  initial vorticity = ',ang

 !FFT qq in longitude (semi-spectral):
call forfft(ng,nt,qq,trig,factors) 

 !Iterate to find diffused solution:
do k=1,nsteps
  ss=qq
  call laplinv(ss,qq,dd)
  qq=two*qq-ss
enddo

 !Get qq in physical space:
call revfft(ng,nt,qq,trig,factors)

 !Compute angular momentum of diffused data:
do i=1,nt
  dd(:,i)=mu(:)*rdt(:)*qq(:,i)
enddo
ang=zero
do i=1,nt
  ang=ang+f1112*(dd(1,i)+dd(ng,i))
  do j=2,ngm1
    ang=ang+dd(j,i)
  enddo
enddo
ang=ang*dsumi
write(*,'(a,f12.9)') ' Angular momentum of diffused vorticity = ',ang

 !Write diffused vorticity field:
open(20,file='dqq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

write(*,*)
write(*,*) ' The diffused vorticity field is available in dqq_init.r8'

!===================================================================

 !Internal subroutine definitions (inherit global variables):

contains 

!===================================================================

subroutine laplinv(src,sol,der)
! Inverts Laplace's operator on src to give sol and tau times its
! latitudinal derivative der.  That is Lap(sol) = src
! and der = tau*d(sol)/dlat.
! *** All fields are in spectral space ***

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: src(ng,nt),sol(ng,nt),der(ng,nt)
 !Local arrays:
double precision:: rhs(ng),pfg(ng)

!-----------------------------------------------------------------
 !Loop over azonal longitudinal wavenumbers and solve the 2x2 
 !block tridiagonal problem:
do m=1,nt
   !Multiply source by -2/(dt*tau^2):
  do j=1,ng
    pfg(j)=-src(j,m)*tsqi(j)
  enddo

  pm=psig(ksig(m))

  rhs(1)=pm*pfg(1)+0.1d0*pfg(2)
  do j=2,ngm1
    rhs(j)=pfg(j)+0.1d0*(pfg(j-1)+pfg(j+1))
  enddo
  rhs(ng)=0.1d0*pfg(ngm1)+pm*pfg(ng)

  sol(1,m)=-rhs(1)*bla0(1,m)
  der(1,m)= rhs(1)*ala0(1,m)

  do j=2,ng
    utdb=      -sol(j-1,m)*amla     -der(j-1,m)*0.2d0
    vtdb=rhs(j)-sol(j-1,m)*cmla(j,m)-der(j-1,m)*dmla(j)
    sol(j,m)=utdb*dla0(j,m)-vtdb*bla0(j,m)
    der(j,m)=vtdb*ala0(j,m)-utdb*cla0(j,m)
  enddo

  do j=ngm1,1,-1
    sol(j,m)=ela1(j,m)*sol(j+1,m)+ela2(j,m)*der(j+1,m)+sol(j,m)
    der(j,m)=fla1(j,m)*sol(j+1,m)+fla2(j,m)*der(j+1,m)+der(j,m)
  enddo

   !Multiply derivative by tau so that tau*d(sol)/dlat is returned:
  do j=1,ng
    der(j,m)=der(j,m)*tau(j)
  enddo

enddo

return
end subroutine

!==========================================================================

 !End main program
end program
