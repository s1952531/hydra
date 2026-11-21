module spectral

! Module containing subroutines for spectral operations, inversion, etc.

use force
use constants
use stafft
use deriv1d

implicit none

 !FFT commons:
integer:: factors(5)
double precision:: trig(2*numLongGridPts),rk(numLongGridPts),wzon(numLongGridPts)

 !Filters and damping operators:
double precision:: d2lon(numLatGridPts,numLongGridPts),diss(numLatGridPts,numLongGridPts),rdis(numLatGridPts,numLongGridPts)
double precision:: plon(numLatGridPts,numLongGridPts),glon(numLatGridPts,numLongGridPts)
double precision:: flonlo(numLatGridPts,numLongGridPts),flonhi(numLatGridPts,numLongGridPts)
double precision:: flatlo(numLongGridPts),flathi(numLongGridPts),hyplat(numLongGridPts)
double precision:: plat(numLongGridPts),glat(numLongGridPts)
double precision:: alm(numLatGridPts)

 !Tridiagonal commons:

 !For taking a derivative in latitude along great circles (periodic):
double precision:: etd1(numLongGridPts),htd1(numLongGridPts),ptd1(numLongGridPts),apd1,denod1
double precision:: etd2(numLongGridPts),htd2(numLongGridPts),ptd2(numLongGridPts),apd2,denod2
double precision:: psig(2),amla
integer:: ksig(numLongGridPts)

 !For inverting the Laplacian:
double precision:: cmla(numLatGridPts,numLongGridPts),dmla(numLatGridPts)
double precision:: ala0(numLatGridPts,numLongGridPts),bla0(numLatGridPts,numLongGridPts)
double precision:: cla0(numLatGridPts,numLongGridPts),dla0(numLatGridPts,numLongGridPts)
double precision:: ela1(numLatGridPts,numLongGridPts),ela2(numLatGridPts,numLongGridPts)
double precision:: fla1(numLatGridPts,numLongGridPts),fla2(numLatGridPts,numLongGridPts)
double precision:: etd0(numLatGridPts),htd0(numLatGridPts)
double precision:: rsumi,dsumi,dl24

 !For inverting the Helmholtzian (mode 1):
double precision:: cm1he(numLatGridPts,numLongGridPts),dm1he(numLatGridPts)
double precision:: a1he0(numLatGridPts,numLongGridPts),b1he0(numLatGridPts,numLongGridPts)
double precision:: c1he0(numLatGridPts,numLongGridPts),d1he0(numLatGridPts,numLongGridPts)
double precision:: e1he1(numLatGridPts,numLongGridPts),e1he2(numLatGridPts,numLongGridPts)
double precision:: f1he1(numLatGridPts,numLongGridPts),f1he2(numLatGridPts,numLongGridPts)

 !For inverting the Helmholtzian (mode 2):
double precision:: cm2he(numLatGridPts,numLongGridPts),dm2he(numLatGridPts)
double precision:: a2he0(numLatGridPts,numLongGridPts),b2he0(numLatGridPts,numLongGridPts)
double precision:: c2he0(numLatGridPts,numLongGridPts),d2he0(numLatGridPts,numLongGridPts)
double precision:: e2he1(numLatGridPts,numLongGridPts),e2he2(numLatGridPts,numLongGridPts)
double precision:: f2he1(numLatGridPts,numLongGridPts),f2he2(numLatGridPts,numLongGridPts)

 !For semi-implicit evolution of delta and gamma (mode 1):
double precision:: cm1si(numLatGridPts,numLongGridPts),dm1si(numLatGridPts)
double precision:: a1si0(numLatGridPts,numLongGridPts),b1si0(numLatGridPts,numLongGridPts)
double precision:: c1si0(numLatGridPts,numLongGridPts),d1si0(numLatGridPts,numLongGridPts)
double precision:: e1si1(numLatGridPts,numLongGridPts),e1si2(numLatGridPts,numLongGridPts)
double precision:: f1si1(numLatGridPts,numLongGridPts),f1si2(numLatGridPts,numLongGridPts)

 !For semi-implicit evolution of delta and gamma (mode 2):
double precision:: cm2si(numLatGridPts,numLongGridPts),dm2si(numLatGridPts)
double precision:: a2si0(numLatGridPts,numLongGridPts),b2si0(numLatGridPts,numLongGridPts)
double precision:: c2si0(numLatGridPts,numLongGridPts),d2si0(numLatGridPts,numLongGridPts)
double precision:: e2si1(numLatGridPts,numLongGridPts),e2si2(numLatGridPts,numLongGridPts)
double precision:: f2si1(numLatGridPts,numLongGridPts),f2si2(numLatGridPts,numLongGridPts)

 !Other arrays:
double precision:: tanLat(numLatGridPts),sinLat(numLatGridPts),cosLat(numLatGridPts),secantLat(numLatGridPts)
double precision:: cof(numLatGridPts),beta(numLatGridPts),fsq(numLatGridPts),dfb(numLatGridPts)

! Weights for evolving the PV forcing (if used):
double precision:: wold,wnew

contains

!=======================================================================

subroutine init_spectral
! Initialises spectral module

implicit none

 !Local variables:
double precision:: norm(minSpherHarmDegree:maxSpherHarmDegree),kbi,wrat
double precision:: alpha(numLatGridPts),sigm(2)
double precision:: radianLat,wm,wmi,fac,rkc,rkci,f283,weff,dfac,rsum
double precision:: a0d1,a0d2,pc12,pc24,a0,b0,c0,d0,deti
double precision:: apla,cpla,dpla,sm,pm
double precision:: cpsi,dpsi,cphe,dphe
integer, dimension(:), allocatable :: seed
integer:: i,j,k,m,order

!----------------------------------------------------------
 !Latitude functions:
do j=1,numLatGridPts
  radianLat=(dble(j)-f12)*dl-hpi
  sinLat(j)=sin(radianLat)
  cosLat(j)=cos(radianLat)
enddo

secantLat=one/cosLat
tanLat=sinLat/cosLat
 !Coriolis frequency (f):
cof=fpole*sinLat
fsq=cof**2
 !beta = df/dlat: = d(fpole*sinLat)/dlat = fpole*d(sinLat)/dlat = fpole*cosLat
beta=fpole*cosLat
 !2*f*beta needed for acceleration divergence tendency:
dfb=two*cof*beta

!-------------------------------------------
 !Initialise the FFT module:
call initfft(numLongGridPts,factors,trig)
 !Initialise the spectral derivative module:
call init_deriv(numLongGridPts,twoPi,rk)

 !Define longitudinal wavenumbers:
wzon(1)=zero
do m=2,numLatGridPts
  wm=rk(2*(m-1))
  wzon(     m)=wm
  wzon(ntp2-m)=wm
enddo
wzon(ngp1)=rk(numLongGridPts)
wmi=one/wzon(numLatGridPts+1)

 !For saving longitudinal spectra (excluding global average):
do m=1,numLatGridPts
  alm(m)=log10(wzon(m+1))
enddo

 !m^2/r^2 factor which appears in Laplace's operator:
do m=1,numLongGridPts
  fac=wzon(m)**2
  d2lon(:,m)=fac*secantLat**2
enddo

!----------------------------------------------------------------
 !Define 2nd-order Butterworth filter (used in subroutines lopass
 !and hipass) and include partial de-aliasing:
rkc=f13*dble(numLatGridPts)
 !rkc: filter wavenumber
rkci=one/rkc
f283=28.d0/3.d0

 !Latitude filters:
flatlo=one/(one+(rkci*wzon)**4)
flathi=(one-flatlo)*f12*(one-erf(16.d0*wzon*wmi-f283))

 !Longitude filters:
do m=1,numLongGridPts
  do j=1,numLatGridPts
    weff=wzon(m)*secantLat(j)
    flonlo(j,m)=one/(one+(rkci*weff)**4)
    flonhi(j,m)=(one-flonlo(j,m))*f12*(one-erf(16.d0*weff*wmi-f283))
  enddo
enddo

!--------------------------------------------------------------
 !Latitudinal hyperviscous damping used on residual PV (qd),
 !divergence (ds), and acceleration divergence (gs):
dfac=hyperviscDampRate*fpole
fac=dt*dfac
hyplat=one/(one+fac*(wmi*wzon)**(2*hyperviscPower))

 !Longitudinal hyperviscous damping rate on wavenumber m:
do m=1,numLongGridPts
  fac=wmi*wzon(m)
  do j=1,numLatGridPts
    diss(j,m)=dfac*(fac*secantLat(j))**(2*hyperviscPower)
  enddo
enddo
rdis=dt2i+diss

!-----------------------------------------------------------------
 ! Initialise PV forcing (if present):

if (forcing) then
   !Initialise normalisation for spherical harmonics.  In general,
   !norm(n) = sqrt(S(n)/(2*n+1)) where S(n) is the spectrum and
   !n is the order.

  if (minSpherHarmDegree > 1) then
     !Here we choose a normalisation corresponding to a flat spectrum
     !S(n) = 1 (usually narrow-band), where n = order:
    do order = minSpherHarmDegree, maxSpherHarmDegree
      norm(order) = one/sqrt(dble(1+2*order))
    enddo
  else
     !We consider the spectrum S(n) = wrat*exp(-wrat) where 
     !wrat = (n/kb)^3 and kb = 4*maxSpherHarmDegree/9:
    kbi = 2.25d0/dble(maxSpherHarmDegree)
    do order = minSpherHarmDegree, maxSpherHarmDegree          !Note minSpherHarmDegree = 1 here
      wrat = (kbi*dble(order))**3
      norm(order) = wrat*exp(-wrat)/sqrt(dble(1+2*order))
    enddo
     !Note, maxSpherHarmDegree is chosen so that norm < 10^{-4} when n = maxSpherHarmDegree+1.
     !That is, 2.25 factor in the definition of kbi ensures this, or
     !equivalently the relation maxSpherHarmDegree = 9*kb/4 ensures this.
  endif

   !Initialise forcing module (and spherical module):
  call init_forcing(norm)

   !Weights needed to have time-correlated forcing:
  wold=one-dt/tPvForcingDecorrelation
  wnew=sqrt(one-wold**2)

endif

!-------------------------------------
 !Initialise tridiagonal coefficients:

!---------------------------
 !First latitude derivative:
a0d1=dl*f43
apd1=dl*f13

htd1(1)=one/a0d1
ptd1(1)=-apd1*htd1(1)
etd1(1)=ptd1(1)

do j=2,numLongGridPts
  htd1(j)=one/(a0d1+apd1*etd1(j-1))
  ptd1(j)=-apd1*ptd1(j-1)*htd1(j)
  etd1(j)=-apd1*htd1(j)
enddo

ptd1(ntm1)=etd1(ntm1)+ptd1(ntm1)

do j=ntm2,1,-1
  ptd1(j)=etd1(j)*ptd1(j+1)+ptd1(j)
enddo

denod1=one/(one-etd1(numLongGridPts)*ptd1(1)-ptd1(numLongGridPts))

!----------------------------
 !Second latitude derivative:
a0d2=dlsq*f56
apd2=dlsq*f112

htd2(1)=one/a0d2
ptd2(1)=-apd2*htd2(1)
etd2(1)=ptd2(1)

do j=2,numLongGridPts
  htd2(j)=one/(a0d2+apd2*etd2(j-1))
  ptd2(j)=-apd2*ptd2(j-1)*htd2(j)
  etd2(j)=-apd2*htd2(j)
enddo

ptd2(ntm1)=etd2(ntm1)+ptd2(ntm1)

do j=ntm2,1,-1
  ptd2(j)=etd2(j)*ptd2(j+1)+ptd2(j)
enddo

denod2=one/(one-etd2(numLongGridPts)*ptd2(1)-ptd2(numLongGridPts))

!--------------------------------------------------------
 !Sign changes for longitudinal wavenumbers in inversion:

 !Odd physical wavenumbers (1, 3, ...):
do m=2,numLongGridPts,2
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

 !Constants used immediately below (some are used in subroutines):
pc12=1.2d0/dl**2
pc24=two*pc12
dl24=dl/24.d0

amla=0.6d0/dl
apla=-amla

!------------------------------------------------
 !Initialise Laplace block-tridiagonal inversion:

 !zonally-symmetric part (physical wavenumber 0):
htd0(1)=one/f74
etd0(1)=-f14*htd0(1)

do j=2,ngm1
  htd0(j)=one/(f32+f14*etd0(j-1))
  etd0(j)=-f14*htd0(j)
enddo

htd0(numLatGridPts)=one/(f74+f14*etd0(ngm1))

rsum=f1112*(cosLat(1)+cosLat(numLatGridPts))+sum(cosLat(2:ngm1))
 !Note, rsum = 1/sin(dl/2) - sin(dl/2)/6 exactly.
rsumi=one/rsum
dsumi=rsumi/dble(numLongGridPts)

 !azonal part:
do m=2,numLongGridPts
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  alpha=d2lon(:,m)

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(1)
  d0=-pm*tanLat(1)

  deti=one/(a0*d0-b0*c0)
  ala0(1,m)=a0*deti
  bla0(1,m)=b0*deti
  cla0(1,m)=c0*deti
  dla0(1,m)=d0*deti

  cpla=pc12-0.1d0*alpha(2)
  dpla=-0.1d0*tanLat(2)
  ela1(1,m)= cpla*bla0(1,m)- apla*dla0(1,m)
  ela2(1,m)= dpla*bla0(1,m)-0.2d0*dla0(1,m)
  fla1(1,m)= apla*cla0(1,m)- cpla*ala0(1,m)
  fla2(1,m)=0.2d0*cla0(1,m)- dpla*ala0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-alpha(j)
    d0=-tanLat(j)
    cmla(j,m)=pc12-0.1d0*alpha(j-1)
    dmla(j)=-0.1d0*tanLat(j-1)

    ala0(j,m)=   ela1(j-1,m)*amla     +fla1(j-1,m)*0.2d0
    bla0(j,m)=b0+ela2(j-1,m)*amla     +fla2(j-1,m)*0.2d0
    cla0(j,m)=c0+ela1(j-1,m)*cmla(j,m)+fla1(j-1,m)*dmla(j)
    dla0(j,m)=d0+ela2(j-1,m)*cmla(j,m)+fla2(j-1,m)*dmla(j)
    deti=one/(ala0(j,m)*dla0(j,m)-bla0(j,m)*cla0(j,m))
    ala0(j,m)=ala0(j,m)*deti
    bla0(j,m)=bla0(j,m)*deti
    cla0(j,m)=cla0(j,m)*deti
    dla0(j,m)=dla0(j,m)*deti

    cpla=pc12-0.1d0*alpha(j+1)
    dpla=-0.1d0*tanLat(j+1)
    ela1(j,m)= cpla*bla0(j,m)- apla*dla0(j,m)
    ela2(j,m)= dpla*bla0(j,m)-0.2d0*dla0(j,m)
    fla1(j,m)= apla*cla0(j,m)- cpla*ala0(j,m)
    fla2(j,m)=0.2d0*cla0(j,m)- dpla*ala0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(numLatGridPts)
  d0=-pm*tanLat(numLatGridPts)

  cmla(numLatGridPts,m)=pc12-0.1d0*alpha(ngm1)
  dmla(numLatGridPts)=-0.1d0*tanLat(ngm1)

  ala0(numLatGridPts,m)=a0+ela1(ngm1,m)*amla      +fla1(ngm1,m)*0.2d0
  bla0(numLatGridPts,m)=b0+ela2(ngm1,m)*amla      +fla2(ngm1,m)*0.2d0
  cla0(numLatGridPts,m)=c0+ela1(ngm1,m)*cmla(numLatGridPts,m)+fla1(ngm1,m)*dmla(numLatGridPts)
  dla0(numLatGridPts,m)=d0+ela2(ngm1,m)*cmla(numLatGridPts,m)+fla2(ngm1,m)*dmla(numLatGridPts)

  deti=one/(ala0(numLatGridPts,m)*dla0(numLatGridPts,m)-bla0(numLatGridPts,m)*cla0(numLatGridPts,m))
  ala0(numLatGridPts,m)=ala0(numLatGridPts,m)*deti
  bla0(numLatGridPts,m)=bla0(numLatGridPts,m)*deti
  cla0(numLatGridPts,m)=cla0(numLatGridPts,m)*deti
  dla0(numLatGridPts,m)=dla0(numLatGridPts,m)*deti
enddo

!-----------------------------------------------------------------------
 !Initialise tri-diagonal arrays for semi-implicit time-stepping of
 !divergence and acceleration divergence for *mode 1*
 !(include hyperviscous damping):
do m=1,numLongGridPts
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  alpha=d2lon(:,m)+(fsq+rdis(:,m)**2)*ctsq1i

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(1)
  d0=-pm*tanLat(1)

  deti=one/(a0*d0-b0*c0)
  a1si0(1,m)=a0*deti
  b1si0(1,m)=b0*deti
  c1si0(1,m)=c0*deti
  d1si0(1,m)=d0*deti

  cpsi=pc12-0.1d0*alpha(2)
  dpsi=-0.1d0*tanLat(2)
  e1si1(1,m)= cpsi*b1si0(1,m)- apla*d1si0(1,m)
  e1si2(1,m)= dpsi*b1si0(1,m)-0.2d0*d1si0(1,m)
  f1si1(1,m)= apla*c1si0(1,m)- cpsi*a1si0(1,m)
  f1si2(1,m)=0.2d0*c1si0(1,m)- dpsi*a1si0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-alpha(j)
    d0=-tanLat(j)
    cm1si(j,m)=pc12-0.1d0*alpha(j-1)
    dm1si(j)=-0.1d0*tanLat(j-1)

    a1si0(j,m)=   e1si1(j-1,m)*amla      +f1si1(j-1,m)*0.2d0
    b1si0(j,m)=b0+e1si2(j-1,m)*amla      +f1si2(j-1,m)*0.2d0
    c1si0(j,m)=c0+e1si1(j-1,m)*cm1si(j,m)+f1si1(j-1,m)*dm1si(j)
    d1si0(j,m)=d0+e1si2(j-1,m)*cm1si(j,m)+f1si2(j-1,m)*dm1si(j)
    deti=one/(a1si0(j,m)*d1si0(j,m)-b1si0(j,m)*c1si0(j,m))
    a1si0(j,m)=a1si0(j,m)*deti
    b1si0(j,m)=b1si0(j,m)*deti
    c1si0(j,m)=c1si0(j,m)*deti
    d1si0(j,m)=d1si0(j,m)*deti

    cpsi=pc12-0.1d0*alpha(j+1)
    dpsi=-0.1d0*tanLat(j+1)
    e1si1(j,m)= cpsi*b1si0(j,m)- apla*d1si0(j,m)
    e1si2(j,m)= dpsi*b1si0(j,m)-0.2d0*d1si0(j,m)
    f1si1(j,m)= apla*c1si0(j,m)- cpsi*a1si0(j,m)
    f1si2(j,m)=0.2d0*c1si0(j,m)- dpsi*a1si0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(numLatGridPts)
  d0=-pm*tanLat(numLatGridPts)

  cm1si(numLatGridPts,m)=pc12-0.1d0*alpha(ngm1)
  dm1si(numLatGridPts)=-0.1d0*tanLat(ngm1)

  a1si0(numLatGridPts,m)=a0+e1si1(ngm1,m)*amla       +f1si1(ngm1,m)*0.2d0
  b1si0(numLatGridPts,m)=b0+e1si2(ngm1,m)*amla       +f1si2(ngm1,m)*0.2d0
  c1si0(numLatGridPts,m)=c0+e1si1(ngm1,m)*cm1si(numLatGridPts,m)+f1si1(ngm1,m)*dm1si(numLatGridPts)
  d1si0(numLatGridPts,m)=d0+e1si2(ngm1,m)*cm1si(numLatGridPts,m)+f1si2(ngm1,m)*dm1si(numLatGridPts)

  deti=one/(a1si0(numLatGridPts,m)*d1si0(numLatGridPts,m)-b1si0(numLatGridPts,m)*c1si0(numLatGridPts,m))
  a1si0(numLatGridPts,m)=a1si0(numLatGridPts,m)*deti
  b1si0(numLatGridPts,m)=b1si0(numLatGridPts,m)*deti
  c1si0(numLatGridPts,m)=c1si0(numLatGridPts,m)*deti
  d1si0(numLatGridPts,m)=d1si0(numLatGridPts,m)*deti
enddo

!-----------------------------------------------------------------------
 !Initialise tri-diagonal arrays for semi-implicit time-stepping of
 !divergence and acceleration divergence for *mode 2*
 !(include hyperviscous damping):
do m=1,numLongGridPts
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  alpha=d2lon(:,m)+(fsq+rdis(:,m)**2)*ctsq2i

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(1)
  d0=-pm*tanLat(1)

  deti=one/(a0*d0-b0*c0)
  a2si0(1,m)=a0*deti
  b2si0(1,m)=b0*deti
  c2si0(1,m)=c0*deti
  d2si0(1,m)=d0*deti

  cpsi=pc12-0.1d0*alpha(2)
  dpsi=-0.1d0*tanLat(2)
  e2si1(1,m)= cpsi*b2si0(1,m)- apla*d2si0(1,m)
  e2si2(1,m)= dpsi*b2si0(1,m)-0.2d0*d2si0(1,m)
  f2si1(1,m)= apla*c2si0(1,m)- cpsi*a2si0(1,m)
  f2si2(1,m)=0.2d0*c2si0(1,m)- dpsi*a2si0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-alpha(j)
    d0=-tanLat(j)
    cm2si(j,m)=pc12-0.1d0*alpha(j-1)
    dm2si(j)=-0.1d0*tanLat(j-1)

    a2si0(j,m)=   e2si1(j-1,m)*amla      +f2si1(j-1,m)*0.2d0
    b2si0(j,m)=b0+e2si2(j-1,m)*amla      +f2si2(j-1,m)*0.2d0
    c2si0(j,m)=c0+e2si1(j-1,m)*cm2si(j,m)+f2si1(j-1,m)*dm2si(j)
    d2si0(j,m)=d0+e2si2(j-1,m)*cm2si(j,m)+f2si2(j-1,m)*dm2si(j)
    deti=one/(a2si0(j,m)*d2si0(j,m)-b2si0(j,m)*c2si0(j,m))
    a2si0(j,m)=a2si0(j,m)*deti
    b2si0(j,m)=b2si0(j,m)*deti
    c2si0(j,m)=c2si0(j,m)*deti
    d2si0(j,m)=d2si0(j,m)*deti

    cpsi=pc12-0.1d0*alpha(j+1)
    dpsi=-0.1d0*tanLat(j+1)
    e2si1(j,m)= cpsi*b2si0(j,m)- apla*d2si0(j,m)
    e2si2(j,m)= dpsi*b2si0(j,m)-0.2d0*d2si0(j,m)
    f2si1(j,m)= apla*c2si0(j,m)- cpsi*a2si0(j,m)
    f2si2(j,m)=0.2d0*c2si0(j,m)- dpsi*a2si0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(numLatGridPts)
  d0=-pm*tanLat(numLatGridPts)

  cm2si(numLatGridPts,m)=pc12-0.1d0*alpha(ngm1)
  dm2si(numLatGridPts)=-0.1d0*tanLat(ngm1)

  a2si0(numLatGridPts,m)=a0+e2si1(ngm1,m)*amla       +f2si1(ngm1,m)*0.2d0
  b2si0(numLatGridPts,m)=b0+e2si2(ngm1,m)*amla       +f2si2(ngm1,m)*0.2d0
  c2si0(numLatGridPts,m)=c0+e2si1(ngm1,m)*cm2si(numLatGridPts,m)+f2si1(ngm1,m)*dm2si(numLatGridPts)
  d2si0(numLatGridPts,m)=d0+e2si2(ngm1,m)*cm2si(numLatGridPts,m)+f2si2(ngm1,m)*dm2si(numLatGridPts)

  deti=one/(a2si0(numLatGridPts,m)*d2si0(numLatGridPts,m)-b2si0(numLatGridPts,m)*c2si0(numLatGridPts,m))
  a2si0(numLatGridPts,m)=a2si0(numLatGridPts,m)*deti
  b2si0(numLatGridPts,m)=b2si0(numLatGridPts,m)*deti
  c2si0(numLatGridPts,m)=c2si0(numLatGridPts,m)*deti
  d2si0(numLatGridPts,m)=d2si0(numLatGridPts,m)*deti
enddo

!--------------------------------------------------------------
 !Initialise tri-diagonal arrays needed to invert dimensionless
 !height anomaly for *mode 1*
do m=1,numLongGridPts
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  alpha=d2lon(:,m)+fsq*ctsq1i

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(1)
  d0=-pm*tanLat(1)

  deti=one/(a0*d0-b0*c0)
  a1he0(1,m)=a0*deti
  b1he0(1,m)=b0*deti
  c1he0(1,m)=c0*deti
  d1he0(1,m)=d0*deti

  cphe=pc12-0.1d0*alpha(2)
  dphe=-0.1d0*tanLat(2)
  e1he1(1,m)= cphe*b1he0(1,m)- apla*d1he0(1,m)
  e1he2(1,m)= dphe*b1he0(1,m)-0.2d0*d1he0(1,m)
  f1he1(1,m)= apla*c1he0(1,m)- cphe*a1he0(1,m)
  f1he2(1,m)=0.2d0*c1he0(1,m)- dphe*a1he0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-alpha(j)
    d0=-tanLat(j)
    cm1he(j,m)=pc12-0.1d0*alpha(j-1)
    dm1he(j)=-0.1d0*tanLat(j-1)

    a1he0(j,m)=   e1he1(j-1,m)*amla      +f1he1(j-1,m)*0.2d0
    b1he0(j,m)=b0+e1he2(j-1,m)*amla      +f1he2(j-1,m)*0.2d0
    c1he0(j,m)=c0+e1he1(j-1,m)*cm1he(j,m)+f1he1(j-1,m)*dm1he(j)
    d1he0(j,m)=d0+e1he2(j-1,m)*cm1he(j,m)+f1he2(j-1,m)*dm1he(j)
    deti=one/(a1he0(j,m)*d1he0(j,m)-b1he0(j,m)*c1he0(j,m))
    a1he0(j,m)=a1he0(j,m)*deti
    b1he0(j,m)=b1he0(j,m)*deti
    c1he0(j,m)=c1he0(j,m)*deti
    d1he0(j,m)=d1he0(j,m)*deti

    cphe=pc12-0.1d0*alpha(j+1)
    dphe=-0.1d0*tanLat(j+1)
    e1he1(j,m)= cphe*b1he0(j,m)- apla*d1he0(j,m)
    e1he2(j,m)= dphe*b1he0(j,m)-0.2d0*d1he0(j,m)
    f1he1(j,m)= apla*c1he0(j,m)- cphe*a1he0(j,m)
    f1he2(j,m)=0.2d0*c1he0(j,m)- dphe*a1he0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(numLatGridPts)
  d0=-pm*tanLat(numLatGridPts)

  cm1he(numLatGridPts,m)=pc12-0.1d0*alpha(ngm1)
  dm1he(numLatGridPts)=-0.1d0*tanLat(ngm1)

  a1he0(numLatGridPts,m)=a0+e1he1(ngm1,m)*amla       +f1he1(ngm1,m)*0.2d0
  b1he0(numLatGridPts,m)=b0+e1he2(ngm1,m)*amla       +f1he2(ngm1,m)*0.2d0
  c1he0(numLatGridPts,m)=c0+e1he1(ngm1,m)*cm1he(numLatGridPts,m)+f1he1(ngm1,m)*dm1he(numLatGridPts)
  d1he0(numLatGridPts,m)=d0+e1he2(ngm1,m)*cm1he(numLatGridPts,m)+f1he2(ngm1,m)*dm1he(numLatGridPts)

  deti=one/(a1he0(numLatGridPts,m)*d1he0(numLatGridPts,m)-b1he0(numLatGridPts,m)*c1he0(numLatGridPts,m))
  a1he0(numLatGridPts,m)=a1he0(numLatGridPts,m)*deti
  b1he0(numLatGridPts,m)=b1he0(numLatGridPts,m)*deti
  c1he0(numLatGridPts,m)=c1he0(numLatGridPts,m)*deti
  d1he0(numLatGridPts,m)=d1he0(numLatGridPts,m)*deti
enddo

!--------------------------------------------------------------
 !Initialise tri-diagonal arrays needed to invert dimensionless
 !height anomaly for *mode 2*
do m=1,numLongGridPts
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  alpha=d2lon(:,m)+fsq*ctsq2i

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(1)
  d0=-pm*tanLat(1)

  deti=one/(a0*d0-b0*c0)
  a2he0(1,m)=a0*deti
  b2he0(1,m)=b0*deti
  c2he0(1,m)=c0*deti
  d2he0(1,m)=d0*deti

  cphe=pc12-0.1d0*alpha(2)
  dphe=-0.1d0*tanLat(2)
  e2he1(1,m)= cphe*b2he0(1,m)- apla*d2he0(1,m)
  e2he2(1,m)= dphe*b2he0(1,m)-0.2d0*d2he0(1,m)
  f2he1(1,m)= apla*c2he0(1,m)- cphe*a2he0(1,m)
  f2he2(1,m)=0.2d0*c2he0(1,m)- dphe*a2he0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-alpha(j)
    d0=-tanLat(j)
    cm2he(j,m)=pc12-0.1d0*alpha(j-1)
    dm2he(j)=-0.1d0*tanLat(j-1)

    a2he0(j,m)=   e2he1(j-1,m)*amla      +f2he1(j-1,m)*0.2d0
    b2he0(j,m)=b0+e2he2(j-1,m)*amla      +f2he2(j-1,m)*0.2d0
    c2he0(j,m)=c0+e2he1(j-1,m)*cm2he(j,m)+f2he1(j-1,m)*dm2he(j)
    d2he0(j,m)=d0+e2he2(j-1,m)*cm2he(j,m)+f2he2(j-1,m)*dm2he(j)
    deti=one/(a2he0(j,m)*d2he0(j,m)-b2he0(j,m)*c2he0(j,m))
    a2he0(j,m)=a2he0(j,m)*deti
    b2he0(j,m)=b2he0(j,m)*deti
    c2he0(j,m)=c2he0(j,m)*deti
    d2he0(j,m)=d2he0(j,m)*deti

    cphe=pc12-0.1d0*alpha(j+1)
    dphe=-0.1d0*tanLat(j+1)
    e2he1(j,m)= cphe*b2he0(j,m)- apla*d2he0(j,m)
    e2he2(j,m)= dphe*b2he0(j,m)-0.2d0*d2he0(j,m)
    f2he1(j,m)= apla*c2he0(j,m)- cphe*a2he0(j,m)
    f2he2(j,m)=0.2d0*c2he0(j,m)- dphe*a2he0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(numLatGridPts)
  d0=-pm*tanLat(numLatGridPts)

  cm2he(numLatGridPts,m)=pc12-0.1d0*alpha(ngm1)
  dm2he(numLatGridPts)=-0.1d0*tanLat(ngm1)

  a2he0(numLatGridPts,m)=a0+e2he1(ngm1,m)*amla       +f2he1(ngm1,m)*0.2d0
  b2he0(numLatGridPts,m)=b0+e2he2(ngm1,m)*amla       +f2he2(ngm1,m)*0.2d0
  c2he0(numLatGridPts,m)=c0+e2he1(ngm1,m)*cm2he(numLatGridPts,m)+f2he1(ngm1,m)*dm2he(numLatGridPts)
  d2he0(numLatGridPts,m)=d0+e2he2(ngm1,m)*cm2he(numLatGridPts,m)+f2he2(ngm1,m)*dm2he(numLatGridPts)

  deti=one/(a2he0(numLatGridPts,m)*d2he0(numLatGridPts,m)-b2he0(numLatGridPts,m)*c2he0(numLatGridPts,m))
  a2he0(numLatGridPts,m)=a2he0(numLatGridPts,m)*deti
  b2he0(numLatGridPts,m)=b2he0(numLatGridPts,m)*deti
  c2he0(numLatGridPts,m)=c2he0(numLatGridPts,m)*deti
  d2he0(numLatGridPts,m)=d2he0(numLatGridPts,m)*deti
enddo

!--------------------------------------------------------
 !Re-define damping operator for use in qd time stepping:
diss=two/(one+dt2*diss)

return 
end subroutine init_spectral

!======================================================================

subroutine main_invert(qSpec1,qSpec2,deltaSpec1,deltaSpec2,gammaSpec1,gammaSpec2,h1,h2,u1,u2,v1,v2,qAnom1,qAnom2,z1,z2)
! Given the PV anomaly qs, divergence ds and acceleration divergence gs
! (all in spectral space and in both layers), this routine computes the
! dimensionless depth anomaly h and the velocity field (u,v) in physical
! space.  It also returns the corrected PV anomaly (q) and relative
! vorticity (z) in physical space (in both layers). The PV anomaly is
! adjusted to ensure that the domain integral of the relative vorticity
! vanishes (a requirement by Stokes' theorem).

! On input, h, u & v contain a guess for the solution; these are zero
! at program start since they are unknown (see caps.f90). Subsequent
! calls start with the existing values of h, u & v at the previous
! time (sub)step.  

implicit none

 !Passed variables:
double precision:: qSpec1(numLatGridPts,numLongGridPts),qSpec2(numLatGridPts,numLongGridPts)
double precision:: deltaSpec1(numLatGridPts,numLongGridPts),deltaSpec2(numLatGridPts,numLongGridPts)
double precision:: gammaSpec1(numLatGridPts,numLongGridPts),gammaSpec2(numLatGridPts,numLongGridPts)

double precision:: h1(numLatGridPts,numLongGridPts),h2(numLatGridPts,numLongGridPts)
double precision:: u1(numLatGridPts,numLongGridPts),u2(numLatGridPts,numLongGridPts)
double precision:: v1(numLatGridPts,numLongGridPts),v2(numLatGridPts,numLongGridPts)
double precision:: qAnom1(numLatGridPts,numLongGridPts),qAnom2(numLatGridPts,numLongGridPts)
double precision:: z1(numLatGridPts,numLongGridPts),z2(numLatGridPts,numLongGridPts)

 !Local parameters and variables:
double precision,parameter:: tole=1.d-12
 !tole: relative energy norm error in successive iterates when finding
 !      hj, uj & vj from qj, dj & gj (for j = 1 & 2).  The energy norm is
 !      <du1^2+dv1^2+c_1^2*dh1^2> + (rho_2/rho_1)*(H_2/H_1)*
 !      <du2^2+dv2^2+c_2^2*dh2^2>
 !      where <:> means a horizontal domain average and (duj,dvj,dhj)
 !      is either the current guess for (uj,vj,hj) or the difference
 !      from the previous guess.

double precision,parameter:: wgt=0.5d0, wgtc=1.d0-wgt
 !wgt: relaxation parameter used to speed convergence (wgt near 1 fails)

double precision:: h1o(numLatGridPts,numLongGridPts),h2o(numLatGridPts,numLongGridPts)
double precision:: u1o(numLatGridPts,numLongGridPts),u2o(numLatGridPts,numLongGridPts)
double precision:: v1o(numLatGridPts,numLongGridPts),v2o(numLatGridPts,numLongGridPts)

double precision:: ud1(numLatGridPts,numLongGridPts),ud2(numLatGridPts,numLongGridPts)
double precision:: vd1(numLatGridPts,numLongGridPts),vd2(numLatGridPts,numLongGridPts)
double precision:: hm1(numLatGridPts,numLongGridPts),hm2(numLatGridPts,numLongGridPts)

double precision:: wka(numLatGridPts,numLongGridPts),wkb(numLatGridPts,numLongGridPts),wkc(numLatGridPts,numLongGridPts)
double precision:: enorm
integer:: i,m

!------------------------------------------------------------------
 !First compute divergent part of the velocity field (this is fixed
 !in the iteration below):

 !Invert Laplace's operator on the divergence (ds) in both layers:
call laplinv(deltaSpec1,wka,vd1)
call laplinv(deltaSpec2,wkb,vd2)
 !Here the divergence potential is wka or wkb, and its latitudinal
 !derivative is vd1 or vd2, the divergent meridional velocity
 !(all in semi-spectral space)

 !Compute divergent zonal velocity ud in both layers:
call deriv(numLatGridPts,numLongGridPts,rk,wka,ud1) 
call deriv(numLatGridPts,numLongGridPts,rk,wkb,ud2) 
do m=1,numLongGridPts
  ud1(:,m)=secantLat*ud1(:,m)
  ud2(:,m)=secantLat*ud2(:,m)
enddo
 !clati = 1/r, where r = cos(lat)

 !Compute total PV in each layer in physical space:
qAnom1=qSpec1
call revfft(numLatGridPts,numLongGridPts,qAnom1,trig,factors) 
qAnom2=qSpec2
call revfft(numLatGridPts,numLongGridPts,qAnom2,trig,factors) 

do i=1,numLongGridPts
  qAnom1(:,i)=qAnom1(:,i)+cof
  qAnom2(:,i)=qAnom2(:,i)+cof
enddo
 !Note, qAnom1 & qAnom2 will be adjusted by a constant offset once the
 !height anomalies in each layer are determined; the offset
 !guarantees that the domain integral of the relative vorticity
 !is zero, as required by Stokes' theorem.

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Iteratively solve for h, u & v in both layers:

 !Energy norm error (must be > tole to start):
enorm=one
do while (enorm > tole)
   !Save previous fields to be able to compute error below:
  h1o=h1
  h2o=h2
  u1o=u1
  u2o=u2
  v1o=v1
  v2o=v2

   !Correct average PV by enforcing zero average vorticity:
  wka=qAnom1*(one+h1)
  qAnom1=qAnom1-average(wka)

  wka=qAnom2*(one+h2)
  qAnom2=qAnom2-average(wka)

   !Compute relative vorticity (now guaranteed to have zero average):
  do i=1,numLongGridPts
    z1(:,i)=qAnom1(:,i)*(one+h1(:,i))-cof
    z2(:,i)=qAnom2(:,i)*(one+h2(:,i))-cof
  enddo

   !-------------------------------------------------------------
   !Complete calculation of zonal and meridional velocity, u & v,
   !and set up rhs for inversion of dimensionless height anomaly:

   !::::::::::::::::::::::::::::::::::::::::::::::::::::::
   !Layer 1:
  
   !Convert the relative vorticity to semi-spectral space and de-alias:
  call dealias(z1)

   !Invert Laplace's operator on relative vorticity:
  call laplinv(z1,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Compute d(wka)/dlon = wkc:
  call deriv(numLatGridPts,numLongGridPts,rk,wka,wkc)

   !Convert h1 to semi-spectral space:
  call forfft(numLatGridPts,numLongGridPts,h1,trig,factors)

  do m=1,numLongGridPts
    u1(:,m)=ud1(:,m)-wkb(:,m)
    v1(:,m)=vd1(:,m)+secantLat*wkc(:,m)
    h1(:,m)=gammaSpec1(:,m)+beta*u1(:,m)-cof*(z1(:,m)-cof*h1(:,m))
  enddo
   !All the terms in this loop are semi-spectral. Note, we will solve
   !f^2*h_1 - Lap(c_1^2*h_1+(rho_2/rho_1)*c_2^2*h_2) = h1, which follows
   !from the definition of gamma_1.

   !::::::::::::::::::::::::::::::::::::::::::::::::::::::
   !Layer 2:
  
   !Convert the relative vorticity to semi-spectral space and de-alias:
  call dealias(z2)
  
   !Invert Laplace's operator on relative vorticity:
  call laplinv(z2,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Compute d(wka)/dlon = wkc:
  call deriv(numLatGridPts,numLongGridPts,rk,wka,wkc)

   !Convert h2 to semi-spectral space:
  call forfft(numLatGridPts,numLongGridPts,h2,trig,factors)

  do m=1,numLongGridPts
    u2(:,m)=ud2(:,m)-wkb(:,m)
    v2(:,m)=vd2(:,m)+secantLat*wkc(:,m)
    h2(:,m)=gammaSpec2(:,m)+beta*u2(:,m)-cof*(z2(:,m)-cof*h2(:,m))
  enddo
   !All the terms in this loop are semi-spectral. Note, we will solve
   !f^2*h_2 - Lap(c_1^2*h_1+c_2^2*h_2) = h2, which follows from
   !the definition of gamma_2.

   !-------------------------------------------------------------------
   !Solve for h1 & h2 in mode space, where the equations are decoupled.

   !Mode 1:
  wka=h1*vm11+h2*vm12
  call helm1(wka,hm1,wkb)
   !Solves f^2*h_1 - \tilde{c}_1^2*Lap(h_1) = wka, where h_1 is h for mode 1
   !The solution is stored in hm1.
  
   !Mode 2:
  wka=h1*vm21+h2*vm22
  call helm2(wka,hm2,wkb)
   !Solves f^2*h_2 - \tilde{c}_2^2*Lap(h_2) = wka, where h_2 is h for mode 2
   !The solution is stored in hm2.

   !Return h1 & h2 to layers as semi-spectral arrays, then physical ones:
  h1=hm1*vt11+hm2*vt12
  call revfft(numLatGridPts,numLongGridPts,h1,trig,factors) 
  h2=hm1*vt21+hm2*vt22
  call revfft(numLatGridPts,numLongGridPts,h2,trig,factors) 

   !Return semi-spectral velocity field to physical space:
  call revfft(numLatGridPts,numLongGridPts,u1,trig,factors)
  call revfft(numLatGridPts,numLongGridPts,u2,trig,factors)

  call revfft(numLatGridPts,numLongGridPts,v1,trig,factors)
  call revfft(numLatGridPts,numLongGridPts,v2,trig,factors)

   !Compute relative error (measured by the energy norm):
  wka=(u1-u1o)**2+(v1-v1o)**2+csq1*(h1-h1o)**2 &
       +mubar*((u2-u2o)**2+(v2-v2o)**2+csq2*(h2-h2o)**2)
  wkb=u1**2+v1**2+csq1*h1**2+mubar*(u2**2+v2**2+csq2*h2**2)
  enorm=average(wka)/(average(wkb)+small)
   !Small = 1.d-12 is added to avoid dividing by zero for a rest state

   !Employ relaxation to improve convergence:
  h1=wgt*h1+wgtc*h1o
  h2=wgt*h2+wgtc*h2o

  u1=wgt*u1+wgtc*u1o
  u2=wgt*u2+wgtc*u2o

  v1=wgt*v1+wgtc*v1o
  v2=wgt*v2+wgtc*v2o

enddo
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !Passing this point, we have converged!

!------------------------------------------------------------------
 !Return semi-spectral relative vorticity field to physical space:
call revfft(numLatGridPts,numLongGridPts,z1,trig,factors)
call revfft(numLatGridPts,numLongGridPts,z2,trig,factors)

 !Return the PV anomaly in qAnom1 & qAnom2:
do i=1,numLongGridPts
  qAnom1(:,i)=qAnom1(:,i)-cof
  qAnom2(:,i)=qAnom2(:,i)-cof
enddo

return
end subroutine main_invert

!===================================================================

subroutine laplinv(src,sol,der)
! Inverts Laplace's operator on src to give sol and its
! latitudinal derivative der.  That is Lap(sol) = src
! and der = d(sol)/dphi.
! *** All fields are in semi-spectral space ***

implicit none

 !Passed variables:
double precision:: src(numLatGridPts,numLongGridPts),sol(numLatGridPts,numLongGridPts),der(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: rhs(numLatGridPts),pfg(numLatGridPts),const,pm,utdb,vtdb
integer:: j,m

!---------------------------------------------------------------
 !Solve for the zonal part of der (ensure src has zero average):
rhs=src(:,1)*cosLat
const=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
rhs=rhs-const*cosLat

pfg(1)=dl24*(21.d0*rhs(1)+rhs(2))
do j=2,ngm1
  pfg(j)=pfg(j-1)+dl24*(rhs(j-1)+22.d0*rhs(j)+rhs(j+1))
enddo

rhs(1)=pfg(1)
do j=2,ngm1
  rhs(j)=pfg(j)+pfg(j-1)
enddo
rhs(numLatGridPts)=pfg(ngm1)

der(1,1)=rhs(1)*htd0(1)

do j=2,numLatGridPts
  der(j,1)=(rhs(j)-f14*der(j-1,1))*htd0(j)
enddo

do j=ngm1,1,-1
  der(j,1)=etd0(j)*der(j+1,1)+der(j,1)
enddo

der(:,1)=der(:,1)*secantLat

!-------------------------------------------------------
 !Solve for the zonal part of sol:
pfg(1)=dl24*(21.d0*der(1,1)+der(2,1))
do j=2,ngm1
  pfg(j)=pfg(j-1)+dl24*(der(j-1,1)+22.d0*der(j,1)+der(j+1,1))
enddo
pfg(numLatGridPts)=pfg(ngm1)+dl24*(der(ngm1,1)+21.d0*der(numLatGridPts,1))

rhs(1)=pfg(1)
do j=2,numLatGridPts
  rhs(j)=pfg(j)+pfg(j-1)
enddo

sol(1,1)=rhs(1)*htd0(1)

do j=2,numLatGridPts
  sol(j,1)=(rhs(j)-f14*sol(j-1,1))*htd0(j)
enddo

do j=ngm1,1,-1
  sol(j,1)=etd0(j)*sol(j+1,1)+sol(j,1)
enddo

 !Remove global mean value of sol:
rhs=sol(:,1)*cosLat
const=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sol(:,1)=sol(:,1)-const

!------------------------------------------------------------------
 !Loop over azonal longitudinal wavenumbers and solve the 2x2 block
 !tridiagonal problem:
do m=2,numLongGridPts
  pm=psig(ksig(m))

  rhs(1)=pm*src(1,m)+0.1d0*src(2,m)
  do j=2,ngm1
    rhs(j)=src(j,m)+0.1d0*(src(j-1,m)+src(j+1,m))
  enddo
  rhs(numLatGridPts)=0.1d0*src(ngm1,m)+pm*src(numLatGridPts,m)

  sol(1,m)=-rhs(1)*bla0(1,m)
  der(1,m)= rhs(1)*ala0(1,m)

  do j=2,numLatGridPts
    utdb=      -sol(j-1,m)*amla     -der(j-1,m)*0.2d0
    vtdb=rhs(j)-sol(j-1,m)*cmla(j,m)-der(j-1,m)*dmla(j)
    sol(j,m)=utdb*dla0(j,m)-vtdb*bla0(j,m)
    der(j,m)=vtdb*ala0(j,m)-utdb*cla0(j,m)
  enddo

  do j=ngm1,1,-1
    sol(j,m)=ela1(j,m)*sol(j+1,m)+ela2(j,m)*der(j+1,m)+sol(j,m)
    der(j,m)=fla1(j,m)*sol(j+1,m)+fla2(j,m)*der(j+1,m)+der(j,m)
  enddo
enddo

return
end subroutine laplinv

!==========================================================================

subroutine simp1(src,sol,der)
! Inverts the operator L = R^2 + f^2 - \tilde{c}_1^2*Lap appearing in
! the semi-implicit time-stepping method on src to give sol and its
! latitudinal derivative der for *mode 1*. Here, R = 2/dt + D, where
! D is the longitudinal part of the hyperviscosity.
!     Explicitly:  L(sol) = src and der = d(sol)/dphi.
! *** All fields are in semi-spectral space ***

implicit none

 !Passed variables:
double precision:: src(numLatGridPts,numLongGridPts),sol(numLatGridPts,numLongGridPts),der(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: rhs(numLatGridPts),pm,utdb,vtdb,avgsol
integer:: j,m

!-----------------------------------------------------------
 !Loop over longitudinal wavenumbers and solve the 2x2 block
 !tridiagonal problem:
do m=1,numLongGridPts
  pm=psig(ksig(m))

  rhs(1)=pm*src(1,m)+0.1d0*src(2,m)
  do j=2,ngm1
    rhs(j)=src(j,m)+0.1d0*(src(j-1,m)+src(j+1,m))
  enddo
  rhs(numLatGridPts)=0.1d0*src(ngm1,m)+pm*src(numLatGridPts,m)

   !This is needed to be able to use the stored coefficients:
  rhs=-ctsq1i*rhs
  sol(1,m)=-rhs(1)*b1si0(1,m)
  der(1,m)= rhs(1)*a1si0(1,m)

  do j=2,numLatGridPts
    utdb=      -sol(j-1,m)*amla      -der(j-1,m)*0.2d0
    vtdb=rhs(j)-sol(j-1,m)*cm1si(j,m)-der(j-1,m)*dm1si(j)
    sol(j,m)=utdb*d1si0(j,m)-vtdb*b1si0(j,m)
    der(j,m)=vtdb*a1si0(j,m)-utdb*c1si0(j,m)
  enddo

  do j=ngm1,1,-1
    sol(j,m)=e1si1(j,m)*sol(j+1,m)+e1si2(j,m)*der(j+1,m)+sol(j,m)
    der(j,m)=f1si1(j,m)*sol(j+1,m)+f1si2(j,m)*der(j+1,m)+der(j,m)
  enddo
enddo

 !Remove global mean value of sol:
rhs=sol(:,1)*cosLat
avgsol=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sol(:,1)=sol(:,1)-avgsol

return
end subroutine simp1

!==========================================================================

subroutine simp2(src,sol,der)
! Inverts the operator L = R^2 + f^2 - \tilde{c}_2^2*Lap appearing in
! the semi-implicit time-stepping method on src to give sol and its
! latitudinal derivative der for *mode 1*. Here, R = 2/dt + D, where
! D is the longitudinal part of the hyperviscosity.
!     Explicitly:  L(sol) = src and der = d(sol)/dphi.
! *** All fields are in semi-spectral space ***

implicit none

 !Passed variables:
double precision:: src(numLatGridPts,numLongGridPts),sol(numLatGridPts,numLongGridPts),der(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: rhs(numLatGridPts),pm,utdb,vtdb,avgsol
integer:: j,m

!-----------------------------------------------------------
 !Loop over longitudinal wavenumbers and solve the 2x2 block
 !tridiagonal problem:
do m=1,numLongGridPts
  pm=psig(ksig(m))

  rhs(1)=pm*src(1,m)+0.1d0*src(2,m)
  do j=2,ngm1
    rhs(j)=src(j,m)+0.1d0*(src(j-1,m)+src(j+1,m))
  enddo
  rhs(numLatGridPts)=0.1d0*src(ngm1,m)+pm*src(numLatGridPts,m)

   !This is needed to be able to use the stored coefficients:
  rhs=-ctsq2i*rhs
  sol(1,m)=-rhs(1)*b2si0(1,m)
  der(1,m)= rhs(1)*a2si0(1,m)

  do j=2,numLatGridPts
    utdb=      -sol(j-1,m)*amla      -der(j-1,m)*0.2d0
    vtdb=rhs(j)-sol(j-1,m)*cm2si(j,m)-der(j-1,m)*dm2si(j)
    sol(j,m)=utdb*d2si0(j,m)-vtdb*b2si0(j,m)
    der(j,m)=vtdb*a2si0(j,m)-utdb*c2si0(j,m)
  enddo

  do j=ngm1,1,-1
    sol(j,m)=e2si1(j,m)*sol(j+1,m)+e2si2(j,m)*der(j+1,m)+sol(j,m)
    der(j,m)=f2si1(j,m)*sol(j+1,m)+f2si2(j,m)*der(j+1,m)+der(j,m)
  enddo
enddo

 !Remove global mean value of sol:
rhs=sol(:,1)*cosLat
avgsol=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sol(:,1)=sol(:,1)-avgsol

return
end subroutine simp2

!==========================================================================

subroutine helm1(src,sol,der)
! Inverts the operator L = f^2 - \tilde{c}_1^2*Lap appearing in the
! equation determining the height anomaly for *mode 1* on src to give
! sol and its latitudinal derivative der.
! Explicitly, L(sol) = src and der = d(sol)/dphi.
! *** All fields are in semi-spectral space ***

implicit none

 !Passed variables:
double precision:: src(numLatGridPts,numLongGridPts),sol(numLatGridPts,numLongGridPts),der(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: rhs(numLatGridPts),pm,utdb,vtdb,avgsol
integer:: j,m

!-----------------------------------------------------------
 !Loop over longitudinal wavenumbers and solve the 2x2 block
 !tridiagonal problem:
do m=1,numLongGridPts
  pm=psig(ksig(m))

  rhs(1)=pm*src(1,m)+0.1d0*src(2,m)
  do j=2,ngm1
    rhs(j)=src(j,m)+0.1d0*(src(j-1,m)+src(j+1,m))
  enddo
  rhs(numLatGridPts)=0.1d0*src(ngm1,m)+pm*src(numLatGridPts,m)

   !This is needed to be able to use the stored coefficients:
  rhs=-ctsq1i*rhs
  sol(1,m)=-rhs(1)*b1he0(1,m)
  der(1,m)= rhs(1)*a1he0(1,m)

  do j=2,numLatGridPts
    utdb=      -sol(j-1,m)*amla      -der(j-1,m)*0.2d0
    vtdb=rhs(j)-sol(j-1,m)*cm1he(j,m)-der(j-1,m)*dm1he(j)
    sol(j,m)=utdb*d1he0(j,m)-vtdb*b1he0(j,m)
    der(j,m)=vtdb*a1he0(j,m)-utdb*c1he0(j,m)
  enddo

  do j=ngm1,1,-1
    sol(j,m)=e1he1(j,m)*sol(j+1,m)+e1he2(j,m)*der(j+1,m)+sol(j,m)
    der(j,m)=f1he1(j,m)*sol(j+1,m)+f1he2(j,m)*der(j+1,m)+der(j,m)
  enddo
enddo

 !Remove global mean value of sol:
rhs=sol(:,1)*cosLat
avgsol=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sol(:,1)=sol(:,1)-avgsol

return
end subroutine helm1

!==========================================================================

subroutine helm2(src,sol,der)
! Inverts the operator L = f^2 - \tilde{c}_2^2*Lap appearing in the
! equation determining the height anomaly for *mode 2* on src to give
! sol and its latitudinal derivative der.
! Explicitly, L(sol) = src and der = d(sol)/dphi.
! *** All fields are in semi-spectral space ***

implicit none

 !Passed variables:
double precision:: src(numLatGridPts,numLongGridPts),sol(numLatGridPts,numLongGridPts),der(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: rhs(numLatGridPts),pm,utdb,vtdb,avgsol
integer:: j,m

!-----------------------------------------------------------
 !Loop over longitudinal wavenumbers and solve the 2x2 block
 !tridiagonal problem:
do m=1,numLongGridPts
  pm=psig(ksig(m))

  rhs(1)=pm*src(1,m)+0.1d0*src(2,m)
  do j=2,ngm1
    rhs(j)=src(j,m)+0.1d0*(src(j-1,m)+src(j+1,m))
  enddo
  rhs(numLatGridPts)=0.1d0*src(ngm1,m)+pm*src(numLatGridPts,m)

   !This is needed to be able to use the stored coefficients:
  rhs=-ctsq2i*rhs
  sol(1,m)=-rhs(1)*b2he0(1,m)
  der(1,m)= rhs(1)*a2he0(1,m)

  do j=2,numLatGridPts
    utdb=      -sol(j-1,m)*amla      -der(j-1,m)*0.2d0
    vtdb=rhs(j)-sol(j-1,m)*cm2he(j,m)-der(j-1,m)*dm2he(j)
    sol(j,m)=utdb*d2he0(j,m)-vtdb*b2he0(j,m)
    der(j,m)=vtdb*a2he0(j,m)-utdb*c2he0(j,m)
  enddo

  do j=ngm1,1,-1
    sol(j,m)=e2he1(j,m)*sol(j+1,m)+e2he2(j,m)*der(j+1,m)+sol(j,m)
    der(j,m)=f2he1(j,m)*sol(j+1,m)+f2he2(j,m)*der(j+1,m)+der(j,m)
  enddo
enddo

 !Remove global mean value of sol:
rhs=sol(:,1)*cosLat
avgsol=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sol(:,1)=sol(:,1)-avgsol

return
end subroutine helm2

!==========================================================================

subroutine latder(var,der)
! Calculates the first latitudinal derivative of var and stores
! the result in der.  That is, der = d(var)/dphi.
! *** Both var and der are in physical space ***

implicit none

 !Passed variables:
double precision:: var(numLatGridPts,numLongGridPts),der(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: src(numLongGridPts),sol(numLongGridPts),solend
integer:: i,ic,j

!----------------------------------------------------------------
 !Loop over great circles in latitude (only half the longitudes):
do i=1,numLatGridPts
  ic=i+numLatGridPts

   !Source vector:
  src(1)=var(2,i)-var(1,ic)
  do j=2,ngm1
    src(j)=var(j+1,i)-var(j-1,i)
  enddo
  src(numLatGridPts)  =var(numLatGridPts,ic)-var(ngm1,i)
  src(ngp1)=var(ngm1,ic)-var(numLatGridPts,i)
  do j=ngp2,ntm1
    src(j)=var(numLongGridPts-j,ic)-var(ntp2-j,ic)
  enddo
  src(numLongGridPts)=var(1,i)-var(2,ic)

   !Get solution by periodic tridiagonal solve:
  sol(1)=src(1)*htd1(1)

  do j=2,numLongGridPts
    sol(j)=(src(j)-apd1*sol(j-1))*htd1(j)
  enddo

  do j=ntm2,1,-1
    sol(j)=etd1(j)*sol(j+1)+sol(j)
  enddo

  sol(numLongGridPts)=(etd1(numLongGridPts)*sol(1)+sol(numLongGridPts))*denod1
  solend=sol(numLongGridPts)

  do j=1,ntm1
    sol(j)=ptd1(j)*solend+sol(j)
  enddo

   !Return sol in der:
  do j=1,numLatGridPts
    der(j,i)=sol(j)
  enddo
  do j=ngp1,numLongGridPts
    der(ntp1-j,ic)=-sol(j)
  enddo
enddo

return
end subroutine latder

!==========================================================================

subroutine laplace(var,sol)
! Calculates the Laplacian of var and stores the result in sol.
! That is, sol = Lap(var).
! *** Both var and sol are in semi-spectral space ***

implicit none

 !Passed variables:
double precision:: var(numLatGridPts,numLongGridPts),sol(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: wka(numLatGridPts,numLongGridPts),der(numLatGridPts,numLongGridPts)
double precision:: src1(numLongGridPts),src2(numLongGridPts),sol1(numLongGridPts),sol2(numLongGridPts)
double precision:: sol1end,sol2end
integer:: i,ic,j

!----------------------------------------------------------------------
 !Convert var to physical space as wka to take latitudinal derivatives:
wka=var
call revfft(numLatGridPts,numLongGridPts,wka,trig,factors)

!----------------------------------------------------------------
 !Loop over great circles in latitude (only half the longitudes):
do i=1,numLatGridPts
  ic=i+numLatGridPts

   !Source vectors:
  src1(1)=wka(2,i)-wka(1,ic)
  src2(1)=wka(2,i)-two*wka(1,i)+wka(1,ic)
  do j=2,ngm1
    src1(j)=wka(j+1,i)-wka(j-1,i)
    src2(j)=wka(j+1,i)-two*wka(j,i)+wka(j-1,i)
  enddo
  src1(numLatGridPts)  =wka(numLatGridPts,ic)-wka(ngm1,i)
  src2(numLatGridPts)  =wka(numLatGridPts,ic) -two*wka(numLatGridPts,i)+wka(ngm1,i)
  src1(ngp1)=wka(ngm1,ic)-wka(numLatGridPts,i)
  src2(ngp1)=wka(ngm1,ic)-two*wka(numLatGridPts,ic)+wka(numLatGridPts,i)
  do j=ngp2,ntm1
    src1(j)=wka(numLongGridPts-j,ic)-wka(ntp2-j,ic)
    src2(j)=wka(numLongGridPts-j,ic)-two*wka(ntp1-j,ic)+wka(ntp2-j,ic)
  enddo
  src1(numLongGridPts)=wka(1,i)-wka(2,ic)
  src2(numLongGridPts)=wka(1,i)-two*wka(1,ic)+wka(2,ic)

   !Get 1st and 2nd derivatives by periodic tridiagonal solve:
  sol1(1)=src1(1)*htd1(1)
  sol2(1)=src2(1)*htd2(1)

  do j=2,numLongGridPts
    sol1(j)=(src1(j)-apd1*sol1(j-1))*htd1(j)
    sol2(j)=(src2(j)-apd2*sol2(j-1))*htd2(j)
  enddo

  do j=ntm2,1,-1
    sol1(j)=etd1(j)*sol1(j+1)+sol1(j)
    sol2(j)=etd2(j)*sol2(j+1)+sol2(j)
  enddo

  sol1(numLongGridPts)=(etd1(numLongGridPts)*sol1(1)+sol1(numLongGridPts))*denod1
  sol2(numLongGridPts)=(etd2(numLongGridPts)*sol2(1)+sol2(numLongGridPts))*denod2

  sol1end=sol1(numLongGridPts)
  sol2end=sol2(numLongGridPts)

  do j=1,ntm1
    sol1(j)=ptd1(j)*sol1end+sol1(j)
    sol2(j)=ptd2(j)*sol2end+sol2(j)
  enddo

   !Combine derivatives (tlat = tan(lat) below):
  do j=1,numLatGridPts
    der(j,i)=sol2(j)-tanLat(j)*sol1(j)
  enddo
  do j=ngp1,numLongGridPts
    der(ntp1-j,ic)=sol2(j)+tanLat(ntp1-j)*sol1(j)
  enddo
enddo

 !Return der to spectral space and add longitudinal part of operator:
call forfft(numLatGridPts,numLongGridPts,der,trig,factors)
sol=der-d2lon*var

return
end subroutine laplace

!==========================================================================

subroutine latdamp(var)
! Applies a latitudinal hyperviscous damping to the field var.
! *** var is assumed to be semi-spectral 
! (physical latitude, spectral longitude)

implicit none

 !Passed variable:
double precision:: var(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: tra(numLatGridPts,numLongGridPts)
integer:: i,ic,j,k

!---------------------------------------------------
 !Inverse FFT in longitude:
call revfft(numLatGridPts,numLongGridPts,var,trig,factors) 

 !Create great circles:
do i=1,numLatGridPts
  ic=i+numLatGridPts
  do j=1,numLatGridPts
    tra(i,j)     =var(j,i)
    tra(i,ntp1-j)=var(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(numLatGridPts,numLongGridPts,tra,trig,factors) 

 !Apply spectral damping operator in latitude (defined in adapt):
do k=1,numLongGridPts
  tra(:,k)=tra(:,k)*hyplat(k)
enddo

 !Inverse FFT in latitude:
call revfft(numLatGridPts,numLongGridPts,tra,trig,factors) 

 !Unpack array:
do i=1,numLatGridPts
  ic=i+numLatGridPts
  do j=1,numLatGridPts
    var(j,i) =tra(i,j)
    var(j,ic)=tra(i,ntp1-j)
  enddo
enddo

 !FFT in longitude:
call forfft(numLatGridPts,numLongGridPts,var,trig,factors) 

return

end subroutine latdamp

!==========================================================================

subroutine lopass(var)
! Applies a low-pass 2nd-order Butterworth filter to the field var
! (individually to latitudinal modes - around great circles -
!  then to longitudinal modes).

! var is in semi-spectral space (physical latitude, spectral longitude)

implicit none

 !Passed variable:
double precision:: var(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: tra(numLatGridPts,numLongGridPts)
integer:: i,ic,j,k

!---------------------------------------------------
call revfft(numLatGridPts,numLongGridPts,var,trig,factors) 

 !Create great circles:
do i=1,numLatGridPts
  ic=i+numLatGridPts
  do j=1,numLatGridPts
    tra(i,j)     =var(j,i)
    tra(i,ntp1-j)=var(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(numLatGridPts,numLongGridPts,tra,trig,factors) 

 !Apply latitudinal filter:
do k=1,numLongGridPts
  tra(:,k)=tra(:,k)*flatlo(k)
enddo

 !Inverse FFT in latitude:
call revfft(numLatGridPts,numLongGridPts,tra,trig,factors) 

 !Unpack array:
do i=1,numLatGridPts
  ic=i+numLatGridPts
  do j=1,numLatGridPts
    var(j,i) =tra(i,j)
    var(j,ic)=tra(i,ntp1-j)
  enddo
enddo

 !FFT in longitude:
call forfft(numLatGridPts,numLongGridPts,var,trig,factors) 

 !Apply longitudinal filter:
var=var*flonlo

return
end subroutine

!==========================================================================

subroutine hipass(var)
! Applies a high-pass 2nd-order Butterworth filter to the field var
! (individually to latitudinal modes - around great circles -
!  then to longitudinal modes).

! var is in semi-spectral space (physical latitude, spectral longitude)

implicit none

 !Passed variable:
double precision:: var(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: tra(numLatGridPts,numLongGridPts)
integer:: i,ic,j,k

!---------------------------------------------------
call revfft(numLatGridPts,numLongGridPts,var,trig,factors) 

 !Create great circles:
do i=1,numLatGridPts
  ic=i+numLatGridPts
  do j=1,numLatGridPts
    tra(i,j)     =var(j,i)
    tra(i,ntp1-j)=var(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(numLatGridPts,numLongGridPts,tra,trig,factors) 

 !Apply latitudinal filter:
do k=1,numLongGridPts
  tra(:,k)=tra(:,k)*flathi(k)
enddo

 !Inverse FFT in latitude:
call revfft(numLatGridPts,numLongGridPts,tra,trig,factors) 

 !Unpack array:
do i=1,numLatGridPts
  ic=i+numLatGridPts
  do j=1,numLatGridPts
    var(j,i) =tra(i,j)
    var(j,ic)=tra(i,ntp1-j)
  enddo
enddo

 !FFT in longitude:
call forfft(numLatGridPts,numLongGridPts,var,trig,factors) 

 !Apply longitudinal filter:
var=var*flonhi

return
end subroutine hipass

!==========================================================================

subroutine dealias(var)
! De-aliases a field var by removing the upper 1/3 of wavenumber in
! longitude and in latitude (around great circles)

! ===> on input, var is in physical space; 
!     on output, var is in semi-spectral space
!     (physical latitude, spectral longitude)

implicit none

 !Passed variable:
double precision:: var(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: tra(numLatGridPts,numLongGridPts)
integer:: i,ic,j

!---------------------------------------------------
 !Create great circles:
do i=1,numLatGridPts
  ic=i+numLatGridPts
  do j=1,numLatGridPts
    tra(i,j)     =var(j,i)
    tra(i,ntp1-j)=var(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(numLatGridPts,numLongGridPts,tra,trig,factors) 

 !Truncate spectrally in latitude:
tra(:,kda1:kda2)=zero

 !Inverse FFT in latitude:
call revfft(numLatGridPts,numLongGridPts,tra,trig,factors) 

 !Unpack array:
do i=1,numLatGridPts
  ic=i+numLatGridPts
  do j=1,numLatGridPts
    var(j,i) =tra(i,j)
    var(j,ic)=tra(i,ntp1-j)
  enddo
enddo

 !FFT in longitude:
call forfft(numLatGridPts,numLongGridPts,var,trig,factors) 

 !Truncate spectrally in longitude:
var(:,kda1:kda2)=zero

return
end subroutine dealias

!=======================================================================

subroutine getrms(var,vrms)
! Computes the rms value of a gridded field var.
! *** Note wka is used as a work array ***

implicit none

 !Passed variables:
double precision:: var(numLatGridPts,numLongGridPts),vrms

 !Local variables:
double precision:: wka(numLatGridPts,numLongGridPts)
integer:: i

 !-------------------------------------------------
do i=1,numLongGridPts
  wka(:,i)=cosLat*var(:,i)**2
enddo

vrms=zero
do i=1,numLongGridPts
  vrms=vrms+f1112*(wka(1,i)+wka(numLatGridPts,i))+sum(wka(2:ngm1,i))
enddo
vrms=sqrt(vrms*dsumi)

return
end subroutine getrms

!==========================================================================

subroutine zeroavg(var)
! Removes the mean of a gridded field var.
! *** Note wka is used as a work array ***

implicit none

 !Passed variables:
double precision:: var(numLatGridPts,numLongGridPts),vsum

 !Local variables:
double precision:: wka(numLatGridPts,numLongGridPts)
integer:: i

 !-----------------------------------------------------
do i=1,numLongGridPts
  wka(:,i)=cosLat*var(:,i)
enddo

 !Compute 4th-order average:
vsum=zero
do i=1,numLongGridPts
  vsum=vsum+f1112*(wka(1,i)+wka(numLatGridPts,i))+sum(wka(2:ngm1,i))
enddo
vsum=vsum*dsumi

 !Remove mean:
var=var-vsum

return
end subroutine zeroavg

!==========================================================================

double precision function average(var)
! Returns the average value of a gridded field var.
! *** Note wka is used as a work array ***

implicit none

 !Passed array:
double precision:: var(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: wka(numLatGridPts,numLongGridPts),vsum
integer:: i

 !---------------------------------------------------
do i=1,numLongGridPts
  wka(:,i)=cosLat*var(:,i)
enddo

 !Compute 4th-order average:
vsum=zero
do i=1,numLongGridPts
  vsum=vsum+f1112*(wka(1,i)+wka(numLatGridPts,i))+sum(wka(2:ngm1,i))
enddo
average=vsum*dsumi

return
end function average

!==========================================================================

subroutine longspec(var,pow)
! Computes the longitudinal spectrum of a field var and creates
! pow, containing the power spectrum for each longitudinal wavenumber
! m from 1 to numLatGridPts+1

! *** var must be in semi-spectral space

implicit none

 !Passed variables:
double precision:: var(numLatGridPts,numLongGridPts),pow(numLatGridPts)

 !Local variables:
integer:: m,mc

!-------------------------------------------------------------
 !Compute power spectrum excluding global average in var(:,1):
do m=2,numLatGridPts
  mc=ntp2-m
  pow(m-1)=sum(cosLat*(var(:,m)**2+var(:,mc)**2))
enddo

pow(numLatGridPts)=sum(cosLat*var(:,ngp1)**2)

return
end subroutine longspec

!==========================================================================

end module spectral
