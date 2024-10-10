module spectral

! Module containing subroutines for spectral operations, inversion, etc.

use constants
use stafft
use deriv1d

implicit none

 !FFT commons:
integer:: factors(5)
double precision:: trig(2*nt),rk(nt),wave(nt)

 !Filters and damping operators:
double precision:: d2lon(ng,nt),diss(ng,nt)
double precision:: rdis(ng,nt),rdisi(ng,nt)
double precision:: plon(ng,nt),glon(ng,nt)
double precision:: flonlo(ng,nt),flonhi(ng,nt)
double precision:: flatlo(nt),flathi(nt),hyplat(nt),plat(nt),glat(nt)

 !Tridiagonal commons:
double precision:: etd1(nt),htd1(nt),ptd1(nt),apd1,denod1
double precision:: etd2(nt),htd2(nt),ptd2(nt),apd2,denod2
double precision:: psig(2),amla
integer:: ksig(nt)

double precision:: cmla(ng,nt),dmla(ng)
double precision:: ala0(ng,nt),bla0(ng,nt)
double precision:: cla0(ng,nt),dla0(ng,nt)
double precision:: ela1(ng,nt),ela2(ng,nt)
double precision:: fla1(ng,nt),fla2(ng,nt)
double precision:: etd0(ng),htd0(ng)
double precision:: rsumi,dsumi,dl24

double precision:: cmhe(ng,nt),dmhe(ng)
double precision:: ahe0(ng,nt),bhe0(ng,nt)
double precision:: che0(ng,nt),dhe0(ng,nt)
double precision:: ehe1(ng,nt),ehe2(ng,nt)
double precision:: fhe1(ng,nt),fhe2(ng,nt)

double precision:: cmpo(ng,nt),dmpo(ng)
double precision:: apo0(ng,nt),bpo0(ng,nt)
double precision:: cpo0(ng,nt),dpo0(ng,nt)
double precision:: epo1(ng,nt),epo2(ng,nt)
double precision:: fpo1(ng,nt),fpo2(ng,nt)

complex*16:: cmsi(ng,ngp1)
complex*16:: asi0(ng,ngp1),bsi0(ng,ngp1)
complex*16:: csi0(ng,ngp1),dsi0(ng,ngp1)
complex*16:: esi1(ng,ngp1),esi2(ng,ngp1)
complex*16:: fsi1(ng,ngp1),fsi2(ng,ngp1)
double precision:: sipf(ng,ngp1),dmsi(ng)

 !Other arrays:
double precision:: tlat(ng),slat(ng),clat(ng),clati(ng),clatsqi(ng)
double precision:: cof(ng),bet(ng),fsq(ng),fri(ng),dfb(ng)

contains

!=======================================================================

subroutine init_spectral
! Initialises spectral module

implicit none

 !Local variables:
double precision:: alpha(ng),sigm(2)
complex*16:: alpcom(ng)
double precision:: rlat,wm,wmi,dnu,fac,rkc,rkci,f283,weff,dfac,rsum
double precision:: a0d1,a0d2,pc12,pc24,a0,b0,c0,d0,deti
double precision:: apla,cpla,dpla,sm,pm
double precision:: cphe,dphe,cppo,dppo,dpsi,xifac
complex*16:: c1,cdeti,cpsi
integer:: i,j,k,m

!----------------------------------------------------------
 !Latitude functions:
do j=1,ng
  rlat=(dble(j)-f12)*dl-hpi
  slat(j)=sin(rlat)         !r
  clat(j)=cos(rlat)         !z
enddo

clati=one/clat              !1/r
clatsqi=clati**2            !1/r^2
tlat=slat/clat              !z/r
cof=fpole*slat              !f = 2*Omega*z
fsq=cof**2                  !f^2
fri=cof*clati               !f/r
bet=fpole*clat              !beta = df/d(lat)
dfb=two*cof*bet             !2*f*beta

!-------------------------------------------
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
wave(ngp1)=rk(nt)
wmi=one/wave(ngp1)          !1/M where M = maximum wavenumber

 !m^2/r^2 factor which appears in Laplace's operator:
do m=1,nt
  fac=wave(m)**2
  d2lon(:,m)=fac*clati**2
enddo

!----------------------------------------------------------------
 !Define 2nd-order Butterworth filter (used in subroutines lopass
 !and hipass) and include partial de-aliasing:
rkc=f13*dble(ng)
 !rkc: filter wavenumber
rkci=one/rkc
f283=28.d0/3.d0

 !Latitude filters:
flatlo=one/(one+(rkci*wave)**4)
flathi=(one-flatlo)*f12*(one-erf(16.d0*wave*wmi-f283))

 !Longitude filters:
do m=1,nt
  do j=1,ng
    weff=wave(m)*clati(j)
    flonlo(j,m)=one/(one+(rkci*weff)**4)
    flonhi(j,m)=(one-flonlo(j,m))*f12*(one-erf(16.d0*weff*wmi-f283))
  enddo
enddo

!--------------------------------------------------------------
 !Latitudinal hyperviscous damping used on residual PV (qd),
 !divergence (ds), and acceleration divergence (gs):
dfac=drate*fpole
fac=dt*dfac
dnu=dble(2*nnu)
hyplat=one/(one+fac*(wmi*wave)**dnu)

 !Longitudinal hyperviscous damping rate on wavenumber m:
wmi=one/(wave(ng+1)-one)
do m=3,nt-1
  fac=wmi*(wave(m)-one)
  do j=1,ng
    diss(j,m)=dfac*(fac*clati(j))**dnu
  enddo
enddo
 !The zonal, cos(lon) and sin(lon) modes are not damped:
diss(:,1)=zero
diss(:,2)=zero
diss(:,nt)=zero

 !The R operator in the notes:
rdis=dt2i+diss
 !R^{-1}:
rdisi=one/rdis

!-------------------------------------
 !Initialise tridiagonal coefficients:

!---------------------------
 !First latitude derivative:
a0d1=dl*f43
apd1=dl*f13

htd1(1)=one/a0d1
ptd1(1)=-apd1*htd1(1)
etd1(1)=ptd1(1)

do j=2,nt
  htd1(j)=one/(a0d1+apd1*etd1(j-1))
  ptd1(j)=-apd1*ptd1(j-1)*htd1(j)
  etd1(j)=-apd1*htd1(j)
enddo

ptd1(ntm1)=etd1(ntm1)+ptd1(ntm1)

do j=ntm2,1,-1
  ptd1(j)=etd1(j)*ptd1(j+1)+ptd1(j)
enddo

denod1=one/(one-etd1(nt)*ptd1(1)-ptd1(nt))

!----------------------------
 !Second latitude derivative:
a0d2=dlsq*f56
apd2=dlsq*f112

htd2(1)=one/a0d2
ptd2(1)=-apd2*htd2(1)
etd2(1)=ptd2(1)

do j=2,nt
  htd2(j)=one/(a0d2+apd2*etd2(j-1))
  ptd2(j)=-apd2*ptd2(j-1)*htd2(j)
  etd2(j)=-apd2*htd2(j)
enddo

ptd2(ntm1)=etd2(ntm1)+ptd2(ntm1)

do j=ntm2,1,-1
  ptd2(j)=etd2(j)*ptd2(j+1)+ptd2(j)
enddo

denod2=one/(one-etd2(nt)*ptd2(1)-ptd2(nt))

!--------------------------------------------------------
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

htd0(ng)=one/(f74+f14*etd0(ngm1))

rsum=f1112*(clat(1)+clat(ng))+sum(clat(2:ngm1))
 !Note, rsum = 1/sin(dl/2) - sin(dl/2)/6 exactly.
rsumi=one/rsum
dsumi=rsumi/dble(nt)

 !azonal part:
do m=2,nt
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  alpha=d2lon(:,m)

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(1)
  d0=-pm*tlat(1)

  deti=one/(a0*d0-b0*c0)
  ala0(1,m)=a0*deti
  bla0(1,m)=b0*deti
  cla0(1,m)=c0*deti
  dla0(1,m)=d0*deti

  cpla=pc12-0.1d0*alpha(2)
  dpla=-0.1d0*tlat(2)
  ela1(1,m)= cpla*bla0(1,m)- apla*dla0(1,m)
  ela2(1,m)= dpla*bla0(1,m)-0.2d0*dla0(1,m)
  fla1(1,m)= apla*cla0(1,m)- cpla*ala0(1,m)
  fla2(1,m)=0.2d0*cla0(1,m)- dpla*ala0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-alpha(j)
    d0=-tlat(j)
    cmla(j,m)=pc12-0.1d0*alpha(j-1)
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

    cpla=pc12-0.1d0*alpha(j+1)
    dpla=-0.1d0*tlat(j+1)
    ela1(j,m)= cpla*bla0(j,m)- apla*dla0(j,m)
    ela2(j,m)= dpla*bla0(j,m)-0.2d0*dla0(j,m)
    fla1(j,m)= apla*cla0(j,m)- cpla*ala0(j,m)
    fla2(j,m)=0.2d0*cla0(j,m)- dpla*ala0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(ng)
  d0=-pm*tlat(ng)

  cmla(ng,m)=pc12-0.1d0*alpha(ngm1)
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

!-----------------------------------------------------------------------
 !Initialise tri-diagonal arrays for semi-implicit time-stepping of
 !divergence and acceleration divergence (include hyperviscous damping):
fac=hbsq3*dt2i
wave(ngp1)=zero   !Needed to treat the Nyquist frequency correctly
do m=1,ngp1
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  sipf(:,m)=-one/(csq*rdisi(:,m)+fac)
  alpcom=d2lon(:,m)-(fsq*rdisi(:,m)+dt2i*cmplx(one,cxi*wave(m)*rdisi(:,m)))*sipf(:,m)

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c1=(sm-two)*pc12-pm*alpcom(1)
  d0=-pm*tlat(1)

  cdeti=one/(a0*d0-b0*c1)
  asi0(1,m)=a0*cdeti
  bsi0(1,m)=b0*cdeti
  csi0(1,m)=c1*cdeti
  dsi0(1,m)=d0*cdeti

  cpsi=pc12-0.1d0*alpcom(2)
  dpsi=-0.1d0*tlat(2)
  esi1(1,m)= cpsi*bsi0(1,m)- apla*dsi0(1,m)
  esi2(1,m)= dpsi*bsi0(1,m)-0.2d0*dsi0(1,m)
  fsi1(1,m)= apla*csi0(1,m)- cpsi*asi0(1,m)
  fsi2(1,m)=0.2d0*csi0(1,m)- dpsi*asi0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c1=-pc24-alpcom(j)
    d0=-tlat(j)
    cmsi(j,m)=pc12-0.1d0*alpcom(j-1)
    dmsi(j)=-0.1d0*tlat(j-1)

    asi0(j,m)=   esi1(j-1,m)*amla      +fsi1(j-1,m)*0.2d0
    bsi0(j,m)=b0+esi2(j-1,m)*amla      +fsi2(j-1,m)*0.2d0
    csi0(j,m)=c1+esi1(j-1,m)*cmsi(j,m)+fsi1(j-1,m)*dmsi(j)
    dsi0(j,m)=d0+esi2(j-1,m)*cmsi(j,m)+fsi2(j-1,m)*dmsi(j)
    cdeti=one/(asi0(j,m)*dsi0(j,m)-bsi0(j,m)*csi0(j,m))
    asi0(j,m)=asi0(j,m)*cdeti
    bsi0(j,m)=bsi0(j,m)*cdeti
    csi0(j,m)=csi0(j,m)*cdeti
    dsi0(j,m)=dsi0(j,m)*cdeti

    cpsi=pc12-0.1d0*alpcom(j+1)
    dpsi=-0.1d0*tlat(j+1)
    esi1(j,m)= cpsi*bsi0(j,m)- apla*dsi0(j,m)
    esi2(j,m)= dpsi*bsi0(j,m)-0.2d0*dsi0(j,m)
    fsi1(j,m)= apla*csi0(j,m)- cpsi*asi0(j,m)
    fsi2(j,m)=0.2d0*csi0(j,m)- dpsi*asi0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c1=(sm-two)*pc12-pm*alpcom(ng)
  d0=-pm*tlat(ng)

  cmsi(ng,m)=pc12-0.1d0*alpcom(ngm1)
  dmsi(ng)=-0.1d0*tlat(ngm1)

  asi0(ng,m)=a0+esi1(ngm1,m)*amla       +fsi1(ngm1,m)*0.2d0
  bsi0(ng,m)=b0+esi2(ngm1,m)*amla       +fsi2(ngm1,m)*0.2d0
  csi0(ng,m)=c1+esi1(ngm1,m)*cmsi(ng,m)+fsi1(ngm1,m)*dmsi(ng)
  dsi0(ng,m)=d0+esi2(ngm1,m)*cmsi(ng,m)+fsi2(ngm1,m)*dmsi(ng)

  cdeti=one/(asi0(ng,m)*dsi0(ng,m)-bsi0(ng,m)*csi0(ng,m))
  asi0(ng,m)=asi0(ng,m)*cdeti
  bsi0(ng,m)=bsi0(ng,m)*cdeti
  csi0(ng,m)=csi0(ng,m)*cdeti
  dsi0(ng,m)=dsi0(ng,m)*cdeti
enddo
wave(ngp1)=rk(nt)

!--------------------------------------------------------
 !Initialise tri-diagonal arrays needed to invert height:
do m=1,nt
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  alpha=d2lon(:,m)+fsq*csqi

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(1)
  d0=-pm*tlat(1)

  deti=one/(a0*d0-b0*c0)
  ahe0(1,m)=a0*deti
  bhe0(1,m)=b0*deti
  che0(1,m)=c0*deti
  dhe0(1,m)=d0*deti

  cphe=pc12-0.1d0*alpha(2)
  dphe=-0.1d0*tlat(2)
  ehe1(1,m)= cphe*bhe0(1,m)- apla*dhe0(1,m)
  ehe2(1,m)= dphe*bhe0(1,m)-0.2d0*dhe0(1,m)
  fhe1(1,m)= apla*che0(1,m)- cphe*ahe0(1,m)
  fhe2(1,m)=0.2d0*che0(1,m)- dphe*ahe0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-alpha(j)
    d0=-tlat(j)
    cmhe(j,m)=pc12-0.1d0*alpha(j-1)
    dmhe(j)=-0.1d0*tlat(j-1)

    ahe0(j,m)=   ehe1(j-1,m)*amla      +fhe1(j-1,m)*0.2d0
    bhe0(j,m)=b0+ehe2(j-1,m)*amla      +fhe2(j-1,m)*0.2d0
    che0(j,m)=c0+ehe1(j-1,m)*cmhe(j,m)+fhe1(j-1,m)*dmhe(j)
    dhe0(j,m)=d0+ehe2(j-1,m)*cmhe(j,m)+fhe2(j-1,m)*dmhe(j)
    deti=one/(ahe0(j,m)*dhe0(j,m)-bhe0(j,m)*che0(j,m))
    ahe0(j,m)=ahe0(j,m)*deti
    bhe0(j,m)=bhe0(j,m)*deti
    che0(j,m)=che0(j,m)*deti
    dhe0(j,m)=dhe0(j,m)*deti

    cphe=pc12-0.1d0*alpha(j+1)
    dphe=-0.1d0*tlat(j+1)
    ehe1(j,m)= cphe*bhe0(j,m)- apla*dhe0(j,m)
    ehe2(j,m)= dphe*bhe0(j,m)-0.2d0*dhe0(j,m)
    fhe1(j,m)= apla*che0(j,m)- cphe*ahe0(j,m)
    fhe2(j,m)=0.2d0*che0(j,m)- dphe*ahe0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(ng)
  d0=-pm*tlat(ng)

  cmhe(ng,m)=pc12-0.1d0*alpha(ngm1)
  dmhe(ng)=-0.1d0*tlat(ngm1)

  ahe0(ng,m)=a0+ehe1(ngm1,m)*amla       +fhe1(ngm1,m)*0.2d0
  bhe0(ng,m)=b0+ehe2(ngm1,m)*amla       +fhe2(ngm1,m)*0.2d0
  che0(ng,m)=c0+ehe1(ngm1,m)*cmhe(ng,m)+fhe1(ngm1,m)*dmhe(ng)
  dhe0(ng,m)=d0+ehe2(ngm1,m)*cmhe(ng,m)+fhe2(ngm1,m)*dmhe(ng)

  deti=one/(ahe0(ng,m)*dhe0(ng,m)-bhe0(ng,m)*che0(ng,m))
  ahe0(ng,m)=ahe0(ng,m)*deti
  bhe0(ng,m)=bhe0(ng,m)*deti
  che0(ng,m)=che0(ng,m)*deti
  dhe0(ng,m)=dhe0(ng,m)*deti
enddo

!--------------------------------------------------------------------
 !Initialise tri-diagonal arrays needed for non-hydrostatic pressure:
do m=1,nt
  k=ksig(m)
  sm=sigm(k)
  pm=psig(k)

  alpha=d2lon(:,m)+hbsq3i

  a0=amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(1)
  d0=-pm*tlat(1)

  deti=one/(a0*d0-b0*c0)
  apo0(1,m)=a0*deti
  bpo0(1,m)=b0*deti
  cpo0(1,m)=c0*deti
  dpo0(1,m)=d0*deti

  cppo=pc12-0.1d0*alpha(2)
  dppo=-0.1d0*tlat(2)
  epo1(1,m)= cppo*bpo0(1,m)- apla*dpo0(1,m)
  epo2(1,m)= dppo*bpo0(1,m)-0.2d0*dpo0(1,m)
  fpo1(1,m)= apla*cpo0(1,m)- cppo*apo0(1,m)
  fpo2(1,m)=0.2d0*cpo0(1,m)- dppo*apo0(1,m)

  b0=0.8d0
  do j=2,ngm1
    c0=-pc24-alpha(j)
    d0=-tlat(j)
    cmpo(j,m)=pc12-0.1d0*alpha(j-1)
    dmpo(j)=-0.1d0*tlat(j-1)

    apo0(j,m)=   epo1(j-1,m)*amla      +fpo1(j-1,m)*0.2d0
    bpo0(j,m)=b0+epo2(j-1,m)*amla      +fpo2(j-1,m)*0.2d0
    cpo0(j,m)=c0+epo1(j-1,m)*cmpo(j,m)+fpo1(j-1,m)*dmpo(j)
    dpo0(j,m)=d0+epo2(j-1,m)*cmpo(j,m)+fpo2(j-1,m)*dmpo(j)
    deti=one/(apo0(j,m)*dpo0(j,m)-bpo0(j,m)*cpo0(j,m))
    apo0(j,m)=apo0(j,m)*deti
    bpo0(j,m)=bpo0(j,m)*deti
    cpo0(j,m)=cpo0(j,m)*deti
    dpo0(j,m)=dpo0(j,m)*deti

    cppo=pc12-0.1d0*alpha(j+1)
    dppo=-0.1d0*tlat(j+1)
    epo1(j,m)= cppo*bpo0(j,m)- apla*dpo0(j,m)
    epo2(j,m)= dppo*bpo0(j,m)-0.2d0*dpo0(j,m)
    fpo1(j,m)= apla*cpo0(j,m)- cppo*apo0(j,m)
    fpo2(j,m)=0.2d0*cpo0(j,m)- dppo*apo0(j,m)
  enddo

  a0=-amla*sm
  b0=0.8d0-0.2d0*sm
  c0=(sm-two)*pc12-pm*alpha(ng)
  d0=-pm*tlat(ng)

  cmpo(ng,m)=pc12-0.1d0*alpha(ngm1)
  dmpo(ng)=-0.1d0*tlat(ngm1)

  apo0(ng,m)=a0+epo1(ngm1,m)*amla       +fpo1(ngm1,m)*0.2d0
  bpo0(ng,m)=b0+epo2(ngm1,m)*amla       +fpo2(ngm1,m)*0.2d0
  cpo0(ng,m)=c0+epo1(ngm1,m)*cmpo(ng,m)+fpo1(ngm1,m)*dmpo(ng)
  dpo0(ng,m)=d0+epo2(ngm1,m)*cmpo(ng,m)+fpo2(ngm1,m)*dmpo(ng)

  deti=one/(apo0(ng,m)*dpo0(ng,m)-bpo0(ng,m)*cpo0(ng,m))
  apo0(ng,m)=apo0(ng,m)*deti
  bpo0(ng,m)=bpo0(ng,m)*deti
  cpo0(ng,m)=cpo0(ng,m)*deti
  dpo0(ng,m)=dpo0(ng,m)*deti
enddo

!--------------------------------------------------------
 !Re-define damping operator for use in qd time stepping:
diss=one/(one+dt2*diss)
 !Also used to damp divergence.

return 
end subroutine init_spectral

!======================================================================

subroutine main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
! Given the PV anomaly qs, divergence ds and acceleration divergence gs
! (all in semi-spectral space), this routine computes the dimensionless
! depth anomaly hh and the velocity field (uu,vv) in physical space.
! It also returns the corrected PV (qq, including f) and the relative
! vorticity (zz) in physical space.

implicit none

 !Passed variables:
double precision:: qs(ng,nt),ds(ng,nt),gs(ng,nt) !Spectral
double precision:: hh(ng,nt),uu(ng,nt),vv(ng,nt) !Physical
double precision:: qq(ng,nt),zz(ng,nt),dd(ng,nt) !Physical

 !Local variables:
double precision,parameter:: tolh=1.d-9
 !tolh: relative error in successive iterates when finding height field

 !Physical work arrays:
double precision:: htot(ng,nt),hs(ng,nt)
double precision:: dda(ng,nt),ddb(ng,nt)
double precision:: wkp(ng,nt),wkq(ng,nt)

 !Spectral work arrays:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt)
double precision:: uds(ng,nt),vds(ng,nt)

 !Other constants:
double precision:: hnorm
integer:: i,m

!------------------------------------------------------------------
 !First compute divergent part of the velocity field (this is fixed
 !in the iteration below):

 !Invert Laplace's operator on the divergence (ds):
call laplinv(ds,wka,vds)
 !Here the divergence potential is wka while d(wka)/dlat = vds,
 !the divergent meridional velocity (all in semi-spectral space)

 !Compute divergent zonal velocity and store in uds:
call deriv(ng,nt,rk,wka,uds) 
do m=1,nt
  uds(:,m)=clati*uds(:,m)
enddo
 !uds = (1/r)*d(wka)/dlon where r = cos(lat)

 !Store scaled divergence derivatives to calculate relative vorticity:
call deriv(ng,nt,rk,ds,dda) 
call revfft(ng,nt,dda,trig,factors)  !dda = d(delta)/d(lon)
wka=ds
call revfft(ng,nt,wka,trig,factors)
call latder(wka,ddb)                 !ddb = d(delta)/d(lat)
do i=1,nt
  dda(:,i)=hbsq3*clati*dda(:,i)
  ddb(:,i)=hbsq3*clati*ddb(:,i)      !hbsq3*clati = (H^2/3)/cos(lat)
enddo

 !Obtain a semi-spectral copy of hh for use in the iteration below:
hs=hh
call forfft(ng,nt,hs,trig,factors) 

 !Obtain a physical space copy of qs:
qq=qs
call revfft(ng,nt,qq,trig,factors) 

 !Add f so that qq contains total PV:
do i=1,nt
  qq(:,i)=qq(:,i)+cof
enddo

 !Define total dimensionless height:
htot=one+hh

!-------------------------------------------------------
 !Iteratively solve for hh, uu & vv:

 !Energy norm error (must be > tolh to start):
hnorm=f12
do while (hnorm .gt. tolh)
   !Correct average PV by enforcing zero average vorticity:
  wka=qq*htot
  qq=qq-average(wka)

   !Calculate (H^2/3)*J(delta,h_tilde):
  call deriv(ng,nt,rk,hs,wkp) 
  call revfft(ng,nt,wkp,trig,factors)  !wkp = d(h_tilde)/d(lon)
  call latder(hh,wkq)                  !wkq = d(h_tilde)/d(lat)
  wkp=dda*wkq-ddb*wkp

   !Compute relative vorticity (now guaranteed to have zero average):
  do i=1,nt
    zz(:,i)=htot(:,i)*qq(:,i)-wkp(:,i)-cof
  enddo
   !Convert to semi-spectral space and dealias to find non-divergent velocity:
  call dealias(zz)

   !Invert Laplace's operator:
  call laplinv(zz,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Compute d(wka)/dlon = wkc:
  call deriv(ng,nt,rk,wka,wkc)

   !Complete calculation of zonal and meridional velocity, uu & vv,
   !and set up rhs for inversion of dimensionless height anomaly:
  do m=1,nt
    uu(:,m)=uds(:,m)-wkb(:,m)
    vv(:,m)=vds(:,m)+clati*wkc(:,m)
    wkb(:,m)=csqi*(cof*(zz(:,m)-cof*hs(:,m))-bet*uu(:,m)-gs(:,m))
  enddo
   !These are in semi-spectral space here in this iterative loop.  
  
   !Invert [f^2-c^2*grad^2]h = gamma + beta*u - f*(zeta - f*h):
  call helminv(wkb,hs,wkc)
   !hs: corrected height field (to be hh below)
   !wkc is not needed but contains d(hh)/dlat.

   !Convert hs to physical space as wka:
  wka=hs
  call revfft(ng,nt,wka,trig,factors) 

   !Compute rms error in hh and re-assign hh & htot:
  wkc=(hh-wka)**2
  hh=wka
  htot=one+hh
  hnorm=sqrt(average(wkc))
enddo
 !Passing this, we have converged!

!------------------------------------------------------------------
 !Get physical space velocity:
call revfft(ng,nt,uu,trig,factors)
call revfft(ng,nt,vv,trig,factors)

 !Return relative vorticity back to physical space:
call revfft(ng,nt,zz,trig,factors)

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
double precision:: src(ng,nt),sol(ng,nt),der(ng,nt)

 !Local variables:
double precision:: rhs(ng),pfg(ng),const,pm,utdb,vtdb
integer:: j,m

!---------------------------------------------------------------
 !Solve for the zonal part of der (ensure src has zero average):
rhs=src(:,1)*clat
const=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
rhs=rhs-const*clat

pfg(1)=dl24*(21.d0*rhs(1)+rhs(2))
do j=2,ngm1
  pfg(j)=pfg(j-1)+dl24*(rhs(j-1)+22.d0*rhs(j)+rhs(j+1))
enddo

rhs(1)=pfg(1)
do j=2,ngm1
  rhs(j)=pfg(j)+pfg(j-1)
enddo
rhs(ng)=pfg(ngm1)

der(1,1)=rhs(1)*htd0(1)

do j=2,ng
  der(j,1)=(rhs(j)-f14*der(j-1,1))*htd0(j)
enddo

do j=ngm1,1,-1
  der(j,1)=etd0(j)*der(j+1,1)+der(j,1)
enddo

der(:,1)=der(:,1)*clati

!-------------------------------------------------------
 !Solve for the zonal part of sol:
pfg(1)=dl24*(21.d0*der(1,1)+der(2,1))
do j=2,ngm1
  pfg(j)=pfg(j-1)+dl24*(der(j-1,1)+22.d0*der(j,1)+der(j+1,1))
enddo
pfg(ng)=pfg(ngm1)+dl24*(der(ngm1,1)+21.d0*der(ng,1))

rhs(1)=pfg(1)
do j=2,ng
  rhs(j)=pfg(j)+pfg(j-1)
enddo

sol(1,1)=rhs(1)*htd0(1)

do j=2,ng
  sol(j,1)=(rhs(j)-f14*sol(j-1,1))*htd0(j)
enddo

do j=ngm1,1,-1
  sol(j,1)=etd0(j)*sol(j+1,1)+sol(j,1)
enddo

 !Remove global mean value of sol:
rhs=sol(:,1)*clat
const=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
sol(:,1)=sol(:,1)-const

!------------------------------------------------------------------
 !Loop over azonal longitudinal wavenumbers and solve the 2x2 block
 !tridiagonal problem:
do m=2,nt
  pm=psig(ksig(m))

  rhs(1)=pm*src(1,m)+0.1d0*src(2,m)
  do j=2,ngm1
    rhs(j)=src(j,m)+0.1d0*(src(j-1,m)+src(j+1,m))
  enddo
  rhs(ng)=0.1d0*src(ngm1,m)+pm*src(ng,m)

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
enddo

return
end subroutine laplinv

!==========================================================================

subroutine simp(dsrc,gsrc,dsol,gsol)
! Solves the semi-implicit system

!      (omega_0*Z - R^{-1}*G)[dsol] = R^{-1}*gsrc + Z[dsrc]
!           gsol = P[omega_0*dsol - dsrc]

! where omega_0 = 2/dt,  G = c^2*Lap - f^2,  Lap = Laplace's operator,
! Z = P + xi*R^{-1}*d/d(lon),  R = omega_0 + D,  P = 1 - (H^2/3)*Lap,
! and xi = 2*Omega*H^2/3.
! *** All fields are in semi-spectral space and are real ***

implicit none

 !Passed variables:
double precision:: dsrc(ng,nt),gsrc(ng,nt),dsol(ng,nt),gsol(ng,nt)

 !Local variables:
complex*16:: src(ng,ngp1),csol(ng,ngp1),cder(ng,ngp1)
complex*16:: rhs(ng),utdb,vtdb
double precision:: wka(ng,nt)
double precision:: pm,avgsol
integer:: j,m

!------------------------------------------------------------------
 !Compute gsrc + Z[dsrc] (use dsol & gsol temporarily):
call project(dsrc,dsol)          !dsol = P[dsrc]
call deriv(ng,nt,rk,dsrc,gsol)   !gsol = d(dsrc)/d(lon)
dsol=dsol+rdisi*(gsrc+cxi*gsol)  !dsol = R^{-1}*gsrc + Z[dsrc]

!------------------------------------------------------------------
 !Pack R^{-1}*gsrc + Z[dsrc] into a complex array (src):
src(:,1)=dsol(:,1)
src(:,ngp1)=dsol(:,ngp1)
do m=2,ng
  src(:,m)=cmplx(dsol(:,m),dsol(:,ntp2-m))
enddo

!------------------------------------------------------------------
 !Multiply src by sipf = -1/(c^2*R^{-1}+(H^2/3)*omega_0) to be able
 !to use the stored tridiagonal coefficients:
src=sipf*src

!------------------------------------------------------------------
 !Loop over longitudinal wavenumbers and solve the 2x2 block
 !tridiagonal problem:
do m=1,ngp1
  pm=psig(ksig(m))

  rhs(1)=pm*src(1,m)+0.1d0*src(2,m)
  do j=2,ngm1
    rhs(j)=src(j,m)+0.1d0*(src(j-1,m)+src(j+1,m))
  enddo
  rhs(ng)=0.1d0*src(ngm1,m)+pm*src(ng,m)

  csol(1,m)=-rhs(1)*bsi0(1,m)
  cder(1,m)= rhs(1)*asi0(1,m)

  do j=2,ng
    utdb=      -csol(j-1,m)*amla     -cder(j-1,m)*0.2d0
    vtdb=rhs(j)-csol(j-1,m)*cmsi(j,m)-cder(j-1,m)*dmsi(j)
    csol(j,m)=utdb*dsi0(j,m)-vtdb*bsi0(j,m)
    cder(j,m)=vtdb*asi0(j,m)-utdb*csi0(j,m)
  enddo

  do j=ngm1,1,-1
    csol(j,m)=esi1(j,m)*csol(j+1,m)+esi2(j,m)*cder(j+1,m)+csol(j,m)
    cder(j,m)=fsi1(j,m)*csol(j+1,m)+fsi2(j,m)*cder(j+1,m)+cder(j,m)
  enddo
enddo
 !cder = d(csol)/d(lat) is not needed further.

!------------------------------------------------------------------
 !Unpack complex solution array (csol) into a real one (dsol):
dsol(:,1)=real(csol(:,1))
dsol(:,ngp1)=real(csol(:,ngp1))
do m=2,ng
  dsol(:,m)=real(csol(:,m))
  dsol(:,ntp2-m)=aimag(csol(:,m))
enddo

!------------------------------------------------------------------
 !Remove global mean value of dsol:
rhs=dsol(:,1)*clat
avgsol=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
dsol(:,1)=dsol(:,1)-avgsol

!------------------------------------------------------------------
 !Compute gsol = P[omega_0*dsol - dsrc]:
wka=dt2i*dsol-dsrc
call project(wka,gsol)

return
end subroutine simp

!==========================================================================

subroutine helminv(src,sol,der)
! Inverts the operator L = Lap - f^2/c^2 appearing in the equation
! determining the height anomaly on src to give sol and its 
! latitudinal derivative der.  
! Explicitly, L(sol) = src and der = d(sol)/dphi.
! *** All fields are in semi-spectral space ***

implicit none

 !Passed variables:
double precision:: src(ng,nt),sol(ng,nt),der(ng,nt)

 !Local variables:
double precision:: rhs(ng),pm,utdb,vtdb,avgsol
integer:: j,m

!-----------------------------------------------------------
 !Loop over longitudinal wavenumbers and solve the 2x2 
 !block tridiagonal problem:
do m=1,nt
  pm=psig(ksig(m))

  rhs(1)=pm*src(1,m)+0.1d0*src(2,m)
  do j=2,ngm1
    rhs(j)=src(j,m)+0.1d0*(src(j-1,m)+src(j+1,m))
  enddo
  rhs(ng)=0.1d0*src(ngm1,m)+pm*src(ng,m)

  sol(1,m)=-rhs(1)*bhe0(1,m)
  der(1,m)= rhs(1)*ahe0(1,m)

  do j=2,ng
    utdb=      -sol(j-1,m)*amla     -der(j-1,m)*0.2d0
    vtdb=rhs(j)-sol(j-1,m)*cmhe(j,m)-der(j-1,m)*dmhe(j)
    sol(j,m)=utdb*dhe0(j,m)-vtdb*bhe0(j,m)
    der(j,m)=vtdb*ahe0(j,m)-utdb*che0(j,m)
  enddo

  do j=ngm1,1,-1
    sol(j,m)=ehe1(j,m)*sol(j+1,m)+ehe2(j,m)*der(j+1,m)+sol(j,m)
    der(j,m)=fhe1(j,m)*sol(j+1,m)+fhe2(j,m)*der(j+1,m)+der(j,m)
  enddo
enddo

 !Remove global mean value of sol:
rhs=sol(:,1)*clat
avgsol=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
sol(:,1)=sol(:,1)-avgsol

return
end subroutine helminv

!==========================================================================

subroutine pope(src,sol,der)
! Inverts the operator L = Lap - 3/H^2 needed for the non-hydrostatic
! pressure.  Explicitly, L(sol) = src and der = d(sol)/dphi.
! *** All fields are in semi-spectral space ***

implicit none

 !Passed variables:
double precision:: src(ng,nt),sol(ng,nt),der(ng,nt)

 !Local variables:
double precision:: rhs(ng),pm,utdb,vtdb,avgsol
integer:: j,m

!-----------------------------------------------------------
 !Loop over longitudinal wavenumbers and solve the 2x2 
 !block tridiagonal problem:
do m=1,nt
  pm=psig(ksig(m))

  rhs(1)=pm*src(1,m)+0.1d0*src(2,m)
  do j=2,ngm1
    rhs(j)=src(j,m)+0.1d0*(src(j-1,m)+src(j+1,m))
  enddo
  rhs(ng)=0.1d0*src(ngm1,m)+pm*src(ng,m)

  sol(1,m)=-rhs(1)*bpo0(1,m)
  der(1,m)= rhs(1)*apo0(1,m)

  do j=2,ng
    utdb=      -sol(j-1,m)*amla     -der(j-1,m)*0.2d0
    vtdb=rhs(j)-sol(j-1,m)*cmpo(j,m)-der(j-1,m)*dmpo(j)
    sol(j,m)=utdb*dpo0(j,m)-vtdb*bpo0(j,m)
    der(j,m)=vtdb*apo0(j,m)-utdb*cpo0(j,m)
  enddo

  do j=ngm1,1,-1
    sol(j,m)=epo1(j,m)*sol(j+1,m)+epo2(j,m)*der(j+1,m)+sol(j,m)
    der(j,m)=fpo1(j,m)*sol(j+1,m)+fpo2(j,m)*der(j+1,m)+der(j,m)
  enddo
enddo

 !Remove global mean value of sol:
rhs=sol(:,1)*clat
avgsol=(f1112*(rhs(1)+rhs(ng))+sum(rhs(2:ngm1)))*rsumi
sol(:,1)=sol(:,1)-avgsol

return
end subroutine pope

!==========================================================================

subroutine latder(var,der)
! Calculates the first latitudinal derivative of var and stores
! the result in der.  That is, der = d(var)/dphi.
! *** Both var and der are in physical space ***

implicit none

 !Passed variables:
double precision:: var(ng,nt),der(ng,nt)

 !Local variables:
double precision:: src(nt),sol(nt),solend
integer:: i,ic,j

!----------------------------------------------------------------
 !Loop over great circles in latitude (only half the longitudes):
do i=1,ng
  ic=i+ng

   !Source vector:
  src(1)=var(2,i)-var(1,ic)
  do j=2,ngm1
    src(j)=var(j+1,i)-var(j-1,i)
  enddo
  src(ng)  =var(ng,ic)-var(ngm1,i)
  src(ngp1)=var(ngm1,ic)-var(ng,i)
  do j=ngp2,ntm1
    src(j)=var(nt-j,ic)-var(ntp2-j,ic)
  enddo
  src(nt)=var(1,i)-var(2,ic)

   !Get solution by periodic tridiagonal solve:
  sol(1)=src(1)*htd1(1)

  do j=2,nt
    sol(j)=(src(j)-apd1*sol(j-1))*htd1(j)
  enddo

  do j=ntm2,1,-1
    sol(j)=etd1(j)*sol(j+1)+sol(j)
  enddo

  sol(nt)=(etd1(nt)*sol(1)+sol(nt))*denod1
  solend=sol(nt)

  do j=1,ntm1
    sol(j)=ptd1(j)*solend+sol(j)
  enddo

   !Return sol in der:
  do j=1,ng
    der(j,i)=sol(j)
  enddo
  do j=ngp1,nt
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
double precision:: var(ng,nt),sol(ng,nt)

 !Local variables:
double precision:: wka(ng,nt),der(ng,nt)
double precision:: src1(nt),src2(nt),sol1(nt),sol2(nt)
double precision:: sol1end,sol2end
integer:: i,ic,j

!----------------------------------------------------------------------
 !Convert var to physical space as wka to take latitudinal derivatives:
wka=var
call revfft(ng,nt,wka,trig,factors)

!----------------------------------------------------------------
 !Loop over great circles in latitude (only half the longitudes):
do i=1,ng
  ic=i+ng

   !Source vectors:
  src1(1)=wka(2,i)-wka(1,ic)
  src2(1)=wka(2,i)-two*wka(1,i)+wka(1,ic)
  do j=2,ngm1
    src1(j)=wka(j+1,i)-wka(j-1,i)
    src2(j)=wka(j+1,i)-two*wka(j,i)+wka(j-1,i)
  enddo
  src1(ng)  =wka(ng,ic)-wka(ngm1,i)
  src2(ng)  =wka(ng,ic) -two*wka(ng,i)+wka(ngm1,i)
  src1(ngp1)=wka(ngm1,ic)-wka(ng,i)
  src2(ngp1)=wka(ngm1,ic)-two*wka(ng,ic)+wka(ng,i)
  do j=ngp2,ntm1
    src1(j)=wka(nt-j,ic)-wka(ntp2-j,ic)
    src2(j)=wka(nt-j,ic)-two*wka(ntp1-j,ic)+wka(ntp2-j,ic)
  enddo
  src1(nt)=wka(1,i)-wka(2,ic)
  src2(nt)=wka(1,i)-two*wka(1,ic)+wka(2,ic)

   !Get 1st and 2nd derivatives by periodic tridiagonal solve:
  sol1(1)=src1(1)*htd1(1)
  sol2(1)=src2(1)*htd2(1)

  do j=2,nt
    sol1(j)=(src1(j)-apd1*sol1(j-1))*htd1(j)
    sol2(j)=(src2(j)-apd2*sol2(j-1))*htd2(j)
  enddo

  do j=ntm2,1,-1
    sol1(j)=etd1(j)*sol1(j+1)+sol1(j)
    sol2(j)=etd2(j)*sol2(j+1)+sol2(j)
  enddo

  sol1(nt)=(etd1(nt)*sol1(1)+sol1(nt))*denod1
  sol2(nt)=(etd2(nt)*sol2(1)+sol2(nt))*denod2

  sol1end=sol1(nt)
  sol2end=sol2(nt)

  do j=1,ntm1
    sol1(j)=ptd1(j)*sol1end+sol1(j)
    sol2(j)=ptd2(j)*sol2end+sol2(j)
  enddo

   !Combine derivatives (tlat = tan(lat) below):
  do j=1,ng
    der(j,i)=sol2(j)-tlat(j)*sol1(j)
  enddo
  do j=ngp1,nt
    der(ntp1-j,ic)=sol2(j)+tlat(ntp1-j)*sol1(j)
  enddo
enddo

 !Return der to spectral space and add longitudinal part of operator:
call forfft(ng,nt,der,trig,factors)
sol=der-d2lon*var

return
end subroutine laplace

!==========================================================================

subroutine project(var,sol)
! Calculates (1 - (H^2/3)*Lap)(var) and stores the result in sol.
! *** Both var and sol are in semi-spectral space ***

implicit none

 !Passed variables:
double precision:: var(ng,nt),sol(ng,nt)

!----------------------------------------------------------------------
 !Calculate Laplacian:
call laplace(var,sol)

 !Finish calculation:
sol=var-hbsq3*sol

return
end subroutine project

!==========================================================================

subroutine latdamp(var)
! Applies a latitudinal hyperviscous damping to the field var.
! *** var is assumed to be semi-spectral 
! (physical latitude, spectral longitude)

implicit none

 !Passed variable:
double precision:: var(ng,nt)

 !Local variables:
double precision:: tra(ng,nt)
integer:: i,ic,j,k

!---------------------------------------------------
 !Inverse FFT in longitude:
call revfft(ng,nt,var,trig,factors) 

 !Create great circles:
do i=1,ng
  ic=i+ng
  do j=1,ng
    tra(i,j)     =var(j,i)
    tra(i,ntp1-j)=var(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(ng,nt,tra,trig,factors) 

 !Apply spectral damping operator in latitude (defined in adapt):
do k=1,nt
  tra(:,k)=tra(:,k)*hyplat(k)
enddo

 !Inverse FFT in latitude:
call revfft(ng,nt,tra,trig,factors) 

 !Unpack array:
do i=1,ng
  ic=i+ng
  do j=1,ng
    var(j,i) =tra(i,j)
    var(j,ic)=tra(i,ntp1-j)
  enddo
enddo

 !FFT in longitude:
call forfft(ng,nt,var,trig,factors) 

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
double precision:: var(ng,nt)

 !Local variables:
double precision:: tra(ng,nt)
integer:: i,ic,j,k

!---------------------------------------------------
call revfft(ng,nt,var,trig,factors) 

 !Create great circles:
do i=1,ng
  ic=i+ng
  do j=1,ng
    tra(i,j)     =var(j,i)
    tra(i,ntp1-j)=var(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(ng,nt,tra,trig,factors) 

 !Apply latitudinal filter:
do k=1,nt
  tra(:,k)=tra(:,k)*flatlo(k)
enddo

 !Inverse FFT in latitude:
call revfft(ng,nt,tra,trig,factors) 

 !Unpack array:
do i=1,ng
  ic=i+ng
  do j=1,ng
    var(j,i) =tra(i,j)
    var(j,ic)=tra(i,ntp1-j)
  enddo
enddo

 !FFT in longitude:
call forfft(ng,nt,var,trig,factors) 

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
double precision:: var(ng,nt)

 !Local variables:
double precision:: tra(ng,nt)
integer:: i,ic,j,k

!---------------------------------------------------
call revfft(ng,nt,var,trig,factors) 

 !Create great circles:
do i=1,ng
  ic=i+ng
  do j=1,ng
    tra(i,j)     =var(j,i)
    tra(i,ntp1-j)=var(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(ng,nt,tra,trig,factors) 

 !Apply latitudinal filter:
do k=1,nt
  tra(:,k)=tra(:,k)*flathi(k)
enddo

 !Inverse FFT in latitude:
call revfft(ng,nt,tra,trig,factors) 

 !Unpack array:
do i=1,ng
  ic=i+ng
  do j=1,ng
    var(j,i) =tra(i,j)
    var(j,ic)=tra(i,ntp1-j)
  enddo
enddo

 !FFT in longitude:
call forfft(ng,nt,var,trig,factors) 

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
double precision:: var(ng,nt)

 !Local variables:
double precision:: tra(ng,nt)
integer:: i,ic,j

!---------------------------------------------------
 !Create great circles:
do i=1,ng
  ic=i+ng
  do j=1,ng
    tra(i,j)     =var(j,i)
    tra(i,ntp1-j)=var(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(ng,nt,tra,trig,factors) 

 !Truncate spectrally in latitude:
tra(:,kda1:kda2)=zero

 !Inverse FFT in latitude:
call revfft(ng,nt,tra,trig,factors) 

 !Unpack array:
do i=1,ng
  ic=i+ng
  do j=1,ng
    var(j,i) =tra(i,j)
    var(j,ic)=tra(i,ntp1-j)
  enddo
enddo

 !FFT in longitude:
call forfft(ng,nt,var,trig,factors) 

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
double precision:: var(ng,nt),vrms

 !Local variables:
double precision:: wka(ng,nt)
integer:: i

 !-------------------------------------------------
do i=1,nt
  wka(:,i)=clat*var(:,i)**2
enddo

vrms=zero
do i=1,nt
  vrms=vrms+f1112*(wka(1,i)+wka(ng,i))+sum(wka(2:ngm1,i))
enddo
vrms=sqrt(vrms*dsumi)

return
end subroutine getrms

!=======================================================================

double precision function rand(i)
! Returns a random number: i is any integer

implicit none

double precision:: r
integer:: i

call random_number(r)
rand=r

return
end function rand

!==========================================================================

subroutine zeroavg(var)
! Removes the mean of a gridded field var.
! *** Note wka is used as a work array ***

implicit none

 !Passed variables:
double precision:: var(ng,nt),vsum

 !Local variables:
double precision:: wka(ng,nt)
integer:: i

 !-----------------------------------------------------
do i=1,nt
  wka(:,i)=clat*var(:,i)
enddo

 !Compute 4th-order average:
vsum=zero
do i=1,nt
  vsum=vsum+f1112*(wka(1,i)+wka(ng,i))+sum(wka(2:ngm1,i))
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
double precision:: var(ng,nt)

 !Local variables:
double precision:: wka(ng,nt),vsum
integer:: i

 !---------------------------------------------------
do i=1,nt
  wka(:,i)=clat*var(:,i)
enddo

 !Compute 4th-order average:
vsum=zero
do i=1,nt
  vsum=vsum+f1112*(wka(1,i)+wka(ng,i))+sum(wka(2:ngm1,i))
enddo
average=vsum*dsumi

return
end function average

!==========================================================================

subroutine longspec(var,pow)
! Computes the longitudinal spectrum of a field var and creates
! pow, containing the power spectrum for each longitudinal wavenumber
! m from 1 to ng+1

! *** var must be in semi-spectral space

implicit none

 !Passed variables:
double precision:: var(ng,nt),pow(ngp1)

 !Local variables:
integer:: m,mc

!---------------------------------------------------
 !Compute power spectrum:
pow(1)=sum(clat*var(:,1)**2)

do m=2,ng
  mc=ntp2-m
  pow(m)=sum(clat*(var(:,m)**2+var(:,mc)**2))
enddo

pow(ngp1)=sum(clat*var(:,ngp1)**2)

return
end subroutine longspec

!==========================================================================

subroutine ranspec(var,vrms)
! Generates a random field var with variance spectrum 
! proportional to k^5*exp(-2k^2/ksr^2), with zero global
! average and an rms value equal to vrms.  ksr is specified
! in params.dat.

implicit none

 !Passed variables:
double precision:: var(ng,nt),vrms

 !Local variables:
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt),wkd(ng,nt)
double precision:: fac,avrms
integer:: i,j,k,ic

!------------------------------------------------------------------
 !Generate random values between -1 and +1 at all grid points:
do i=1,nt
  do j=1,ng
    wkc(j,i)=two*rand(0)-one
  enddo
enddo

 !Create great circles:
do i=1,ng
  ic=i+ng
  do j=1,ng
    wka(i,j)     =wkc(j,i)
    wka(i,ntp1-j)=wkc(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(ng,nt,wka,trig,factors)

 !Apply latitudinal spectrum:
do k=1,nt
  do i=1,ng
    wka(i,k)=wka(i,k)*glat(k)
    wkb(i,k)=wka(i,k)*plat(k)
  enddo
enddo

 !Inverse FFT in latitude:
call revfft(ng,nt,wka,trig,factors)
call revfft(ng,nt,wkb,trig,factors)

 !Unpack arrays:
do i=1,ng
  ic=i+ng
  do j=1,ng
    wkc(j,i) =wka(i,j)
    wkc(j,ic)=wka(i,ntp1-j)
    wkd(j,i) =wkb(i,j)
    wkd(j,ic)=wkb(i,ntp1-j)
  enddo
enddo

 !FFT in longitude:
call forfft(ng,nt,wkc,trig,factors)
call forfft(ng,nt,wkd,trig,factors)

 !Apply longitudinal spectrum and define var:
var=glon*(plon*wkc+wkd)

 !Return var to physical space:
call revfft(ng,nt,var,trig,factors)

 !Remove global mean:
call zeroavg(var)

 !Find rms value:
call getrms(var,avrms)

 !Normalise:
fac=vrms/avrms
var=fac*var

return
end subroutine ranspec

!==========================================================================

end module spectral
