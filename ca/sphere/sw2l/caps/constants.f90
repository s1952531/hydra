module constants

 !Module containing all non-modifiable parameters.

use parameters

!-------------------------------------------------------------------------
 !Integer kind parameters used in contours.f90 and congen.f90:
 !  =>  dbleint represents up to 10^18
 !  =>  halfint represents up to 127
integer,parameter:: dbleint=selected_int_kind(16)
integer,parameter:: halfint=selected_int_kind(-1)

!-------------------------------------------------------------------------
 !Other integer parameters associated with resolution, etc:
integer,parameter:: ngm1=numLatGridPts-1, ngp1=numLatGridPts+1, ngp2=numLatGridPts+2
integer,parameter:: ntm2=numLongGridPts-2, ntm1=numLongGridPts-1, ntp1=numLongGridPts+1, ntp2=numLongGridPts+2

integer,parameter:: ngridp=numLatGridPts*numLongGridPts,nbytes=4*(ngridp+1)

integer,parameter:: kda1=ngp1-numLatGridPts/3, kda2=ngp1+numLatGridPts/3
!      Wavenumbers between kda1 and kda2 are removed for de-aliasing;
!      this is done in longitude and in latitude (using great circles)

integer,parameter:: mgf=4, ngf=numLatGridPts*mgf, ntf=numLongGridPts*mgf
!      mgf:      Fine conversion grid ratio
!      ntf,ngf:  number of grid boxes in the x & y directions (fine grid)
integer,parameter:: mgu=16,ntu=mgu*numLongGridPts, ngu=mgu*numLatGridPts
!      mgu:      Ultra-fine conversion grid ratio (used in module congen.f90)
!      ntu,ngu:  number of grid boxes in the x & y directions (ultra-fine grid)

 !Maximum number of contour node points:
integer,parameter:: npm=200*ngridp
 !Maximum number of contours:
integer,parameter:: nm=npm/20+npm/200
 !Maximum number of nodes on any single contour:
integer,parameter:: nprm=npm/10
 !Maximum number of nodes in any contour level:
integer,parameter:: nplm=npm/4

!-------------------------------------------------------------------------
 !Numerical constants:
double precision,parameter:: zero=0.d0,  one=1.d0, two=2.d0,  three=3.d0
double precision,parameter:: four=4.d0, five=5.d0, six=6.d0, twelve=12.d0
double precision,parameter:: f12=one/two,  f13=one/three,  f23=two/three
double precision,parameter:: f14=one/four, f32=three/two,  f43=four/three
double precision,parameter:: f53=five/three, f56=five/six, f74=7.d0/four
double precision,parameter:: f16=one/six, f112=one/twelve, f1112=11.d0/twelve
double precision,parameter:: hpi=pi/two, twoPi=two*pi, fourpi=four*pi
double precision,parameter:: small=1.d-12, small3=small*small*small

!-------------------------------------------------------------------------
 !Grid spacing constants:
double precision,parameter:: hlx=pi, hlxi=one/(pi+small), oms=one-small
 !The "x" range is  -pi  <= x <=  pi   (longitude)
 !The "y" range is -pi/2 <= y <= pi/2  (latitude)
double precision,parameter:: glxu=twoPi/dble(ntu), glyu=pi/dble(ngu)

double precision,parameter:: dl=twoPi/dble(numLongGridPts), dli=dble(numLongGridPts)/(twoPi+small)
double precision,parameter:: dlsq=dl**2, hpidl=(pi+dl)/two

 !Time step related parameters:
double precision,parameter:: dt2=dt*f12, dt4=dt*f14
double precision,parameter:: dt2i=one/dt2, dt4i=one/dt4

!-------------------------------------------------------------------------
 !Derived physical constants: 

 !Planetary rotation rate to ensure 1 unit of time is 1 "day":
double precision,parameter:: omega=two*pi, fpole=two*omega

 !Squared gravity wave speeds and their inverses in each layer:
double precision,parameter:: csq1=gravity*hAvgLayer1, csq1i=one/csq1
double precision,parameter:: csq2=gravity*hAvgLayer2, csq2i=one/csq2

 !Mode -> layer coefficients needed in time stepping and height inversion:
double precision,parameter:: vmdet=vm11*vm22-vm12*vm21
double precision,parameter:: vt11= vm22/vmdet, vt12=-vm12/vmdet
double precision,parameter:: vt21=-vm21/vmdet, vt22= vm11/vmdet

 !Inversed squared modal gravity wave speeds:
double precision,parameter:: ctsq1i=one/ctsq1, ctsq2i=one/ctsq2

 !Constant offset to remove for defining available potential energy:
double precision,parameter:: dape=hAvgLayer1**2+densityRatio*hAvgLayer2*(two*hAvgLayer1+hAvgLayer2)

 !Mean layer depth & mass ratio, and related constants:
double precision,parameter:: hrat=hAvgLayer2/hAvgLayer1,mubar=densityRatio*hrat
double precision,parameter:: mrat1=one/(one+mubar),mrat2=mubar*mrat1

!-------------------------------------------------------------------------
 !For optional stochastic forcing; number of spherical harmonics needed:
integer,parameter:: ntot=maxSpherHarmDegree+1-minSpherHarmDegree+(maxSpherHarmDegree*(maxSpherHarmDegree+1)-minSpherHarmDegree*(minSpherHarmDegree-1))/2

!-------------------------------------------------------------------------
 !PV jump and its inverse:
double precision,parameter:: dq=two*fpole/dble(numContoursTot), dqi=one/dq

!-------------------------------------------------------------------------
 !Logicals:
logical,parameter:: forcing=(rmsPvForcing > zero)
logical,parameter:: thermal=(thermDampRate > zero)

end module constants
