module constants

 !Module containing all non-modifiable parameters.

use parameters

 !Integer kind parameters used in contours.f90 and congen.f90:
 !  =>  dbleint represents up to 10^18
 !  =>  halfint represents up to 127
integer,parameter:: dbleint=selected_int_kind(16)
integer,parameter:: halfint=selected_int_kind(-1)

 !Size of record used in unformatted writes of real*4 data:
integer,parameter:: nbytes=4*(ng*ng+1)

 !Fine grid used normally in contour -> grid conversion: 
integer,parameter:: mgf=4,ngf=mgf*ng
 !mgf:  fine grid/coarse grid ratio (4 is the standard value)

 !Ultra-fine grid used in contouring: 
integer,parameter:: mgu=16,ngu=mgu*ng
 !mgu:  ultra-fine grid/coarse grid ratio (16 is the standard value)

 !Maximum number of contour levels (used in surgery and congen):
integer,parameter:: nlevm=2000
 !nlevm: up to 2*nlevm contour levels are allowed (very generous)

 !Maximum number of contour nodes:
integer,parameter:: npm=400*ng*ng
 !Maximum number of contours:
integer,parameter:: nm=npm/20+npm/200
 !Maximum number of nodes on any single contour:
integer,parameter:: nprm=npm/10
 !Maximum number of nodes in any contour level:
integer,parameter:: nplm=npm/4

 !Generic double precision numerical constants:
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0,three=3.d0
double precision,parameter:: four=4.d0,six=6.d0,f12=one/two,f13=one/three
double precision,parameter:: f14=one/four,f16=one/six,f23=two/three
double precision,parameter:: twopi=two*pi,thrpi=three*pi,pinv=one/pi
double precision,parameter:: small=1.d-12,small3=small*small*small
double precision,parameter:: oms=one-small

 !Grid constants:
double precision,parameter:: domarea=twopi*twopi
double precision,parameter:: gl=twopi/dble(ng),gli=dble(ng)/twopi
double precision,parameter:: garea=gl*gl,dsumi=one/dble(ng*ng)

! Time step related parameters:
double precision,parameter:: dt2=dt*f12,dt4=dt*f14
double precision,parameter:: dt2i=one/dt2,dt4i=one/dt4

! Squared gravity wave speeds in each layer:
double precision,parameter:: csq=cgw**2,hbar=hbar1+hbar2,gravity=csq/hbar
double precision,parameter:: csq1=gravity*hbar1,csq2=gravity*hbar2
double precision,parameter:: dape=hbar1**2+alpha*hbar2*(two*hbar1+hbar2)
double precision,parameter:: ekmf=f12*garea/hbar,epmf=gravity*ekmf

! Mean layer thickness and mass constants:
double precision,parameter:: hrat=hbar2/hbar1,hrati=hbar1/hbar2
double precision,parameter:: mubar=alpha*hrat,ttbar=six/(four+three*mubar)
double precision,parameter:: ahrsq=mubar*hrat,hhrati=f12*hrati
double precision,parameter:: wmass=one/(one+mubar),cio=wmass*dsumi
double precision,parameter:: hbsq1=hbar1**2,hbsq2=hbar2**2,hb1hb2=hbar1*hbar2
double precision,parameter:: cona1=three*ttbar,cona2=cona1*hrati
double precision,parameter:: cona3=six*(two-ttbar),mubar3=three*mubar

end module
