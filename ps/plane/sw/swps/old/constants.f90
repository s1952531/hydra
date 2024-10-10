module constants

 !Module containing all non-modifiable parameters.

use parameters

 !Size of record used in unformatted writes of real*4 data:
integer,parameter:: nbytes=4*(ng*ng+1)

 !Generic double precision numerical constants:
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0,four=4.d0
double precision,parameter:: f12=one/two,f14=one/four
double precision,parameter:: twopi=two*pi

 !Grid constants:
double precision,parameter:: domarea=twopi*twopi,dsumi=one/dble(ng*ng)
double precision,parameter:: gl=twopi/dble(ng),garea=gl*gl

! Time step related parameters:
double precision,parameter:: dt2=dt*f12,dt4=dt*f14
double precision,parameter:: hdti=one/dt2,hdtisq=hdti**2
double precision,parameter:: cfac0=one/(one+(dt2*cof)**2),cfac1=dt2*cof*cfac0

! Squared gravity wave speed:
double precision,parameter:: csq=cgw**2

end module
