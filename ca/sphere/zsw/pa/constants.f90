module constants

 !Module containing all non-modifiable parameters.

use parameters

integer, parameter:: n1=n-1, n2=n-2
integer, parameter:: nbytes=8+8*n !for saving data

 !Commonly used numerical constants:
double precision, parameter:: zero=0.d0, one=1.d0, two=2.d0, four=4.d0
double precision, parameter:: f12=one/two, f14=one/four
double precision, parameter:: hpi=pi/two

 !Physical constants:
double precision, parameter:: cgw=rrad*fpole, csq=cgw**2, omega=fpole/two

 !Time step (fixed); based on c*dt/dphi = cfl0 initially:
double precision, parameter:: dt=cfl0/(four*rrad*dble(n))

 !For 4th-order Runge-Kutta integration:
double precision, parameter:: dt2=dt/2.d0, dt3=dt/3.d0, dt6=dt/6.d0

end module constants
