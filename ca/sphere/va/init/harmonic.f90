program harmonic

! Sets up a PV anomaly field consisting of the superposition of 
! three spherical harmonics, Y_1^1, Y_2^1 & Y_3^3.

! Note, Y_1^1 = -(3/8*pi)^{1/2} r exp(i lambda)
! Note, Y_2^1 = -(15/8*pi)^{1/2} zr exp(i lambda)
!  and  Y_3^3 = -(5/16*pi)^{1/2} r^3 exp(3i lambda)

! where r = (1-z^2)^{1/2} and z = sin(latitude).  Here, lambda
! is longitude.

use constants

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: ay11(ng),ay21(ng),ay33(ng),cof(ng)
double precision:: clona(nt),clonb(nt),c3lon(nt)
double precision:: hh(ng,nt),qq(ng,nt)

!------------------------------------------------------------
 !Define latitude functions:
c11=-sqrt(3.d0/(8.d0*pi))
c21=-sqrt(15.d0/(8.d0*pi))
c33=-sqrt(5.d0/(16.d0*pi))
do j=1,ng
  rlat=dl*(dble(j)-f12)-hpi
  z=sin(rlat)
  r=cos(rlat)
  ay11(j)=c11*r
  ay21(j)=c21*z*r
  ay33(j)=c33*r**3
  cof(j)=fpole*z
enddo

write(*,*) '    We consider a PV anomaly (divided by f_pole) of the form '
write(*,*) ' q'' = A_1*Re(Y_1^1)+A_2*Re(Y_2^1*exp(2i*p_2))+A_3*Re(Y_3^3*exp(3i*p_3))'
write(*,*)
write(*,*) ' Enter A_1, A_2 & A_3:'
read(*,*) a1,a2,a3
write(*,*) ' Enter p_2 & p_3 (degrees):'
read(*,*) p2,p3

p2=p2*pi/180.d0
p3=p3*pi/180.d0

 !Define longitude functions:
do i=1,nt
  rlon=dl*dble(i-1)-pi
  clona(i)=cos(rlon)
  clonb(i)=cos(rlon+p_2)
  c3lon(i)=cos(three*(rlon+p_3))
enddo

 !Define hh (used for writing h, d & h_equil) and PV:
do i=1,nt
  do j=1,ng
    hh(j,i)=zero
    qq(j,i)=cof(j)+fpole*(a1*ay11(j)*clona(i)+a2*ay21(j)*clonb(i)+ &
                        & a3*ay33(j)*c3lon(i))
  enddo
enddo

 !Write equilibrium height field:
open(20,file='hequil.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
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
