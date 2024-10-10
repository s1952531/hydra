program quadtopo
! Initialises a quadrupolar topographic pattern h. 
! Note that the pattern is a combination of 
! spherical harmonics of degree l=2. Hence the Laplacian
! of h is just -l*(l+1)*h. This is needed in evolution.f90.

use constants

implicit double precision(a-h,o-z)

double precision:: clat(ng),slat(ng)
double precision:: clon(nt),slon(nt)
double precision:: hh(ng,nt),hb(ng,nt),qq(ng,nt)

!--------------------------------------------------------------
 !Read in angles:
write(*,*) ' Enter the angles alpha & theta (degrees):'
read(*,*) alp,the
alp=alp*pi/180.d0
the=the*pi/180.d0
ca=cos(alp)
sa=sin(alp)
ct=cos(the)
st=sin(the)

!------------------------------------------------------------
 !Define latitude and longitude arrays:
do j=1,ng
  rlat=dl*(dble(j)-f12)-hpi
  clat(j)=cos(rlat)
  slat(j)=sin(rlat)
enddo

do i=1,nt
  rlon=dl*dble(i-1)-pi
  clon(i)=cos(rlon)
  slon(i)=sin(rlon)
enddo

!----------------------------------------------------------
do i=1,nt
  do j=1,ng
    x=clat(j)*clon(i)
    y=clat(j)*slon(i)
    z=slat(j)
    hb(j,i)=-six*((x**2-(y*ct-z*st)**2)*ca+2.d0*x*(y*ct-z*st)*sa)
    hh(j,i)=zero
    qq(j,i)=fpole*slat(j)
  enddo
enddo

open(20,file='topogr.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hb
close(20)

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
