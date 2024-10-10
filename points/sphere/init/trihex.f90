program trihex
!===============================================================
! Initialises triangular clusters of point vortices in hexagons
! on a sphere.  

! The hexagons are prevented from overlapping in order to study 
! the inverse energy cascade.

! Adapted from ihex.F on 25 July 2012 by D G Dritschel @ NY

! Output files:
! -------------------------------------------------
! points.dat:    x,y,z positions of all point vortices
! strengths.dat: strengths s of all point vortices
! energy.dat:    interaction energy
!==============================================================

 ! Import parameters and constants:
use parameters
use constants
use variables

implicit double precision (a-h,o-z)

 ! Routine-specific parameters:
integer,parameter:: nclum=2048, nvorm=1260
 ! nclum: the maximum number of hexagonal clusters
 ! nvorm: the maximum number of point vortices per cluster

 ! Constants:
double precision,parameter:: pi=3.1415926535897932385d0, twopi=two*pi

 ! Routine-specific arrays:
double precision:: sv(n)
double precision:: radv(nvorm), phiv(nvorm), strv(nvorm)
double precision:: sc(nclum), xc(nclum), yc(nclum), zc(nclum)
double precision:: phioff(nclum), phit(0:6)

 !--------------------------------------------------------------------
write(*,*) ' Number of hexagons, m?'
read(*,*) nclu

write(*,*) ' Each hexagon consists of 6 triangular clusters of'
write(*,*) ' alternating signed point vortices.  Each triangle'
write(*,*) ' consists of n(n+1)/2 vortices.  Enter n:'
read(*,*) ndiv

 ! Number of vortices per triangle:
nvpt=ndiv*(ndiv+1)/2
 ! Number of vortices per hexagon:
nvpc=6*nvpt
 ! Number of vortices altogether:
nvor=nclu*nvpc

 ! Safety check:
if (nvor .ne. n) then
  write(*,*)
  write(*,*) ' *** Stopping ***'
  write(*,*)
  write(*,'(a,i6,a)') '  Edit parameters.f90 and make sure   n = ',nvor, &
                    & '   then remake.'
  write(*,*)
  stop
endif

 ! Choose vortex strengths as +/-sqrt(nclu)/n:
str=sqrt(dble(nclu))/dble(nvor)

write(*,*) ' Area fraction covered by all triangles?'
read(*,*) afrac

write(*,*) ' Separation constant c (3.6 recommended)?'
read(*,*) cdis

 ! Area of each triangle:
area=4.d0*pi*afrac/dble(6*nclu)

sq3=sqrt(3.d0)
r=two*sqrt(area/(three*sq3))
radt=sq3*r
wid=radt/dble(2*ndiv)
hgt=r/dble(2*ndiv)

 ! Fill one triangle:
k=0
do j=1,ndiv
  xoff=wid*dble(j)
  yv=hgt*dble(3*j-2)
  do i=1,ndiv+1-j
    k=k+1
    xv=xoff+wid*dble(2*i-2)
    radv(k)=sqrt(xv**2+yv**2)
    phiv(k)=atan2(yv,xv)
    strv(k)=str
  enddo
enddo

 ! Fill remaining 5 triangles by rotation:
phit(0)=zero
sstr=str
do l=1,5
  sstr=-sstr
  phit(l)=pi*dble(l)/three
  koff=nvpt*l
  do k=1,nvpt
    radv(koff+k)=radv(k)
    phiv(koff+k)=phiv(k)+phit(l)
    strv(koff+k)=sstr
  enddo
enddo
phit(6)=twopi

 !----------------------------------------------------------------
 ! Now one hexagon is specified.  Next randomly place hexagon
 ! centres, with sufficient spacing:
write(*,*) ' Enter a random seed (integer):'
read(*,*) iseed
 ! Initialize random number generator:
do i=1,iseed
  uni=rand(0)
enddo

 ! Place first hexagon:
zc(1)=two*rand(0)-one
rc=sqrt(one-zc(1)**2)
rlonc=twopi*rand(0)-pi
xc(1)=rc*cos(rlonc)
yc(1)=rc*sin(rlonc)
write(*,'(a,i3)') ' Placed hexagon ',1

 ! Place remaining ones, ensuring adequate spacing:
do j=2,nclu
1 zc(j)=two*rand(0)-one
  rc=sqrt(one-zc(j)**2)
  rlonc=twopi*rand(0)-pi
  xc(j)=rc*cos(rlonc)
  yc(j)=rc*sin(rlonc)
  dotm=one-cdis/dble(j)
  do i=1,j-1
    if (xc(i)*xc(j)+yc(i)*yc(j)+zc(i)*zc(j) .gt. dotm) goto 1
  enddo
  write(*,'(a,i3)') ' Placed hexagon ',j
enddo

 ! Check minimum separation:
dmax=zero
do j=2,nclu
  do i=1,j-1
    dmax=max(dmax,xc(i)*xc(j)+yc(i)*yc(j)+zc(i)*zc(j))
  enddo
enddo
rmin=sqrt(two*(one-dmax))
write(*,*)
write(*,'(a,f12.8)') ' min|X_i-X_j|/R = ',rmin/radt

 ! Now place clusters of vortices within each hexagon:
do j=1,nclu
  slatc=zc(j)
  clatc=sqrt(one-slatc**2)
  clonc=xc(j)/clatc
  slonc=yc(j)/clatc

  ioff=(j-1)*nvpc
   ! Randomly rotate each cluster:
  phioff(j)=twopi*rand(0)-pi
  do k=1,nvpc
    zv=cos(radv(k))
    rv=sin(radv(k))
    xv=rv*cos(phiv(k)+phioff(j))
    yv=rv*sin(phiv(k)+phioff(j))
    xm=zv*clatc+xv*slatc
    zm=zv*slatc-xv*clatc
    x(ioff+k)=xm*clonc-yv*slonc
    y(ioff+k)=yv*clonc+xm*slonc
    z(ioff+k)=zm
    sv(ioff+k)=strv(k)
  enddo
enddo

 ! Compute energy:
ene=zero
do j=2,nvor
  do i=1,j-1
    ene=ene-sv(i)*sv(j)*log(one-x(i)*x(j)-y(i)*y(j)-z(i)*z(j))
  enddo
enddo
ene=two*ene

write(*,'(a,1p,e14.7)') ' Point vortex interaction energy = ',ene

 ! Write vortex strengths:
open(20,file='strengths.dat',status='unknown')
do i=1,nvor
  write(20,'(1x,f14.11)') sv(i)
enddo
close(20)

 ! Write vortex positions:
open(11,file='points.dat',status='unknown')
write(11,'(f12.5)') zero
do i=1,nvor
  write(11,'(3(1x,f12.9))') x(i),y(i),z(i)
enddo
close(11)

 ! Write energy to a file:
open(44,file='energy.dat',status='unknown')
write(44,'(1x,f12.5,1x,1p,e14.7)') zero,ene
close(44)


 ! Internal subroutine definitions:
contains 

!=======================================================================

function rand(i)
! returns a random number: i is any integer

implicit none
double precision:: rand,r
integer:: i

call random_number(r)
rand=r

return 
end function


 !End main program
end program
!=======================================================================
