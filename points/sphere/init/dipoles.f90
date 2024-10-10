program dipoles
!===============================================================
! Initialises dipole pairs of point vortices (+kappa,-kappa),
! with the strength kappa fixed, on a sphere.  

! The strengths of the vortices are chosen so that a dipole
! propagates a mean inter-vortex spacing, L = sqrt(4*pi/N), 
! in unit time.  Here, N is the total number of vortices.
! Below, the dipole spacing a divided by L is specified,
! while parameters.f90 contains N.

! Written on 30 April 2013 by D G Dritschel @ NY

! Output files:
! -------------------------------------------------
! points.dat:    x,y,z positions of all point vortices
! strengths.dat: strengths s of all point vortices
! energy.dat:    interaction energy
!==============================================================

 ! Import parameters and constants:
use parameters
use constants

implicit double precision (a-h,o-z)

 ! Routine-specific parameters:
integer,parameter:: ndip=n/2
 ! ndip: the number of dipoles (n is set in parameters.f90)

 ! Constants:
double precision,parameter:: pi=3.1415926535897932385d0, twopi=two*pi

 ! Routine-specific arrays:
double precision:: sv(n), x(n), y(n), z(n)

 !----------------------------------------------------------------
 ! Specify dipole spacing (between + and - vortices):
write(*,'(a,i4)') ' Note, the number of dipoles = ',ndip
ell=sqrt(4.d0*pi/dble(n))
write(*,'(a,f9.7)') ' and the mean inter-vortex spacing L = ',ell

write(*,*)
write(*,*) ' Enter alpha, the dipole separation divided by L:'
read(*,*) alp

 ! Choose vortex strengths as +/-strv
strv=alp*ell**2
do j=1,ndip
  sv(2*j-1)=strv
  sv(2*j)=-strv
enddo
sum=zero
do i=1,n
  sum=sum+sv(i)**2
enddo
eps=one/sum

 ! Distance from centre of dipole to one of its two vortices:
radv=alp*ell/two
zv=cos(radv)
rv=sin(radv)

 !----------------------------------------------------------------
 ! Next randomly place dipole centres and vortices in each dipole:
write(*,*) ' Energy range to accept initial state?'
read(*,*) emin,emax

write(*,*)
write(*,*) ' Searching for an initial state ...'

tene=emax+one
do while ((tene-emin)*(tene-emax) .gt. zero)

do j=1,ndip
   ! Place dipole centre at (clatc*clonc, clatc*slonc, slatc):
  slatc=two*rand(0)-one
  clatc=sqrt(one-slatc**2)
  rlonc=twopi*rand(0)-pi
  clonc=cos(rlonc)
  slonc=sin(rlonc)

   ! Randomly rotate each dipole:
  phioff=twopi*rand(0)-pi

   ! Place +ve vortex:
  xv=rv*cos(phioff)
  yv=rv*sin(phioff)
  xm=zv*clatc+xv*slatc
  zm=zv*slatc-xv*clatc
  x(2*j-1)=xm*clonc-yv*slonc
  y(2*j-1)=yv*clonc+xm*slonc
  z(2*j-1)=zm

   ! Place -ve vortex:
  xv=-xv
  yv=-yv
  xm=zv*clatc+xv*slatc
  zm=zv*slatc-xv*clatc
  x(2*j)=xm*clonc-yv*slonc
  y(2*j)=yv*clonc+xm*slonc
  z(2*j)=zm
enddo

   ! Slightly adjust vortex positions to ensure zero angular momentum:
angm=one
do while (angm .gt. 1.d-14)
  ax=zero
  ay=zero
  az=zero
  do i=1,n
    ax=ax+sv(i)*x(i)
    ay=ay+sv(i)*y(i)
    az=az+sv(i)*z(i)
  enddo

  do i=1,n
    fac=eps*sv(i)
    xx=x(i)-fac*ax
    yy=y(i)-fac*ay
    zz=z(i)-fac*az
    fac=one/sqrt(xx**2+yy**2+zz**2)
    x(i)=fac*xx
    y(i)=fac*yy
    z(i)=fac*zz
  enddo

  angm=sqrt(ax**2+ay**2+az**2)
enddo

 ! Compute energy:
ene=zero
do j=2,n
  do i=1,j-1
    ene=ene-sv(i)*sv(j)*log(one-x(i)*x(j)-y(i)*y(j)-z(i)*z(j))
  enddo
enddo
tene=two*ene

enddo
write(*,'(a,1p,e14.7)') ' Point vortex interaction energy = ',tene

 ! Write vortex strengths:
open(20,file='strengths.dat',status='unknown')
do i=1,n
  write(20,'(1x,f14.11)') sv(i)
enddo
close(20)

 ! Write vortex positions:
open(11,file='points.dat',status='unknown')
write(11,'(f12.5)') zero
do i=1,n
  write(11,'(3(1x,f12.9))') x(i),y(i),z(i)
enddo
close(11)

 ! Write energy to a file:
open(44,file='energy.dat',status='unknown')
write(44,'(1x,f12.5,1x,1p,e14.7)') zero,tene
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
