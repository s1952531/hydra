program powhex
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Initialises variable strength point vortices in local hexagons on 
! a sphere.  The number density of vortices of strength s is taken 
! to be proportional to 1/s.

! Here, one chooses (1) the number of hexagons, (2) the max to min
! strength ratio, and (3) the radius of the hexagon containing the 
! strongest vortices.  The radius r of other hexagons is proportional 
! to the strength s of vortices in that hexagon.  Without loss of 
! generality, we take s = r, effectively giving a unit propagation 
! speed of all dipoles that spread away from each hexagon.

! Written 27 Oct 2012 by D G Dritschel @ St Andrews

! Output files:
! -------------
!   points.dat:    Contains x,y,z positions of all vortices
!   strengths.dat: Contains the strengths s of all vortices
!   energy.dat:    Contains the energy of the configuration
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 ! Import parameters and constants:
use parameters
use constants
use variables

implicit double precision (a-h,o-z)

integer,parameter:: m=n/6
 ! m: number of hexagonal clusters

 ! Constants:
double precision,parameter:: pi=3.1415926535897932385d0,twopi=two*pi

 ! Routine-specific arrays:
double precision:: sc(m),radc(m),phiv(6)

!--------------------------------------------------------------------
dphi=pi/three
do k=1,6
  phiv(k)=dphi*dble(k-1)
enddo

write(*,'(a,i5,a)') ' There are ',m,' hexagonal clusters initially.'

write(*,*) ' Max to min vortex strength ratio?'
read(*,*) rho

davg=sqrt(4.d0*pi/dble(n))
write(*,'(a,f11.9)') & 
  &' Note, the average vortex separation d_avg = sqrt(4*pi/n) = ',davg

write(*,*) ' Enter the radius of the largest hexagon r_max / d_avg:'
read(*,*) rmax
rmax=rmax*davg

 ! Take a mean dipole propagation speed of O(1):
smax=rmax
 ! smax is the maximum vortex strength
sc(1)=smax
radc(1)=rmax
alpha=rho**(-one/dble(m-1))
do j=2,m
  sc(j)=sc(j-1)*alpha
  radc(j)=radc(j-1)*alpha
enddo

 !----------------------------------------------------------------
 ! Next randomly place hexagon centres:

write(*,*) ' Enter a random seed (integer):'
read(*,*) iseed
 ! Initialize random # generator:
do i=1,iseed
  uni=rand(0)
enddo

write(*,*) ' Range of energy (E_min,E_max) to accept state?'
read(*,*) emin,emax

ene=emax+one
 ! Search for a configuration with this energy:
do while ((ene-emin)*(emax-ene) .lt. zero)

   ! Place clusters of vortices within each hexagon:
  do j=1,m
    zc=two*rand(0)-one
    rc=sqrt(one-zc**2)
    rlonc=twopi*rand(0)-pi
    clonc=cos(rlonc)
    slonc=sin(rlonc)
    xc=rc*clonc
    yc=rc*slonc
    zv=cos(radc(j))
    rv=sin(radc(j))

    ioff=(j-1)*6
     ! Randomly rotate each cluster:
    phioff=twopi*rand(0)-pi
    ss=one
    do k=1,6
      phi=phioff+phiv(k)
      xv=rv*cos(phi)
      yv=rv*sin(phi)
      xm=zv*rc+xv*zc
      zm=zv*zc-xv*rc
      x(ioff+k)=xm*clonc-yv*slonc
      y(ioff+k)=yv*clonc+xm*slonc
      z(ioff+k)=zm
      s(ioff+k)=sc(j)*ss
      ss=-ss
    enddo
  enddo

   ! Compute energy:
  ene=zero
  do j=2,n
    do i=1,j-1
      ene=ene-s(i)*s(j)*log(one-x(i)*x(j)-y(i)*y(j)-z(i)*z(j))
    enddo
  enddo
  ene=two*ene

  write(*,'(a,1p,e14.7)') ' E = ',ene

enddo

 ! Write vortex strengths:
open(20,file='strengths.dat',status='unknown')
do i=1,n
  write(20,'(1x,f14.11)') s(i)
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
write(44,'(1x,f12.5,1x,1p,e14.7)') zero,ene
close(44)

write(*,'(a,f12.5,a,f14.8)') ' t = ',zero,' Energy = ',ene


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
