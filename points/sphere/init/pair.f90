program pair
!===============================================================
! Initialises point vortices of strengths kappa_1 and kappa_2, 
! on a sphere. Here N=2 is the total number of vortices.

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

 ! Constants:
double precision,parameter:: pi=3.1415926535897932385d0, twopi=two*pi

 ! Routine-specific array:
double precision:: sv(n)

 !----------------------------------------------------------------
 ! Read in the strengths of the two vortices:
write(*,*) ' Strength of vortex 1 & vortex 2:'
read(*,*) sv(1),sv(2)

 !----------------------------------------------------------------
 ! Next place vortex centres:
write(*,*) ' Initial lat/long of vortex 1:'
read(*,*) rlat1,rlon1

write(*,*) ' Initial lat/long of vortex 2:'
read(*,*) rlat2,rlon2

rlatc=rlat1
rlonc=rlon1
clatc=cos(rlatc)
z(1)=sin(rlatc)
x(1)=clatc*cos(rlonc)
y(1)=clatc*sin(rlonc)
rlatc=rlat2
rlonc=rlon2
clatc=cos(rlatc)
z(2)=sin(rlatc)
x(2)=clatc*cos(rlonc)
y(2)=clatc*sin(rlonc)

 ! Compute energy:
ene=zero
do j=2,n
  do i=1,j-1
    ene=ene-sv(i)*sv(j)*log(one-x(i)*x(j)-y(i)*y(j)-z(i)*z(j))
  enddo
enddo
tene=two*ene

write(*,*)
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

 !End main program
end program
!=======================================================================
