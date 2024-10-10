program dipoles
!===============================================================
! Forms the energy PDF based on random initialisations of 
! point vortex dipoles on the sphere.  

! Uses E_rms/10 as the energy interval for forming the PDF.

! The strengths of the vortices are chosen so that a dipole
! propagates a mean inter-vortex spacing, L = sqrt(4*pi/N), 
! in unit time.  Here, N is the total number of vortices.
! Below, the dipole spacing a divided by L is specified,
! while parameters.f90 contains N.

! Written on 2 May 2013 by D G Dritschel @ NY

! Output file:
! -------------------------------------------------
! energy-pdf.dat:    The PDF of energy
!==============================================================

 ! Import parameters and constants:
use parameters
use constants

implicit double precision (a-h,o-z)

 ! Routine-specific parameters:
integer,parameter:: nbm=2000
 ! nbm: the maximum number of bins to form the PDF

integer,parameter:: ndip=n/2
 ! ndip: the number of dipoles (n is set in parameters.f90)

 ! Constants:
double precision,parameter:: pi=3.1415926535897932385d0, twopi=two*pi

 ! Routine-specific arrays:
double precision:: sv(n), x(n), y(n), z(n)
double precision:: pe(nbm)

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

 ! Distance from centre of dipole to one of its two vortices:
radv=alp*ell/two
zv=cos(radv)
rv=sin(radv)

write(*,*) ' Number of samples to form energy PDF?'
read(*,*) nsamp

open(10,file='energy-pdf.dat',status='unknown')
emin=1.d14
emax=-emin
se1=zero
se2=zero

 !----------------------------------------------------------------
 ! Generate vortex positions to ensure zero angular momentum:
do isamp=1,nsamp

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
sum=zero
do i=1,n
  sum=sum+sv(i)**2
enddo
eps=one/sum

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
ene=two*ene
emin=min(emin,ene)
emax=max(emax,ene)
se1=se1+ene
se2=se2+ene**2

 ! Write energy to a file:
write(10,'(1x,1p,e14.7)') ene

enddo
 ! Above ends loop over samples

close(10)

 ! Get mean and rms energy:
ebar=se1/dble(nsamp)
erms=sqrt(se2/dble(nsamp)-ebar**2)

 ! Use erms/10 as the energy interval for forming the PDF:
de=erms/10.d0
 ! Number of bins based on range of energies found:
nbe=nint((emax-emin)/de)+1

write(*,*)
write(*,'(a,f12.8)') ' Note, the PDF of E uses equally-spaced bins in E,'// &
                     ' with dE = E_rms/10 = ',de

write(*,*)
write(*,'(2(a,f14.8))') ' min(E) = ',emin,' & max(E) = ',emax

 ! Initialise PDF:
do ibe=1,nbe
  pe(ibe)=zero
enddo

 ! Re-open file and read energy to form PDF:
open(10,file='energy-pdf.dat',status='old')
rewind 10
do isamp=1,nsamp
  read(10,*) ene
  ibe=1+nint((ene-emin)/de)
  pe(ibe)=pe(ibe)+one
enddo
close(10)

 ! Normalise PDF so that sum{p(E)*dE} = 1:
sum=half*(pe(1)+pe(nbe))
do ibe=2,nbe-1
  sum=sum+pe(ibe)
enddo
fac=one/(de*sum)
do ibe=1,nbe
  pe(ibe)=pe(ibe)*fac
enddo

 ! Overwrite file with the PDF and its integral:
open(10,file='energy-pdf.dat',status='unknown')
sum=zero
hde=half*de
write(10,'(1x,f14.8,1x,e14.7,1x,f10.8)') emin,pe(1),sum
do ibe=2,nbe
  sum=sum+hde*(pe(ibe-1)+pe(ibe))
  e=emin+de*dble(ibe-1)
  write(10,'(1x,f14.9,1x,e14.7,1x,f10.8)') e,pe(ibe),sum
enddo

write(*,*) ' See energy-pdf.dat for E vs p(E) & int{p(E)dE}.'

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
