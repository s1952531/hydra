program axisym_top

!=============================================================
!  Sets up an axisymmetric hill or dip at a specified point
!=============================================================

! Import constants & parameters:
use constants

! Declarations:
implicit none

double precision:: botTopo(numLatGridPts,numLongGridPts), amp, wid, fac, lon, lat
double precision:: phic, lamc, xc, yc, zc, x, y, z, zz, arg
integer:: i, j

!--------------------------------------------------------------------------
write(*,*) ' We take eta_b/H_1 = A*exp(-s), where s = 4*(1-Z)/(w^2*(1+Z)),'
write(*,*) ' and where Z is the axial coordinate pointing in the direction'
write(*,*) ' of latitude phi_top and longitude lambda_top.'
write(*,*)

write(*,*) ' Enter the amplitude A and width w (degrees):'
read(*,*) amp, wid
! Convert to radians:
fac = pi/180.d0
wid = fac*wid

write(*,*) ' Enter phi_top and lambda_top (degrees):'
read(*,*) phic, lamc
! Convert to radians:
phic = fac*phic
lamc = fac*lamc

! Cartesian centre of the topography:
xc = cos(phic)*cos(lamc)
yc = cos(phic)*sin(lamc)
zc = sin(phic)

! Define scaled topography botTopo:
fac = four/wid**2
do i = 1,numLongGridPts
   lon = dl*dble(i-1)-pi
   do j = 1,numLatGridPts
      lat = (dble(j)-f12)*dl-hpi
      x = cos(phi)*cos(lam)
      y = cos(phi)*sin(lam)
      z = sin(phi)
      zz = x*xc+y*yc+z*zc
      arg = fac*(one-zz)/(one+zz)
      ! Avoid possible underflow:
      if (arg < 40.d0) then
         botTopo(j,i) = amp*exp(-arg)
      else
         botTopo(j,i) = zero
      endif
   enddo
enddo

!Write data:
open(11,file = 'bb_init.r8',form = 'unformatted', &
     access = 'direct',status = 'replace',recl = 2*nbytes)
write(11,rec=1) zero,botTopo
close(11)

end program axisym_top
