program vort

use constants

! This routine sets up initial buoyancy and (zero) vorticity fields
! in a general domain.

! We consider a Gaussian buoyancy anomaly centred at a specified location.

! Revised by dgd & hjd on 17 Jan 2019 to allow for a conformal domain.

implicit none

 !Local vorticity and buoyancy arrays:
double precision:: zz(0:ny,0:nx),bb(0:ny,0:nx)

 !Arrays related to the conformal map (if used):
double precision:: xori(0:ny,0:nx),yori(0:ny,0:nx)

 !Other variables:
double precision:: xorimin,xorimax
double precision:: yorimin,yorimax
double precision:: xfrac,yfrac
double precision:: x0,y0
double precision:: wfrac,wtran
double precision:: bamp,zamp
double precision:: gfac,wfac
integer:: ix,iy

!-----------------------------------------------------------------------
! Read coords.r8 to get the original coordinates as a function of the
! conformal coordinates:
open(11,file='coords.r8',form='unformatted', &
     access='direct',status='old',recl=2*nbytes)
read(11,rec=1) xori
read(11,rec=2) yori
close(11)

xorimin=minval(xori)
xorimax=maxval(xori)

yorimin=minval(yori)
yorimax=maxval(yori)

write(*,*) ' Enter the x location of the vortex centre as a fraction'
write(*,*) ' of the domain width:'
read(*,*) xfrac
x0=xorimin+xfrac*(xorimax-xorimin)

write(*,*) ' Enter the y location of the vortex centre as a fraction'
write(*,*) ' of the domain height:'
read(*,*) yfrac
y0=yorimin+yfrac*(yorimax-yorimin)

write(*,*) ' Enter the width of the vortex as a fraction'
write(*,*) ' of the domain width:'
read(*,*) wfrac
wtran=wfrac*(xorimax-xorimin)
wfac=one/wtran**2

write(*,*) ' Enter the maximum buoyancy anomaly in the vortex:'
read(*,*) bamp

write(*,*) ' Enter the central vorticity value in the vortex:'
read(*,*) zamp

 !Set up buoyancy and vorticity distributions:
do ix=0,nx
  do iy=0,ny
    gfac=exp(-wfac*((xori(iy,ix)-x0)**2+(yori(iy,ix)-y0)**2))
    bb(iy,ix)=bamp*gfac
    zz(iy,ix)=zamp*gfac
  enddo
enddo

 !Write vorticity distribution to file:
open(20,file='zz_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,zz
close(20)

 !Write buoyancy distribution to file:
open(20,file='bb_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,bb
close(20)

end program
