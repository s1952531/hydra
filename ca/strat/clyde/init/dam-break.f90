program dambreak

use constants

! This routine sets up initial buoyancy and (zero) vorticity fields
! in a general domain.

! We consider two uniform buoyancy regions with a transition centred
! at a specific value of x.

! Revised by dgd & hjd on 17 Jan 2019 to allow for a conformal domain.
implicit none

 !Local vorticity and buoyancy arrays:
double precision:: zz(0:ny,0:nx),bb(0:ny,0:nx)

 !Arrays related to the conformal map (if used):
double precision:: xori(0:ny,0:nx),yori(0:ny,0:nx)

 !Other variables:
double precision:: xorimin,xorimax
double precision:: xfrac,xtran
double precision:: wfrac,wtran
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

write(*,*) ' Enter the x location of the density transition zone'
write(*,*) ' as a fraction of the domain width:'
read(*,*) xfrac
xtran=xorimin+xfrac*(xorimax-xorimin)

write(*,*) ' Enter the width of the transition also as a fraction'
write(*,*) ' of the domain width:'
read(*,*) wfrac
wtran=wfrac*(xorimax-xorimin)

 !Set up buoyancy distribution:
do ix=0,nx
  do iy=0,ny
    zz(iy,ix)=zero
    bb(iy,ix)=-f12*erfc((xori(iy,ix)-xtran)/wtran)
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

!==========================================================================
end program
