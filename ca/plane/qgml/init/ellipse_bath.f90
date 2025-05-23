program ellipse_bath
! --------------------------------------------------------------|
! |   This routine sets up an elliptically-shaped bathymetry    |
!---------------------------------------------------------------|

 !Import constants & parameters:
use constants

implicit none

double precision:: bb(0:ny,0:nxm1)
double precision:: amp,ael,bel,phi,cop,sip
double precision:: xg,yg,xx,yy
integer:: ix,iy

write(*,*) ' The PV due to bathymetry is taken to be'
write(*,*) '   q_b = A*exp(-(X/a)^2-(Y/b)^2)'
write(*,*) ' where X & Y are rotated CCW by an angle phi from the x & y axes,'
write(*,*) ' taken to cross through the domain centre.'
write(*,*)
write(*,*) ' Enter A/pi:'
read(*,*) amp
amp=pi*amp
write(*,*) ' Enter a:'
read(*,*) ael
write(*,*) ' Enter b:'
read(*,*) bel
write(*,*) ' Enter phi:'
read(*,*) phi
phi=phi*pi/180.d0

 !---------------------------------------------------------------------------
 !Generate and write the PV due to the bathymetry:
cop=cos(phi)
sip=sin(phi)
do ix=0,nxm1
  xg=xmin+glx*dble(ix)
  do iy=0,ny
    yg=ymin+gly*dble(iy)
    xx=xg*cop+yg*sip
    yy=yg*cop-xg*sip
    bb(iy,ix)=amp*exp(-(xx/ael)**2-(yy/bel)**2)
  enddo
enddo

open(11,file='bath.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nhbytes)
write(11,rec=1) zero,bb
close(11)

end program ellipse_bath
