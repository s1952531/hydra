program bump
! Initialises a bump in h centred in the middle of the upper boundary,
! together with uniform PV (so zeta = f*h_tilde) and divergence chosen
! to be consistent with the linearised boundary conditions
! d(delta)/dy = -f*d(h_tilde)/dx.

use constants

implicit none

double precision:: hh(0:ny,0:nxm1),gg(0:ny,0:nxm1),dd(0:ny,0:nxm1)
double precision:: uub(0:nxm1)
double precision:: eps,rad,radsqi,rsq,fac,x,y
integer:: ix,iy

!-----------------------------------------------------------------
write(*,*) ' We take h_tilde = eps*exp(-r^2/R^2) + C, where'
write(*,*) ' r^2 = x^2 + (y-y_max)^2, together with uniform PV'
write(*,*) ' and zero divergence. C ensures <h_tilde> = 0.'
write(*,*)
write(*,*) ' Enter eps and R:'
read(*,*) eps,rad
radsqi=one/rad**2

!--------------------------------------------------
 !Define fields:
fac=four*csq*radsqi
do ix=0,nxm1
  x=xmin+glx*dble(ix)
  do iy=0,ny
    y=ymin+gly*dble(iy)
    rsq=radsqi*(x**2+(y-ymax)**2)
    hh(iy,ix)=eps*exp(-rsq)
    gg(iy,ix)=(cof**2-fac*(rsq-one))*hh(iy,ix)
  enddo
enddo

dd(ny,:)=zero
fac=cof*radsqi*gly
do ix=0,nxm1
  x=xmin+glx*dble(ix)
  do iy=ny-1,0,-1
    dd(iy,ix)=dd(iy+1,ix)-fac*x*(hh(iy,ix)+hh(iy+1,ix))
  enddo
enddo

 !Remove averages:
fac=dsumi*(f12*sum(hh(0,:)+hh(ny,:))+sum(hh(1:nym1,:)))
hh=hh-fac

fac=dsumi*(f12*sum(gg(0,:)+gg(ny,:))+sum(gg(1:nym1,:)))
gg=gg-fac

fac=dsumi*(f12*sum(dd(0,:)+dd(ny,:))+sum(dd(1:nym1,:)))
dd=dd-fac

!--------------------------------------------------
 !Write data:
open(11,file='hh_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,hh
close(11)

open(11,file='dd_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,dd
close(11)

open(11,file='gg_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,gg
close(11)

hh=cof
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,hh
close(11)

uub=zero
open(11,file='uum_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nxbytes)
write(11,rec=1) uub
close(11)

open(11,file='uup_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nxbytes)
write(11,rec=1) uub
close(11)

end program bump
