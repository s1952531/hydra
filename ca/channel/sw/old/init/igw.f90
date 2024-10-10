program igw
! Initialises an inertia--gravity wave

use constants

implicit none

double precision:: hh(0:ny,0:nxm1),uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
double precision:: qq(0:ny,0:nxm1),dd(0:ny,0:nxm1),gg(0:ny,0:nxm1)
double precision:: zz(0:ny,0:nxm1)
double precision:: eps,ss,fac,om,x,y,cokx,sikx,coly,sily,hamp
double precision:: k,l
integer:: ix,iy

!-----------------------------------------------------------------
write(*,*) ' We consider an inertia--gravity wave with y velocity'
write(*,*) ' v = eps*cos(k*x-omega*t)*cos(l*y), where k = 2*pi/L_x,'
write(*,*) ' l = pi/L_y, and omega = s*sqrt(f^2 + c^2*(k^2+l^2)).'
write(*,*)
write(*,*) ' Enter eps << 1 and s (+/-1):'
read(*,*) eps,ss

!--------------------------------------------------
 !Define fields:
k=twopi/ellx
l=   pi/elly
om=ss*sqrt(cof**2+csq*(k**2+l**2))
fac=eps/(cof**2+csq*l**2)
do ix=0,nxm1
  x=xmin+glx*dble(ix)
  cokx=cos(k*x)
  sikx=sin(k*x)
  do iy=0,ny
    y=ymin+gly*dble(iy)
    coly=cos(l*y) 
    sily=sin(l*y)
    hamp=-fac*(om*l*sily+cof*k*coly)
    hh(iy,ix)=hamp*sikx
    uu(iy,ix)=-fac*sikx*(om*cof*coly+csq*k*l*sily)
    vv(iy,ix)=eps*cokx*coly
    zz(iy,ix)=cof*hh(iy,ix)
    qq(iy,ix)=cof
    dd(iy,ix)=om*hamp*cokx
    gg(iy,ix)=om**2*hh(iy,ix)
  enddo
enddo

!--------------------------------------------------
 !Write data:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,qq
close(11)

open(11,file='dd_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,dd
close(11)

open(11,file='gg_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,gg
close(11)

open(11,file='uum_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nxbytes)
write(11,rec=1) uu(0,:)
close(11)

open(11,file='uup_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nxbytes)
write(11,rec=1) uu(ny,:)
close(11)

open(11,file='zz_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,zz
close(11)

open(11,file='hh_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,hh
close(11)

open(11,file='uu_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,uu
close(11)

open(11,file='vv_init.r8',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
write(11,rec=1) zero,vv
close(11)

end program igw
