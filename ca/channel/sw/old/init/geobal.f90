program geobal
! Initialises a height field in near geostrophic balance.

use spectral

implicit none

double precision:: hh(0:ny,0:nxm1),uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
double precision:: qq(0:ny,0:nxm1),dd(0:ny,0:nxm1),gg(0:ny,0:nxm1)
double precision:: zz(0:ny,0:nxm1),aa(0:ny,0:nxm1),bb(0:ny,0:nxm1)
double precision:: d10(0:nxm1),d21(0:nxm1),d32(0:nxm1)
double precision:: a,b,c,eps,fac,x,bcox,esix,y,yhat
double precision:: k,l,phase
integer:: ix,iy

!-----------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!-----------------------------------------------------------------
write(*,*) ' We consider a dimensionless height anomaly of the form'
write(*,*) ' h = h_0 + eps*h_1, where'
write(*,*) ' h_0 = A*cos(ly)*log(1 + B*y*cos(kx)) + C, and'
write(*,*) ' h_1 = (4*eta^3 - 3*eta)*sin(kx) with eta = y/L_y.'
write(*,*) ' Here k = 2*pi/L_x and l = pi/L_y, and C is found by'
write(*,*) ' requiring <h> = 0.'
write(*,*)
write(*,*) ' Enter A and B:'
read(*,*) a,b
if (abs(b)*hly .ge. one) then
  write(*,*) ' |B|*L_y/2 is too large!  *** Stopping!'
  stop
endif
write(*,*) ' When eps = 0, the ageostrophic vorticity = 0.  Enter eps:'
read(*,*) eps
write(*,*) ' Enter phase offset in x:'
read(*,*) phase

 !Define h_0 & eps*h_1 separately:
k=twopi/ellx
l=   pi/elly
do ix=0,nxm1
  x=xmin+glx*dble(ix)
  bcox=  b*cos(k*x+phase)
  esix=eps*sin(k*x+phase)
  do iy=0,ny
    y=ymin+gly*dble(iy)
    aa(iy,ix)=a*cos(l*y)*log(one+y*bcox)
    yhat=y/elly
    bb(iy,ix)=yhat*(four*yhat**2-three)*esix
  enddo
enddo

 !Remove domain mean from h_0 in aa:
call restore(aa,zero)

 !Complete height anomaly:
hh=aa+bb

 !Calculate dh_0/dy:
call yderiv(aa,uu)
 !Use cubic interpolation to correct boundary values:
d10=aa(1,:)-aa(0,:)
d21=aa(2,:)-aa(1,:)
d32=aa(3,:)-aa(2,:)
uu(0,:)=f13*hglyi*(11.d0*d10-7.d0*d21+2.d0*d32)
d10=aa(ny,:)-aa(ny-1,:)
d21=aa(ny-1,:)-aa(ny-2,:)
d32=aa(ny-2,:)-aa(ny-3,:)
uu(ny,:)=f13*hglyi*(11.d0*d10-7.d0*d21+2.d0*d32)
 !Calculate dh_0/dx:
call forfft(nyp1,nx,aa,xtrig,xfactors)
call xderiv(aa,vv)
call revfft(nyp1,nx,vv,xtrig,xfactors)

write(*,*) ' The flow is made divergent by multiplying v by'
write(*,*) ' a factor F not equal to 1.  Enter F:'
read(*,*) fac
uu=   -geo*uu
vv=fac*geo*vv

 !Calculate the relative vorticity:
aa=vv
call forfft(nyp1,nx,aa,xtrig,xfactors)
call xderiv(aa,zz)
call revfft(nyp1,nx,zz,xtrig,xfactors)
call yderiv(uu,bb)
 !Correct boundaries:
d10=uu(1,:)-uu(0,:)
d21=uu(2,:)-uu(1,:)
d32=uu(3,:)-uu(2,:)
bb(0,:)=f13*hglyi*(11.d0*d10-7.d0*d21+2.d0*d32)
d10=uu(ny,:)-uu(ny-1,:)
d21=uu(ny-1,:)-uu(ny-2,:)
d32=uu(ny-2,:)-uu(ny-3,:)
bb(ny,:)=f13*hglyi*(11.d0*d10-7.d0*d21+2.d0*d32)
zz=zz-bb

 !Calculate the potential vorticity:
qq=(zz+cof)/(one+hh)

 !Calculate the divergence:
aa=uu
call forfft(nyp1,nx,aa,xtrig,xfactors)
call xderiv(aa,dd)
call revfft(nyp1,nx,dd,xtrig,xfactors)
call yderiv(vv,bb)
 !Correct boundaries:
d10=vv(1,:)
d21=vv(2,:)-vv(1,:)
d32=vv(3,:)-vv(2,:)
bb(0,:)=f13*hglyi*(11.d0*d10-7.d0*d21+2.d0*d32)
d10=-vv(ny-1,:)
d21=vv(ny-1,:)-vv(ny-2,:)
d32=vv(ny-2,:)-vv(ny-3,:)
bb(ny,:)=f13*hglyi*(11.d0*d10-7.d0*d21+2.d0*d32)
dd=dd+bb

 !Calculate the acceleration divergence:
aa=hh
call forfft(nyp1,nx,aa,xtrig,xfactors)
call xderiv(aa,bb)
call xderiv(bb,gg)
call revfft(nyp1,nx,gg,xtrig,xfactors)
call yderiv(hh,aa)
 !Correct boundaries:
aa(0,:)=-uu(0,:)/geo
aa(ny,:)=-uu(ny,:)/geo
call yderiv(aa,bb)
 !Correct boundaries:
d10=aa(1,:)-aa(0,:)
d21=aa(2,:)-aa(1,:)
d32=aa(3,:)-aa(2,:)
bb(0,:)=f13*hglyi*(11.d0*d10-7.d0*d21+2.d0*d32)
d10=aa(ny,:)-aa(ny-1,:)
d21=aa(ny-1,:)-aa(ny-2,:)
d32=aa(ny-2,:)-aa(ny-3,:)
bb(ny,:)=f13*hglyi*(11.d0*d10-7.d0*d21+2.d0*d32)
gg=cof*zz-csq*(gg+bb)

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

end program geobal
