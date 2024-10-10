program r4toc2
!  ---------------------------------------------------------------------
!  |   Creates character maps for producing images of the vortices     |
!  |   either at a fixed time or over a selected range of times.       |
!  |                                                                   |
!  |   Images are created in orthographic projection.                  |
!  |                                                                   |
!  |   The input data file (points.r4) is assumed to be unfortmatted   |
!  |   and direct access.  The ".r4" extension indicates a single-     |
!  |   precision real.                                                 |
!  |                                                                   |
!  |   Output files (ending in .c2) are also unformatted and direct    |
!  |   access.  Here, they contain a pair of characters (hence "c2")   |
!  |   per pixel to enable practically continuous true colour;         |
!  |   images can be converted to various formats using the script     |
!  |   c2image, and they can be made into movies using c2movie.        |
!  |                                                                   |
!  |         Completed 16 July 2013 by D G Dritschel @ NYC             |
!  ---------------------------------------------------------------------

use constants

implicit none

integer,parameter:: nbpc=1
! npbc: number of bytes used to represent a single character;
!       normally 1 but can be 4 (as when using Intel Fortran)

real,parameter:: cfrac=4096./(256.*256)
 !This parameter is used to access the fraction of the colourmap 
 !extending from cfrac to 1-cfrac (to avoid oversaturation)

integer,parameter:: np=20, npt=n*np
 ! np:  the number of nodes used to image each vortex as a small patch
 ! npt: the total number of contour nodes

real,parameter:: hpi=0.5*pi

real:: sv(n),xv(n),yv(n),zv(n),xd(npt),yd(npt),zd(npt)
real:: xc(0:np-1),yc(0:np-1)
integer:: next(npt)

real:: dth,th,svmax,qqmin,qqmax,cfac,crange,coff,rchar,rmax,frad,t
integer:: i,j,k,sopt,npix,kbeg,kend,kint,period,nbread,nbprec
character(len=4):: suffix

!---------------------------------------------------------
 !Quantities required for converting data into characters:
crange=float(256**2)
coff=cfrac*crange

qqmin=-1.0
qqmax= 1.0
cfac=(crange-2.*coff)/(qqmax-qqmin)

!--------------------------------------------------------
 !For defining circular patches around each point vortex:
dth=twopi/np
do i=0,np-1
  th=dth*i
  xc(i)=cos(th)
  yc(i)=sin(th)
enddo

!-------------------------------------------------------------
 !For use in p2c, define next(i), which gives the contour node
 !following node i:
do i=1,npt-1
  next(i)=i+1
enddo
do j=1,n
  next(j*np)=(j-1)*np+1
enddo

!---------------------------------------------------
 !Read vortex strengths:
open(20,file='strengths.dat',status='old')
do k=1,n
  read(20,*) sv(k)
enddo
close(20)

svmax=0.
do k=1,n
  svmax=max(svmax,abs(sv(k)))
enddo

rchar=2./sqrt(float(n))
write(*,'(a,f9.7)') & 
  &' Note, if we take pi*R^2 = 4*pi/n, we get R = 2/sqrt(n) = ',rchar

write(*,*) & 
  &' Radius of the maximum strength vortex, relative to R (i.e. R_max/R)?'
read(*,*) rmax

rmax=rmax*rchar
 ! taking vortex radius rv(j) = sqrt(frad*sv(j))
frad=rmax**2/svmax

!------------------------------------------------------
 !Select time evolution or fixed time:
write(*,*) ' Here we show an orthographic projection of the flow.'
write(*,*)
write(*,*) ' Image sequence type, (1) time evolution or (2) fixed time?'
read(*,*) sopt

if (sopt .lt. 1 .or. sopt .gt. 2) then
  write(*,*) ' Option not recognised.  *** Stopping ***'
  stop
endif

!------------------------------------------------------
 !Select image size and size of imaged vortices:
write(*,*) ' Image size (width = height) in pixels?'
read(*,*) npix

!------------------------------------------------------
 !Define input and output record lengths (in bytes):
nbread=4*(3*n+1)
nbprec=2*npix*npix/nbpc

!-----------------------------------------------------------------------
 !Open input file to image along with the associated output file:
open(11,file='points.r4',form='unformatted',access='direct', &
          & status='old',recl=nbread)

if (sopt .eq. 1) then
   !Image a time sequence:
  write(*,*) & 
   & ' Enter the beginning & ending frames, and the interval (e.g. 1):'
  read(*,*) kbeg,kend,kint
  kbeg=kbeg+1
  kend=kend+1
   !Open output data file:
  open(22,file='evolution.c2',form='unformatted',access='direct', & 
           & status='unknown',recl=nbprec)
else
   !Image a fixed time:
  write(*,*) ' Time frame (choose 0 for the initial one)?'
  read(*,*) period
  write(suffix(1:4),'(i4.4)') period
  kbeg=period+1
  kend=kbeg
  kint=1
   !Open output data file:
  open(22,file='frame'//suffix//'.c2',form='unformatted',access='direct', &
                   & status='unknown',recl=nbprec)
endif

!-------------------------------------------------------------
 !Convert data to character maps:
call opgrid(npix,kbeg,kend,kint,sopt)

close(11)
close(22)

!============================================================================

 ! Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine p2c
! Converts points (xv,yv,zv) & sv to contours (xd,yd,zd)

implicit none

real:: rv,fac,u1x,u1y,u2x,u2y,u2z,rp,xx,yy,zz,cc
integer:: j,idir,ibeg,i,ii

do j=1,n
  rv =sqrt(xv(j)**2+yv(j)**2)
  fac=1./rv
  u1x=-fac*yv(j)
  u1y= fac*xv(j)
  u2x= fac*xv(j)*zv(j)
  u2y= fac*yv(j)*zv(j)
  u2z=-rv

  if (sv(j) .gt. 0.) then
    idir=-1
    ibeg=j*np
  else
    idir=1
    ibeg=(j-1)*np+1
  endif

  rp=sqrt(frad*abs(sv(j)))

  do ii=0,np-1
    xx=xv(j)+rp*(xc(ii)*u1x+yc(ii)*u2x)
    yy=yv(j)+rp*(xc(ii)*u1y+yc(ii)*u2y)
    zz=zv(j)+rp*yc(ii)*u2z
    cc=1./sqrt(xx**2+yy**2+zz**2)
    i=ibeg+ii*idir
    xd(i)=cc*xx
    yd(i)=cc*yy
    zd(i)=cc*zz
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine c2g(npix,qq)
! Converts contours (xd,yd,zd) to gridded values (qq):

implicit none

 ! Local work arrays:
real:: qq(0:2*npix+1,4*npix)
real:: clon(4*npix),slon(4*npix)
real:: clat(2*npix)
real:: a(npt),c(npt),d(npt)
integer:: ntc(npt),ilm1(npt)

real:: dl,dli,hpidl,rlon,rlat,f1112,rsum,dsumi,deni,sig,rlatc,p,qbar
integer:: npix,nlat,nlon,i,j,k,ka,jump,ioff,ncr,ic

nlat=2*npix
nlon=4*npix

 !----------------------------------------------------
dl=twopi/float(nlon)
dli=float(nlon)/(twopi+1.e-7)
hpidl=(pi+dl)/2.

do i=1,nlon
  rlon=dl*float(i-1)-pi
  clon(i)=cos(rlon)
  slon(i)=sin(rlon)
enddo

do j=1,nlat
  rlat=(float(j)-0.5)*dl-hpi
  clat(j)=cos(rlat)
enddo

f1112=11./12.
rsum=f1112*(clat(1)+clat(nlat))
do j=2,nlat-1
  rsum=rsum+clat(j)
enddo
dsumi=1./(rsum*float(nlon))

 !------------------------------------------------------
do k=1,npt
  ilm1(k)=int(dli*(pi+atan2(yd(k),xd(k))))
enddo

do k=1,npt
  ka=next(k)
  a(k)=xd(k)*yd(ka)-yd(k)*xd(ka)
  c(k)=zd(k)*yd(ka)-yd(k)*zd(ka)
  d(k)=xd(k)*zd(ka)-zd(k)*xd(ka)
  ntc(k)=ilm1(ka)-ilm1(k)
enddo

do k=1,npt
  deni=1./a(k)
  c(k)=c(k)*deni
  d(k)=d(k)*deni
  sig=sign(1.,a(k))
  a(k)=sig
  ntc(k)=ntc(k)-nlon*((2*ntc(k))/nlon)
  if (sig*float(ntc(k)) .lt. 0.) ntc(k)=-ntc(k)
enddo

do i=1,nlon
  do j=1,nlat+1
    qq(j,i)=0.
  enddo
enddo

do k=1,npt
  if (ntc(k) .ne. 0) then
    jump=sign(1,ntc(k))
    ioff=nlon+ilm1(k)+(1+jump)/2
    ncr=0
175   i=1+mod(ioff+ncr,nlon)
      rlatc=dli*(hpi+atan(c(k)*clon(i)+d(k)*slon(i)))
      j=int(rlatc)+1
      p=rlatc-float(j-1)
      qq(j,i)=  qq(j,i)+(1.-p)*a(k)
      qq(j+1,i)=qq(j+1,i)+   p*a(k)
      ncr=ncr+jump
      if (ncr .ne. ntc(k)) goto 175
  endif
enddo

do i=1,nlon
  do j=2,nlat
    qq(j,i)=qq(j,i)+qq(j-1,i)
  enddo
enddo
 ! Here, qq(j,i) stands for the vorticity at latitude j-1/2,
 ! from j = 1, ..., nlat.

 ! Next compute and remove the global mean value of qq:
qbar=0.
do i=1,nlon
  qbar=qbar+f1112*(clat(1)*qq(1,i)+clat(nlat)*qq(nlat,i))
  do j=2,nlat-1
    qbar=qbar+clat(j)*qq(j,i)
  enddo
enddo
qbar=qbar*dsumi

do i=1,nlon
  do j=1,nlat
    qq(j,i)=qq(j,i)-qbar
  enddo
enddo

 ! Copy latitudes adjacent to poles (j = 1 and nlat) with a pi
 ! shift in longitude to simplify interpolation below:
do i=1,nlat
  ic=i+nlat
  qq(0,i)=qq(1,ic)
  qq(0,ic)=qq(1,i)
  qq(nlat+1,i)=qq(nlat,ic)
  qq(nlat+1,ic)=qq(nlat,i)
enddo

return
end subroutine

!=======================================================================

subroutine opgrid(npix,kbeg,kend,kint,sopt)
! Converts data to a character map in an orthographic projection

implicit none

integer:: npix,ii

integer:: kbeg,kend,kint,sopt,i,j,k,loop,ich,ich1,ich2
integer:: nlat,nlon,iy,iz,ic,ip1,jp1,nalp

real:: qq(0:2*npix+1,4*npix)
real:: xg(npix,npix),yg(npix),zg(npix)
real:: dl,dli,hpidl
real:: dy,dz,rlatc,clatc,slatc,rlonc,clonc,slonc
real:: xp,yp,zp,xm,xt,yt,zt,ri,rj,aa,bb,cc,dd,qqt
real:: t11,t12,t21,t22,t31,alp,dalp,calp,salp
real:: a11,a12,a13,a21,a22,a23,a31,a32,a33

character(len=2):: bqq(npix*npix),white

logical:: inside(npix,npix)

!------------------------------------------------------------------
 !For creating a white surround of each image:
white=char(0)//char(0)

 !Work out points inside disk of view:
dy=2./float(npix)
dz=dy
do iy=1,npix
  yg(iy)=dy*(float(iy)-0.5)-1.
  do iz=1,npix
    zg(iz)=dz*(float(iz)-0.5)-1.
    inside(iz,iy)=(yg(iy)**2+zg(iz)**2 .lt. 1.)
    if (inside(iz,iy)) xg(iz,iy)=sqrt(1.-yg(iy)**2-zg(iz)**2)
  enddo
enddo

 !For converting gridded values from a lat/lon grid:
nlat=2*npix
nlon=4*npix

dl=twopi/float(nlon)
dli=float(nlon)/(twopi+1.e-7)
hpidl=(pi+dl)/2.

!------------------------------------------------------------------
 !Either show a time sequence from a fixed perspective or a fixed 
 !time rotated about a specified axis:
if (sopt .eq. 1) then
   !Show a time sequence from a fixed perspective:

  write(*,*) ' Latitude & longitude of the direction of view (degrees)?'
  read(*,*) rlatc,rlonc
  rlatc=rlatc*pi/180.
  rlonc=rlonc*pi/180.

  clonc=cos(rlonc)
  slonc=sin(rlonc)
  clatc=cos(rlatc)
  slatc=sin(rlatc)

  loop=0
  do k=kbeg,kend,kint
    loop=loop+1
     !Read each frame of the data:
    read(11,rec=k) t,xv,yv,zv

     !Generate contour nodes:
    call p2c

     !Convert contours to gridded values:
    call c2g(npix,qq)

     !Do interpolation in chosen perspective and convert each data value
     !to a pair of characters:
    do iy=1,npix
      yp=yg(iy)
      do iz=1,npix
        if (inside(iz,iy)) then
          xp=xg(iz,iy)
          zp=zg(iz)

          xm=xp*clatc-zp*slatc
          zt=zp*clatc+xp*slatc
          yt=yp*clonc+xm*slonc
          xt=xm*clonc-yp*slonc

           !Find lat & lon then bi-linearly interpolate qq:
          ri=dli*(pi+atan2(yt,xt))
          i=1+int(ri)
          ip1=1+mod(i,nlon)
          bb=float(i)-ri
          aa=1.-bb

          rj=dli*(hpidl+asin(zt))
          j=int(rj)
          jp1=j+1
          cc=rj-float(j)
          dd=1.-cc

          qqt=bb*(dd*qq(j,i)+cc*qq(jp1,i))+aa*(dd*qq(j,ip1)+cc*qq(jp1,ip1))
          ich=nint(coff+cfac*(min(qqmax,max(qqmin,qqt))-qqmin))
          ich1=ich/256+1
          ich2=ich-(ich1-1)*256+1
          bqq(npix*(iz-1)+iy)=char(ich1)//char(ich2)
        else
          bqq(npix*(iz-1)+iy)=white
        endif
      enddo
    enddo

     !Write character image:
    write(22,rec=loop) bqq
  enddo

else
   !Show a fixed time rotated about a specified axis:
  write(*,*) ' Latitude & longitude of the axis of rotation (degrees)?'
  read(*,*) rlatc,rlonc
  rlatc=rlatc*pi/180.
  rlonc=rlonc*pi/180.

   !Prepare rotation matrix:
  clonc=cos(rlonc)
  slonc=sin(rlonc)
  clatc=cos(rlatc)
  slatc=sin(rlatc)
  t11=clonc*slatc
  t12=-slonc
  a13=clonc*clatc
  t21=slonc*slatc
  t22=clonc
  a23=slonc*clatc
  t31=-clatc
  a33=slatc

  write(*,*) ' Number of frames to display in one full rotation?'
  read(*,*) nalp
  dalp=twopi/float(nalp)

   !Read a single frame of data:
  read(11,rec=kbeg) t,xv,yv,zv

   !Generate contour nodes:
  call p2c

   !Convert contours to gridded values:
  call c2g(npix,qq)

   !Loop over rotation angle alpha
  do loop=1,nalp
    alp=dalp*float(loop-1)
    calp=cos(alp)
    salp=sin(alp)

    a11=t11*calp-t12*salp
    a12=t12*calp+t11*salp
    a21=t21*calp-t22*salp
    a22=t22*calp+t21*salp
    a31=t31*calp
    a32=t31*salp

     !Do interpolation in chosen perspective and convert each data value
     !to a pair of characters:
    do iy=1,npix
      yp=yg(iy)
      do iz=1,npix
        if (inside(iz,iy)) then
          xp=xg(iz,iy)
          zp=zg(iz)

          xt=a11*xp+a12*yp+a13*zp
          yt=a21*xp+a22*yp+a23*zp
          zt=a31*xp+a32*yp+a33*zp

           !Find lat & lon then bi-linearly interpolate qq:
          ri=dli*(pi+atan2(yt,xt))
          i=1+int(ri)
          ip1=1+mod(i,nlon)
          bb=float(i)-ri
          aa=1.-bb

          rj=dli*(hpidl+asin(zt))
          j=int(rj)
          jp1=j+1
          cc=rj-float(j)
          dd=1.-cc

          qqt=bb*(dd*qq(j,i)+cc*qq(jp1,i))+aa*(dd*qq(j,ip1)+cc*qq(jp1,ip1))
          ich=nint(coff+cfac*(min(qqmax,max(qqmin,qqt))-qqmin))
          ich1=ich/256+1
          ich2=ich-(ich1-1)*256+1
          bqq(npix*(iz-1)+iy)=char(ich1)//char(ich2)
        else
          bqq(npix*(iz-1)+iy)=white
        endif
      enddo
    enddo

     !Write character image:
    write(22,rec=loop) bqq
  enddo
endif

write(*,*)
write(*,'(2(a,i5))') ' Note, the image has dimensions ',npix,' by ',npix

end subroutine

!=========================================================

end program
