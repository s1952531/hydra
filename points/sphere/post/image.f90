program image
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Converts point vortices in points.dat to gridded values on a
! longitude-latitude grid of dimensions (2*ng) x ng and creates
! a sequence of images in orthographic perspective.

! One can either choose to image a fixed time, viewed rotated about
! a fixed axis, or a range of times, viewed from one or several 
! directions.

! The image dimensions are ng x ng.  The images can be viewed using
! xvidi, or converted to a gif movie using the script "b2gif".

! Note: the value of ng is specified by the user in parameters.f90

!     Written 30 October 2012 by D G Dritschel @ St Andrews
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 ! Import parameters and constants:
use parameters
use constants

implicit double precision (a-h,o-z)

integer,parameter:: nt=2*ng
 ! nt: number of longitudes
integer,parameter:: ntp1=nt+1, ngp1=ng+1

integer,parameter:: np=20, npt=n*np
 ! np:  the number of nodes used to image each vortex as a small patch
 ! npt: the total number of contour nodes

integer,parameter:: ny=ng, nz=ng, ngp=nz*ny
 ! ny, nz: image dimensions (to keep subroutine g2i general)

real,parameter:: cfrac=1024./(256.*256)
 !This parameter is used to access the fraction of the colourmap 
 !extending from cfrac to 1-cfrac (to avoid oversaturation)
real:: crange,coff,cfac

double precision,parameter:: pi=3.1415926535897932385d0, twopi=two*pi
double precision,parameter:: hpi=half*pi
double precision,parameter:: f1112=11.d0/12.d0

 ! Angular spacing of grid, etc:
double precision,parameter:: dl=twopi/dble(nt), dli=dble(nt)/(twopi+1.d-12)
double precision,parameter:: hpidl=(pi+dl)/two

 ! point vortex arrays:
double precision:: sv(n),xv(n),yv(n),zv(n)
double precision:: xc(0:np-1),yc(0:np-1)

 ! vorticity contour arrays:
double precision:: xd(npt),yd(npt),zd(npt)
integer:: next(npt)

 ! Contour -> Grid arrays:
double precision:: qq(0:ngp1,ntp1)
double precision:: clon(nt),slon(nt)
double precision:: clat(ng)

 ! Various global constants:
double precision:: frad,dsumi

 ! For creating character movie image:
double precision:: xg(nz,ny),yg(ny),zg(nz)
logical:: inside(nz,ny)
character:: xvfile*11,cfr1*3,cfr2*3

 !-------------------------------------------
 ! Initialise all fixed constants and arrays:
call init

 !Quantities required for converting data into characters:
crange=float(256**2)
coff=cfrac*crange
cfac=crange/2.0-coff

 !---------------------------------------------------------------------------
 ! Select type of image sequence:
write(*,*) ' Two types of imaging are possible; choose one of the following:'
write(*,*) ' (1) fixed time, rotated about an axis, or'
write(*,*) ' (2) range of times, from a chosen viewpoint?'
read(*,*) iopt

if (iopt .eq. 1) then
  write(*,*) ' Time frame to image (this is 0 for t = 0)?'
  read(*,*) kfr1
  kfr2=kfr1
  write(*,*) '  Enter the latitude and longitude of the'
  write(*,*) 'axis of rotation (in degrees):'
else
  write(*,*) ' Range of time frames to image (e.g. 0 100)?'
  read(*,*) kfr1,kfr2
  write(*,*) '  Enter the latitude and longitude of the'
  write(*,*) 'direction of view (in degrees):'
endif
read(*,*) rlatc,rlonc
rlatc=rlatc*pi/180.d0
rlonc=rlonc*pi/180.d0

clonc=cos(rlonc)
slonc=sin(rlonc)
clatc=cos(rlatc)
slatc=sin(rlatc)

if (iopt .eq. 1) then
   ! Prepare rotation matrix for a fixed time:
  t11=clonc*slatc
  t12=-slonc
  a13=clonc*clatc
  t21=slonc*slatc
  t22=clonc
  a23=slonc*clatc
  t31=-clatc
  a33=slatc

  write(*,*) '  Enter the number of frames to display'
  write(*,*) 'in one full rotation:'
  read(*,*) nalp
  dalp=twopi/dble(nalp)

else
   ! Define fixed rotation matrix for all times:
  a11=clatc*clonc
  a12=-slonc
  a13=-slatc*clonc
  a21=clatc*slonc
  a22=clonc
  a23=-slatc*slonc
  a31=slatc
  a32=zero
  a33=clatc
endif

 !--------------------------------------
 ! Create and open character movie file:
write(cfr1(1:3),'(i3.3)') kfr1
write(cfr2(1:3),'(i3.3)') kfr2
xvfile='z'//cfr1//'-'//cfr2//'.c2'
open(21,file=xvfile,form='unformatted',access='direct', &
               &  status='unknown',recl=2*ngp)
 ! recl = ngp generally, but sometimes ngp/4.
loop=0

 !----------------------------------------
 ! Read vortex positions and process data:
open(11,file='points.dat',status='old')

if (kfr1 .gt. 0) then
   ! skip data:
  do k=1,(n+1)*kfr1
    read(11,*)
  enddo
endif

 ! Loop over selected time frames:
do kfr=kfr1,kfr2
   !Read each frame of the data with error handler:
  call readframe(ierr)
  if (ierr .eq. 1) exit

  write(*,'(a,f12.5)') ' *** Imaging t = ',t

   ! Generate contour nodes:
  call p2c

   ! Convert contours to gridded values:
  call c2g

  if (iopt .eq. 1) then
     ! Loop over rotation angle alpha
    do loop=1,nalp
      alp=-dalp*dble(loop-1)
      calp=cos(alp)
      salp=sin(alp)

      a11=t11*calp+t12*salp
      a12=t12*calp-t11*salp
      a21=t21*calp+t22*salp
      a22=t22*calp-t21*salp
      a31=t31*calp
      a32=-t31*salp

       ! Convert gridded values to a character map (movie image):
      call g2i(loop)
    enddo
  else
     ! Here, we are looping over time with a fixed viewpoint;
     ! Convert gridded values to a character map (movie image):
    loop=loop+1
    call g2i(loop)
  endif

enddo

 ! Close files:
close(11)
close(21)

 ! Finish up:
write(*,*)
write(*,*) ' All done.  To display the movie, type'
write(*,*)
write(*,'(a,2(i4),a)') ' dataview -ndim ',ny,nz,' '//xvfile//' &'

!============================================================================

 ! Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine init
! Converts points (xv,yv,zv) & sv to contours (xd,yd,zd)

implicit double precision (a-h,o-z)

 !--------------------------------------------------------
 ! For defining circular patches around each point vortex:
dth=twopi/dble(np)
do i=0,np-1
  th=dth*dble(i)
  xc(i)=cos(th)
  yc(i)=sin(th)
enddo

 !---------------------------------------------------
 ! Read vortex strengths:
open(20,file='strengths.dat',status='old')
do k=1,n
  read(20,*) sv(k)
enddo
close(20)

svmax=zero
do k=1,n
  svmax=max(svmax,abs(sv(k)))
enddo

write(*,'(2(a,i4))') ' *** Note, the image dimensions will be ',ny,' x ',nz
write(*,*)

rchar=two/sqrt(dble(n))
write(*,'(a,f11.9)') & 
  &' Note, if we take pi*R^2 = 4*pi/n, we get R = 2/sqrt(n) = ',rchar

write(*,*) & 
  &' Radius of the maximum strength vortex, relative to R (i.e. R_max/R)?'
read(*,*) rmaxnd

rmax=rmaxnd*rchar
 ! taking vortex radius rv(j) = sqrt(frad*sv(j))
frad=rmax**2/svmax

 !---------------------------------------------
 ! Quantities used in contour->grid conversion:
do i=1,nt
  rlon=dl*dble(i-1)-pi
  clon(i)=cos(rlon)
  slon(i)=sin(rlon)
enddo

do j=1,ng
  rlat=(dble(j)-half)*dl-hpi
  clat(j)=cos(rlat)
enddo

rsum=f1112*(clat(1)+clat(ng))
do j=2,ng-1
  rsum=rsum+clat(j)
enddo
dsumi=one/(rsum*dble(nt))

 !--------------------------------------------
 ! Initialise grid points inside disk of view:
dy=two/dble(ny)
dz=dy
do iy=1,ny
yg(iy)=dy*(dble(iy)-half)-one
  do iz=1,nz
    zg(iz)=dz*(dble(iz)-half)-one
    inside(iz,iy)=(yg(iy)**2+zg(iz)**2 .lt. one)
    if (inside(iz,iy)) xg(iz,iy)=sqrt(one-yg(iy)**2-zg(iz)**2)
  enddo
enddo

 !-------------------------------------------------------------
 ! For use in p2c, define next(i), which gives the contour node
 ! following node i:
do i=1,npt-1
  next(i)=i+1
enddo
do j=1,n
  next(j*np)=(j-1)*np+1
enddo

return
end subroutine

!=======================================================================

subroutine p2c
! Converts points (xv,yv,zv) & sv to contours (xd,yd,zd)

implicit double precision (a-h,o-z)

do j=1,n
  rv =sqrt(xv(j)**2+yv(j)**2)
  fac=one/rv
  u1x=-fac*yv(j)
  u1y= fac*xv(j)
  u2x= fac*xv(j)*zv(j)
  u2y= fac*yv(j)*zv(j)
  u2z=-rv

  if (sv(j) .gt. zero) then
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
    cc=one/sqrt(xx**2+yy**2+zz**2)
    i=ibeg+ii*idir
    xd(i)=cc*xx
    yd(i)=cc*yy
    zd(i)=cc*zz
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine c2g
! Converts contours (xd,yd,zd) to gridded values (qq):

implicit double precision (a-h,o-z)

 ! Local work arrays:
double precision:: a(npt), c(npt), d(npt)
integer:: ntc(npt),ilm1(npt)

do k=1,npt
  ilm1(k)=int(dli*(pi+atan2(yd(k),xd(k)+1.d-10)))
enddo

do k=1,npt
  ka=next(k)
  a(k)=xd(k)*yd(ka)-yd(k)*xd(ka)
  c(k)=zd(k)*yd(ka)-yd(k)*zd(ka)
  d(k)=xd(k)*zd(ka)-zd(k)*xd(ka)
  ntc(k)=ilm1(ka)-ilm1(k)
enddo

do k=1,npt
  deni=one/a(k)
  c(k)=c(k)*deni
  d(k)=d(k)*deni
  sig=sign(one,a(k))
  a(k)=sig
  ntc(k)=ntc(k)-nt*((2*ntc(k))/nt)
  if (sig*dble(ntc(k)) .lt. zero) ntc(k)=-ntc(k)
enddo

do i=1,nt
  do j=1,ngp1
    qq(j,i)=zero
  enddo
enddo

do k=1,npt
  if (ntc(k) .ne. 0) then
    jump=sign(1,ntc(k))
    ioff=nt+ilm1(k)+(1+jump)/2
    ncr=0
175   i=1+mod(ioff+ncr,nt)
      rlatc=dli*(hpi+atan(c(k)*clon(i)+d(k)*slon(i)))
      j=int(rlatc)+1
      p=rlatc-dble(j-1)
      qq(j,i)=  qq(j,i)+(one-p)*a(k)
      qq(j+1,i)=qq(j+1,i)+    p*a(k)
      ncr=ncr+jump
      if (ncr .ne. ntc(k)) goto 175
  endif
enddo

do i=1,nt
  do j=2,ng
    qq(j,i)=qq(j,i)+qq(j-1,i)
  enddo
enddo
 ! Here, qq(j,i) stands for the vorticity at latitude j-1/2,
 ! from j = 1, ..., ng.

 ! Next compute and remove the global mean value of qq:
qbar=zero
do i=1,nt
  qbar=qbar+f1112*(clat(1)*qq(1,i)+clat(ng)*qq(ng,i))
  do j=2,ng-1
    qbar=qbar+clat(j)*qq(j,i)
  enddo
enddo
qbar=qbar*dsumi

do i=1,nt
  do j=1,ng
    qq(j,i)=qq(j,i)-qbar
  enddo
enddo

 ! Copy latitudes adjacent to poles (j = 1 and ng) with a pi shift
 ! in longitude to simplify interpolation in subroutine g2i:
do i=1,ng
  ic=i+ng
  qq(0,i)=qq(1,ic)
  qq(0,ic)=qq(1,i)
  qq(ngp1,i)=qq(ng,ic)
  qq(ngp1,ic)=qq(ng,i)
enddo

 ! Add a periodic column at i = ntp1:
i=ntp1
do j=0,ngp1
  qq(j,i)=qq(j,1)
enddo

return
end subroutine

!=======================================================================

subroutine g2i(loop)
! Converts points (xv,yv,zv) & sv to contours (xd,yd,zd)

implicit double precision (a-h,o-z)

character(len=2):: bqq(nz,ny), white

 ! For creating "white" colour:
white=char(0)//char(0)

do iy=1,ny
  y=yg(iy)
  do iz=1,nz
    if (inside(iz,iy)) then
      x=xg(iz,iy)
      z=zg(iz)

       !Apply rotation matrix:
      xt=a11*x+a12*y+a13*z
      yt=a21*x+a22*y+a23*z
      zt=a31*x+a32*y+a33*z

       !Find latitude & longitude then bi-linearly interpolate qq:
      ri=dli*(pi+atan2(yt,xt))
      i=1+int(ri)
      ip1=i+1
      bb=dble(i)-ri
      aa=one-bb

      rj=dli*(hpidl+asin(zt))
      j=int(rj)
      jp1=j+1
      cc=rj-dble(j)
      dd=one-cc

      qqt=bb*(dd*qq(j,i)+cc*qq(jp1,i))+aa*(dd*qq(j,ip1)+cc*qq(jp1,ip1))

      ich=nint(coff+cfac*(min(1.0,max(-1.0,qq(j,i)))+1.0))
      ich1=ich/256+1
      ich2=ich-(ich1-1)*256+1
      bqq(iz,iy)=char(ich1)//char(ich2)
    else
      bqq(iz,iy)=white
    endif
  enddo
enddo

 ! Write character pixel map for this frame:
write(21,rec=loop) bqq

return
end subroutine

!=======================================================================

subroutine readframe(ierr)

! Subroutine reads a frame of the file. Exits with ierr=0 if 
! the frame was read successfully, 1 otherwise.

implicit double precision(a-h,o-z)
implicit integer(i-n)

ierr=0
read(11,*,iostat=itmp) t
if (itmp .ne. 0) then
  ierr=1
  return
endif

do k=1,n
  read(11,*,iostat=itmp) xv(k),yv(k),zv(k)
  if (itmp .ne. 0) then 
    ierr=1
    return 
  endif
enddo

end subroutine

!=======================================================================

 !End main program
end program
!=======================================================================
