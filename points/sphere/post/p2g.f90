program p2g
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Converts point vortices in points.dat to gridded values on a
! longitude-latitude grid of dimensions (2*ng) x ng and writes 
! the data to fgz???, where ??? = 000, 001 ... is the chosen frame.

! Note: the value of ng is specified by the user in parameters.f90
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 ! Import parameters and constants:
use parameters
use constants

implicit double precision (a-h,o-z)

integer,parameter:: nt=2*ng
 ! nt: number of longitudes
integer,parameter:: ntp1=nt+1,ngp1=ng+1

integer,parameter:: np=20, npt=n*np
 ! np:  the number of nodes used to image each vortex as a small patch
 ! npt: the total number of contour nodes

double precision,parameter:: pi=3.1415926535897932385d0,twopi=two*pi
double precision,parameter:: hpi=half*pi

 ! point vortex arrays:
double precision:: sv(n),xv(n),yv(n),zv(n)
double precision:: xc(0:np-1),yc(0:np-1)

 ! vorticity contour arrays:
double precision:: xd(npt),yd(npt),zd(npt)
double precision:: dx(npt),dy(npt),dz(npt)
double precision:: a(npt), c(npt), d(npt)
integer:: next(npt)

 ! Contour -> Grid arrays:
double precision:: qa(0:ngp1,ntp1)
double precision:: clon(nt),slon(nt)
double precision:: clat(ng)
integer:: ntc(npt),ilm1(npt)

character:: cper*3

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
  read(20,*) sv(n)
enddo
close(20)

svmax=zero
do k=1,n
  svmax=max(svmax,abs(sv(k)))
enddo

 !-------------------------------------------------
 ! Constants used in contour->grid conversion:
dl =twopi/dble(nt)
dli=dble(nt)/(twopi+small)

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

 !----------------------------------------------
 ! Select period to process and open files:
write(*,*) ' Time frame to image (this is 0 for t = 0)?'
read(*,*) kfr
write(*,*) ' Radius of the maximum strength vortex?'
read(*,*) rmax

frad=rmax/svmax

 ! Read vortex positions:
open(11,file='points.dat',status='old')

if (kfr .gt. 0) then
   ! skip data:
  do k=1,(n+1)*kfr
    read(11,*)
  enddo
endif

read(11,*) t
do k=1,n
  read(11,*) xv(k),yv(k),zv(k)
enddo
close(11)

write(*,*)
write(*,'(a,f12.5)') ' *** Imaging t = ',t

 ! Generate contour nodes:
do j=1,n
  rv =sqrt(xv(j)**2+yv(j)**2)
  fac=one/rv
  u1x=-fac*yv(j)
  u1y= fac*xv(j)
  u2x= fac*xv(j)*zv(j)
  u2y= fac*yv(j)*zv(j)
  u2z=-rv

  if (sv(j) .gt. zero) then
    idir=1
    ibeg=(j-1)*np+1
  else
    idir=-1
    ibeg=j*np
  endif

  rp=frad*abs(sv(j))

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

 ! Define next(i), which gives the node following node i:
do i=1,npt-1
  next(i)=i+1
enddo
do j=1,n
  next(j*np)=(j-1)*np+1
enddo

 ! Convert contours to gridded values:
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
    qa(j,i)=zero
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
      qa(j,i)=  qa(j,i)+(one-p)*a(k)
      qa(j+1,i)=qa(j+1,i)+    p*a(k)
      ncr=ncr+jump
      if (ncr .ne. ntc(k)) goto 175
  endif
enddo

do i=1,nt
  do j=2,ng
    qa(j,i)=qa(j,i)+qa(j-1,i)
  enddo
enddo
 ! Here, qa(j,i) stands for the vorticity at latitude j-1/2,
 ! from j = 1, ..., ng.

 ! Next compute and remove the global mean value of qa:
qbar=zero
do i=1,nt
  qbar=qbar+f1112*(clat(1)*qa(1,i)+clat(ng)*qa(ng,i))
  do j=2,ng-1
    qbar=qbar+clat(j)*qa(j,i)
  enddo
enddo
qbar=qbar*dsumi

do i=1,nt
  do j=1,ng
    qa(j,i)=qa(j,i)-qbar
  enddo
enddo

 ! Copy latitudes adjacent to poles (j = 1 and ng) with a pi
 ! shift in longitude to simplify interpolation in zfgrid.F:
do i=1,ng
  ic=i+ng
  qa(0,i)=qa(1,ic)
  qa(0,ic)=qa(1,i)
  qa(ngp1,i)=qa(ng,ic)
  qa(ngp1,ic)=qa(ng,i)
enddo

 ! Add a periodic column at i = ntp1:
i=ntp1
do j=0,ngp1
  qa(j,i)=qa(j,1)
enddo

 ! Write data:
write(cper(1:3),'(i3.3)') kfr
open(18,file='fgz'//cper,form='unformatted',status='unknown')
write(18) qa
close(18)

write(*,*) 
write(*,*) ' Fine grid vorticity written to fgz'//cper

 !End main program
end program
!=======================================================================
