module generic

! Module contains generic subroutines for casl in a doubly periodic geometry. 
!          *** These should not be modified ***

use constants

implicit none

contains 

!=======================================================================

subroutine interpol(fp,fe,xp,yp)

! Interpolates the field fp(iy,ix) at a given set of points
! xp(iy,ix), yp(iy,ix) (given in grid units) using 
! bi-cubic Lagrange interpolation.

! The interpolated field is written into the array fe. 
! Note, the domain average value of fe is removed.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: fp(ny,nx),fe(ny,nx)
double precision:: xp(ny,nx),yp(ny,nx)
 !Local arrays:
double precision:: phix(-1:2),phiy(-1:2)
integer:: ixper(-1:nx+1),iyper(-1:ny+1)
!---------------------------------------------------------------
 !Grid box reference index used below for x periodicity:
ixper(-1)=nx
do ix=0,nx-1
  ixper(ix)=ix+1
enddo
ixper(nx)=1
ixper(nx+1)=2

iyper(-1)=ny
do iy=0,ny-1
  iyper(iy)=iy+1
enddo
iyper(ny)=1
iyper(ny+1)=2

do ix=1,nx
  do iy=1,ny
    ix0=int(xp(iy,ix))
     !ix0: the x grid box containing the point
    px0=xp(iy,ix)-dble(ix0)

    pxm=one+px0
    px1=one-px0
    px2=two-px0
    phix(-1)=-f16*px0*px1*px2
    phix( 0)= f12*pxm*px1*px2
    phix( 1)= f12*pxm*px0*px2
    phix( 2)=-f16*pxm*px0*px1

    iy0=int(yp(iy,ix))
     !iy0: the y grid box containing the point
    py0=yp(iy,ix)-dble(iy0)

    pym=one+py0
    py1=one-py0
    py2=two-py0
    phiy(-1)=-f16*py0*py1*py2
    phiy( 0)= f12*pym*py1*py2
    phiy( 1)= f12*pym*py0*py2
    phiy( 2)=-f16*pym*py0*py1

    fe(iy,ix)=zero
    do jx=-1,2
      ixi=ixper(ix0+jx)
      do jy=-1,2
        iyi=iyper(iy0+jy)
        fe(iy,ix)=fe(iy,ix)+fp(iyi,ixi)*phix(jx)*phiy(jy)
      enddo
    enddo
     !Clip function to min/max values at corners of grid box:
    ix0p1=ixper(ix0+1)
    ix0  =ixper(ix0)
    iy0p1=iyper(iy0+1)
    iy0  =iyper(iy0)
    femin=min(fp(iy0,ix0),fp(iy0p1,ix0),fp(iy0,ix0p1),fp(iy0p1,ix0p1))
    femax=max(fp(iy0,ix0),fp(iy0p1,ix0),fp(iy0,ix0p1),fp(iy0p1,ix0p1))
    fe(iy,ix)=min(femax,max(femin,fe(iy,ix)))
  enddo
enddo

 !Remove domain average:
call average(fe,favg)
do ix=1,nx
  do iy=1,ny
    fe(iy,ix)=fe(iy,ix)-favg
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine combine(vqw,vqc,vqs,vqd)

! Combines contour (vqc), large scale (vqs), and residual (vqd) fields into 
! the full PV field (vqw).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: vqw(ny,nx),vqc(ny,nx),vqs(ny,nx),vqd(ny,nx)
 !Define local array:
double precision:: wka(ny,nx)

!Define vq = F[vqs-vqc]+vqc+vqd
!where F is a low pass filter (see subroutine filter)

do ix=1,nx
  do iy=1,ny
    wka(iy,ix)=vqs(iy,ix)-vqc(iy,ix)
  enddo
enddo

call filter(wka,2)

do ix=1,nx
  do iy=1,ny
    vqw(iy,ix)=wka(iy,ix)+vqc(iy,ix)+vqd(iy,ix)
  enddo
enddo

return
end subroutine
!=======================================================================

subroutine reset(vqc,vqs,vqd)

! Resets the gridded fields vqs & vqd

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: vqc(ny,nx),vqs(ny,nx),vqd(ny,nx)
 !Define local array:
double precision:: wka(ny,nx)

 !------------------------------------------------------------
 !Reset vqs = vq = F[vqs-vqc]+vqc+vqd, where F is a low pass filter:
do ix=1,nx
  do iy=1,ny
    wka(iy,ix)=vqs(iy,ix)-vqc(iy,ix)
  enddo
enddo

call filter(wka,2)

 !Reset vqd = vq-vqc-F[vq-vqc]
do ix=1,nx
  do iy=1,ny
    vqs(iy,ix)=wka(iy,ix)+vqc(iy,ix)+vqd(iy,ix)
    vqd(iy,ix)=vqs(iy,ix)-vqc(iy,ix)
    wka(iy,ix)=vqd(iy,ix)
  enddo
enddo

call filter(wka,2)

do ix=1,nx
  do iy=1,ny
    vqd(iy,ix)=vqd(iy,ix)-wka(iy,ix)
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine filter(qw,nrep)

! Performs nrep 1-2-1 filters in each direction to return a 
! low pass filtered version of the original array.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qw(ny,nx)
 !Define local array:
double precision:: wka(ny,nx)

 !Perform 1-2-1 average nrep times utilising intermediate work array:
do j=1,nrep
  do iy=1,ny
    wka(iy, 1)=f12*qw(iy, 1)+f14*(qw(iy,  nx)+qw(iy,2))
    wka(iy,nx)=f12*qw(iy,nx)+f14*(qw(iy,nxm1)+qw(iy,1))
  enddo
  do ix=2,nxm1
    do iy=1,ny
      wka(iy,ix)=f12*qw(iy,ix)+f14*(qw(iy,ix-1)+qw(iy,ix+1))
    enddo
  enddo

  do ix=1,nx
    qw(1, ix)=f12*wka(1 ,ix)+f14*(wka(ny  ,ix)+wka(2,ix))
    qw(ny,ix)=f12*wka(ny,ix)+f14*(wka(nym1,ix)+wka(1,ix))
  enddo
  do ix=1,nx 
    do iy=2,nym1
      qw(iy,ix)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
    enddo
  enddo
enddo

return
end subroutine

!=======================================================================

subroutine average(qw,qavg)

! Computes the average value of a field qw and returns the result in qavg

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qw(ny,nx)

qavg=zero
do ix=1,nx
  do iy=1,ny
    qavg=qavg+qw(iy,ix)
  enddo
enddo
qavg=qavg/dble(nx*ny)

return
end subroutine

!=======================================================================

subroutine l1norm(qw,ql1)

! Computes the L1 norm of a field qw and returns the result in ql1:
!          ql1 = int_xmin^xmax{int_ymin^ymax{|qw| dxdy}}

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qw(ny,nx)

ql1=zero
do ix=1,nx
  do iy=1,ny
    ql1=ql1+abs(qw(iy,ix))
  enddo
enddo

ql1=garea*ql1
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

subroutine l2norm(qw,ql2)

! Computes the L2 norm of a field qw and returns the result in ql2:
!          ql2 = int_xmin^xmax{int_ymin^ymax{qw^2 dxdy}}

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qw(ny,nx)

ql2=zero
do ix=1,nx
  do iy=1,ny
    ql2=ql2+qw(iy,ix)**2
  enddo
enddo

ql2=garea*ql2
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

subroutine binorm(qw1,qw2,qwbi)

! Computes the integral of qw1*qw2 over the domain

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qw1(ny,nx),qw2(ny,nx)

qwbi=zero
do ix=1,nx
  do iy=1,ny
    qwbi=qwbi+qw1(iy,ix)*qw2(iy,ix)
  enddo
enddo

qwbi=garea*qwbi
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

 !Main end module
end module
