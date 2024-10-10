module generic

! Module contains generic subroutines for the casl algorithm 
! in a multi-layer aperiodic geometry. 

!          *** These should not be modified ***

use constants

implicit none

contains 

!=======================================================================

subroutine interpol(fp,fe,xp,yp)

! Interpolates the field fp(iy,ix,iz) at a given set of points
! xp(iy,ix,iz), yp(iy,ix,iz) (given in grid units) using 
! bi-cubic Lagrange interpolation.

! The interpolated field is written into the array fe.

implicit none

 !Passed arrays:
double precision:: fp(0:ny,0:nxm1,nz),fe(0:ny,0:nxm1,nz)
double precision:: xp(0:ny,0:nxm1,nz),yp(0:ny,0:nxm1,nz)
 !Local arrays:
double precision:: phix(-1:2),phiy(-1:2)
integer:: ixper(-1:nx+2)

 !Other local variables:
double precision:: px0,pxm,px1,px2,py0,pym,py1,py2
integer:: ix,iy,iz,ix0,iy0,ixi,iyi,jx,jy

!---------------------------------------------------------------
 !Grid box reference index used below for x periodicity:
ixper(-1)=nxm1
do ix=0,nxm1
   ixper(ix)=ix
enddo
ixper(nx  )=0
ixper(nx+1)=1
ixper(nx+2)=2

do iz=1,nz
   do ix=0,nxm1
      do iy=0,ny
         ix0=int(xp(iy,ix,iz))
         !ix0: the x grid box containing the point
         px0=xp(iy,ix,iz)-dble(ix0)

         pxm=one+px0
         px1=one-px0
         px2=two-px0
         phix(-1)=-f16*px0*px1*px2
         phix( 0)= f12*pxm*px1*px2
         phix( 1)= f12*pxm*px0*px2
         phix( 2)=-f16*pxm*px0*px1

         iy0=min(int(yp(iy,ix,iz)),nym1)
         !iy0: the y grid box containing the point
         py0=yp(iy,ix,iz)-dble(iy0)

         pym=one+py0
         py1=one-py0
         py2=two-py0
         phiy(-1)=-f16*py0*py1*py2
         phiy( 0)= f12*pym*py1*py2
         phiy( 1)= f12*pym*py0*py2
         phiy( 2)=-f16*pym*py0*py1

         fe(iy,ix,iz)=zero
         do jx=-1,2
            ixi=ixper(ix0+jx)
            do jy=-1,2
               iyi=ny-abs(ny-abs(iy0+jy))
               fe(iy,ix,iz)=fe(iy,ix,iz)+fp(iyi,ixi,iz)*phix(jx)*phiy(jy)
            enddo
         enddo
      enddo
   enddo
enddo

return
end subroutine interpol

!=======================================================================

subroutine combine(qq,qc,qs,qd,qavg)

! Combines contour (qc), large scale (qs), and residual (qd) fields 
! into qq, ensuring the domain average of qq = qavg in each layer.

implicit none

 !Define passed arrays:
double precision:: qq(0:ny,0:nxm1,nz),qc(0:ny,0:nxm1,nz)
double precision:: qs(0:ny,0:nxm1,nz),qd(0:ny,0:nxm1,nz)
double precision:: qavg(nz)
 !Define local arrays:
double precision:: wka(0:ny,0:nxm1,nz)
 !Local variables:
double precision:: qavg0,qqadd
integer:: iz

!--------------------------------------------------------------
 !Define q = F[qs-qc]+qc+qd where F is a low pass filter
 !(see subroutine filter):
wka=qs-qc
call filter(wka,0,2)
qq=wka+qc+qd

 !Restore domain average:
do iz=1,nz
   call average(qq(0,0,iz),qavg0)
   qqadd=qavg(iz)-qavg0
   qq(:,:,iz)=qq(:,:,iz)+qqadd
enddo

return
end subroutine combine

!=======================================================================

subroutine reset(qc,qs,qd,qavg)

! Resets the gridded fields qs & qd and ensures that <qs> = qavg
! in each layer

implicit none

 !Define passed arrays:
double precision:: qc(0:ny,0:nxm1,nz),qs(0:ny,0:nxm1,nz),qd(0:ny,0:nxm1,nz)
double precision:: qavg(nz)
 !Define local array:
double precision:: wka(0:ny,0:nxm1,nz)
 !Local variables:
double precision:: qavg0,qsadd
integer:: iz

 !------------------------------------------------------------
 !Reset qs = q = F[qs-qc]+qc+qd, where F is a low pass filter:
wka=qs-qc
call filter(wka,0,2)
qs=wka+qc+qd

 !------------------------------------------------------------
 !Reset qd = q-qc-F[q-qc] and restore domain average for qs:
do iz=1,nz
   call average(qs(0,0,iz),qavg0)
   qsadd=qavg(iz)-qavg0
   qs(:,:,iz)=qs(:,:,iz)+qsadd
enddo
qd=qs-qc
wka=qd
 !Note: qd has a zero average by construction since the domain
 !      average of qc, like qs, is qavg.

call filter(wka,0,2)
qd=qd-wka

return
end subroutine reset

!=======================================================================

subroutine filter(qq,isym,nrep)

! Performs nrep 1-2-1 filters in each direction to return a 
! low pass filtered version of the original array.
! If isym = 0, var is assumed to be symmetric across the boundaries, 
! while if isym = 1, var is assumed to be anti-symmetric.

implicit none

 !Define passed variables:
double precision:: qq(0:ny,0:nxm1,nz)
integer:: isym,nrep
 !Define local array:
double precision:: wka(0:ny,0:nxm1)
 !Other local variables:
integer:: j,ix,iy,iz

!--------------------------------------------------------------------
 !Perform 1-2-1 average nrep times utilising intermediate work array:
if (isym .eq. 0) then
   !The function is symmetric across y boundaries:

   do j=1,nrep
      do iz=1,nz
         do iy=0,ny
            wka(iy,   0)=f12*qq(iy,   0,iz)+f14*(qq(iy,nxm1,iz)+qq(iy,1,iz))
            wka(iy,nxm1)=f12*qq(iy,nxm1,iz)+f14*(qq(iy,nxm2,iz)+qq(iy,0,iz))
         enddo
         do ix=1,nxm2
            do iy=0,ny
               wka(iy,ix)=f12*qq(iy,ix,iz)+f14*(qq(iy,ix-1,iz)+qq(iy,ix+1,iz))
            enddo
         enddo

         do ix=0,nxm1
            qq(0, ix,iz)=f12*(wka(0, ix)+wka(1,   ix))
            do iy=1,nym1
               qq(iy,ix,iz)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
            enddo
            qq(ny,ix,iz)=f12*(wka(ny,ix)+wka(nym1,ix))
         enddo
      enddo
   enddo

else
   !Function is anti-symmetric across y boundaries:

   do j=1,nrep
      do iz=1,nz
         do iy=0,ny
            wka(iy,   0)=f12*qq(iy,   0,iz)+f14*(qq(iy,nxm1,iz)+qq(iy,1,iz))
            wka(iy,nxm1)=f12*qq(iy,nxm1,iz)+f14*(qq(iy,nxm2,iz)+qq(iy,0,iz))
         enddo
         do ix=1,nxm2
            do iy=0,ny
               wka(iy,ix)=f12*qq(iy,ix,iz)+f14*(qq(iy,ix-1,iz)+qq(iy,ix+1,iz))
            enddo
         enddo

         do ix=0,nxm1
            qq(0, ix,iz)=f12*wka(0, ix)
            do iy=1,nym1
               qq(iy,ix,iz)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
            enddo
            qq(ny,ix,iz)=f12*wka(ny,ix)
         enddo
      enddo
   enddo

endif

return
end subroutine filter

!=======================================================================

subroutine average(var,vavg)

! Computes the average value of a 2D field var and returns 
! the result in vavg

implicit none

 !Passed variables:
double precision:: var(0:ny,0:nxm1),vavg

 !Use trapezoidal rule in both directions:
vavg=dsumi*(f12*sum(var(0,:)+var(ny,:))+sum(var(1:nym1,:)))

return
end subroutine average

!=======================================================================

subroutine l1norm(var,vl1)

! Computes the L1 norm of a 2D field var and returns the result in vl1:
!          vl1 = int_xmin^xmax{int_ymin^ymax{|var| dxdy}}

implicit none

 !Define passed variables:
double precision:: var(0:ny,0:nxm1),vl1

 !Use trapezoidal rule in both directions:
vl1=garea*(f12*sum(abs(var(0,:)+var(ny,:)))+sum(abs(var(1:nym1,:))))
 !Note: garea is the grid box area, glx*gly

return
end subroutine l1norm

!=======================================================================

subroutine l2norm(var,vl2)

! Computes the L2 norm of a 2D field var and returns the result in vl2:
!          vl2 = int_xmin^xmax{int_ymin^ymax{var^2 dxdy}}

implicit none

 !Define passed variables:
double precision:: var(0:ny,0:nxm1),vl2

 !Use trapezoidal rule in both directions:
vl2=garea*(f12*sum(var(0,:)**2+var(ny,:)**2)+sum(var(1:nym1,:)**2))

return
end subroutine l2norm

!=======================================================================

subroutine contint(qq,nqjumps,dq,qqmin,qqmax)

! Computes a contour interval for a field qq from (qq_max-qq_min)/n_qjumps

implicit none

 !Passed variables:
double precision:: qq(0:ny,0:nxm1,nz)
double precision:: dq(nz),qqmin(nz),qqmax(nz)
integer:: nqjumps

 !Local variable:
integer:: iz

!-----------------------------------------------------------
do iz=1,nz
  qqmax(iz)=maxval(qq(:,:,iz))
  qqmin(iz)=minval(qq(:,:,iz))
  dq(iz)=(qqmax(iz)-qqmin(iz))/dble(nqjumps)
enddo

return
end subroutine contint

!=======================================================================

 !Main end module
end module generic
