module common

 !Import contants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral

 !Define quantities which need to be preserved between recontouring and evolution:

 !Gridded buoyancy & vorticity fields:
double precision:: bb(0:ny,0:nx),zs(0:ny,0:nx),zd(0:ny,0:nx)

 !Time stepping parameters:
double precision:: dtmax,tfin,tgrid

 !Domain area & reference potential energy:
double precision:: domsumi,domarea,eperef

contains 


!=======================================================================

subroutine combine(qq,qc,qs,qd,qavg)

! Combine contour (qc), large scale (qs), and residual (qd) fields 
! into qq, ensuring the domain average of qq = qavg.

implicit none

 !Passed arrays and variable:
double precision:: qq(0:ny,0:nx),qc(0:ny,0:nx),qs(0:ny,0:nx),qd(0:ny,0:nx)
double precision:: qavg
 !Local array and variables:
double precision:: wka(0:ny,0:nx)
double precision:: qavg0,qadd

!----------------------------------------------------------------------
! Define q = F[qs-qc]+qc+qd
! where F is a low pass filter (see subroutine filter)
wka=qs-qc
call filter(wka,0,2)
qq=wka+qc

 !Restore domain average:
call average(qq,qavg0)
qadd=qavg-qavg0
qq=qq+qadd+qd

return
end subroutine
!=======================================================================

subroutine reset(qc,qs,qd,qavg)

! Resets the gridded fields qs & qd and ensures that <qs> = qavg

implicit none

 !Passed arrays and variable:
double precision:: qc(0:ny,0:nx),qs(0:ny,0:nx),qd(0:ny,0:nx)
double precision:: qavg
 !Local array and variables:
double precision:: wka(0:ny,0:nx)
double precision:: qavg0,qsadd

!----------------------------------------------------------------------
 !Reset qs = q = F[qs-qc]+qc+qd, where F is a low pass filter:
wka=qs-qc
call filter(wka,0,2)
qs=wka+qc

 !Restore domain average:
call average(qs,qavg0)
qsadd=qavg-qavg0

!------------------------------------------------------------
 !Reset qd = q-qc-F[q-qc]
qs=qs+qsadd+qd
qd=qs-qc
wka=qd
call filter(wka,0,2)
qd=qd-wka

!Recompute domain average (which evolves because of qd):
call average(qs,qavg)

return
end subroutine

!=======================================================================

subroutine filter(qq,isym,nrep)

! Performs nrep 1-2-1 filters in each direction to return a 
! low pass filtered version of the original array.
! If isym = 0, var is assumed to be symmetric across the boundaries, 
! while if isym = 1, var is assumed to be anti-symmetric.

implicit none

 !Passed array and variables:
double precision:: qq(0:ny,0:nx)
integer:: isym,nrep

 !Local array and variables:
double precision:: wka(0:ny,0:nx)
integer:: ix,iy,j

 !Perform 1-2-1 average nrep times utilising intermediate work array:
if (isym .eq. 0) then
   !The function is symmetric across boundaries:

  do j=1,nrep
    do iy=0,ny
      wka(iy, 0)=f12*(qq(iy, 0)+qq(iy,   1))
      wka(iy,nx)=f12*(qq(iy,nx)+qq(iy,nxm1))
    enddo
    do ix=1,nxm1
      do iy=0,ny
        wka(iy,ix)=f12*qq(iy,ix)+f14*(qq(iy,ix-1)+qq(iy,ix+1))
      enddo
    enddo

    do ix=0,nx
      qq(0, ix)=f12*(wka(0, ix)+wka(1,   ix))
      do iy=1,nym1
        qq(iy,ix)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
      enddo
      qq(ny,ix)=f12*(wka(ny,ix)+wka(nym1,ix))
    enddo
  enddo

else
   !Function is anti-symmetric across boundaries:

  do j=1,nrep
    do iy=0,ny
      wka(iy, 0)=f12*qq(iy, 0)
      wka(iy,nx)=f12*qq(iy,nx)
    enddo
    do ix=1,nxm1
      do iy=0,ny
        wka(iy,ix)=f12*qq(iy,ix)+f14*(qq(iy,ix-1)+qq(iy,ix+1))
      enddo
    enddo

    do ix=0,nx
      qq(0, ix)=f12*wka(0, ix)
      do iy=1,nym1
        qq(iy,ix)=f12*wka(iy,ix)+f14*(wka(iy-1,ix)+wka(iy+1,ix))
      enddo
      qq(ny,ix)=f12*wka(ny,ix)
    enddo
  enddo

endif

return
end subroutine

!=======================================================================

subroutine average(qq,qavg)

! Computes the average value of a field qq and returns the result in qavg

implicit none

 !Passed array and variable:
double precision:: qq(0:ny,0:nx),qavg
 !Local array and variables:
double precision:: ss(0:ny,0:nx)
integer:: ix,iy

!----------------------------------------------------------------------
 !Include conformal factor in integration:
ss=confac*qq

 !Use trapezoidal rule in both directions:
qavg=zero
do ix=1,nxm1
  qavg=qavg+ss(0,ix)+ss(ny,ix)
enddo
do iy=1,nym1
  qavg=qavg+ss(iy,0)+ss(iy,nx)
enddo
qavg=f12*qavg+f14*(ss(0,0)+ss(ny,0)+ss(0,nx)+ss(ny,nx))

do ix=1,nxm1
  do iy=1,nym1
    qavg=qavg+ss(iy,ix)
  enddo
enddo

qavg=qavg*domsumi
 !domsumi = glx*gly/domarea where domarea is the domain area.

return
end subroutine

!=======================================================================

subroutine l1norm(qq,ql1)

! Computes the L1 norm of a field qq and returns the result in ql1:
!          ql1 = int_xmin^xmax{int_ymin^ymax{|qq| dxdy}}

implicit none

 !Passed array and variable:
double precision:: qq(0:ny,0:nx),ql1
 !Local array and variables:
double precision:: ss(0:ny,0:nx)
integer:: ix,iy

!----------------------------------------------------------------------
 !Include conformal factor in integration:
ss=confac*abs(qq)

 !Use trapezoidal rule in both directions:
ql1=zero
do ix=1,nxm1
  ql1=ql1+ss(0,ix)+ss(ny,ix)
enddo
do iy=1,nym1
  ql1=ql1+ss(iy,0)+ss(iy,nx)
enddo
ql1=f12*ql1+f14*(ss(0,0)+ss(ny,0)+ss(0,nx)+ss(ny,nx))

do ix=1,nxm1
  do iy=1,nym1
    ql1=ql1+ss(iy,ix)
  enddo
enddo

ql1=garea*ql1
 !Note: garea is the conformal domain grid box area, glx*gly.

return
end subroutine

!=======================================================================

subroutine l2norm(qq,ql2)

! Computes the L2 norm of a field qq and returns the result in ql2:
!          ql2 = int_xmin^xmax{int_ymin^ymax{qq^2 dxdy}}

implicit none

 !Passed array and variable:
double precision:: qq(0:ny,0:nx),ql2
 !Local array and variables:
double precision:: ss(0:ny,0:nx)
integer:: ix,iy

!----------------------------------------------------------------------
 !Include conformal factor in integration:
ss=confac*qq**2

 !Use trapezoidal rule in both directions:
ql2=zero
do ix=1,nxm1
  ql2=ql2+ss(0,ix)+ss(ny,ix)
enddo
do iy=1,nym1
  ql2=ql2+ss(iy,0)+ss(iy,nx)
enddo
ql2=f12*ql2+f14*(ss(0,0)+ss(ny,0)+ss(0,nx)+ss(ny,nx))

do ix=1,nxm1
  do iy=1,nym1
    ql2=ql2+ss(iy,ix)
  enddo
enddo

ql2=garea*ql2
 !Note: garea is the conformal domain grid box area, glx*gly.

return
end subroutine

!=======================================================================

subroutine binorm(qq1,qq2,qqbi)

! Computes the integral of qq1*qq2 over the domain

implicit none

 !Passed arrays and variable:
double precision:: qq1(0:ny,0:nx),qq2(0:ny,0:nx),qqbi
 !Local array and variables:
double precision:: ss(0:ny,0:nx)
integer:: ix,iy

!----------------------------------------------------------------------
 !Include conformal factor in integration:
ss=confac*qq1*qq2

 !Use trapezoidal rule in both directions:
qqbi=zero
do ix=1,nxm1
  qqbi=qqbi+ss(0,ix)+ss(ny,ix)
enddo
do iy=1,nym1
  qqbi=qqbi+ss(iy,0)+ss(iy,nx)
enddo
qqbi=f12*qqbi+f14*(ss(0,0)+ss(ny,0)+ss(0,nx)+ss(ny,nx))

do ix=1,nxm1
  do iy=1,nym1
    qqbi=qqbi+ss(iy,ix)
  enddo
enddo

qqbi=garea*qqbi
 !Note: garea is the conformal domain grid box area, glx*gly.

return
end subroutine

!=======================================================================

subroutine kinetic(uu,vv,eke)

! Computes the kinetic energy  

implicit none

 !Passed arrays and variable:
double precision:: uu(0:ny,0:nx),vv(0:ny,0:nx),eke
 !Local array and variables:
double precision:: ss(0:ny,0:nx)
integer:: ix,iy

!----------------------------------------------------------------------
 !Include conformal factor in integration:
ss=(uu**2+vv**2)*confac**2

 !Use trapezoidal rule in both directions:
eke=zero
do ix=1,nxm1
  eke=eke+ss(0,ix)+ss(ny,ix)
enddo
do iy=1,nym1
  eke=eke+ss(iy,0)+ss(iy,nx)
enddo
eke=f12*eke+f14*(ss(0,0)+ss(ny,0)+ss(0,nx)+ss(ny,nx))

do ix=1,nxm1
  do iy=1,nym1
    eke=eke+ss(iy,ix)
  enddo
enddo

eke=f12*garea*eke
 !Note: garea is the conformal domain grid box area, glx*gly.

return
end subroutine

!=======================================================================

subroutine potential(qq,epe)

! Computes the potential energy, the area integral of -Y*q where
! q = buoyancy

implicit none

 !Passed array and variable:
double precision:: qq(0:ny,0:nx),epe
 !Local array and variables:
double precision:: ss(0:ny,0:nx)
integer:: ix,iy

!----------------------------------------------------------------------
 !Include conformal factor in integration:
ss=qq*yo*confac

 !Use trapezoidal rule in both directions:
epe=zero
do ix=1,nxm1
  epe=epe+ss(0,ix)+ss(ny,ix)
enddo
do iy=1,nym1
  epe=epe+ss(iy,0)+ss(iy,nx)
enddo
epe=f12*epe+f14*(ss(0,0)+ss(ny,0)+ss(0,nx)+ss(ny,nx))

do ix=1,nxm1
  do iy=1,nym1
    epe=epe+ss(iy,ix)
  enddo
enddo

epe=-garea*epe-eperef
 !Note: garea is the conformal domain grid box area, glx*gly, while
 !      eperef is a reference potential energy (computed in caps.f90).

return
end subroutine

!=======================================================================

subroutine restore(qq,qavg)

! Restores the average of qq to the value qavg

implicit none

 !Passed array and variable:
double precision:: qq(0:ny,0:nx),qavg
 !Local variables:
double precision:: qavg0,qadd

 !Restore average (qavg):
call average(qq,qavg0)

qadd=qavg-qavg0
qq=qq+qadd
 !Now qq has the correct average

return
end subroutine

!=======================================================================

end module
