module generic

! Module contains generic subroutines for casl in a doubly periodic geometry. 
!          *** These should not be modified ***

use constants

implicit none

contains 

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
