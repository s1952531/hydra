module generic

! Module contains generic subroutines for casl in a doubly periodic geometry. 
!          *** These should not be modified ***

use constants

implicit none

contains 

!=======================================================================

subroutine average(qq,qavg)

! Computes the average value of a field qq and returns the result in qavg

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(ny,nx)

qavg=zero
do ix=1,nx
  do iy=1,ny
    qavg=qavg+qq(iy,ix)
  enddo
enddo
qavg=qavg/dble(nx*ny)

return
end subroutine

!=======================================================================

subroutine l1norm(qq,ql1)

! Computes the L1 norm of a field qq and returns the result in ql1:
!          ql1 = int_xmin^xmax{int_ymin^ymax{|qq| dxdy}}

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(ny,nx)

ql1=zero
do ix=1,nx
  do iy=1,ny
    ql1=ql1+abs(qq(iy,ix))
  enddo
enddo

ql1=garea*ql1
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

subroutine l2norm(qq,ql2)

! Computes the L2 norm of a field qq and returns the result in ql2:
!          ql2 = int_xmin^xmax{int_ymin^ymax{qq^2 dxdy}}

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq(ny,nx)

ql2=zero
do ix=1,nx
  do iy=1,ny
    ql2=ql2+qq(iy,ix)**2
  enddo
enddo

ql2=garea*ql2
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

subroutine binorm(qq1,qq2,qqbi)

! Computes the integral of qq1*qq2 over the domain

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define passed arrays:
double precision:: qq1(ny,nx),qq2(ny,nx)

qqbi=zero
do ix=1,nx
  do iy=1,ny
    qqbi=qqbi+qq1(iy,ix)*qq2(iy,ix)
  enddo
enddo

qqbi=garea*qqbi
 !Note: garea is the grid box area, glx*gly

return
end subroutine

!=======================================================================

 !Main end module
end module
