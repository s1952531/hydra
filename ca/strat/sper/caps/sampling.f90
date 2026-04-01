module sampling
! Provides infrastructure for sampling contour-to-grid and ugrid-to-contour
! operations at specified times during the simulation.

use constants

implicit none

 !Toggle baseline logging for con2grid input/output:
logical,parameter:: log_con2grid=.true.
integer,parameter:: n_con2grid_samples=100

 !Toggle baseline logging for ugrid2con input/output:
logical,parameter:: log_ugrid2con=.true.
integer,parameter:: n_ugrid2con_samples=100

contains

!=======================================================================

logical function con2grid_save_time(tnow, tsim)
! Determines if current time tnow corresponds to a con2grid sample time.
! Sample times are evenly spaced over [0, tsim].

implicit none

double precision, intent(in):: tnow, tsim
double precision:: dts, tsamp, tol
integer:: isamp

if (n_con2grid_samples .le. 1) then
  con2grid_save_time=.false.
  return
endif

dts=tsim/dble(n_con2grid_samples-1)
isamp=nint(tnow/dts)
isamp=max(0,min(n_con2grid_samples-1,isamp))
tsamp=dble(isamp)*dts
tol=1.d-10*max(one,abs(tsim))

con2grid_save_time=(abs(tnow-tsamp) .le. tol)

end function

!=======================================================================

logical function ugrid2con_save_time(tnow, tsim)
! Determines if current time tnow corresponds to a ugrid2con sample time.
! Sample times are evenly spaced over [0, tsim].

implicit none

double precision, intent(in):: tnow, tsim
double precision:: dts, tsamp, tol
integer:: isamp

if (n_ugrid2con_samples .le. 1) then
  ugrid2con_save_time=.false.
  return
endif

dts=tsim/dble(n_ugrid2con_samples-1)
isamp=nint(tnow/dts)
isamp=max(0,min(n_ugrid2con_samples-1,isamp))
tsamp=dble(isamp)*dts
tol=1.d-10*max(one,abs(tsim))

ugrid2con_save_time=(abs(tnow-tsamp) .le. tol)

end function

!=======================================================================

subroutine write_con2grid_input(xq,yq,dqin,qavgin,nextq,nptqin,ioptin)

implicit none

 !Passed arrays:
double precision, intent(in):: xq(:),yq(:)
integer, intent(in):: nextq(:)

 !Passed scalars:
double precision, intent(in):: dqin,qavgin
integer, intent(in):: nptqin,ioptin

 !Local:
integer:: iu

open(newunit=iu,file='c2g_inputs.dat',status='unknown',position='append', &
 & action='write',access='stream',form='unformatted')
write(iu) nptqin,dqin,qavgin,ioptin
write(iu) xq(1:nptqin)
write(iu) yq(1:nptqin)
write(iu) nextq(1:nptqin)
close(iu)
end subroutine

!=======================================================================

subroutine write_con2grid_output(qq)

implicit none

 !Passed arrays:
double precision, intent(in):: qq(:,:)

 !Local:
integer:: iu

open(newunit=iu,file='c2g_outputs.dat',status='unknown',position='append', &
 & action='write',access='stream',form='unformatted')
write(iu) qq
close(iu)

end subroutine

!=======================================================================

subroutine write_ugrid2con_input(qa,dq)

implicit none

 !Passed arrays:
double precision, intent(in):: qa(:,:)

 !Passed scalars:
double precision, intent(in):: dq

 !Local:
integer:: iu

open(newunit=iu,file='g2c_inputs.dat',status='unknown',position='append', &
 & action='write',access='stream',form='unformatted')
write(iu) dq
write(iu) qa
close(iu)
end subroutine

!=======================================================================

subroutine write_ugrid2con_output(xa,ya,nextq,npta)

implicit none

 !Passed arrays:
double precision, intent(in):: xa(:),ya(:)
integer, intent(in):: nextq(:)

 !Passed scalars:
integer, intent(in):: npta

 !Local:
integer:: iu

open(newunit=iu,file='g2c_outputs.dat',status='unknown',position='append', &
 & action='write',access='stream',form='unformatted')
write(iu) npta
write(iu) xa(1:npta)
write(iu) ya(1:npta)
write(iu) nextq(1:npta)
close(iu)

end subroutine

!=======================================================================

end module sampling
