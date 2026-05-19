module sampling
! Provides infrastructure for sampling contour-to-grid and ugrid-to-contour
! operations at specified times during the simulation.

use constants

implicit none

 !Toggle baseline logging for con2grid input/output:
logical,parameter:: log_con2grid=.false.
integer,parameter:: n_con2grid_samples=100

 !Toggle baseline logging for ugrid2con input/output:
logical,parameter:: log_ugrid2con=.false.
integer,parameter:: n_ugrid2con_samples=100

 !Toggle baseline logging for getzzsrc input/output:
logical,parameter:: log_getzzsrc=.false.
integer,parameter:: n_getzzsrc_samples=100

 !Internal sampling cursors to ensure exactly N writes when times advance:
integer,save:: con2grid_next_sample=0
integer,save:: ugrid2con_next_sample=0
integer,save:: getzzsrc_next_sample=0

contains

!=======================================================================

logical function con2grid_save_time(tnow, tsim)
! Determines if current time tnow corresponds to a con2grid sample time.
! Sample times are evenly spaced over [0, tsim], with deterministic
! monotonic triggering so exactly n_con2grid_samples writes are produced
! (provided tnow advances through tsim).

implicit none

double precision, intent(in):: tnow, tsim
double precision:: dts, tsamp, tol

if (n_con2grid_samples .le. 0) then
  con2grid_save_time=.false.
  return
endif

if (n_con2grid_samples .eq. 1) then
  if (con2grid_next_sample .eq. 0) then
    con2grid_save_time=(tnow .ge. 0.d0)
    if (con2grid_save_time) con2grid_next_sample=1
  else
    con2grid_save_time=.false.
  endif
  return
endif

dts=tsim/dble(n_con2grid_samples-1)
tol=1.d-12*max(one,abs(tsim))

if (con2grid_next_sample .ge. n_con2grid_samples) then
  con2grid_save_time=.false.
  return
endif

tsamp=dble(con2grid_next_sample)*dts

if (tnow .ge. tsamp-tol) then
  con2grid_save_time=.true.
  con2grid_next_sample=con2grid_next_sample+1
else
  con2grid_save_time=.false.
endif

end function

!=======================================================================

logical function ugrid2con_save_time(tnow, tsim)
! Determines if current time tnow corresponds to a ugrid2con sample time.
! Sample times are evenly spaced over [0, tsim], with deterministic
! monotonic triggering so exactly n_ugrid2con_samples writes are produced
! (provided tnow advances through tsim).

implicit none

double precision, intent(in):: tnow, tsim
double precision:: dts, tsamp, tol

if (n_ugrid2con_samples .le. 0) then
  ugrid2con_save_time=.false.
  return
endif

if (n_ugrid2con_samples .eq. 1) then
  if (ugrid2con_next_sample .eq. 0) then
    ugrid2con_save_time=(tnow .ge. 0.d0)
    if (ugrid2con_save_time) ugrid2con_next_sample=1
  else
    ugrid2con_save_time=.false.
  endif
  return
endif

dts=tsim/dble(n_ugrid2con_samples-1)
tol=1.d-12*max(one,abs(tsim))

if (ugrid2con_next_sample .ge. n_ugrid2con_samples) then
  ugrid2con_save_time=.false.
  return
endif

tsamp=dble(ugrid2con_next_sample)*dts

if (tnow .ge. tsamp-tol) then
  ugrid2con_save_time=.true.
  ugrid2con_next_sample=ugrid2con_next_sample+1
else
  ugrid2con_save_time=.false.
endif

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
character(len=7):: fstatus

if (con2grid_next_sample .eq. 1) then
  fstatus='replace'
else
  fstatus='unknown'
endif

open(newunit=iu,file='c2g_inputs.dat',status=fstatus,position='append', &
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
character(len=7):: fstatus

if (con2grid_next_sample .eq. 1) then
  fstatus='replace'
else
  fstatus='unknown'
endif

open(newunit=iu,file='c2g_outputs.dat',status=fstatus,position='append', &
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
character(len=7):: fstatus

if (ugrid2con_next_sample .eq. 1) then
  fstatus='replace'
else
  fstatus='unknown'
endif

open(newunit=iu,file='ug2c_inputs.dat',status=fstatus,position='append', &
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
character(len=7):: fstatus

if (ugrid2con_next_sample .eq. 1) then
  fstatus='replace'
else
  fstatus='unknown'
endif

open(newunit=iu,file='ug2c_outputs.dat',status=fstatus,position='append', &
 & action='write',access='stream',form='unformatted')
write(iu) npta
write(iu) xa(1:npta)
write(iu) ya(1:npta)
write(iu) nextq(1:npta)
close(iu)

end subroutine

!=======================================================================

logical function getzzsrc_save_time(tnow, tsim)
! Determines if current time tnow corresponds to a getzzsrc sample time.
! Sample times are evenly spaced over [0, tsim], with deterministic
! monotonic triggering so exactly n_getzzsrc_samples writes are produced
! (provided tnow advances through tsim).

implicit none

double precision, intent(in):: tnow, tsim
double precision:: dts, tsamp, tol

if (n_getzzsrc_samples .le. 0) then
  getzzsrc_save_time=.false.
  return
endif

if (n_getzzsrc_samples .eq. 1) then
  if (getzzsrc_next_sample .eq. 0) then
    getzzsrc_save_time=(tnow .ge. 0.d0)
    if (getzzsrc_save_time) getzzsrc_next_sample=1
  else
    getzzsrc_save_time=.false.
  endif
  return
endif

dts=tsim/dble(n_getzzsrc_samples-1)
tol=1.d-12*max(one,abs(tsim))

if (getzzsrc_next_sample .ge. n_getzzsrc_samples) then
  getzzsrc_save_time=.false.
  return
endif

tsamp=dble(getzzsrc_next_sample)*dts

if (tnow .ge. tsamp-tol) then
  getzzsrc_save_time=.true.
  getzzsrc_next_sample=getzzsrc_next_sample+1
else
  getzzsrc_save_time=.false.
endif

end function

!=======================================================================

subroutine write_getzzsrc_input(xb, yb, nextb, i1b, i2b, nptb, nb, bjumpin)

implicit none

 !Passed arrays:
double precision, intent(in):: xb(:), yb(:)
integer, intent(in):: nextb(:), i1b(:), i2b(:)

 !Passed scalars:
integer, intent(in):: nptb, nb
double precision, intent(in):: bjumpin

 !Local:
integer:: iu
character(len=7):: fstatus

if (getzzsrc_next_sample .eq. 1) then
  fstatus='replace'
else
  fstatus='unknown'
endif

open(newunit=iu,file='getzzsrc_inputs.dat',status=fstatus,position='append', &
 & action='write',access='stream',form='unformatted')
write(iu) nptb, nb, bjumpin
write(iu) xb(1:nptb)
write(iu) yb(1:nptb)
write(iu) nextb(1:nptb)
write(iu) i1b(1:nb)
write(iu) i2b(1:nb)
close(iu)

end subroutine

!=======================================================================

subroutine write_getzzsrc_output(dzdt)

implicit none

 !Passed arrays:
double precision, intent(in):: dzdt(:,:)

 !Local:
integer:: iu
character(len=7):: fstatus

if (getzzsrc_next_sample .eq. 1) then
  fstatus='replace'
else
  fstatus='unknown'
endif

open(newunit=iu,file='getzzsrc_outputs.dat',status=fstatus,position='append', &
 & action='write',access='stream',form='unformatted')
write(iu) dzdt
close(iu)

end subroutine

!=======================================================================

end module sampling
