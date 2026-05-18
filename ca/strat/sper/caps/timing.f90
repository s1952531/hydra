module timing
! Simple wall-clock timing helpers for ugrid2con profiling.

use omp_lib
implicit none

logical, save :: timing_on = .false.

real(8), save :: l0TotTime = 0.d0  ! overall level loop total
real(8), save :: l1TotTime = 0.d0  ! q-range scan
real(8), save :: l2TotTime = 0.d0  ! x-crossings
real(8), save :: l3TotTime = 0.d0  ! bottom edge crossings
real(8), save :: l4TotTime = 0.d0  ! top edge crossings
real(8), save :: l5TotTime = 0.d0  ! interior y crossings
real(8), save :: l6TotTime = 0.d0  ! open contour assembly
real(8), save :: l7TotTime = 0.d0  ! open contour renoding
real(8), save :: l8TotTime = 0.d0  ! closed contour assembly
real(8), save :: l9TotTime = 0.d0  ! closed contour renoding
real(8), save :: u2cParTotTime = 0.d0      ! ugrid2con parallel region
real(8), save :: u2cMasterTotTime = 0.d0   ! ugrid2con master combine
real(8), save :: u2cRebuildTotTime = 0.d0  ! ugrid2con rebuild section
real(8), save :: g0TotTime = 0.d0  ! getzzsrc overall
real(8), save :: g1TotTime = 0.d0  ! getzzsrc init dzdtf
real(8), save :: g2TotTime = 0.d0  ! getzzsrc contour processing
real(8), save :: g2OpenTotTime = 0.d0  ! getzzsrc open contour coeffs
real(8), save :: g2ClosedTotTime = 0.d0  ! getzzsrc closed contour coeffs
real(8), save :: g2AccumTotTime = 0.d0  ! getzzsrc source accumulation
real(8), save :: g3TotTime = 0.d0  ! getzzsrc edge doubling
real(8), save :: g4TotTime = 0.d0  ! getzzsrc coarsen

contains

subroutine timer_start(tstart)
implicit none
real(8), intent(out) :: tstart
tstart = omp_get_wtime()
end subroutine timer_start

subroutine timer_stop(tstart, tend, ttime, ttot)
implicit none
real(8), intent(in) :: tstart
real(8), intent(out) :: tend, ttime
real(8), intent(inout) :: ttot
tend = omp_get_wtime()
ttime = tend - tstart
ttot = ttot + ttime
end subroutine timer_stop

subroutine timing_reset()
implicit none
l0TotTime = 0.d0
l1TotTime = 0.d0
l2TotTime = 0.d0
l3TotTime = 0.d0
l4TotTime = 0.d0
l5TotTime = 0.d0
l6TotTime = 0.d0
l7TotTime = 0.d0
l8TotTime = 0.d0
l9TotTime = 0.d0
u2cParTotTime = 0.d0
u2cMasterTotTime = 0.d0
u2cRebuildTotTime = 0.d0
g0TotTime = 0.d0
g1TotTime = 0.d0
g2TotTime = 0.d0
g2OpenTotTime = 0.d0
g2ClosedTotTime = 0.d0
g2AccumTotTime = 0.d0
g3TotTime = 0.d0
g4TotTime = 0.d0
end subroutine timing_reset

subroutine timing_report()
implicit none
real(8) :: totalTime
real(8) :: getzzsrcTotalTime

if (.not. timing_on) return

totalTime = u2cParTotTime + u2cRebuildTotTime

if (u2cParTotTime .gt. 0.d0 .or. u2cMasterTotTime .gt. 0.d0 .or. u2cRebuildTotTime .gt. 0.d0) then
  write(*,'(a)') '=========================================='
  write(*,'(a)') 'UGrid2Con timing totals'
  write(*,'(a,f12.6)') '  parallel region total:      ', u2cParTotTime
  write(*,'(a,f12.6)') '    master combine total:     ', u2cMasterTotTime
  write(*,'(a,f12.6)') '  rebuild section total:      ', u2cRebuildTotTime
  write(*,'(a,f12.6)') '  total (excl. nested):       ', totalTime
  write(*,'(a)') '=========================================='
endif

getzzsrcTotalTime = g1TotTime + g2TotTime + g3TotTime + g4TotTime

if (g0TotTime .gt. 0.d0) then
  write(*,'(a)') '=========================================='
  write(*,'(a)') 'GetZZSrc timing totals'
  write(*,'(a,f12.6)') '  OUTER total:               ', g0TotTime
  write(*,'(a,f12.6)') '  init dzdtf total:          ', g1TotTime
  write(*,'(a,f12.6)') '  contour processing total:  ', g2TotTime
  write(*,'(a,f12.6)') '    open coeff total:        ', g2OpenTotTime
  write(*,'(a,f12.6)') '    closed coeff total:      ', g2ClosedTotTime
  write(*,'(a,f12.6)') '    source accum total:      ', g2AccumTotTime
  write(*,'(a,f12.6)') '  edge doubling total:       ', g3TotTime
  write(*,'(a,f12.6)') '  coarsen total:             ', g4TotTime
  write(*,'(a,f12.6)') '  total (excl. nested):      ', getzzsrcTotalTime
  write(*,'(a)') '=========================================='
endif
end subroutine timing_report

end module timing
