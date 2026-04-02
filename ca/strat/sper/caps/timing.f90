module timing
! Simple wall-clock timing helpers for ugrid2con profiling.

use omp_lib
implicit none

logical, save :: timing_on = .false.

real(8), save :: l1TotTime = 0.d0  ! q-range scan
real(8), save :: l2TotTime = 0.d0  ! x-crossings
real(8), save :: l3TotTime = 0.d0  ! bottom edge crossings
real(8), save :: l4TotTime = 0.d0  ! top edge crossings
real(8), save :: l5TotTime = 0.d0  ! interior y crossings
real(8), save :: l6TotTime = 0.d0  ! open contour assembly
real(8), save :: l7TotTime = 0.d0  ! open contour renoding
real(8), save :: l8TotTime = 0.d0  ! closed contour assembly
real(8), save :: l9TotTime = 0.d0  ! closed contour renoding

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
 l1TotTime = 0.d0
 l2TotTime = 0.d0
 l3TotTime = 0.d0
 l4TotTime = 0.d0
 l5TotTime = 0.d0
 l6TotTime = 0.d0
 l7TotTime = 0.d0
 l8TotTime = 0.d0
 l9TotTime = 0.d0
end subroutine timing_reset

subroutine timing_report()
implicit none
real(8) :: totalTime

if (.not. timing_on) return

totalTime = l1TotTime + l2TotTime + l3TotTime + l4TotTime + l5TotTime + l6TotTime + l7TotTime + l8TotTime + l9TotTime

write(*,'(a)') '=========================================='
write(*,'(a)') 'UGrid2Con timing totals'
write(*,'(a,f12.6)') '  q-range scan total:         ', l1TotTime
write(*,'(a,f12.6)') '  x-crossings total:          ', l2TotTime
write(*,'(a,f12.6)') '  bottom edge crossings total:', l3TotTime
write(*,'(a,f12.6)') '  top edge crossings total:   ', l4TotTime
write(*,'(a,f12.6)') '  interior y crossings total: ', l5TotTime
write(*,'(a,f12.6)') '  open assembly total:        ', l6TotTime
write(*,'(a,f12.6)') '  open renode total:          ', l7TotTime
write(*,'(a,f12.6)') '  closed assembly total:      ', l8TotTime
write(*,'(a,f12.6)') '  closed renode total:        ', l9TotTime
write(*,'(a,f12.6)') '  overall total: ', totalTime
write(*,'(a)') '=========================================='
end subroutine timing_report

end module timing
