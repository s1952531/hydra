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
real(8), save :: l10TotTime = 0.d0  ! level reorder/repack
real(8), save :: g0TotTime = 0.d0  ! getzzsrc overall
real(8), save :: g1TotTime = 0.d0  ! getzzsrc init dzdtf
real(8), save :: g2TotTime = 0.d0  ! getzzsrc contour processing
real(8), save :: g2OpenTotTime = 0.d0  ! getzzsrc open contour coeffs
real(8), save :: g2ClosedTotTime = 0.d0  ! getzzsrc closed contour coeffs
real(8), save :: g2AccumTotTime = 0.d0  ! getzzsrc source accumulation
real(8), save :: g3TotTime = 0.d0  ! getzzsrc edge doubling
real(8), save :: g4TotTime = 0.d0  ! getzzsrc coarsen
real(8), save :: c0TotTime = 0.d0  ! con2grid overall
real(8), save :: c1TotTime = 0.d0  ! con2grid setup / lower boundary
real(8), save :: c2TotTime = 0.d0  ! con2grid accumulation / sweep
real(8), save :: c3TotTime = 0.d0  ! con2grid APE
real(8), save :: c4TotTime = 0.d0  ! con2grid coarsen / restore average
real(8), save :: c1aTotTime = 0.d0  ! con2grid ixc setup
real(8), save :: c1bTotTime = 0.d0  ! con2grid qjx reset
real(8), save :: c1cTotTime = 0.d0  ! con2grid boundary segment scan
real(8), save :: c1dTotTime = 0.d0  ! con2grid qbot accumulation
real(8), save :: c2aTotTime = 0.d0  ! con2grid qa reset
real(8), save :: c2bTotTime = 0.d0  ! con2grid source accumulation
real(8), save :: c2cTotTime = 0.d0  ! con2grid y sweep
real(8), save :: c4aTotTime = 0.d0  ! con2grid coarsen
real(8), save :: c4bTotTime = 0.d0  ! con2grid average
real(8), save :: c4cTotTime = 0.d0  ! con2grid restore average loop

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
l10TotTime = 0.d0
g0TotTime = 0.d0
g1TotTime = 0.d0
g2TotTime = 0.d0
g2OpenTotTime = 0.d0
g2ClosedTotTime = 0.d0
g2AccumTotTime = 0.d0
g3TotTime = 0.d0
g4TotTime = 0.d0
c0TotTime = 0.d0
c1TotTime = 0.d0
c2TotTime = 0.d0
c3TotTime = 0.d0
c4TotTime = 0.d0
c1aTotTime = 0.d0
c1bTotTime = 0.d0
c1cTotTime = 0.d0
c1dTotTime = 0.d0
c2aTotTime = 0.d0
c2bTotTime = 0.d0
c2cTotTime = 0.d0
c4aTotTime = 0.d0
c4bTotTime = 0.d0
c4cTotTime = 0.d0
end subroutine timing_reset

subroutine timing_report()
implicit none
real(8) :: totalTime
real(8) :: getzzsrcTotalTime
real(8) :: con2gridTotalTime

if (.not. timing_on) return

totalTime = l1TotTime + l2TotTime + l3TotTime + l4TotTime + l5TotTime + l6TotTime + l8TotTime + l10TotTime

if (l0TotTime .gt. 0.d0) then
  write(*,'(a)') '=========================================='
  write(*,'(a)') 'UGrid2Con timing totals'
  write(*,'(a,f12.6)') '  OUTER loop total:           ', l0TotTime
  write(*,'(a,f12.6)') '  q-range scan total:         ', l1TotTime
  write(*,'(a,f12.6)') '  x-crossings total:          ', l2TotTime
  write(*,'(a,f12.6)') '  bottom edge crossings total:', l3TotTime
  write(*,'(a,f12.6)') '  top edge crossings total:   ', l4TotTime
  write(*,'(a,f12.6)') '  interior y crossings total: ', l5TotTime
  write(*,'(a,f12.6)') '  open assembly total:        ', l6TotTime
  write(*,'(a,f12.6)') '    of which renoding:        ', l7TotTime
  write(*,'(a,f12.6)') '  closed assembly total:      ', l8TotTime
  write(*,'(a,f12.6)') '    of which renoding:        ', l9TotTime
  write(*,'(a,f12.6)') '  reorder/repack total:       ', l10TotTime
  write(*,'(a,f12.6)') '  total (excl. nested):       ', totalTime
  write(*,'(a)') '=========================================='
endif

getzzsrcTotalTime = g1TotTime + g2TotTime + g3TotTime + g4TotTime
con2gridTotalTime = c1TotTime + c2TotTime + c3TotTime + c4TotTime

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

if (c0TotTime .gt. 0.d0) then
  write(*,'(a)') '==========================================='
  write(*,'(a)') 'Con2Grid timing totals'
  write(*,'(a,f12.6)') '  OUTER total:               ', c0TotTime
  write(*,'(a,f12.6)') '  setup / lower boundary (c1): ', c1TotTime
  write(*,'(a,f12.6)') '    ixc setup (c1a):         ', c1aTotTime
  write(*,'(a,f12.6)') '    qjx reset (c1b):         ', c1bTotTime
  write(*,'(a,f12.6)') '    boundary scan (c1c):     ', c1cTotTime
  write(*,'(a,f12.6)') '    qbot accumulation (c1d): ', c1dTotTime
  write(*,'(a,f12.6)') '  accumulation / sweep (c2): ', c2TotTime
  write(*,'(a,f12.6)') '    qa reset (c2a):          ', c2aTotTime
  write(*,'(a,f12.6)') '    source accumulation (c2b): ', c2bTotTime
  write(*,'(a,f12.6)') '    y sweep (c2c):           ', c2cTotTime
  write(*,'(a,f12.6)') '  APE total (c3):            ', c3TotTime
  write(*,'(a,f12.6)') '  coarsen / restore avg (c4):', c4TotTime
  write(*,'(a,f12.6)') '    coarsen (c4a):           ', c4aTotTime
  write(*,'(a,f12.6)') '    average (c4b):           ', c4bTotTime
  write(*,'(a,f12.6)') '    restore loop (c4c):      ', c4cTotTime
  write(*,'(a,f12.6)') '  total (excl. nested):      ', con2gridTotalTime
  write(*,'(a)') '==========================================='
endif
end subroutine timing_report

end module timing
