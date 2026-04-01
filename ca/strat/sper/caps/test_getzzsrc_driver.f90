program test_getzzsrc_driver
! Driver program to read baseline cases for getzzsrc and verify results
! Compares computed outputs against recorded baselines

use parameters
use constants
use contours
implicit none

integer:: iu_in, iu_out, ios, case_num
integer:: i, j
double precision:: rmsdiff, maxdiff, rel_rms
double precision:: dzdt_computed(0:ny,0:nxm1), dzdt_baseline(0:ny,0:nxm1)
double precision:: global_maxdiff, global_rmsdiff
integer:: total_cases, cases_with_error
integer:: nptb_read, nb_read
logical:: file_exists
character(len=6):: status_str
integer:: threshold_exp

! Initialize
case_num = 0
total_cases = 0
cases_with_error = 0
global_maxdiff = 0.d0
global_rmsdiff = 0.d0

print *, "=========================================="
print *, "Getzzsrc Baseline Verification Driver"
print *, "=========================================="
print *, ""

! Check if baseline files exist
inquire(file='getzzsrc_inputs.dat', exist=file_exists)
if (.not. file_exists) then
   print *, "ERROR: getzzsrc_inputs.dat not found!"
   stop
endif

inquire(file='getzzsrc_outputs.dat', exist=file_exists)
if (.not. file_exists) then
   print *, "ERROR: getzzsrc_outputs.dat not found!"
   stop
endif

! Ensure contour geometry/interpolation tables are initialised
call init_contours

! Open baseline files
open(newunit=iu_in, file='getzzsrc_inputs.dat', form='unformatted', access='stream', &
     action='read', status='old', iostat=ios)
if (ios /= 0) then
   print *, "ERROR opening g2s_inputs.dat: iostat =", ios
   stop
endif

open(newunit=iu_out, file='getzzsrc_outputs.dat', form='unformatted', access='stream', &
     action='read', status='old', iostat=ios)
if (ios /= 0) then
   print *, "ERROR opening g2s_outputs.dat: iostat =", ios
   stop
endif

print *, "Opened baseline files. Beginning verification..."
print *, ""
print *, "Case#  MaxDiff    RMS Diff   Rel RMS   Status"
print *, "-----  ---------  ---------  --------  --------"

! Main loop: read cases and verify
do
   ! Read input metadata
   read(iu_in, iostat=ios) nptb_read, nb_read
   if (ios /= 0) exit  ! End of file or read error
   
   case_num = case_num + 1
   total_cases = total_cases + 1
   
   ! Read input arrays
   read(iu_in, iostat=ios) xb(1:nptb_read)
   if (ios /= 0) then
      print *, "ERROR reading xb for case", case_num
      exit
   endif
   
   read(iu_in, iostat=ios) yb(1:nptb_read)
   if (ios /= 0) then
      print *, "ERROR reading yb for case", case_num
      exit
   endif
   
   read(iu_in, iostat=ios) nextb(1:nptb_read)
   if (ios /= 0) then
      print *, "ERROR reading nextb for case", case_num
      exit
   endif
   
   read(iu_in, iostat=ios) i1b(1:nb_read)
   if (ios /= 0) then
      print *, "ERROR reading i1b for case", case_num
      exit
   endif
   
   read(iu_in, iostat=ios) i2b(1:nb_read)
   if (ios /= 0) then
      print *, "ERROR reading i2b for case", case_num
      exit
   endif
   
   ! Read baseline output
   read(iu_out, iostat=ios) dzdt_baseline
   if (ios /= 0) then
      print *, "ERROR reading baseline dzdt for case", case_num
      exit
   endif
   
   ! Set module variables and call getzzsrc
   nptb = nptb_read
   nb = nb_read
   dzdt_computed = 0.d0
   call getzzsrc(dzdt_computed)
   
   ! Compute differences (all (ny+1) x nx points)
   maxdiff = 0.d0
   rmsdiff = 0.d0
   do i = 0, nxm1
      do j = 0, ny
         maxdiff = max(maxdiff, abs(dzdt_computed(j,i) - dzdt_baseline(j,i)))
         rmsdiff = rmsdiff + (dzdt_computed(j,i) - dzdt_baseline(j,i))**2
      enddo
   enddo
   
   rmsdiff = sqrt(rmsdiff / dble((ny+1)*nx))
   
   ! Compute relative RMS difference
   rel_rms = 0.d0
   do i = 0, nxm1
      do j = 0, ny
         rel_rms = rel_rms + dzdt_baseline(j,i)**2
      enddo
   enddo
   rel_rms = sqrt(rmsdiff**2 / dble((ny+1)*nx)) / sqrt(rel_rms / dble((ny+1)*nx))
   if (isnan(rel_rms)) rel_rms = 0.d0
   
   ! Track global statistics
   global_maxdiff = max(global_maxdiff, maxdiff)
   global_rmsdiff = max(global_rmsdiff, rmsdiff)
   
   ! Determine pass/fail and output
   threshold_exp = 100
   if (maxdiff > threshold_exp * small .or. rmsdiff > threshold_exp * small) then
      status_str = "FAIL  "
      cases_with_error = cases_with_error + 1
   else
      status_str = "PASS  "
   endif
   
   write(*, '(I5,2X,E9.2,2X,E9.2,2X,E8.1,2X,A)') &
      case_num, maxdiff, rmsdiff, rel_rms, status_str

end do

close(iu_in)
close(iu_out)

! Print summary
print *, ""
print *, "=========================================="
print *, "Verification Summary"
print *, "=========================================="
print *, "Total cases processed:", total_cases
print *, "Cases with errors:   ", cases_with_error
print *, "Global max diff:     ", global_maxdiff
print *, "Global max RMS diff: ", global_rmsdiff
print *, ""

if (cases_with_error == 0) then
   print *, "RESULT: ALL TESTS PASSED ✓"
else
   print *, "RESULT: SOME TESTS FAILED ✗"
endif

end program test_getzzsrc_driver
