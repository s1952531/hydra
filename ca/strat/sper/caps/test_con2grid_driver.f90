program test_con2grid_driver
! Driver program to read baseline cases for con2grid and verify results
! Compares computed outputs against recorded baselines

use parameters
use constants
use contours
implicit none

integer:: iu_in, iu_out, ios, case_num
integer:: nptq, iopt, i, j
double precision:: dq, qavg, rmsdiff, maxdiff, rel_rms, qq_baseline_val
double precision:: xq_work(npm), yq_work(npm)
double precision:: nextq_int(npm)
integer:: nextq_work(npm)
double precision:: qq_computed(0:ny,0:nxm1), qq_baseline(0:ny,0:nxm1)
double precision:: global_maxdiff, global_rmsdiff
integer:: total_cases, cases_with_error
logical:: file_exists

! Initialize
case_num = 0
total_cases = 0
cases_with_error = 0
global_maxdiff = 0.d0
global_rmsdiff = 0.d0

print *, "=========================================="
print *, "Con2Grid Baseline Verification Driver"
print *, "=========================================="
print *, ""

! Check if baseline files exist
inquire(file='c2g_inputs.dat', exist=file_exists)
if (.not. file_exists) then
   print *, "ERROR: c2g_inputs.dat not found!"
   stop
endif

inquire(file='c2g_outputs.dat', exist=file_exists)
if (.not. file_exists) then
   print *, "ERROR: c2g_outputs.dat not found!"
   stop
endif

! Open baseline files
open(newunit=iu_in, file='c2g_inputs.dat', form='unformatted', access='stream', &
     action='read', status='old', iostat=ios)
if (ios /= 0) then
   print *, "ERROR opening c2g_inputs.dat: iostat =", ios
   stop
endif

open(newunit=iu_out, file='c2g_outputs.dat', form='unformatted', access='stream', &
     action='read', status='old', iostat=ios)
if (ios /= 0) then
   print *, "ERROR opening c2g_outputs.dat: iostat =", ios
   stop
endif

print *, "Opened baseline files. Beginning verification..."
print *, ""
print *, "Case#  Npts  MaxDiff    RMS Diff   Rel RMS   Status"
print *, "-----  ----  ---------  ---------  --------  --------"

! Main loop: read cases and verify
do
   ! Read input metadata
   read(iu_in, iostat=ios) nptq, dq, qavg, iopt
   if (ios /= 0) exit  ! End of file or read error
   
   case_num = case_num + 1
   total_cases = total_cases + 1
   
   ! Read input arrays (only first nptq elements)
   read(iu_in, iostat=ios) xq_work(1:nptq), yq_work(1:nptq), nextq_int(1:nptq)
   if (ios /= 0) then
      print *, "ERROR reading input arrays for case", case_num
      exit
   endif
   
   ! Convert nextq to integer (in case there's numerical type mismatch)
   nextq_work(1:nptq) = nint(nextq_int(1:nptq))
   
   ! Read baseline output
   read(iu_out, iostat=ios) qq_baseline
   if (ios /= 0) then
      print *, "ERROR reading baseline output for case", case_num
      exit
   endif
   
   ! Reset computed field
   qq_computed = 0.d0
   
   ! Call con2grid with baseline inputs
   ! Temporarily disable logging to avoid overwriting baseline files
   call con2grid(qq_computed, xq_work, yq_work, dq, qavg, nextq_work, nptq, iopt)
   
   ! Compute L2 norm and max difference
   rmsdiff = 0.d0
   maxdiff = 0.d0
   
   do i = 0, nxm1
      do j = 0, ny
         rmsdiff = rmsdiff + (qq_computed(j,i) - qq_baseline(j,i))**2
         maxdiff = max(maxdiff, abs(qq_computed(j,i) - qq_baseline(j,i)))
      enddo
   enddo
   
   rmsdiff = sqrt(rmsdiff / dble((nx)*(ny+1)))
   
   ! Compute relative RMS (avoid division by zero)
   rel_rms = 0.d0
   if (sum(abs(qq_baseline)) > small) then
      rel_rms = rmsdiff / (sum(abs(qq_baseline)) / dble((nx)*(ny+1)))
   endif
   
   ! Update global statistics
   global_maxdiff = max(global_maxdiff, maxdiff)
   global_rmsdiff = sqrt(global_rmsdiff**2 + rmsdiff**2)
   
   ! Print case result
   if (maxdiff > small*100.d0 .or. rmsdiff > small*100.d0) then
      cases_with_error = cases_with_error + 1
      write(*, '(I5,2X,I4,2X,E9.2,2X,E9.2,2X,E8.2,2X,A)') &
         case_num, nptq, maxdiff, rmsdiff, rel_rms, "FAIL"
   else
      write(*, '(I5,2X,I4,2X,E9.2,2X,E9.2,2X,E8.2,2X,A)') &
         case_num, nptq, maxdiff, rmsdiff, rel_rms, "PASS"
   endif
   
   if (case_num >= 10000) exit  ! Safety limit to prevent infinite loops
   
enddo

! Summary
print *, ""
print *, "=========================================="
print *, "Summary:"
print *, "  Total cases verified:", total_cases
print *, "  Cases with errors:", cases_with_error
print *, "  Global max difference:", global_maxdiff
print *, "  Global RMS difference:", global_rmsdiff / sqrt(dble(max(total_cases, 1)))
print *, "=========================================="

if (cases_with_error == 0 .and. total_cases > 0) then
   print *, "SUCCESS: All cases passed verification!"
else if (total_cases == 0) then
   print *, "WARNING: No baseline cases found!"
else
   print *, "FAILURE: Some cases did not pass verification."
endif

! Close files
close(iu_in)
close(iu_out)

end program test_con2grid_driver
