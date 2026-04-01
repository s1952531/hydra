program test_ugrid2con_driver
! Driver program to verify ugrid2con against recorded baselines.
! Reads ug2c_inputs.dat and ug2c_outputs.dat, re-runs ugrid2con, and compares.

use parameters
use constants
use contours
use common
use congen
implicit none

integer:: iu_in, iu_out, ios, case_num
integer:: i, npta_base, npta_work
double precision:: dq, rmsdiff, maxdiff
double precision:: qa_input(0:nyu,0:nxum1)
double precision:: xa_base(npm), ya_base(npm)
integer:: nextq_base(npm), nextq_work(npm)
double precision:: global_maxdiff, global_rmsdiff
integer:: total_cases, cases_with_error
logical:: file_exists

case_num = 0
total_cases = 0
cases_with_error = 0
global_maxdiff = 0.d0
global_rmsdiff = 0.d0

print *, "=========================================="
print *, "UGrid2Con Baseline Verification Driver"
print *, "=========================================="
print *, ""

inquire(file='ug2c_inputs.dat', exist=file_exists)
if (.not. file_exists) then
   print *, "ERROR: ug2c_inputs.dat not found!"
   stop 1
endif

inquire(file='ug2c_outputs.dat', exist=file_exists)
if (.not. file_exists) then
   print *, "ERROR: ug2c_outputs.dat not found!"
   stop 1
endif

! Ensure contour geometry/interpolation tables are initialised.
call init_contours

! Keep ugrid2con sampling disabled during test runs (t outside [0,tsim]).
t = -1.d0

open(newunit=iu_in, file='ug2c_inputs.dat', form='unformatted', access='stream', &
     action='read', status='old', iostat=ios)
if (ios /= 0) then
   print *, "ERROR opening ug2c_inputs.dat: iostat =", ios
   stop 1
endif

open(newunit=iu_out, file='ug2c_outputs.dat', form='unformatted', access='stream', &
     action='read', status='old', iostat=ios)
if (ios /= 0) then
   print *, "ERROR opening ug2c_outputs.dat: iostat =", ios
   stop 1
endif

print *, "Opened baseline files. Beginning verification..."
print *, ""
print *, "Case#  Npts  MaxDiff    RMS Diff   Topology  Status"
print *, "-----  ----  ---------  ---------  --------  --------"

do
   read(iu_in, iostat=ios) dq
   if (ios /= 0) exit

   read(iu_in, iostat=ios) qa_input
   if (ios /= 0) then
      print *, "ERROR reading input qa for case", case_num+1
      exit
   endif

   read(iu_out, iostat=ios) npta_base
   if (ios /= 0) then
      print *, "ERROR reading output npta for case", case_num+1
      exit
   endif

   case_num = case_num + 1
   total_cases = total_cases + 1

   if (npta_base > 0) then
      read(iu_out, iostat=ios) xa_base(1:npta_base)
      if (ios /= 0) exit
      read(iu_out, iostat=ios) ya_base(1:npta_base)
      if (ios /= 0) exit
      read(iu_out, iostat=ios) nextq_base(1:npta_base)
      if (ios /= 0) exit
   endif

   qa = qa_input
   nextq_work = 0
   call ugrid2con(dq, nextq_work)
   npta_work = npta

   if (npta_work /= npta_base) then
      cases_with_error = cases_with_error + 1
      write(*,'(I5,2X,I8,2X,A,2X,A,2X,A)') case_num, npta_work, "   n/a   ", "   n/a   ", "MISMATCH", "FAIL"
      cycle
   endif

   rmsdiff = 0.d0
   maxdiff = 0.d0

   if (npta_base > 0) then
      do i=1,npta_base
         rmsdiff = rmsdiff + (xa(i)-xa_base(i))**2 + (ya(i)-ya_base(i))**2
         maxdiff = max(maxdiff, abs(xa(i)-xa_base(i)))
         maxdiff = max(maxdiff, abs(ya(i)-ya_base(i)))
      enddo
      rmsdiff = sqrt(rmsdiff / dble(2*npta_base))
   endif

   global_maxdiff = max(global_maxdiff, maxdiff)
   global_rmsdiff = global_rmsdiff + rmsdiff**2

   if (npta_base > 0) then
      if (any(nextq_work(1:npta_base) /= nextq_base(1:npta_base))) then
         cases_with_error = cases_with_error + 1
         write(*,'(I5,2X,I8,2X,E9.2,2X,E9.2,2X,A,2X,A)') case_num, npta_base, maxdiff, rmsdiff, "MISMATCH", "FAIL"
      else if (maxdiff > small*100.d0 .or. rmsdiff > small*100.d0) then
         cases_with_error = cases_with_error + 1
         write(*,'(I5,2X,I8,2X,E9.2,2X,E9.2,2X,A,2X,A)') case_num, npta_base, maxdiff, rmsdiff, "OK", "FAIL"
      else
         write(*,'(I5,2X,I8,2X,E9.2,2X,E9.2,2X,A,2X,A)') case_num, npta_base, maxdiff, rmsdiff, "OK", "PASS"
      endif
   else
      write(*,'(I5,2X,I8,2X,E9.2,2X,E9.2,2X,A,2X,A)') case_num, npta_base, maxdiff, rmsdiff, "OK", "PASS"
   endif
enddo

print *, ""
print *, "=========================================="
print *, "Summary:"
print *, "  Total cases verified:", total_cases
print *, "  Cases with errors:", cases_with_error
print *, "  Global max difference:", global_maxdiff
if (total_cases > 0) then
   print *, "  Global RMS difference:", sqrt(global_rmsdiff/dble(total_cases))
else
   print *, "  Global RMS difference: n/a"
endif
print *, "=========================================="

if (cases_with_error == 0 .and. total_cases > 0) then
   print *, "SUCCESS: All cases passed verification!"
else if (total_cases == 0) then
   print *, "WARNING: No baseline cases found!"
else
   print *, "FAILURE: Some cases did not pass verification."
endif

close(iu_in)
close(iu_out)

end program test_ugrid2con_driver
