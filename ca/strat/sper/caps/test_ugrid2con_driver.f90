program test_ugrid2con_driver

use parameters
use constants
use contours
use common
use congen
use timing
implicit none

type contour_t
   integer :: field_level
   integer :: num_nodes
   logical :: is_open
   double precision :: first_x, first_y
   double precision :: last_x,  last_y
   double precision, allocatable :: x(:), y(:)
end type contour_t

integer:: iu_in, iu_out, ios, case_num
integer:: npta_base, npta_work
integer:: na_base, na_work
double precision:: dq, rmsdiff, maxdiff
double precision:: qa_input(0:nyup1,0:nxum1)
double precision:: xa_base(npm), ya_base(npm)
integer:: nextq_base(npm), nextq_work(npm)
integer:: inda_base(nm), npa_base(nm), i1_base(nm), i2_base(nm)
integer:: inda_work(nm), npa_work(nm), i1_work(nm), i2_work(nm)
double precision:: global_maxdiff, global_rmsdiff
integer:: total_cases, cases_with_error
logical:: file_exists, contour_match_ok
type(contour_t), allocatable :: base_contours(:), work_contours(:)
double precision, parameter:: tol = small*100.d0

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

call init_contours
timing_on = .true.
call timing_reset()

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
print *, "Case#  Npts  Nctrs  MaxDiff    RMS Diff   Topology  Status"
print *, "-----  ----  -----  ---------  ---------  --------  --------"

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

   read(iu_out, iostat=ios) na_base
   if (ios /= 0) then
      print *, "ERROR reading output na for case", case_num+1
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

   if (na_base > 0) then
      read(iu_out, iostat=ios) inda_base(1:na_base)
      if (ios /= 0) exit
      read(iu_out, iostat=ios) npa_base(1:na_base)
      if (ios /= 0) exit
      read(iu_out, iostat=ios) i1_base(1:na_base)
      if (ios /= 0) exit
      read(iu_out, iostat=ios) i2_base(1:na_base)
      if (ios /= 0) exit
   endif

   qa = qa_input
   npta = 0
   na = 0
   nextq_work = 0
   call ugrid2con(dq, nextq_work)
   npta_work = npta
   na_work = na

   if (na_work > 0) then
      inda_work(1:na_work) = inda(1:na_work)
      npa_work(1:na_work) = npa(1:na_work)
      i1_work(1:na_work) = i1a(1:na_work)
      i2_work(1:na_work) = i2a(1:na_work)
   endif

   rmsdiff = 0.d0
   maxdiff = 0.d0
   contour_match_ok = .false.

   if (npta_work /= npta_base .or. na_work /= na_base) then
      cases_with_error = cases_with_error + 1
      write(*,'(I5,2X,I8,2X,I8,2X,A,2X,A,2X,A,2X,A)') case_num, npta_work, na_work, "   n/a   ", "   n/a   ", "MISMATCH", "FAIL"
      cycle
   endif

   call build_contour_list(na_base, npta_base, xa_base, ya_base, nextq_base, inda_base, npa_base, i1_base, i2_base, base_contours, contour_match_ok)
   if (.not. contour_match_ok) then
      cases_with_error = cases_with_error + 1
      write(*,'(I5,2X,I8,2X,I8,2X,A,2X,A,2X,A,2X,A)') case_num, npta_base, na_base, "   n/a   ", "   n/a   ", "INVALID ", "FAIL"
      cycle
   endif

   call build_contour_list(na_work, npta_work, xa, ya, nextq_work, inda_work, npa_work, i1_work, i2_work, work_contours, contour_match_ok)
   if (.not. contour_match_ok) then
      cases_with_error = cases_with_error + 1
      call free_contour_list(base_contours)
      write(*,'(I5,2X,I8,2X,I8,2X,A,2X,A,2X,A,2X,A)') case_num, npta_work, na_work, "   n/a   ", "   n/a   ", "INVALID ", "FAIL"
      cycle
   endif

   call compare_contour_lists(base_contours, na_base, work_contours, na_work, tol, contour_match_ok, maxdiff, rmsdiff)

   global_maxdiff = max(global_maxdiff, maxdiff)
   global_rmsdiff = global_rmsdiff + rmsdiff**2

   if (.not. contour_match_ok) then
      cases_with_error = cases_with_error + 1
      write(*,'(I5,2X,I8,2X,I8,2X,E9.2,2X,E9.2,2X,A,2X,A)') case_num, npta_base, na_base, maxdiff, rmsdiff, "MISMATCH", "FAIL"
   else
      write(*,'(I5,2X,I8,2X,I8,2X,E9.2,2X,E9.2,2X,A,2X,A)') case_num, npta_base, na_base, maxdiff, rmsdiff, "OK", "PASS"
   endif

   call free_contour_list(base_contours)
   call free_contour_list(work_contours)

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

call timing_report()

contains

subroutine build_contour_list(na_in, npta_in, xa_in, ya_in, nextq_in, inda_in, npa_in, i1_in, i2_in, clist, ok)
implicit none
integer, intent(in):: na_in, npta_in
double precision, intent(in):: xa_in(npm), ya_in(npm)
integer, intent(in):: nextq_in(npm)
integer, intent(in):: inda_in(nm), npa_in(nm), i1_in(nm), i2_in(nm)
type(contour_t), allocatable, intent(out):: clist(:)
logical, intent(out):: ok

integer:: j, k, np, idx(nprm)
logical:: chain_ok, chain_open

ok = .false.

if (allocated(clist)) deallocate(clist)
if (na_in < 0 .or. na_in > nm) return

allocate(clist(na_in))
if (na_in == 0) then
   ok = .true.
   return
endif

do j = 1, na_in
   if (npa_in(j) <= 0 .or. npa_in(j) > nprm) return
   np = npa_in(j)

   call trace_contour(np, i1_in(j), i2_in(j), nextq_in, npta_in, idx, chain_ok, chain_open)
   if (.not. chain_ok) return

   clist(j)%field_level = inda_in(j)
   clist(j)%num_nodes = np
   clist(j)%is_open = chain_open
   clist(j)%first_x = xa_in(idx(1))
   clist(j)%first_y = ya_in(idx(1))
   clist(j)%last_x = xa_in(idx(np))
   clist(j)%last_y = ya_in(idx(np))

   allocate(clist(j)%x(np), clist(j)%y(np))
   do k = 1, np
      clist(j)%x(k) = xa_in(idx(k))
      clist(j)%y(k) = ya_in(idx(k))
   enddo
enddo

ok = .true.

end subroutine build_contour_list

subroutine compare_contour_lists(base_list, na_b, work_list, na_w, tol_loc, ok, maxdiff_out, rmsdiff_out)
implicit none
integer, intent(in):: na_b, na_w
type(contour_t), intent(in):: base_list(:), work_list(:)
double precision, intent(in):: tol_loc
logical, intent(out):: ok
double precision, intent(out):: maxdiff_out, rmsdiff_out

logical:: used_w(nm), found, pair_ok
integer:: jb, jw, npts_pair, total_pts
double precision:: pair_maxdiff, pair_sumsq, sumsq

ok = .true.
maxdiff_out = 0.d0
rmsdiff_out = 0.d0
sumsq = 0.d0
total_pts = 0
used_w = .false.

if (na_b /= na_w) then
   ok = .false.
   return
endif

if (na_b == 0) return

do jb = 1, na_b
   found = .false.
   do jw = 1, na_w
      if (used_w(jw)) cycle
      call compare_contour_pair(base_list(jb), work_list(jw), tol_loc, pair_ok, pair_maxdiff, pair_sumsq, npts_pair)
      if (pair_ok) then
         used_w(jw) = .true.
         found = .true.
         maxdiff_out = max(maxdiff_out, pair_maxdiff)
         sumsq = sumsq + pair_sumsq
         total_pts = total_pts + npts_pair
         exit
      endif
   enddo
   if (.not. found) then
      ok = .false.
      return
   endif
enddo

if (total_pts > 0) rmsdiff_out = sqrt(sumsq / dble(2*total_pts))

end subroutine compare_contour_lists

subroutine compare_contour_pair(base_c, work_c, tol_loc, pair_ok, pair_maxdiff, pair_sumsq, npts_pair)
implicit none
type(contour_t), intent(in):: base_c, work_c
double precision, intent(in):: tol_loc
logical, intent(out):: pair_ok
double precision, intent(out):: pair_maxdiff, pair_sumsq
integer, intent(out):: npts_pair

integer:: k
double precision:: dx, dy

pair_ok = .false.
pair_maxdiff = 0.d0
pair_sumsq = 0.d0
npts_pair = 0

if (base_c%field_level /= work_c%field_level) return
if (base_c%num_nodes /= work_c%num_nodes) return
if (base_c%is_open .neqv. work_c%is_open) return

if (abs(base_c%first_x - work_c%first_x) > tol_loc) return
if (abs(base_c%first_y - work_c%first_y) > tol_loc) return
if (abs(base_c%last_x  - work_c%last_x ) > tol_loc) return
if (abs(base_c%last_y  - work_c%last_y ) > tol_loc) return

if (.not. allocated(base_c%x) .or. .not. allocated(base_c%y)) return
if (.not. allocated(work_c%x) .or. .not. allocated(work_c%y)) return

npts_pair = base_c%num_nodes

do k = 1, npts_pair
   dx = abs(base_c%x(k) - work_c%x(k))
   dy = abs(base_c%y(k) - work_c%y(k))
   if (dx > tol_loc .or. dy > tol_loc) return
   pair_maxdiff = max(pair_maxdiff, dx, dy)
   pair_sumsq = pair_sumsq + dx*dx + dy*dy
enddo

pair_ok = .true.

end subroutine compare_contour_pair

subroutine trace_contour(np, istart, iend, nextq_arr, npta_arr, idx, ok, is_open)
implicit none
integer, intent(in):: np, istart, iend, npta_arr
integer, intent(in):: nextq_arr(npm)
integer, intent(out):: idx(nprm)
logical, intent(out):: ok, is_open
integer:: k

ok = .false.
is_open = .false.

if (np <= 0 .or. np > nprm) return
if (istart < 1 .or. istart > npta_arr) return
if (iend < 1 .or. iend > npta_arr) return

idx(1) = istart

do k = 2, np
   idx(k) = nextq_arr(idx(k-1))
   if (idx(k) < 1 .or. idx(k) > npta_arr) return
enddo

if (idx(np) /= iend) return

if (nextq_arr(iend) == 0) then
   is_open = .true.
else if (nextq_arr(iend) == istart) then
   is_open = .false.
else
   return
endif

ok = .true.

end subroutine trace_contour

subroutine free_contour_list(clist)
implicit none
type(contour_t), allocatable, intent(inout):: clist(:)
integer:: j

if (.not. allocated(clist)) return

do j = 1, size(clist)
   if (allocated(clist(j)%x)) deallocate(clist(j)%x)
   if (allocated(clist(j)%y)) deallocate(clist(j)%y)
enddo

deallocate(clist)

end subroutine free_contour_list

end program test_ugrid2con_driver
