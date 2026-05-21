module congen
! Converts contours (xq,yq) to gridded values on an ultra-fine grid
! of dimensions mgu*nx x mgu*ny (in a closed rectangular domain),
! optionally adds a residual field (interpolated to the ultra-fine
! grid), and then creates new contours.

! Open contours originating and terminating in a boundary added by
! D G Dritschel on 18 June 2012 @ Moscow

use common
use generic
use timing
use omp_lib

implicit none

double precision:: qa(0:(nyup1 + 1) * (nxum1 + 1) - 1)
double precision:: xa(npm),ya(npm)
integer:: inda(nm),npa(nm),i1a(nm),i2a(nm)
integer:: na,npta

!Marching squares definitions:
integer, parameter:: LEFT=1, RIGHT=2, BOTTOM=3, TOP=4

contains

integer function qa_idx(iy,ix)
  implicit none
  integer, intent(in) :: iy, ix

  qa_idx = iy + (nyup1 + 1) * ix
end function qa_idx

!=====================================================================

subroutine recontour(qq,xq,yq,dq,qavg,nextq,indq,npq,i1q,i2q,nq,nptq,iopt)
! Main routine for recontouring (from D & Ambaum, 1996, QJRMS)

! qq           : a gridded field added to that due to contours if iopt = 1
! xq(i),yq(i)  : location of node i in the domain
! dq           : contour interval
! qavg         : average value of field (computed if nptq=0)
! nextq(i)     : index of the node following node i
!                *** this must be zero for an endpoint on the boundary ***
! indq(j)      : field level (integer) of contour j
! npq(j)       : number of nodes on contour j
! i1q(j)       : beginning node index on contour j
! i2q(j)       : ending node index on contour j
! nq           : number of contours
! nptq         : total number of nodes
! iopt         : if 1, always combine a residual gridded field with
!                that due to any contours for recontouring;
!                if 0, and if nq = 0, qq is assumed to contain the
!                full field to be contoured (this is normal at t = 0).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: qq(0:ny,0:nxm1)
double precision:: xq(npm),yq(npm)
integer:: nextq(npm),indq(nm),npq(nm),i1q(nm),i2q(nm)

 !Local quantities:
logical:: contours,residual

 !-------------------------------------------------------------------
 !First check if there are any contours to convert to gridded values:
contours=nq .gt. 0

 !See if there is a residual field to add:
residual=(.not. contours) .or. (iopt .eq. 1)

 !-----------------------------------------------------------------
 !Counters for total number of nodes and contours:
npta=0
na=0

 !Obtain fine grid field qa:
if (contours) then
   !Convert contours to gridded values (in the array qa):
  call con2ugrid(xq,yq,dq,qavg,nextq,nptq)

  if (residual) then
     !Bi-linear interpolate the residual qq to the fine grid and add to qa:
    do ix=0,nxum1
      ixf=ixfw(ix)
      ix0=ix0w(ix)
      ix1=ix1w(ix)

      do iy=0,nyu
        iyf=iyfw(iy)
        iy0=iy0w(iy)
        iy1=iy1w(iy)

        qa(qa_idx(iy,ix))=qa(qa_idx(iy,ix))+w00(iyf,ixf)*qq(iy0,ix0) &
                         & +w10(iyf,ixf)*qq(iy1,ix0) &
                         & +w01(iyf,ixf)*qq(iy0,ix1) &
                         & +w11(iyf,ixf)*qq(iy1,ix1)

      enddo
    enddo
  endif

else

   !Check if field requires contouring by computing l1 norm of qq:
  call l1norm(qq,qql1)
  if (qql1 .lt. small) then
    qavg=zero
    return
  endif

   !Compute average value of the field (qavg):
  call average(qq,qavg)

   !No contours: interpolate qq (which here contains the full field)
   !to the fine grid as qa:
  do ix=0,nxum1
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    do iy=0,nyu
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)

      qa(qa_idx(iy,ix))=w00(iyf,ixf)*qq(iy0,ix0)+w10(iyf,ixf)*qq(iy1,ix0) &
             & +w01(iyf,ixf)*qq(iy0,ix1)+w11(iyf,ixf)*qq(iy1,ix1)

    enddo
  enddo

endif

 !Generate new contours (xa,ya) from qa array:
call ugrid2con(dq,nextq)

 !Copy arrays back to those in the argument of the subroutine:
do i=1,npta
  xq(i)=xa(i)
  yq(i)=ya(i)
enddo

do j=1,na
  i1q(j)=i1a(j)
  i2q(j)=i2a(j)
  npq(j)=npa(j)
  indq(j)=inda(j)
enddo

nq=na
nptq=npta

return
end subroutine

!==========================================================================

! subroutine ugrid2con(dq,nextq)
! ! Generates contours (xa,ya) from the gridded field qa for the levels
! ! +/-dq/2, +/-3*dq/2, ....

! implicit double precision(a-h,o-z)
! implicit integer(i-n)

!  !Passed array:
! integer:: nextq(npm)

!  !Local parameters and arrays:
! integer,parameter:: ncrm=3*nplm/4
!  !ncrm:  max number of contour crossings of a single contour level
!  !nplm:  max number of nodes in any contour level

! integer,parameter:: nxny=nxu*nyu, koff=nxu*(nyu-1)

! double precision:: ycr(ncrm),xcr(ncrm)
! double precision:: qdx(0:nxu),qdy(0:nyu)
! double precision:: xd(nprm),yd(nprm)
! integer:: isx(0:nxu),isy(0:nyu)
! integer:: kib(ncrm),icre(nm)
! integer:: icrtab(nxny,2)
! integer*1:: noctab(nxny)
! logical:: free(ncrm),keep

!  !initialise constants and arrays:
! dqi=one/dq
! qoff=dq*dble(nlevm)
!  !qoff: should be a large integer multiple of the contour interval, dq.
!  !The multiple should exceed the maximum expected number of contour levels.

!  !--------------------------------------------------------
!  !First get the beginning and ending contour levels:
! qamax=qa(0,0)
! qamin=qa(0,0)
! do ix=0,nxum1
!   do iy=0,nyu
!     qamax=max(qamax,qa(iy,ix))
!     qamin=min(qamin,qa(iy,ix))
!   enddo
! enddo

! levbeg=int((qoff+qamin)*dqi+f12)+1
! levend=int((qoff+qamax)*dqi+f12)

! if (levbeg .le. levend) then
!  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  !Loop over contour levels and process:
! do lev=levbeg,levend
!  !Integer index giving contour level:
! indq=lev-nlevm+(lev-1)/nlevm-1

!  !Counter for total number of grid line crossings:
! ncr=0

!  !Counter for total number of open contours originating in an edge:
! npe=0

!  !Contour level being sought:
! qtmp=(dble(lev)-f12)*dq-qoff

!  !Below, kib = grid box into which the contour (containing ncr) is going
!  !       kob =   "   "  out of "    "     "         "       "    " coming
!  !      [kob -> ncr -> kib:  ncr lies at the boundary between kob & kib]

!  !   *** grid boxes are numbered 1 (lower left) to nxu*nyu (upper right) ***

!  !Initialise number of crossings per box:
! do k=1,nxny
!   noctab(k)=0
! enddo

!  !-----------------------------------------------------------
!  !Find x grid line crossings first:
! do ix=0,nxum1
!   xgt=xgu(ix)

!   do iy=0,nyu
!     qdy(iy)=qa(iy,ix)-qtmp
!     isy(iy)=sign(one,qdy(iy))
!   enddo

!   do iy=0,nyu-1
!     if (isy(iy) .ne. isy(iy+1)) then
!       ncr=ncr+1
!       inc=(1-isy(iy))/2
!       kaa=iy*nxu
!       kib(ncr)=kaa+ibx(ix+inc)
!       kob=kaa+ibx(ix+1-inc)
!       noctab(kob)=noctab(kob)+1
!       icrtab(kob,noctab(kob))=ncr
!       xcr(ncr)=xgt
!       ycr(ncr)=ygu(iy)-glyu*qdy(iy)/(qdy(iy+1)-qdy(iy))
!     endif
!   enddo

! enddo

!  !----------------------------------------------------------
!  !Find y grid line crossings next (edge values are special):
!  !Bottom edge:
! iy=0
! ygt=ygu(iy)

! do ix=0,nxum1
!   qdx(ix)=qa(iy,ix)-qtmp
!   isx(ix)=sign(one,qdx(ix))
! enddo
! qdx(nxu)=qdx(0)
! isx(nxu)=isx(0)

! do ix=0,nxum1
!   if (isx(ix) .ne. isx(ix+1)) then
!     ncr=ncr+1
!     if (isx(ix) .gt. 0) then
!        !A contour comes out of the boundary at this point:
!       kib(ncr)=ix+1
!       npe=npe+1
!       icre(npe)=ncr
!     else
!        !A contour goes into the boundary at this point:
!       kib(ncr)=0
!       kob=ix+1
!       noctab(kob)=noctab(kob)+1
!       icrtab(kob,noctab(kob))=ncr
!     endif
!     ycr(ncr)=ymin
!     xx=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
!     xcr(ncr)=oms*(xx-ellx*dble(int(xx*hlxi)))
!   endif
! enddo

!  !Top edge:
! iy=nyu
! ygt=ygu(iy)

! do ix=0,nxum1
!   qdx(ix)=qa(iy,ix)-qtmp
!   isx(ix)=sign(one,qdx(ix))
! enddo
! qdx(nxu)=qdx(0)
! isx(nxu)=isx(0)

! do ix=0,nxum1
!   if (isx(ix) .ne. isx(ix+1)) then
!     ncr=ncr+1
!     if (isx(ix) .lt. 0) then
!        !A contour comes out of the boundary at this point:
!       kib(ncr)=koff+ix+1
!       npe=npe+1
!       icre(npe)=ncr
!     else
!        !A contour goes into the boundary at this point:
!       kib(ncr)=0
!       kob=koff+ix+1
!       noctab(kob)=noctab(kob)+1
!       icrtab(kob,noctab(kob))=ncr
!     endif
!     ycr(ncr)=ymax
!     xx=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
!     xcr(ncr)=oms*(xx-ellx*dble(int(xx*hlxi)))
!   endif
! enddo
!  !koff = nxu*(nyu-1) above

!  !Interior y = constant grid lines:
! do iy=1,nyu-1
!   ygt=ygu(iy)

!   do ix=0,nxum1
!     qdx(ix)=qa(iy,ix)-qtmp
!     isx(ix)=sign(one,qdx(ix))
!   enddo
!   qdx(nxu)=qdx(0)
!   isx(nxu)=isx(0)

!   do ix=0,nxum1
!     if (isx(ix) .ne. isx(ix+1)) then
!       ncr=ncr+1
!       inc=(1-isx(ix))/2
!       kaa=(iy-1)*nxu+ix+1
!       kib(ncr)=kaa+(1-inc)*nxu
!       kob=kaa+inc*nxu
!       noctab(kob)=noctab(kob)+1
!       icrtab(kob,noctab(kob))=ncr
!       ycr(ncr)=ygt
!       xx=xgu(ix)-glxu*qdx(ix)/(qdx(ix+1)-qdx(ix))
!       xcr(ncr)=oms*(xx-ellx*dble(int(xx*hlxi)))
!     endif
!   enddo

! enddo

! !find which boxes have noctab(k) = 2 hence these are case 5/10
! !for these, determine storage/build order by printing:
! ! determine the left and right x coords and the lower and upper y coords of the box:
! !get the crossing coordinates: xcr(irctab(k, 1)), ycr(irctab(k, 1) and xcr(irctab(k, 2)), ycr(irctab(k, 2))
! !determine what edge the crossings lie on (left, right, bottom or top)
! !Alternatively determine what edge the crossings lie on by looking at kib(irctab(k, 1)) and kib(irctab(k, 2))
! !if kib is +1 then crossing is on right edge, if -1 then crossing is on left edge, if +nxu then crossing is on top edge, if -nxu then crossing is on bottom edge
! !if kib is zero then top if k>koff else bottom

! do k=1,nxny
!   if (noctab(k) .eq. 2) then

!     !determine where the edge of the first crossing (exit1)
!       cr1_boxdiff= kib(icrtab(k,1)) - k !this decides what box it goes out to the exit edge
!     if (kib(icrtab(k,1)) == 0) then

!         if (k .gt. koff) then
!             exit1 = TOP
!             exit2 = BOTTOM
!         else
!             exit1 = BOTTOM
!             exit2 = TOP
!         endif

!     else
!       cr1_boxdiff = kib(icrtab(k,1)) - k
!       select case(int(cr1_boxdiff))
!       case(1)
!           exit1 = RIGHT
!           exit2 = LEFT

!       case(-1)
!           exit1 = LEFT
!           exit2 = RIGHT

!       case(nxu)
!           exit1 = TOP
!           exit2 = BOTTOM

!       case(-nxu)
!           exit1 = BOTTOM
!           exit2 = TOP

!       end select
!     endif

!     ! if (exit1 .eq. LEFT .or. exit1 .eq. RIGHT) then
!     !   write(*,*) 'Case 5 at box=', k
!     ! elseif (exit1 .eq. TOP .or. exit1 .eq. BOTTOM) then
!     !   write(*,*) 'Case 10 at box=', k
!     ! endif
!   endif
! enddo

!  !----------------------------------------------------------------
!  !Now re-build contours:
! do icr=1,ncr
!   free(icr)=.true.
! enddo

!  !First deal with any open contours attached to boundaries:
! if (npe .gt. 0) then
!   do ie=1,npe
!      !A new contour (indexed na) starts here:
!     na=na+1
!     inda(na)=indq
!     ibeg=npta+1
!     i1a(na)=ibeg

!      !The starting node on the contour (coming out of a boundary):
!     icr=icre(ie)

!     ! write(*,*) 'DEBUG open contour start:', 'lev=', lev, 'ie=', ie, 'icr=', icr, &
!     !       'x=', xcr(icr), 'y=', ycr(icr), 'kib=', kib(icr)

!      !First point on the contour:
!     npd=1
!     xd(1)=xcr(icr)
!     yd(1)=ycr(icr)

!      !Find remaining points on the contour:
!     k=kib(icr)     !k is the box the contour is entering (0 if going into a boundary)
!     do while (k .ne. 0)
!       noc=noctab(k)
!        !Use last crossing in this box (noc) as the next node:
!       icrn=icrtab(k,noc)
!        !icrn gives the next point after icr (icrn is leaving box k)
!       noctab(k)=noc-1
!        !noctab is usually zero now except for boxes with a
!        !maximum possible 2 crossings
!       npd=npd+1
!        !Coordinates of new node:
!       xd(npd)=xcr(icrn)
!       yd(npd)=ycr(icrn)
!       ! write(*,*) 'DEBUG open contour coord', 'x=', xcr(icrn), 'y=', ycr(icrn)
!       free(icrn)=.false.
!       k=kib(icrn)

!       !if next box is saddle, print the current box
!       ! if (kib(icrn) .eq. 1977934) then
!       !   write(*,*) 'DEBUG ', 'kib(icr)=', kib(icr), ' (open contour)'
!       ! endif

!       !if current box is saddle, print the next box
!       ! if (kib(icr) .eq. 1977934) then
!       !   write(*,*) 'DEBUG ', 'kib(icrn)=', kib(icrn), ' (open contour)'
!       ! endif
!     enddo

!     ! write(*,*) 'DEBUG open contour before re-node:', 'lev=', lev, 'ie=', ie, &
!     !            'npd=', npd

!      !Re-distribute nodes on this contour 3 times to reduce complexity:
!     keep=.false.
!     do
!       call renode_open(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
!        !Delete contour if deemed too small (see renode_open):
!       if (npa(na) .eq. 0) exit
!       call renode_open(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
!        !Delete contour if deemed too small (see renode_open):
!       if (npd .eq. 0) exit
!       call renode_open(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
!        !Delete contour if deemed too small (see renode_open):
!       if (npa(na) .eq. 0) exit
!        !Contour is big enough to keep:
!       keep=.true.
!       exit
!     enddo

!     if (keep) then
!       npta=npta+npa(na)
!       iend=ibeg+npa(na)-1
!       i2a(na)=iend
!       do i=ibeg,iend-1
!         nextq(i)=i+1
!       enddo
!       nextq(iend)=0
!     else
!       na=na-1
!     endif

!     free(icr)=.false.
!   enddo
! endif

!  !Next deal with remaining closed contours:
! do icr=1,ncr
!   if (free(icr)) then

!      !A new contour (indexed na) starts here:
!     na=na+1
!     inda(na)=indq
!     ibeg=npta+1
!     i1a(na)=ibeg

!      !First point on the contour:
!     npd=1
!     xd(1)=xcr(icr)
!     yd(1)=ycr(icr)

!     ! write(*,*) 'DEBUG closed contour start:', 'lev=', lev, 'icr=', icr, &
!     !   'x=', xcr(icr), 'y=', ycr(icr), 'kib=', kib(icr)

!      !Find remaining points on the contour:
!     k=kib(icr)
!      !k is the box the contour is entering
!     noc=noctab(k)
!      !Use last crossing (noc) in this box (k) as the next node:
!     icrn=icrtab(k,noc)
!      !icrn gives the next point after icr (icrn is leaving box k)

!     ! write(*,*) 'DEBUG ', 'icrn=', icrn, 'icrtab(',k,',', noc,')=', icrtab(k,noc)
!           !if start box is saddle, print the next box
!     ! if (kib(icr) .eq. 1977934) then
!     !   write(*,*) 'DEBUG ', 'kib(icrn)=', kib(icrn), ' (closed contour)'
!     ! endif
!     do while (icrn .ne. icr)
!       noctab(k)=noc-1
!        !noctab is usually zero now except for boxes with a
!        !maximum possible 2 crossings
!       npd=npd+1
!       xd(npd)=xcr(icrn)
!       yd(npd)=ycr(icrn)
!       ! write(*,*) 'DEBUG closed contour coord', 'x=', xcr(icrn), 'y=', ycr(icrn)
!       free(icrn)=.false.

!       !if box is saddle, print the current box
!       ! if (kib(icrn) .eq. 1666646) then
!       !   write(*,*) 'DEBUG ', 'kib(icr)=', kib(icr), 'kib(icrn)=', kib(icrn), ' (closed contour)'
!       ! endif

!       ! if (kib(icr) .eq. 2166330) then
!       !   write(*,*) 'DEBUG ', 'kib(icrn)=', kib(icrn), ' (closed contour)'
!       ! endif

!       k=kib(icrn)
!       noc=noctab(k)
!       icrn=icrtab(k,noc)

!     enddo

!     ! write(*,*) 'DEBUG closed contour before re-node:', 'lev=', lev, 'icr=', icr, &
!     !            'npd=', npd

!      !Re-distribute nodes on this contour 3 times to reduce complexity:
!     keep=.false.
!     do
!       call renode_closed(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
!        !Delete contour if deemed too small (see renode_closed):
!       if (npa(na) .eq. 0) exit
!       call renode_closed(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
!        !Delete contour if deemed too small (see renode_closed):
!       if (npd .eq. 0) exit
!       call renode_closed(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
!        !Delete contour if deemed too small (see renode_closed):
!       if (npa(na) .eq. 0) exit
!        !Contour is big enough to keep:
!       keep=.true.
!       exit
!     enddo

!     if (keep) then
!       npta=npta+npa(na)
!       iend=ibeg+npa(na)-1
!       i2a(na)=iend
!       do i=ibeg,iend-1
!         nextq(i)=i+1
!       enddo
!       nextq(iend)=ibeg
!     else
!       na=na-1
!     endif

!     free(icr)=.false.
!   endif
! enddo

! enddo
!  !End of loop over contour levels
!  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! endif

! return
! end subroutine

! !=======================================================================

subroutine init_marching_squares_lookup(nseg, edge1_lookup, edge2_lookup)
  integer, intent(out):: nseg(0:15)
  integer, intent(out):: edge1_lookup(0:15, 0:1, 2)
  integer, intent(out):: edge2_lookup(0:15, 0:1, 2)

  nseg(0)=0
  edge1_lookup(0,0,1)=0
  edge2_lookup(0,0,1)=0

  nseg(1)=1
  edge1_lookup(1,0,1)=BOTTOM
  edge2_lookup(1,0,1)=LEFT

  nseg(2)=1
  edge1_lookup(2,0,1)=RIGHT
  edge2_lookup(2,0,1)=BOTTOM

  nseg(3)=1
  edge1_lookup(3,0,1)=RIGHT
  edge2_lookup(3,0,1)=LEFT

  nseg(4)=1
  edge1_lookup(4,0,1)=TOP
  edge2_lookup(4,0,1)=RIGHT

  nseg(5)=2
  !not separated
  edge1_lookup(5,0,1)=TOP
  edge2_lookup(5,0,1)=LEFT
  edge1_lookup(5,0,2)=BOTTOM
  edge2_lookup(5,0,2)=RIGHT

  !separated
  edge1_lookup(5,1,1)=BOTTOM
  edge2_lookup(5,1,1)=LEFT
  edge1_lookup(5,1,2)=TOP
  edge2_lookup(5,1,2)=RIGHT

  nseg(6)=1
  edge1_lookup(6,0,1)=TOP
  edge2_lookup(6,0,1)=BOTTOM

  nseg(7)=1
  edge1_lookup(7,0,1)=TOP
  edge2_lookup(7,0,1)=LEFT

  nseg(8)=1
  edge1_lookup(8,0,1)=LEFT
  edge2_lookup(8,0,1)=TOP

  nseg(9)=1
  edge1_lookup(9,0,1)=BOTTOM
  edge2_lookup(9,0,1)=TOP

  nseg(10)=2
  !not separated
  edge1_lookup(10,0,1)=LEFT
  edge2_lookup(10,0,1)=BOTTOM
  edge1_lookup(10,0,2)=RIGHT
  edge2_lookup(10,0,2)=TOP

  !separated
  edge1_lookup(10,1,1)=LEFT
  edge2_lookup(10,1,1)=TOP
  edge1_lookup(10,1,2)=RIGHT
  edge2_lookup(10,1,2)=BOTTOM

  nseg(11)=1
  edge1_lookup(11,0,1)=RIGHT
  edge2_lookup(11,0,1)=TOP

  nseg(12)=1
  edge1_lookup(12,0,1)=LEFT
  edge2_lookup(12,0,1)=RIGHT

  nseg(13)=1
  edge1_lookup(13,0,1)=BOTTOM
  edge2_lookup(13,0,1)=RIGHT

  nseg(14)=1
  edge1_lookup(14,0,1)=LEFT
  edge2_lookup(14,0,1)=BOTTOM

  nseg(15)=0
  edge1_lookup(15,0,1)=0
  edge2_lookup(15,0,1)=0

end subroutine

!=====================================================================
! Helper subroutines for ugrid2con - extracted for profiling
!=====================================================================

subroutine ug2c_init_allocate(kib, icre, ncrm, nm)
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  integer, intent(in) :: ncrm, nm
  integer, allocatable, intent(inout) :: kib(:), icre(:)
  
  !deallocate from previous level
  if (allocated(kib)) deallocate(kib)
  if (allocated(icre)) deallocate(icre)

  allocate(kib(ncrm))!thread_max_ncr))
  allocate(icre(nm))!thread_max_npe))
end subroutine

!=====================================================================

subroutine ug2c_get_box_coords(box_ID, iy, ix, ixp1)
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  integer, intent(in) :: box_ID
  integer, intent(out) :: iy, ix, ixp1
  
  !lower left corner has indices (iy,ix) and is at coords xgu(ix), ygu(iy)
  !boxes are numbered by x-slices, so iy advances fastest
  iy = mod(box_ID - 1, nyu)
  ix = (box_ID - 1) / nyu
  ixp1 = mod(ix + 1, nxu)
end subroutine

!=====================================================================

subroutine ug2c_get_corner_vals(iy, ix, ixp1, ll, ul, ur, lr)
  use common
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  integer, intent(in) :: iy, ix, ixp1
  double precision, intent(out) :: ll, ul, ur, lr
  
  !get corner values and note min and max
  ll = qa(qa_idx(iy, ix))
  ul = qa(qa_idx(iy + 1, ix))
  ur = qa(qa_idx(iy + 1, ixp1))
  lr = qa(qa_idx(iy, ixp1))
end subroutine

!=====================================================================

subroutine ug2c_corner_min_max(ll, ul, ur, lr, minVal, maxVal)
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  double precision, intent(in) :: ll, ul, ur, lr
  double precision, intent(out) :: minVal, maxVal
  
  !get min and max of the corners
  minVal = min(ll, ul, ur, lr)
  maxVal = max(ll, ul, ur, lr)
end subroutine

!=====================================================================

subroutine ug2c_get_ms_case(ll, ul, ur, lr, qtmp, ms_case)
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  double precision, intent(in) :: ll, ul, ur, lr, qtmp
  integer, intent(out) :: ms_case
  
  !get marching squares case
  ms_case = 0
  if (ll >= qtmp) ms_case = ms_case + 1
  if (lr >= qtmp) ms_case = ms_case + 2
  if (ur >= qtmp) ms_case = ms_case + 4
  if (ul >= qtmp) ms_case = ms_case + 8
  !Using >= condition: What happens if >?
end subroutine

!=====================================================================

subroutine ug2c_interp(ll, ul, ur, lr, qtmp, iy, ix, ixp1, &
                       x_b_interp, y_b_interp, x_t_interp, y_t_interp, &
                       x_l_interp, y_l_interp, x_r_interp, y_r_interp)
  use common
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  double precision, intent(in) :: ll, ul, ur, lr, qtmp
  integer, intent(in) :: iy, ix, ixp1
  double precision, intent(out) :: x_b_interp, y_b_interp, x_t_interp, y_t_interp
  double precision, intent(out) :: x_l_interp, y_l_interp, x_r_interp, y_r_interp
  
  double precision :: dz, dz_safe, t_interp
  double precision, parameter :: eps = 1.0e-12
  
  !calculate all interpolation points for this cell and level (some may not be needed but this minimises branching):
  ! bottom: ll -> lr
  dz = lr - ll
  dz_safe = sign(max(abs(dz), eps), dz)
  t_interp = (qtmp - ll) / dz_safe
  x_b_interp = xgu(ix) + t_interp * glxu
  y_b_interp = ygu(iy)
  x_b_interp = oms * (x_b_interp - ellx * dble(int(x_b_interp * hlxi)))
  
  ! top: ul -> ur
  dz = ur - ul
  dz_safe = sign(max(abs(dz), eps), dz)
  t_interp = (qtmp - ul) / dz_safe
  x_t_interp = xgu(ix) + t_interp * glxu
  y_t_interp = ygu(iy + 1)
  x_t_interp = oms * (x_t_interp - ellx * dble(int(x_t_interp * hlxi)))
  
  ! left: ll -> ul
  dz = ul - ll
  dz_safe = sign(max(abs(dz), eps), dz)
  t_interp = (qtmp - ll) / dz_safe
  x_l_interp = xgu(ix)
  y_l_interp = ygu(iy) + t_interp * glyu
  
  ! right: lr -> ur
  dz = ur - lr
  dz_safe = sign(max(abs(dz), eps), dz)
  t_interp = (qtmp - lr) / dz_safe
  x_r_interp = xgu(ixp1)
  y_r_interp = ygu(iy) + t_interp * glyu
end subroutine

!=====================================================================

subroutine ug2c_ms_edge_lookup(edge1_lookup, edge2_lookup, ms_case, separated, seg, &
                               edge1, edge2)
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  integer, intent(in) :: edge1_lookup(0:15, 0:1, 2), edge2_lookup(0:15, 0:1, 2)
  integer, intent(in) :: ms_case, separated, seg
  integer, intent(out) :: edge1, edge2
  
  edge1 = edge1_lookup(ms_case, separated, seg)
  edge2 = edge2_lookup(ms_case, separated, seg)
end subroutine

!=====================================================================

subroutine ug2c_store_crossing_coords(edge1, edge2, &
                                      x_b_interp, y_b_interp, &
                                      x_t_interp, y_t_interp, &
                                      x_l_interp, y_l_interp, &
                                      x_r_interp, y_r_interp, &
                                      x1, y1, x2, y2)
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  integer, parameter :: BOTTOM=3, TOP=4, LEFT=1, RIGHT=2
  integer, intent(in) :: edge1, edge2
  double precision, intent(in) :: x_b_interp, y_b_interp, x_t_interp, y_t_interp
  double precision, intent(in) :: x_l_interp, y_l_interp, x_r_interp, y_r_interp
  double precision, intent(out) :: x1, y1, x2, y2
  
  ! Get coordinates for edge1
  select case(edge1)
    case(BOTTOM)
      x1 = x_b_interp
      y1 = y_b_interp
    case(TOP)
      x1 = x_t_interp
      y1 = y_t_interp
    case(LEFT)
      x1 = x_l_interp
      y1 = y_l_interp
    case(RIGHT)
      x1 = x_r_interp
      y1 = y_r_interp
  end select
  
  ! Get coordinates for edge2
  select case(edge2)
    case(BOTTOM)
      x2 = x_b_interp
      y2 = y_b_interp
    case(TOP)
      x2 = x_t_interp
      y2 = y_t_interp
    case(LEFT)
      x2 = x_l_interp
      y2 = y_l_interp
    case(RIGHT)
      x2 = x_r_interp
      y2 = y_r_interp
  end select
end subroutine

!=====================================================================

subroutine ug2c_boundary_edge_storage(edge1, x1, y1, box_ID, iy, ix, &
                                      ncr, npe, xcr, ycr, kib, icre, &
                                      noctab, icrtab_local_ncr, icrtab_threadID, thread_id)
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  integer, intent(in) :: edge1, box_ID, iy, ix, thread_id
  double precision, intent(in) :: x1, y1
  integer, intent(inout) :: ncr, npe
  double precision, intent(inout) :: xcr(:), ycr(:)
  integer, intent(inout) :: kib(:), icre(:)
  integer*1, intent(inout) :: noctab(:)
  integer, intent(inout) :: icrtab_local_ncr(:,:), icrtab_threadID(:,:)
  
  integer :: kob
  
  !if box on top row (iy = nyu-1) need to store entry at top
  if (edge1 .eq. TOP .and. iy == nyu - 1) then
    ncr = ncr + 1
    xcr(ncr) = x1
    ycr(ncr) = y1
    kib(ncr) = box_ID
    npe = npe + 1
    icre(npe) = ncr
  endif
  
  !if box on bottom row (iy = 0) need to store entry at bottom
  if (edge1 .eq. BOTTOM .and. iy == 0) then
    ncr = ncr + 1
    xcr(ncr) = x1
    ycr(ncr) = y1
    kib(ncr) = box_ID
    npe = npe + 1
    icre(npe) = ncr
  endif
  
  !if box on left column (ix = 0) need to store entry at left
  if (edge1 .eq. LEFT .and. ix == 0) then
    ncr = ncr + 1
    kob = box_ID + (nxu - 1) * nyu !contour is wrapping around from right to left
    xcr(ncr) = x1
    ycr(ncr) = y1
    kib(ncr) = box_ID
    
    !$OMP CRITICAL
    noctab(kob) = noctab(kob) + 1
    icrtab_local_ncr(kob, noctab(kob)) = ncr
    icrtab_threadID(kob, noctab(kob)) = thread_id
    !$OMP END CRITICAL
  endif
  
  !if box on right column (ix = nxu-1) need to store entry at right
  if (edge1 .eq. RIGHT .and. ix == nxu - 1) then
    ncr = ncr + 1
    kob = box_ID - (nxu - 1) * nyu !contour is wrapping around from left to right
    xcr(ncr) = x1
    ycr(ncr) = y1
    kib(ncr) = box_ID
    
    !$OMP CRITICAL
    noctab(kob) = noctab(kob) + 1
    icrtab_local_ncr(kob, noctab(kob)) = ncr
    icrtab_threadID(kob, noctab(kob)) = thread_id
    !$OMP END CRITICAL
  endif
end subroutine

!=====================================================================

subroutine ug2c_calc_edge2_kib(edge2, box_ID, iy, ix, kib_val)
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  integer, intent(in) :: edge2, box_ID, iy, ix
  integer, intent(out) :: kib_val
  
  select case(edge2)
    case(BOTTOM)
      kib_val = box_ID - 1
      if (iy == 0) then
        kib_val = 0 !contour is going into the boundary
      endif
    case(TOP)
      kib_val = box_ID + 1
      if (iy == nyu - 1) then
        kib_val = 0 !contour is going into the boundary
      endif
    case(LEFT)
      kib_val = box_ID - nyu
      if (ix == 0) then
        kib_val = box_ID + (nxu - 1) * nyu !contour is wrapping around from left to right
      endif
    case(RIGHT)
      kib_val = box_ID + nyu
      if (ix == nxu - 1) then
        kib_val = box_ID - (nxu - 1) * nyu !contour is wrapping around from right to left
      endif
  end select
end subroutine

!=====================================================================

subroutine ug2c_store_crossing_and_connectivity(x2, y2, ncr, box_ID, ix, &
                                                xcr, ycr, kob, noctab, &
                                                icrtab_local_ncr, icrtab_threadID, &
                                                thread_id)
  implicit double precision(a-h,o-z)
  implicit integer(i-n)
  
  double precision, intent(in) :: x2, y2
  integer, intent(in) :: ncr, box_ID, ix, thread_id
  double precision, intent(inout) :: xcr(:), ycr(:)
  integer, intent(inout) :: kob
  integer*1, intent(inout) :: noctab(:)
  integer, intent(inout) :: icrtab_local_ncr(:,:), icrtab_threadID(:,:)
  
  xcr(ncr) = x2
  ycr(ncr) = y2
  kob = box_ID
  
  !contention is on left and right boundary boxes so if box being processed is there do critical else parallel update of noctab and icrtab
  if (ix == 0 .or. ix == nxu - 1) then
    !$OMP CRITICAL
    noctab(kob) = noctab(kob) + 1

    !original store/build order:
    !icrtab(kob,noctab(kob))=ncr
    icrtab_local_ncr(kob, noctab(kob)) = ncr
    icrtab_threadID(kob, noctab(kob)) = thread_id
    !$OMP END CRITICAL
  else
    noctab(kob) = noctab(kob) + 1

    !original store/build order:
    !icrtab(kob,noctab(kob))=ncr
    icrtab_local_ncr(kob, noctab(kob)) = ncr
    icrtab_threadID(kob, noctab(kob)) = thread_id
  endif
end subroutine

!=====================================================================

subroutine ugrid2con(dq,nextq)
!subroutine marching_squares_ug2c(dq,nextq)
  !original ug2c definitions
  implicit double precision(a-h,o-z)
  implicit integer(i-n)

  !Passed array:
  integer:: nextq(npm)

  !Local parameters and arrays:
  integer,parameter:: ncrm=3*nplm/4
  !ncrm:  max number of contour crossings of a single contour level
  !nplm:  max number of nodes in any contour level

  integer,parameter:: nxny=nxu*nyu, koff=nxu*(nyu-1)

  double precision:: ycr_global(ncrm),xcr_global(ncrm)
  double precision:: qdx(0:nxu),qdy(0:nyu)
  double precision:: xd(nprm),yd(nprm)
  integer:: isx(0:nxu),isy(0:nyu)
  ! integer:: kib(ncrm),icre(nm)
  integer:: icrtab(nxny,2)
  integer*1:: noctab(nxny)
  logical:: free(ncrm),keep

  !end original ug2c definitions
  integer:: icrtab_global(nxny,2)
  integer:: icrtab_threadID(nxny,2)
  integer:: icrtab_local_ncr(nxny,2)
  integer, allocatable:: ncr_offset(:)
  !do i need an npe offset too?

  integer:: kib_global(ncrm),icre_global(nm)
  !integer:: icrtab_global(nxny,2)
  !integer*1:: noctab_global(nxny)

  integer:: ixp1

  !lookup tables
  integer :: nseg(0:15) !the number of segments in each case
  integer :: edge1_lookup(0:15, 0:1, 2) !(case, variant, segment)
  integer :: edge2_lookup(0:15, 0:1, 2) !(case, variant, segment)
  !variant is for the disambiguation of case 5/10's saddle point: separated or not

  double precision:: ll, ul, ur, lr  !corner values
  double precision:: minVal, maxVal
  integer:: lev, levbeg, levend, box_ID, seg, ms_case, separated, edge1, edge2
  integer:: kob !loop and lookup variables

  !output helpers for qa-qtmp write
  double precision :: qout(0:nyup1,0:nxum1)
  character(len=200) :: fname
  character(len=200) :: outdir, mkdir_cmd
  integer :: outunit, ios
  integer :: s_dim, t_dim
  integer :: mkdir_stat
  integer, save :: sample_counter = 0

  !the coordinates of the crossing point
  double precision:: x_b_interp, y_b_interp !at bottom edge
  double precision:: x_t_interp, y_t_interp !at top edge
  double precision:: x_l_interp, y_l_interp !at left edge
  double precision:: x_r_interp, y_r_interp !at right edge

  !interpolation definitions
  double precision :: dz, dz_safe, t_interp
  double precision, parameter :: eps = 1.0e-12

  !omp variables
  integer :: thread_id, num_threads, thread_max_ncr, thread_max_npe
  integer, allocatable:: ncr_thread_list(:)  !thread_id
  integer, allocatable:: npe_thread_list(:)  !thread_id
  integer, allocatable:: icre_thread_list(:, :)  !npe, thread_id
  double precision, allocatable:: xcr_thread_list(:, :)   !ncr, thread_id
  double precision, allocatable:: ycr_thread_list(:, :)   !ncr, thread_id
  integer, allocatable:: kib_thread_list(:, :)   !ncr, thread_id

  !thread local versions of global arrays limit large memory copies
  integer:: ncr, npe
  !integer, allocatable:: noctab(:), icrtab(:,:)
  integer, allocatable:: kib(:),icre(:)
  double precision:: ycr(ncrm),xcr(ncrm)
  double precision:: u2cParStart,u2cParEnd,u2cParTime
  double precision:: u2cMasterStart,u2cMasterEnd,u2cMasterTime
  double precision:: u2cRebuildStart,u2cRebuildEnd,u2cRebuildTime
  !kib, icre

  !get omp_get_max_threads() and allocate
  num_threads = omp_get_max_threads()
  if (num_threads .lt. 1) then
    num_threads = 1
  endif

  ! write(*,*) 'DEBUG omp_get_max_threads()=', omp_get_max_threads()

  allocate(ncr_offset(num_threads))

  !worst case there are 2 crossings per box
  thread_max_ncr = ncrm!/num_threads
  thread_max_npe = nm!/num_threads

  allocate(ncr_thread_list(num_threads))
  allocate(npe_thread_list(num_threads))
  allocate(icre_thread_list(thread_max_npe, num_threads))
  allocate(xcr_thread_list(thread_max_ncr, num_threads))
  allocate(ycr_thread_list(thread_max_ncr, num_threads))
  allocate(kib_thread_list(thread_max_ncr, num_threads))

  !Saddle point ambiguity terminology:
  !marked = point is above the contour level
  !separated = two edges common to a single marked point are connected

  !--------------------------------------------------------
  !initialise constants and arrays:
  dqi=one/dq
  qoff=dq*dble(nlevm)
  !qoff: should be a large integer multiple of the contour interval, dq.
  !The multiple should exceed the maximum expected number of contour levels.

  !First get the beginning and ending contour levels:
  qamax=qa(qa_idx(0,0))
  qamin=qa(qa_idx(0,0))
  do ix=0,nxum1
    do iy=0,nyu
      qamax=max(qamax,qa(qa_idx(iy,ix)))
      qamin=min(qamin,qa(qa_idx(iy,ix)))
    enddo
  enddo

  levbeg=int((qoff+qamin)*dqi+f12)+1
  levend=int((qoff+qamax)*dqi+f12)

  !first call the subroutine to initialise the marching squares lookup tables:
  call init_marching_squares_lookup(nseg, edge1_lookup, edge2_lookup)

  !--------------------------------------------------------

  !TODO:
  !loop over tiles
  !loop over cells in tile
  !loop over levels
  !rebuild in parallel (c.f. segment joining and fragment joining)
  !rebuild inside the level loop avoids having to store all crossings for all levels

  !regardless of loop order,
  !the shared write to storage (ncr, icrtab, noctab) probably causes lots idle threads

  !currently looping over levels then cells since reconstruction is serial over field levels for now
  !this loses the cache benefit of loading a chunk of corner values
  !and calculating the cell max/mins only once per level

  !------ code to store qdiff for debugging -----
  ! Create one output directory per call (sample): Sample_X_qdiff
  ! sample_counter = sample_counter + 1
  ! write(outdir,'(A,I0,A)') 'Sample_', sample_counter, '_qdiff'
  ! mkdir_cmd = 'mkdir -p ' // trim(outdir)
  ! call execute_command_line(trim(mkdir_cmd), wait=.true., exitstat=mkdir_stat)
  ! if (mkdir_stat .ne. 0) then
  !   write(*,*) 'ERROR creating directory ', trim(outdir), ' exitstat=', mkdir_stat
  ! endif
  !----------------------------------------------

  if (levbeg .le. levend) then
  do lev=levbeg, levend

    ! ncr=0 !Counter for total number of grid line crossings
    ! npe=0 !Counter for total number of open contours originating in an edge

    !Initialise number of crossings per box:
    do k=1,nxny
      noctab(k)=0
    enddo

    qtmp=(dble(lev)-f12)*dq-qoff !Contour level being sought
    indq=lev-nlevm+(lev-1)/nlevm-1  !Integer index giving contour level

    !------ code to store qdiff for debugging -----
    ! !--- Write raw binary of qa - qtmp for this level (stream unformatted)
    ! qout = qa - qtmp
    ! s_dim = size(qout,1)
    ! t_dim = size(qout,2)
    ! write(fname,'(A,A,I4.4,A)') trim(outdir), '/qa_minus_qtmp_lev', lev, '.bin'
    ! open(newunit=outunit, file=fname, access='stream', form='unformatted', status='replace', action='write', iostat=ios)
    ! if (ios .ne. 0) then
    !   write(*,*) 'ERROR opening file', trim(fname), 'iostat=', ios
    ! else
    !   write(outunit) s_dim, t_dim
    !   write(outunit) qout
    !   close(outunit)
    ! endif
    ! !--- end write

    ! write(*,*) 'DEBUG level start:', 'lev=', lev, 'qtmp=', qtmp, 'indq=', indq

    !zero icrtab from previous level
    icrtab_global=0
    icrtab_threadID=0
    icrtab_local_ncr=0

    !zero per-thread counters and thread-local storage from previous level
    if (allocated(ncr_thread_list)) then
      ncr_thread_list = 0
    endif
    if (allocated(npe_thread_list)) then
      npe_thread_list = 0
    endif
    if (allocated(ncr_offset)) then
      ncr_offset = 0
    endif
    if (allocated(icre_thread_list)) then
      icre_thread_list = 0
    endif
    if (allocated(kib_thread_list)) then
      kib_thread_list = 0
    endif
    if (allocated(xcr_thread_list)) then
      xcr_thread_list = 0.0d0
    endif
    if (allocated(ycr_thread_list)) then
      ycr_thread_list = 0.0d0
    endif

    if (timing_on) call timer_start(u2cParStart)

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(box_ID, iy, ix, ixp1, &
    !$OMP         ll, ul, ur, lr, minVal, maxVal, ms_case, separated, &
    !$OMP         dz, dz_safe, t_interp, &
    !$OMP         x_b_interp, y_b_interp, &
    !$OMP         x_t_interp, y_t_interp, &
    !$OMP         x_l_interp, y_l_interp, &
    !$OMP         x_r_interp, y_r_interp, &
    !$OMP         ncr, npe, &
    !$OMP         kob, kib, icre, ycr, xcr, &
    !$OMP         seg, edge1, edge2, x1, y1, x2, y2, &
    !$OMP         thread_id) &
    !$OMP SHARED(qtmp, qa, xgu, ygu, &
    !$OMP        nseg, edge1_lookup, edge2_lookup, &
    !$OMP        kib_thread_list,  kib_global, &
    !$OMP        icre_thread_list, icre_global, &
    !$OMP        ycr_thread_list,  ycr_global, &
    !$OMP        xcr_thread_list,  xcr_global, &
    !$OMP        ncr_thread_list,  ncr_global, &
    !$OMP        npe_thread_list,  npe_global, &
    !$OMP        noctab, &
    !$OMP        icrtab_global, icrtab_threadID, icrtab_local_ncr, &
    !$OMP        num_threads, thread_max_ncr, thread_max_npe, ncr_offset, &
    !$OMP        u2cMasterStart,u2cMasterEnd,u2cMasterTime,u2cMasterTotTime, timing_on)

    !allocate(noctab(nxny))!/num_threads)) !if indexed using box_ID then would be out of range for nxny/num_threads
    !allocate(icrtab(nxny), 2)!/num_threads,2))

    !---- begin init_allocate()
    call ug2c_init_allocate(kib, icre, ncrm, nm)
    !---- end init_allocate()

    ncr = 0
    npe = 0

    thread_id = omp_get_thread_num() + 1 !for 1-based indexing of thread_id

    !$OMP DO SCHEDULE(static, nyu)
    do box_ID=1,nxny !grid boxes are numbered 1 (lower left) to nxu*nyu (upper right)

      !---- begin get_box_coords()
      call ug2c_get_box_coords(box_ID, iy, ix, ixp1)
      !---- end get_box_coords()

      !---- begin get_corner_vals()
      call ug2c_get_corner_vals(iy, ix, ixp1, ll, ul, ur, lr)
      !---- end get_corner_vals()

      !---- begin corner_min_max()
      call ug2c_corner_min_max(ll, ul, ur, lr, minVal, maxVal)
      !---- end corner_min_max()

      !loop over levels here when parallel rebuilding (per tile) is implemented

      !could skip this level early if it does not cross the cell
      !should this use lt, gt or le, ge?
      ! if ((qtmp .le. minVal) .or. (qtmp .ge. maxVal)) cycle

      !---- begin get_ms_case()
      call ug2c_get_ms_case(ll, ul, ur, lr, qtmp, ms_case)
      !---- end get_ms_case()

      !if ms_case is 0 or 15 then there are no crossings i.e. nseg=0 so the loop is skipped and we move to the next cell.
      if (ms_case == 0 .or. ms_case == 15) cycle

      ! !$OMP CRITICAL

      !disambiguate saddle cases (ms_case = 5 or 10) by using asymptotic decider
      !note, edge-based re-build cannot actually utilise this at the moment hence commented out
      separated=0 !reset this for the cell
      ! if (ms_case == 5 .or. ms_case == 10) then
      !   ! separated=1
      !   discriminant = (ul-qtmp)*(lr-qtmp) - (ll-qtmp)*(ur-qtmp)
      !   if (discriminant > 0) then !separated
      !       separated=1
      !   endif
      ! endif

      !---- begin interp()
      call ug2c_interp(ll, ul, ur, lr, qtmp, iy, ix, ixp1, &
                       x_b_interp, y_b_interp, x_t_interp, y_t_interp, &
                       x_l_interp, y_l_interp, x_r_interp, y_r_interp)
      !---- end interp()

      !lookup case to get edges that are crossed
      do seg=1,nseg(ms_case)
        !---- begin ms_edge_lookup()
        call ug2c_ms_edge_lookup(edge1_lookup, edge2_lookup, ms_case, separated, seg, &
                                 edge1, edge2)
        !---- end ms_edge_lookup()

        !get the coordinates of the crossing points for edge1 and edge2
        !---- begin store_crossing_coords()
        call ug2c_store_crossing_coords(edge1, edge2, &
                                        x_b_interp, y_b_interp, &
                                        x_t_interp, y_t_interp, &
                                        x_l_interp, y_l_interp, &
                                        x_r_interp, y_r_interp, &
                                        x1, y1, x2, y2)
        !---- end store_crossing_coords()

        !MS stores, for each box, the enter and exit crossing points
        !This gives duplicate crossings for shared edges
        !to store only unique values there are three options
        !1. store only the entry (/exit) point except for boundary edges where exit (/entry) may determine contour type (open or closed)
        !2. store all points but only on left+bottom (/right+top) except for boundary edges where non-saved edges may determine contour type (open or closed)
        !3. store all and filter duplicates in a post-processing step (likely messier and more expensive)

        !if on bottom row (boxID <= nxu) definitely store the bottom edge crossing point
        !if on the top row (boxID > koff) definitely store the top edge crossing point

        !currently storing exits (except for boundary edges) - original stored exits

        ! !-------work on first crossing (edge1)
        ! ncr = ncr + 1

        ! !find the box the contour is coming out of (kob) given that it enters at edge1
        ! select case(edge1)
        !   case(BOTTOM)
        !     kob = box_ID - nxu
        !     if (kob < 1) then
        !       kob = 0 !contour is coming out of the boundary
        !     endif

        !   case(TOP)
        !     kob = box_ID + nxu
        !     if (kob > nxny) then
        !       kob = 0 !contour is coming out of the boundary
        !     endif

        !   case(LEFT)
        !     kob = box_ID - 1 !if mod(box_ID, nxu) == 1 then box_ID is on the left edge
        !     if (mod(box_ID, nxu) == 1) then
        !       kob = box_ID + nxu - 1 !contour is wrapping around from right to left
        !     endif

        !   case(RIGHT)
        !     kob = box_ID + 1
        !     if (mod(box_ID, nxu) == 0) then
        !       kob = box_ID - nxu + 1 !contour is wrapping around from left to right
        !     endif

        ! end select

        ! !store crossing points and connectivity information
        ! xcr(ncr) = x1
        ! ycr(ncr) = y1
        ! kib(ncr) = box_ID

        ! if (kob == 0) then
        !   npe = npe + 1
        !   icre(npe)=ncr
        ! else
        !   !if kob is zero then out of bounds for notcab and icrtab
        !   !& kob wasn't used in original ug2c closed crossing detection
        !   noctab(kob)=noctab(kob)+1
        !   icrtab(kob,noctab(kob))=ncr
        ! endif

        ! ------- boundary box entry storage since generally storing exits only
        !---- begin boundary_edge_storage()
        call ug2c_boundary_edge_storage(edge1, x1, y1, box_ID, iy, ix, &
                                        ncr, npe, xcr, ycr, kib, icre, &
                                        noctab, icrtab_local_ncr, icrtab_threadID, thread_id)
        !---- end boundary_edge_storage()

        ! -------work on second (exit) crossing (edge2)
        ncr = ncr + 1

        
        !find the box the contour is going into (kib) given that it leaves at edge2
        !---- begin calc_edge2_kib()
        call ug2c_calc_edge2_kib(edge2, box_ID, iy, ix, kib(ncr))
        !---- end calc_edge2_kib()

        !store crossing points and connectivity information
        !---- begin store_crossing_and_connectivity()
        call ug2c_store_crossing_and_connectivity(x2, y2, ncr, box_ID, ix, &
                                                  xcr, ycr, kob, noctab, &
                                                  icrtab_local_ncr, icrtab_threadID, &
                                                  thread_id)
        !---- end store_crossing_and_connectivity()


        !opposite (of original) store/build order if box is ambiguous (case 5 or 10)
        ! if (ms_case == 5 .or. ms_case == 10) then
        !   icrtab(kob, mod(noctab(kob),2)+1) = ncr
        ! else
        !   icrtab(kob,noctab(kob))=ncr
        ! endif

        !TODO: figure out how to utilise ambiguity resolved separation in rebuild/icrtab store order
        ! lookup table has all required info, rebuild utilises edges not segments
        ! To swap store (thus build) order use icrtab(kob, mod(noctab(kob),2)+1) = ncr instead

        !local orientation (cw or ccw) doesn't help since case A sep B will always look the same for same A and B
        !for a given separation: the first entry edge encountered during build traversal determines which exit should be first used

      enddo !end loop over segments in the cell
      ! !$OMP END CRITICAL
    enddo !loop over cells
    !$OMP END DO

    !store thread-local crossing data into shared arrays
    ncr_thread_list(thread_id) = ncr
    npe_thread_list(thread_id) = npe
    icre_thread_list(:, thread_id) = icre(:)
    xcr_thread_list(:, thread_id) = xcr(:)
    ycr_thread_list(:, thread_id) = ycr(:)
    kib_thread_list(:, thread_id) = kib(:)

    !ensure all threads have stored their data before MASTER thread gathers
    !$OMP BARRIER

    !$OMP MASTER
    !gather thread data into global storage
    
     !FYI:
      ! ncr 	  (global crossing ID)
      ! icrtab  ({box ID, box local crossing ID} -> global crossing ID. Maximum 2 crossings per box)
      ! noctab  (box ID -> total crossings in box)
      ! kib	    (global crossing ID -> box ID of that being entered after the given global crossing ID)
      ! xcr	    (x coord of crossing)
      ! ycr	    (y coord of crossing)
      ! npe	    (Counter for total number of open contours originating in an edge)
      ! icre	  (open contour ID -> global crossing ID of its first crossing)

    !check if num_threads = omp_get_num_threads()
    if (num_threads /= omp_get_num_threads()) then
      write(*,*) 'ERROR: num_threads=', num_threads, ' does not match omp_get_num_threads()=', omp_get_num_threads()
    endif

    if (timing_on) call timer_start(u2cMasterStart)
    !combine thread-local crossing data into shared arrays (serial section)
    ncr_global = 0
    npe_global = 0

    !calc ncr_offset for each thread to know where its local ncr fits into global
    ncr_offset(1)=0
    do i=2,num_threads
      ncr_offset(i) = ncr_offset(i-1) + ncr_thread_list(i-1)
    enddo

    !for each thread, i
    do i=1,num_threads
      !for each crossing, j, recorded by this thread, copy from thread-local arrays to global arrays
      do j=1,ncr_thread_list(i)
        ncr_global = ncr_global + 1
        xcr_global(ncr_global) = xcr_thread_list(j, i)
        ycr_global(ncr_global) = ycr_thread_list(j, i)
        kib_global(ncr_global) = kib_thread_list(j, i)
      enddo

      do j=1,npe_thread_list(i)
        npe_global = npe_global + 1
        icre_global(npe_global) = icre_thread_list(j, i) + ncr_offset(i)
      enddo
    enddo

    do i=1,nxny
      do j=1,2
        thread_id = icrtab_threadID(i,j)
        if (thread_id .ge. 1) then !if there was a crossing, combine into icrtab
          ! write(*,*) 'DEBUG icrtab entry for box', i, 'slot', j, 'thread_id=', thread_id, &
          !             'local_ncr=', icrtab_local_ncr(i,j), 'ncr_offset=', ncr_offset(thread_id)
          icrtab_global(i,j) = icrtab_local_ncr(i,j) + ncr_offset(thread_id)
        endif
      enddo
    enddo
    !----- Basic bounds checks after MASTER gather to catch corruption early -----

    !print the ncr_offsets
    ! do i=1,num_threads
    !   write(*,*) 'DEBUG ncr_offset for thread', i, '=', ncr_offset(i)
    ! enddo
    
    ! do ie=1,npe_global
    !   if (icre_global(ie) < 1 .or. icre_global(ie) > ncr_global) then
    !     write(*,*) 'ERROR: icre_global out of range after gather:', &
    !                 'ie=', ie, 'icre=', icre_global(ie), 'ncr_global=', ncr_global
    !     stop 1
    !   endif
    ! enddo

    ! do icr=1,ncr_global
    !   if (kib_global(icr) < 0 .or. kib_global(icr) > nxny) then
    !     write(*,*) 'ERROR: kib_global out of range after gather:', &
    !                 'icr=', icr, 'kib=', kib_global(icr), 'nxny=', nxny
    !     stop 1
    !   endif
    ! enddo

    ! do i=1,nxny
    !   do j=1,2
    !     if (icrtab_global(i,j) .ne. 0) then
    !       if (icrtab_global(i,j) < 1 .or. icrtab_global(i,j) > ncr_global) then
    !         write(*,*) 'ERROR: icrtab_global out of range after gather:', &
    !                     'box=', i, 'slot=', j, 'icrtab=', icrtab_global(i,j), 'ncr=', ncr_global
    !         stop 1
    !       endif
    !     endif
    !   enddo
    ! enddo

    if (timing_on) call timer_stop(u2cMasterStart,u2cMasterEnd,u2cMasterTime,u2cMasterTotTime)
    !$OMP END MASTER

    ! Deallocate thread-local arrays before exiting parallel region
    if (allocated(kib)) deallocate(kib)
    if (allocated(icre)) deallocate(icre)

    !$OMP END PARALLEL

    if (timing_on) call timer_stop(u2cParStart,u2cParEnd,u2cParTime,u2cParTotTime)

    !this saves renaming *_global below
    npe = npe_global
    ncr = ncr_global

    ! write(*,*) 'DEBUG level sweep done:', 'lev=', lev, 'ncr=', ncr, 'npe=', npe

    !Now re-build contours:
    if (timing_on) call timer_start(u2cRebuildStart)
    do icr=1,ncr
      free(icr)=.true.
    enddo

    !First deal with any open contours attached to boundaries:
    if (npe .gt. 0) then
      do ie=1,npe
        !A new contour (indexed na) starts here:
        na=na+1
        inda(na)=indq
        ibeg=npta+1
        i1a(na)=ibeg

        !The starting node on the contour (coming out of a boundary):
        icr=icre_global(ie)

        ! write(*,*) 'DEBUG open contour start:', 'lev=', lev, 'ie=', ie, 'icr=', icr, &
        !      'x=', xcr_global(icr), 'y=', ycr_global(icr), 'kib=', kib_global(icr)

        !First point on the contour:
        npd=1
        xd(1)=xcr_global(icr)
        yd(1)=ycr_global(icr)

        !Find remaining points on the contour:
        k=kib_global(icr)
        !k is the box the contour is entering (0 if going into a boundary)
        do while (k .ne. 0)
          noc=noctab(k)
          !Use last crossing in this box (noc) as the next node:
          icrn=icrtab_global(k,noc)
          !icrn gives the next point after icr (icrn is leaving box k)
          noctab(k)=noc-1
          !noctab is usually zero now except for boxes with a
          !maximum possible 2 crossings
          npd=npd+1
          !Coordinates of new node:
          xd(npd)=xcr_global(icrn)
          yd(npd)=ycr_global(icrn)
          ! write(*,*) 'DEBUG open contour coord', 'x=', xcr_global(icrn), 'y=', ycr_global(icrn)
          free(icrn)=.false.
          k=kib_global(icrn)
        enddo

        ! write(*,*) 'DEBUG open contour before re-node:', 'lev=', lev, 'ie=', ie, 'icr=', icr, &
        !            'npd=', npd

        !Re-distribute nodes on this contour 3 times to reduce complexity:
        keep=.false.
        do
          call renode_open(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
          !Delete contour if deemed too small (see renode_open):
          if (npa(na) .eq. 0) exit
          call renode_open(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
          !Delete contour if deemed too small (see renode_open):
          if (npd .eq. 0) exit
          call renode_open(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
          !Delete contour if deemed too small (see renode_open):
          if (npa(na) .eq. 0) exit
          !Contour is big enough to keep:
          keep=.true.
          exit
        enddo

        if (keep) then
          npta=npta+npa(na)
          iend=ibeg+npa(na)-1
          i2a(na)=iend
          do i=ibeg,iend-1
            nextq(i)=i+1
          enddo
          nextq(iend)=0
        else
          na=na-1
        endif

        free(icr)=.false.
      enddo
    endif

    !Next deal with remaining closed contours:
    do icr=1,ncr
      if (free(icr)) then
        !A new contour (indexed na) starts here:
        na=na+1
        inda(na)=indq
        ibeg=npta+1
        i1a(na)=ibeg

        !First point on the contour:
        npd=1
        xd(1)=xcr_global(icr)
        yd(1)=ycr_global(icr)

        ! write(*,*) 'DEBUG closed contour start:', 'lev=', lev, 'icr=', icr, &
        !      'x=', xcr_global(icr), 'y=', ycr_global(icr), 'kib=', kib_global(icr)

        !Find remaining points on the contour:
        k=kib_global(icr)
        !k is the box the contour is entering
        ! if (k < 1 .or. k > nxny) then
        !   write(*,*) 'DEBUG closed start k out of range:', 'lev=', lev, 'icr=', icr, &
        !              'k=', k, 'kib=', kib_global(icr), 'ncr=', ncr
        !   stop 1
        ! endif
        noc=noctab(k)
        ! if (noc < 1 .or. noc > 2) then
        !   write(*,*) 'DEBUG closed noc out of range:', 'lev=', lev, 'icr=', icr, &
        !              'k=', k, 'noc=', int(noc), 'noctab=', int(noctab(k)), &
        !              'icrtab1=', icrtab_global(k,1), 'icrtab2=', icrtab_global(k,2), 'ncr=', ncr
        !   stop 1
        ! endif
        !Use last crossing (noc) in this box (k) as the next node:
        icrn=icrtab_global(k,noc)
        !icrn gives the next point after icr (icrn is leaving box k)
        do while (icrn .ne. icr)
          noctab(k)=noc-1
          !noctab is usually zero now except for boxes with a
          !maximum possible 2 crossings
          npd=npd+1
          xd(npd)=xcr_global(icrn)
          yd(npd)=ycr_global(icrn)
          ! write(*,*) 'DEBUG closed contour coord', 'x=', xcr_global(icrn), 'y=', ycr_global(icrn)

          free(icrn)=.false.
          k=kib_global(icrn)
          ! if (k < 1 .or. k > nxny) then
          !   write(*,*) 'DEBUG closed hop k out of range:', 'lev=', lev, 'icr=', icr, &
          !              'icrn=', icrn, 'k=', k, 'ncr=', ncr
          !   write(*,*) 'DEBUG ', 'noc=', int(noc), 'icrn=', icrn
          !   stop 1
          ! endif

          noc=noctab(k)
          ! if (noc < 1 .or. noc > 2) then
          !   write(*,*) 'DEBUG closed noc out of range:', 'lev=', lev, 'icr=', icr, &
          !             'k=', k, 'noc=', int(noc), 'noctab=', int(noctab(k)), &
          !             'icrtab1=', icrtab_global(k,1), 'icrtab2=', icrtab_global(k,2), 'ncr=', ncr
          !   stop 1
          ! endif
          icrn=icrtab_global(k,noc)
        enddo

        ! write(*,*) 'DEBUG closed contour before re-node:', 'lev=', lev, 'icr=', icr, &
        !            'npd=', npd

        !Re-distribute nodes on this contour 3 times to reduce complexity:
        keep=.false.
        do
          call renode_closed(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
          !Delete contour if deemed too small (see renode_closed):
          if (npa(na) .eq. 0) exit
          call renode_closed(xa(ibeg),ya(ibeg),npa(na),xd,yd,npd)
          !Delete contour if deemed too small (see renode_closed):
          if (npd .eq. 0) exit
          call renode_closed(xd,yd,npd,xa(ibeg),ya(ibeg),npa(na))
          !Delete contour if deemed too small (see renode_closed):
          if (npa(na) .eq. 0) exit
          !Contour is big enough to keep:
          keep=.true.
          exit
        enddo

        if (keep) then
          npta=npta+npa(na)
          iend=ibeg+npa(na)-1
          i2a(na)=iend
          do i=ibeg,iend-1
            nextq(i)=i+1
          enddo
          nextq(iend)=ibeg
        else
          na=na-1
        endif

        free(icr)=.false.
      endif
    enddo !closed contour ncr loop
    if (timing_on) call timer_stop(u2cRebuildStart,u2cRebuildEnd,u2cRebuildTime,u2cRebuildTotTime)
  enddo !loop over levels
endif

  ! Deallocate per-thread arrays
  if (allocated(ncr_offset)) deallocate(ncr_offset)
  if (allocated(ncr_thread_list)) deallocate(ncr_thread_list)
  if (allocated(npe_thread_list)) deallocate(npe_thread_list)
  if (allocated(icre_thread_list)) deallocate(icre_thread_list)
  if (allocated(xcr_thread_list)) deallocate(xcr_thread_list)
  if (allocated(ycr_thread_list)) deallocate(ycr_thread_list)
  if (allocated(kib_thread_list)) deallocate(kib_thread_list)

return
end subroutine

! =======================================================================

subroutine con2ugrid(xq,yq,dq,qavg,nextq,nptq)
! Contour -> grid conversion.  The contours are represented by
! nodes (xq(i),yq(i)), i = 1, ..., nptq, where nextq(i) gives the
! index of the node following i, dq is the jump in q across all
! contours, and qavg is average value of the field.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: xq(npm),yq(npm)
integer:: nextq(npm)

 !Local parameters and arrays:
double precision:: qjx(0:nxum1),qbot(0:nxum1)
double precision:: dx(nptq),dy(nptq)
integer:: ixc(nptq),nxc(nptq)
logical:: crossx(nptq)

!----------------------------------------------------------------
 !Initialise interior x grid line crossing information and fill the
 !q jump array along lower boundary:
do i=1,nptq
  ixc(i)=1+int(glxui*(xq(i)-xmin))
enddo

do ix=0,nxum1
  qjx(ix)=zero
enddo

do i=1,nptq
  ia=nextq(i)
  if (ia .gt. 0) then
     !A node with ia = 0 terminates a contour at a boundary
    xx=xq(ia)-xq(i)
    dx(i)=xx-ellx*dble(int(xx*hlxi))
    dy(i)=yq(ia)-yq(i)
    ixdif=ixc(ia)-ixc(i)
    nxc(i)=ixdif-nxu*((2*ixdif)/nxu)
    crossx(i)=(nxc(i) .ne. 0)
    if ((yq(ia)-ybeg)*(ybeg-yq(i)) .gt. zero) then
       !The contour segment (i,ia) crosses y = ybeg; find x location:
      py0=(ybeg-yq(i))/dy(i)
      xx=xq(i)+py0*dx(i)
      xx=oms*(xx-ellx*dble(int(xx*hlxi)))
      ix=int(glxui*(xx-xmin))
      qjx(ix)=qjx(ix)-dq*sign(one,dy(i))
       !Note: qjx gives the jump going from ix to ix+1
    endif
  else
     !Here, there is no segment (i,next(i)) to consider:
    crossx(i)=.false.
  endif
enddo
 !Above, ybeg is very slightly greater than ymin to detect boundary crossings

 !Sum q jumps to obtain the gridded q along lower boundary:
qbot(0)=zero
 !Corner value cannot be determined a priori; qavg is used for this below
do ix=0,nxum2
  qbot(ix+1)=qbot(ix)+qjx(ix)
enddo

!----------------------------------------------------------------
 !Initialise interior q jump array:
do ix=0,nxum1
  do iy=0,nyup1
    qa(qa_idx(iy,ix))=zero
  enddo
enddo

 !Determine x grid line crossings and accumulate q jumps:
do i=1,nptq
  if (crossx(i)) then
    jump=sign(1,nxc(i))
    ixbeg=ixc(i)+(jump-1)/2+nxu
    sdq=dq*sign(one,dx(i))
    ncr=0
    do while (ncr .ne. nxc(i))
      ix=mod(ixbeg+ncr,nxu)
      xx=xgu(ix)-xq(i)
      px0=(xx-ellx*dble(int(xx*hlxi)))/dx(i)
       !The contour crossed the fine grid line ix at the point
       !   x = xq(i) + px0*dx(i) and y = yq(i) + px0*dy(i):
      iy=int(one+dyyui*(yq(i)+px0*dy(i)-ybeg))
       !Increment q jump between the grid lines iy-1 & iy:
      qa(qa_idx(iy,ix))=qa(qa_idx(iy,ix))+sdq
       !Go on to consider next x grid line (if there is one):
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get q values by sweeping through y:
do ix=0,nxum1
  qa(qa_idx(0,ix))=qbot(ix)
  do iy=1,nyu
    qa(qa_idx(iy,ix))=qa(qa_idx(iy,ix))+qa(qa_idx(iy-1,ix))
  enddo
enddo

 !Restore average (use qjx as temp array):
do ix=0,nxum1
  qjx(ix)=f12*(qa(qa_idx(0,ix))+qa(qa_idx(nyu,ix)))
  do iy=1,nyu-1
    qjx(ix)=qjx(ix)+qa(qa_idx(iy,ix))
  enddo
enddo

qavg0=zero
do ix=0,nxum1
  qavg0=qavg0+qjx(ix)
enddo
qavg0=qavg0/dble(nxu*nyu)

qadd=qavg-qavg0
do ix=0,nxum1
  do iy=0,nyu
    qa(qa_idx(iy,ix))=qa(qa_idx(iy,ix))+qadd
  enddo
enddo

return
end subroutine

!==========================================================================

 !Main end module
end module
