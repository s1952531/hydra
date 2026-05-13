module congen
! Converts contours (xq,yq) to gridded values on an ultra-fine grid
! of dimensions mgu*nx x mgu*ny (in a closed rectangular domain),
! optionally adds a residual field (interpolated to the ultra-fine
! grid), and then creates new contours.

! Open contours originating and terminating in a boundary added by
! D G Dritschel on 18 June 2012 @ Moscow

use common
use generic

implicit none

double precision:: qa(0:nyup1,0:nxum1)
double precision:: xa(npm),ya(npm)
integer:: inda(nm),npa(nm),i1a(nm),i2a(nm)
integer:: na,npta

!Marching squares definitions:
integer, parameter:: LEFT=1, RIGHT=2, BOTTOM=3, TOP=4

contains

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

        qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qq(iy0,ix0) &
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

      qa(iy,ix)=w00(iyf,ixf)*qq(iy0,ix0)+w10(iyf,ixf)*qq(iy1,ix0) &
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

  double precision:: ycr(ncrm),xcr(ncrm)
  double precision:: qdx(0:nxu),qdy(0:nyu)
  double precision:: xd(nprm),yd(nprm)
  integer:: isx(0:nxu),isy(0:nyu)
  integer:: kib(ncrm),icre(nm)
  integer:: icrtab(nxny,2)
  integer*1:: noctab(nxny)
  logical:: free(ncrm),keep
  logical:: is_open

  !end original ug2c definitions

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
  qamax=qa(0,0)
  qamin=qa(0,0)
  do ix=0,nxum1
    do iy=0,nyu
      qamax=max(qamax,qa(iy,ix))
      qamin=min(qamin,qa(iy,ix))
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

    ncr=0 !Counter for total number of grid line crossings
    npe=0 !Counter for total number of open contours originating in an edge

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

    do box_ID=1,nxny !grid boxes are numbered 1 (lower left) to nxu*nyu (upper right)

      !lower left corner has indices (iy,ix) and is at coords xgu(ix), ygu(iy)
      iy=(box_ID-1)/nxu
      ix=mod(box_ID-1,nxu)

      !get corner values and note min and max
      ll=qa(iy,ix)
      ul=qa(iy+1,ix)
      ur=qa(iy+1,ix+1)
      lr=qa(iy,ix+1)

      !get min and max of the corners
      minVal=min(ll,ul,ur,lr)
      maxVal=max(ll,ul,ur,lr)

      !loop over levels here when parallel rebuilding (per tile) is implemented

      !could skip this level early if it does not cross the cell
      !should this use lt, gt or le, ge?
      ! if ((qtmp .le. minVal) .or. (qtmp .ge. maxVal)) cycle

      !get marching squares case
      ms_case = 0
      if (ll >= qtmp) ms_case = ms_case + 1
      if (lr >= qtmp) ms_case = ms_case + 2
      if (ur >= qtmp) ms_case = ms_case + 4
      if (ul >= qtmp) ms_case = ms_case + 8
      !Using >= condition: What happens if >?

      !if ms_case is 0 or 15 then there are no crossings i.e. nseg=0 so the loop is skipped and we move to the next cell.
      if (ms_case == 0 .or. ms_case == 15) cycle

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

      !calculate all interpolation points for this cell and level (some may not be needed but this minimises branching):
      ! bottom: ll -> lr
      dz = lr - ll
      dz_safe = sign(max(abs(dz), eps), dz)
      t_interp = (qtmp - ll) / dz_safe
      x_b_interp = xgu(ix) + t_interp*glxu
      y_b_interp = ygu(iy)
      x_b_interp=oms*(x_b_interp-ellx*dble(int(x_b_interp*hlxi)))

      ! top: ul -> ur
      dz = ur - ul
      dz_safe = sign(max(abs(dz), eps), dz)
      t_interp = (qtmp - ul) / dz_safe
      x_t_interp = xgu(ix) + t_interp*glxu
      y_t_interp = ygu(iy+1)
      x_t_interp=oms*(x_t_interp-ellx*dble(int(x_t_interp*hlxi)))

      ! left: ll -> ul
      dz = ul - ll
      dz_safe = sign(max(abs(dz), eps), dz)
      t_interp = (qtmp - ll) / dz_safe
      x_l_interp = xgu(ix)
      y_l_interp = ygu(iy) + t_interp*glyu

      ! right: lr -> ur
      dz = ur - lr
      dz_safe = sign(max(abs(dz), eps), dz)
      t_interp = (qtmp - lr) / dz_safe
      x_r_interp = xgu(ix+1)
      y_r_interp = ygu(iy) + t_interp*glyu

      !lookup case to get edges that are crossed
      do seg=1,nseg(ms_case)
        edge1 = edge1_lookup(ms_case,separated,seg)
        edge2 = edge2_lookup(ms_case,separated,seg)

        !get the coordinates of the crossing points for edge1 and edge2
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

        ! if ((box_ID .eq. 3837947) .or. (box_ID .eq. 3837948) .or. (box_ID .eq. 3837949) .or. (box_ID .eq. 3829757)) then
        !   write (*,*) 'DEBUG ', 'kib= ', kib(ncr), 'kob= ', kob, 'at crossing coords: ', '(', x1, ',', y1, ')', ' on edge1= ', edge1
        ! endif

        ! ------- boundary box entry storage since generally storing exits only

        !if box on top row (box_ID > koff) need to store entry at top
        if (edge1 .eq. TOP .and. box_ID > koff) then
          ncr = ncr + 1
          kob = 0
          xcr(ncr) = x1
          ycr(ncr) = y1
          kib(ncr) = box_ID
          npe = npe + 1
          icre(npe)=ncr
        endif

        !if box on bottom row (box_ID <= nxu) need to store entry at bottom
        if (edge1 .eq. BOTTOM .and. box_ID <= nxu) then
          ncr = ncr + 1
          kob = 0
          xcr(ncr) = x1
          ycr(ncr) = y1
          kib(ncr) = box_ID
          npe = npe + 1
          icre(npe)=ncr
        endif

        !if box on left column (mod(box_ID, nxu) == 1) need to store entry at left
        if (edge1 .eq. LEFT .and. mod(box_ID, nxu) == 1) then
          ncr = ncr + 1
          kob = box_ID + nxu - 1 !contour is wrapping around from right to left
          xcr(ncr) = x1
          ycr(ncr) = y1
          kib(ncr) = box_ID
          noctab(kob)=noctab(kob)+1
          icrtab(kob,noctab(kob))=ncr
        endif

        !if box on right column (mod(box_ID, nxu) == 0) need to store entry at right
        if (edge1 .eq. RIGHT .and. mod(box_ID, nxu) == 0) then
          ncr = ncr + 1
          kob = box_ID - nxu + 1 !contour is wrapping around from left to right
          xcr(ncr) = x1
          ycr(ncr) = y1
          kib(ncr) = box_ID
          noctab(kob)=noctab(kob)+1
          icrtab(kob,noctab(kob))=ncr
        endif

        ! -------work on second (exit) crossing (edge2)
        ncr = ncr + 1

        !find the box the contour is going into (kib) given that it leaves at edge2
        select case(edge2)
          case(BOTTOM)
            kib(ncr) = box_ID - nxu
            if (kib(ncr) < 1) then
              kib(ncr) = 0 !contour is going into the boundary
            endif

          case(TOP)
            kib(ncr) = box_ID + nxu
            if (kib(ncr) > nxny) then
              kib(ncr) = 0 !contour is going into the boundary
            endif

          case(LEFT)
            kib(ncr) = box_ID - 1
            if (mod(box_ID, nxu) == 1) then
              kib(ncr) = box_ID + nxu - 1 !contour is wrapping around from left to right
            endif

          case(RIGHT)
            kib(ncr) = box_ID + 1
            if (mod(box_ID, nxu) == 0) then
              kib(ncr) = box_ID - nxu + 1 !contour is wrapping around from right to left
            endif

        end select

        !store crossing points and connectivity information
        xcr(ncr) = x2
        ycr(ncr) = y2
        kob = box_ID
        noctab(kob)=noctab(kob)+1

        !original store/build order:
        icrtab(kob,noctab(kob))=ncr

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

    enddo !loop over cells

    ! write(*,*) 'DEBUG level sweep done:', 'lev=', lev, 'ncr=', ncr, 'npe=', npe

    !Now re-build contours:
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
        icr=icre(ie)

        ! write(*,*) 'DEBUG open contour start:', 'lev=', lev, 'ie=', ie, 'icr=', icr, &
        !      'x=', xcr(icr), 'y=', ycr(icr), 'kib=', kib(icr)

        !First point on the contour:
        npd=1
        xd(1)=xcr(icr)
        yd(1)=ycr(icr)

        !Find remaining points on the contour:
        k=kib(icr)
        !k is the box the contour is entering (0 if going into a boundary)
        do while (k .ne. 0)
          noc=noctab(k)
          !Use last crossing in this box (noc) as the next node:
          icrn=icrtab(k,noc)
          !icrn gives the next point after icr (icrn is leaving box k)
          noctab(k)=noc-1
          !noctab is usually zero now except for boxes with a
          !maximum possible 2 crossings
          npd=npd+1
          !Coordinates of new node:
          xd(npd)=xcr(icrn)
          yd(npd)=ycr(icrn)
          ! write(*,*) 'DEBUG open contour coord', 'x=', xcr(icrn), 'y=', ycr(icrn)
          free(icrn)=.false.
          k=kib(icrn)
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
        xd(1)=xcr(icr)
        yd(1)=ycr(icr)

        ! write(*,*) 'DEBUG closed contour start:', 'lev=', lev, 'icr=', icr, &
        !      'x=', xcr(icr), 'y=', ycr(icr), 'kib=', kib(icr)

        !Find remaining points on the contour:
        k=kib(icr)
        !k is the box the contour is entering
        if (k < 1 .or. k > nxny) then
          write(*,*) 'DEBUG closed start k out of range:', 'lev=', lev, 'icr=', icr, &
                     'k=', k, 'kib=', kib(icr), 'ncr=', ncr
          stop 1
        endif
        noc=noctab(k)
        if (noc < 1 .or. noc > 2) then
          write(*,*) 'DEBUG closed noc out of range:', 'lev=', lev, 'icr=', icr, &
                     'k=', k, 'noc=', int(noc), 'noctab=', int(noctab(k)), &
                     'icrtab1=', icrtab(k,1), 'icrtab2=', icrtab(k,2), 'ncr=', ncr
          stop 1
        endif
        !Use last crossing (noc) in this box (k) as the next node:
        icrn=icrtab(k,noc)

        !icrn gives the next point after icr (icrn is leaving box k)
        do while (icrn .ne. icr)
          ! write(*,*) 'DEBUG :', 'lev=', lev, 'icr=', icr, 'k=', k, &
          !            'noc=', int(noc), 'icrn=', icrn, 'next_k=', kib(icrn), &
          !            'noctab_next=', int(noctab(kib(icrn))), 'ncr=', ncr
          noctab(k)=noc-1
          !noctab is usually zero now except for boxes with a
          !maximum possible 2 crossings
          npd=npd+1
          xd(npd)=xcr(icrn)
          yd(npd)=ycr(icrn)
          ! write(*,*) 'DEBUG closed contour coord', 'x=', xcr(icrn), 'y=', ycr(icrn)

          free(icrn)=.false.
          ! if (kib(icrn) == 0) then
          !   write(*,*) 'DEBUG closed boundary transition:', 'lev=', lev, 'icr=', icr, &
          !              'k=', k, 'noc=', int(noc), 'icrn=', icrn, 'x=', xcr(icrn), &
          !              'y=', ycr(icrn), 'box1=', icrtab(k,1), 'box2=', icrtab(k,2)
          !   write(*,*) 'DEBUG ', 'k=', k, 'mod(k, nxu)= ', mod(k, nxu), 'k/nxu=', k/nxu, 'nxu=', nxu, 'nyu=', nyu
          ! endif
          k=kib(icrn)
          if (k < 1 .or. k > nxny) then
            write(*,*) 'DEBUG closed hop k out of range:', 'lev=', lev, 'icr=', icr, &
                       'icrn=', icrn, 'k=', k, 'ncr=', ncr
            write(*,*) 'DEBUG ', 'noc=', int(noc), 'icrn=', icrn
            stop 1
          endif

          noc=noctab(k)
          icrn=icrtab(k,noc)
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
  enddo !loop over levels
endif
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
    qa(iy,ix)=zero
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
      qa(iy,ix)=qa(iy,ix)+sdq
       !Go on to consider next x grid line (if there is one):
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get q values by sweeping through y:
do ix=0,nxum1
  qa(0,ix)=qbot(ix)
  do iy=1,nyu
    qa(iy,ix)=qa(iy,ix)+qa(iy-1,ix)
  enddo
enddo

 !Restore average (use qjx as temp array):
do ix=0,nxum1
  qjx(ix)=f12*(qa(0,ix)+qa(nyu,ix))
  do iy=1,nyu-1
    qjx(ix)=qjx(ix)+qa(iy,ix)
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
    qa(iy,ix)=qa(iy,ix)+qadd
  enddo
enddo

return
end subroutine

!==========================================================================

 !Main end module
end module
