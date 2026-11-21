module congen

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!  Converts PV contours to gridded values on an ultra-fine grid of
!  dimensions mgu*numLongGridPts x mgu*numLatGridPts, adds the residual qResid (interpolated to
!  the ultra-fine grid), then creates new contours.

!  Adapted from ~dgd/cs/spe/sources/mcongen.F on 26/2/13
!  by Stuart King & DG Dritschel @ St Andrews

!  Revised 6/11/2026 by D G Dritschel @ Kyoto to use any number of layers.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

use common

implicit none

 !Array for storing the PV field interpolated to the ultra-fine grid: 
double precision:: qa(0:ngu+1,ntu+1)

 !Temporary arrays for contour storage:
double precision:: xa(npm),ya(npm),za(npm)
integer:: nexta(npm),inda(nm),npa(nm),i1a(nm),i2a(nm)
integer:: na,npta
 !Note: contours are built layer by layer and level by level to keep
 !storage to a minimum.


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains
 
!==========================================================================

subroutine recontour
! Main routine for recontouring (from Dritschel & Ambaum, 1996, QJRMS)
! but adapted to spherical geometry.

! The following quantities are available in commons.f90, included here:
  
! qResid             : a gridded field added to that due to contours;
!                  if there are no contours (as when t = 0),
!                  then qResid is the full PV field
! qSPole, qNPole       : the PV at the south and north poles; these are
!                  needed to uniquely define the combined field
!                  (the array qa in this module)
! x(i),y(i),z(i) : Cartesian position of node i on the sphere
! nextNode(i)        : index of the node following node i 
! contourIdxToLayerIdx(j)         : layer (integer) containing contour j
! contsQIdx(j)         : field level (integer) of contour j
! numNodesOnContourList(j)          : number of nodes on contour j
! firstNodeIdxList(j)          : beginning node index on contour j
! lastNodeIdxList(j)          : ending node index on contour j
! layerFirstContIdx(iz)        : beginning contour index in layer iz
! layerLastContIdx(iz)        : ending contour index in layer iz
! layerFirstNodeIdx(iz)        : beginning node index in layer iz
! layerLastNodeIdx(iz)        : ending node index in layer iz
! numContoursTot              : number of contours
! numNodesTot            : total number of nodes

implicit none

 !Local variables:
double precision:: qq0(0:numLatGridPts,numLongGridPts)
double precision:: xtd(numLongGridPts),utd(numLongGridPts),xend
double precision:: qrsp,qrnp,qinc
integer:: ix,ixf,ix0,ix1
integer:: iy,iyf,iy0,iy1
integer:: ic,i,j,iz

 !-----------------------------------------------------------------
 !Counters for total number of nodes and contours:
npta=0
na=0

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !              Begin a major loop over layers
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do iz=1,numLayers
   !Obtain ultra-fine grid field qa in layer iz:

   !Interpolate the residual qResid to the fine grid and add to qa
   !(note: qResid contains the full field at the start of a simulation);

   !loop over great circles in latitude (only half the longitudes):
  do ix=1,numLatGridPts
    ic=ix+numLatGridPts

    !Interpolate qResid (store in qq0) by 4th-order method (periodic):
    utd(1)=f23*(qResid(1,ix,iz)+qResid(1,ic,iz))
    do j=2,numLatGridPts
      utd(j)=f23*(qResid(j,ix,iz)+qResid(j-1,ix,iz))
    enddo
    utd(ngp1)=f23*(qResid(numLatGridPts,ic,iz)+qResid(numLatGridPts,ix,iz))
    do j=ngp2,numLongGridPts
      utd(j)=f23*(qResid(ntp2-j,ic,iz)+qResid(ntp1-j,ic,iz))
    enddo

    xtd(1)=utd(1)*htd(1)
    do j=2,numLongGridPts
      xtd(j)=(utd(j)-f16*xtd(j-1))*htd(j)
    enddo
    do j=ntm2,1,-1
      xtd(j)=etd(j)*xtd(j+1)+xtd(j)
    enddo
    xtd(numLongGridPts)=(etd(numLongGridPts)*xtd(1)+xtd(numLongGridPts))*xndeno
    xend=xtd(numLongGridPts)

    do j=1,ntm1
      xtd(j)=ptd(j)*xend+xtd(j)
    enddo

    do j=0,numLatGridPts
      qq0(j,ix)=xtd(j+1)
    enddo
    qq0(0,ic)=xtd(1)
    do j=1,numLatGridPts
      qq0(j,ic)=xtd(ntp1-j)
    enddo
  enddo
  !Ends loops over great circles. Interpolation complete.

  !Obtain unique polar values of the residual qq0 for use below:
  qrsp=sum(qq0( 0,:))/dble(numLongGridPts)
  qrnp=sum(qq0(numLatGridPts,:))/dble(numLongGridPts)
  qq0(0 ,:)=qrsp
  qq0(numLatGridPts,:)=qrnp

  if (layerLastContIdx(iz) > 0) then
    !This means that there are also contours in this layer;
    !convert them to gridded values (in the array qa):
    call con2ugrid(iz)

    !Bi-linearly interpolate qq0 to the fine grid and add to qa:
    do ix=1,ntu
      ixf=ixfw(ix)
      ix0=ix0w(ix)
      ix1=ix1w(ix)
  
      qa(0,ix)=qa(0,ix)+qrsp
      do iy=1,ngu-1
        iyf=iyfw(iy)
        iy0=iy0w(iy)
        iy1=iy1w(iy)
        qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qq0(iy0,ix0) &
                           +w10(iyf,ixf)*qq0(iy1,ix0) &
                           +w01(iyf,ixf)*qq0(iy0,ix1) &
                           +w11(iyf,ixf)*qq0(iy1,ix1)
      enddo
      qa(ngu,ix)=qa(ngu,ix)+qrnp
    enddo
    !qrsp & qrnp are the polar qResid values (necessarily uniform).

    !Next adjust qa by an overall constant so that the PV at the south 
    !pole is the same as that in that of the input field, qSPole(iz):
    qinc=qSPole(iz)-qa(0,1)
    !note: qa does not vary with longitude at either pole
    qa(0:ngu,1:ntu)=qa(0:ngu,1:ntu)+qinc

  else
    !There are no contours (the usual situation at t = 0);
    !interpolate qq0 (which here contains the full field)
    !directly to the fine grid as qa:
    do ix=1,ntu
      ixf=ixfw(ix)
      ix0=ix0w(ix)
      ix1=ix1w(ix)

      qa(0,ix)=qSPole(iz)
      do iy=1,ngu-1
        iyf=iyfw(iy)
        iy0=iy0w(iy)
        iy1=iy1w(iy)
        qa(iy,ix)=w00(iyf,ixf)*qq0(iy0,ix0)+w10(iyf,ixf)*qq0(iy1,ix0)+ &
                  w01(iyf,ixf)*qq0(iy0,ix1)+w11(iyf,ixf)*qq0(iy1,ix1)
      enddo
      qa(ngu,ix)=qNPole(iz)
    enddo
    !qSPole & qNPole are the polar PV values (necessarily uniform).

  endif

  !Add a periodic column at ix = ntu+1:
  ix=ntu+1
  qa(0:ngu,ix)=qa(0:ngu,1)

  !Reset starting contour and node indices in each layer:
  layerFirstContIdx(iz)=na+1
  layerFirstNodeIdx(iz)=npta+1

  !Generate new contours (xa,ya) from qa array:
  call ugrid2con

  !Reset ending contour and node indices in each layer:
  layerLastContIdx(iz)=na
  layerLastNodeIdx(iz)=npta

enddo
! Ends the loop over layers

 !Copy contour arrays and indices back to those in the main code:
do i=1,npta
  x(i)=xa(i)
  y(i)=ya(i)
  z(i)=za(i)
  nextNode(i)=nexta(i)
enddo

do j=1,na
  firstNodeIdxList(j)=i1a(j)
  lastNodeIdxList(j)=i2a(j)
  numNodesOnContourList(j)=npa(j)
  contsQIdx(j)=inda(j)
enddo

do iz=1,numLayers
  if (layerLastContIdx(iz) .gt. 0) then
    do j=layerFirstContIdx(iz),layerLastContIdx(iz)
      contourIdxToLayerIdx(j)=iz
    enddo
  endif
enddo

numContoursTot=na
numNodesTot=npta

return
end subroutine recontour

!=======================================================================

subroutine con2ugrid(iz)
! Contour -> ultra-fine grid conversion in layer iz.  

implicit none

 !Passed variable:
integer:: iz

 !Local variables:
double precision:: cx(numNodesTot),cy(numNodesTot),cz(numNodesTot),sq(numNodesTot)
double precision:: sig,rlatc,p,qasp,qanp
integer:: ntc(numNodesTot),ilm1(numNodesTot)
integer:: k,ka,jump,ioff,i,j,ncr

!----------------------------------------------------------------
 !Initialise crossing information:
do k=layerFirstNodeIdx(iz),layerLastNodeIdx(iz)
  ilm1(k)=int(dlui*(pi+atan2(y(k),x(k))))
enddo

do k=layerFirstNodeIdx(iz),layerLastNodeIdx(iz)
  ka=nextNode(k)
  cx(k)=z(k)*y(ka)-y(k)*z(ka)
  cy(k)=x(k)*z(ka)-z(k)*x(ka)
  cz(k)=x(k)*y(ka)-y(k)*x(ka)
  ntc(k)=ilm1(ka)-ilm1(k)
enddo

do k=layerFirstNodeIdx(iz),layerLastNodeIdx(iz)
  sig=sign(one,cz(k))
  sq(k)=dq*sig
  ntc(k)=ntc(k)-ntu*((2*ntc(k))/ntu)
  if (sig*dble(ntc(k)) < zero) ntc(k)=-ntc(k)
  if (abs(cz(k)) > zero) then
    cx(k)=cx(k)/cz(k)
    cy(k)=cy(k)/cz(k)
  endif
enddo

!----------------------------------------------------------------------
 !Initialise PV jump array:
qa=zero

 !Determine crossing indices:
do k=layerFirstNodeIdx(iz),layerLastNodeIdx(iz)
  if (ntc(k) .ne. 0) then
    jump=sign(1,ntc(k))
    ioff=ntu+ilm1(k)+(1+jump)/2
    ncr=0
    do while (ncr .ne. ntc(k))
      i=1+mod(ioff+ncr,ntu)
      rlatc=dlui*(hpi+atan(cx(k)*clonu(i)+cy(k)*slonu(i)))
      j=int(rlatc)+1
      p=rlatc-dble(j-1)
      qa(j,i)=  qa(j,i)+(one-p)*sq(k)
      qa(j+1,i)=qa(j+1,i)+    p*sq(k)
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get PV values, at half latitudes, by sweeping through latitudes:
do i=1,ntu
  do j=2,ngu
    qa(j,i)=qa(j,i)+qa(j-1,i)
  enddo
enddo
 !Here, qa(j,i) stands for the PV at latitude j-1/2,
 !from j = 1, ..., ngu.

 !Determine unique polar values:
qasp=sum(qa(1  ,1:ntu))/dble(ntu)
qanp=sum(qa(ngu,1:ntu))/dble(ntu)

 !Finally, average half-grid PV to full grid where it is needed for
 !contour generation (ugrid2con):
do i=1,ntu
  qa(0,i)=qasp
  do j=1,ngu-1
    qa(j,i)=f12*(qa(j,i)+qa(j+1,i))
  enddo
  qa(ngu,i)=qanp
enddo

return
end subroutine con2ugrid

!=======================================================================
subroutine ugrid2con
! Generates contours (xa,ya,za) from the gridded field qa for the PV
! levels +/-dq/2, +/-3*dq/2, ....

implicit none

integer,parameter:: ncrm=3*(npm/2)
! ncrm: max number of contour crossings of a single field
!       level on the finest grid

 !Local Grid -> Contour arrays:
integer(kind=dbleint),parameter:: ntng=int(ntu,kind=dbleint)*int(ngu,kind=dbleint)
integer(kind=dbleint):: kob,kib(ncrm)
double precision:: ycr(ncrm),xcr(ncrm)
double precision:: qrx(ntu+1),qry(0:ngu)
double precision:: xd(nprm),yd(nprm),zd(nprm)
double precision:: qamax,qamin,qtmp,xgt,ygt,coslat
integer:: isx(ntu+1),isy(0:ngu)
integer:: icrtab(ntng,2)
integer:: levbeg,levend,popQIdx,indq
integer:: ncr,ix,iy,inc,icr,npd,noc,icrn,firstNodeIdx,lastNodeIdx,i
integer(kind=halfint):: noctab(ntng)
logical:: free(ncrm),keep

!--------------------------------------------------------
 !First get the beginning and ending contour levels:
qamin=minval(qa(0:ngu,1:ntu))
levbeg=int((qoff+qamin)*dqi+f12)+1

qamax=maxval(qa(0:ngu,1:ntu))
levend=int((qoff+qamax)*dqi+f12)

 !Return if no levels to process:
if (levbeg > levend) return 

 !Loop over contour levels and process:
do popQIdx=levbeg,levend
   !Integer index giving contour level:
  indq=popQIdx-nlevm+(popQIdx-1)/nlevm-1

   !Counter for total number of grid line crossings:
  ncr=0

   !Contour level being sought (defined in contours.f90):
  qtmp=qlev(popQIdx)

   !Initialise number of crossings per box:
  do kob=1,ntng
    noctab(kob)=0
  enddo

   !Find x grid line crossings first:
  do ix=1,ntu
    xgt=xgu(ix)

    do iy=0,ngu
      qry(iy)=qa(iy,ix)-qtmp
      isy(iy)=sign(one,qry(iy))
    enddo

    do iy=0,ngu-1
      if (isy(iy) .ne. isy(iy+1)) then
        ncr=ncr+1
        inc=(1-isy(iy))/2
        kib(ncr)=iy+1+ibx(ix,inc)
        kob=iy+1+ibx(ix,1-inc)
        noctab(kob)=noctab(kob)+1
        icrtab(kob,noctab(kob))=ncr
        xcr(ncr)=xgt
        ycr(ncr)=ygu(iy)-glyu*qry(iy)/(qry(iy+1)-qry(iy))
      endif
    enddo

  enddo

!   Above, kib = grid box into which the contour (containing icr) is going
!          kob =   "   "  out of "    "     "         "       "    " coming
!     [kob -> icr -> kib:  icr lies at the boundary between kob & kib]

 !Find y grid line crossings next (no crossings can occur at iy=0,ngu):
  do iy=1,ngu-1
    ygt=ygu(iy)

    do ix=1,ntu+1
      qrx(ix)=qa(iy,ix)-qtmp
      isx(ix)=sign(one,qrx(ix))
    enddo

    do ix=1,ntu
      if (isx(ix) .ne. isx(ix+1)) then
        ncr=ncr+1
        inc=(1-isx(ix))/2
        kib(ncr)=ibx(ix,1)+iy+1-inc
        kob=ibx(ix,1)+iy+inc
        noctab(kob)=noctab(kob)+1
        icrtab(kob,noctab(kob))=ncr
        ycr(ncr)=ygt
        xcr(ncr)=xgu(ix)-glxu*qrx(ix)/(qrx(ix+1)-qrx(ix))
      endif
    enddo

  enddo

!------------------------------------------------------------------------
   !Now re-build contours - converting to spherical geometry:
  free(1:ncr)=.true.

  do icr=1,ncr
    if (free(icr)) then
       !A new contour (indexed na) starts here:
      na=na+1
      inda(na)=indq
      firstNodeIdx=npta+1
      i1a(na)=firstNodeIdx

       !First point on the contour:
      npd=1
      coslat=cos(ycr(icr))
      xd(1)=coslat*cos(xcr(icr))
      yd(1)=coslat*sin(xcr(icr))
      zd(1)=sin(ycr(icr))

       !Find remaining points on the contour:
      kob=kib(icr)
       !kib(icr) is the box the contour is entering
      noc=noctab(kob)
       !Use last crossing (noc) in this box (kob) as the next node:
      icrn=icrtab(kob,noc)
       !icrn gives the next point after icr (icrn is leaving box kob)
      do while (icrn .ne. icr)
        noctab(kob)=noc-1
         !noctab is usually zero now except for boxes with a
         !maximum possible 2 crossings
        npd=npd+1
        coslat=cos(ycr(icrn))
        xd(npd)=coslat*cos(xcr(icrn))
        yd(npd)=coslat*sin(xcr(icrn))
        zd(npd)=sin(ycr(icrn))
        free(icrn)=.false.
        kob=kib(icrn)
        noc=noctab(kob)
        icrn=icrtab(kob,noc)
      enddo

       !Re-distribute nodes on this contour 3 times to reduce complexity:
      keep=.false.
      do
        call renode(xd,yd,zd,npd,xa(firstNodeIdx),ya(firstNodeIdx),za(firstNodeIdx),npa(na))
         !Delete contour if deemed too small (see renode):
        if (numNodesOnContourList(numContoursTot) == 0) exit
        call renode(xa(firstNodeIdx),ya(firstNodeIdx),za(firstNodeIdx),npa(na),xd,yd,zd,npd)
         !Delete contour if deemed too small (see renode):
        if (npd == 0) exit
        call renode(xd,yd,zd,npd,xa(firstNodeIdx),ya(firstNodeIdx),za(firstNodeIdx),npa(na))
         !Delete contour if deemed too small (see renode):
        if (npa(na) == 0) exit
         !Contour is big enough to keep:
        keep=.true.
        exit
      enddo

      if (keep) then 
        npta=npta+npa(na)
        lastNodeIdx=firstNodeIdx+npa(na)-1
        i2a(na)=lastNodeIdx
        do i=firstNodeIdx,lastNodeIdx-1
          nexta(i)=i+1
        enddo
        nexta(lastNodeIdx)=firstNodeIdx
      else
        na=na-1
      endif

      free(icr)=.false.
    endif
  enddo

enddo
 !End of loop over contour levels
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

return
end subroutine ugrid2con

!=======================================================================

 !Main end module
end module congen
