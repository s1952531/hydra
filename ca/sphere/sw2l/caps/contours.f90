module contours

! Module contains all subroutines related to contour advection, for the 
! spe suite of f90 codes.

use constants

implicit none

 !Maximum number of PV levels:
integer,parameter:: nlevm=10000

 !Contour arrays:
double precision:: x(npm),y(npm),z(npm)
integer:: firstNodeIdxList(nm),lastNodeIdxList(nm),numNodesOnContourList(nm),contsQIdx(nm),contourToLayerIdx(nm)
integer:: layerFirstNodeIdx(numLayers),layerLastNodeIdx(numLayers),layerFirstContIdx(numLayers),layerLastContIdx(numLayers)
integer:: next(0:npm),numContoursTot,numNodesTot

 !Contour to grid conversion quantities (fine grid used in this module):
double precision:: clonf(ntf),slonf(ntf)
double precision,parameter:: dlf=twoPi/dble(ntf),dlfi=dble(ntf)/(twoPi+small)

 !Contour to grid conversion quantities (ultra-fine grid used in congen):
double precision:: clonu(ntu),slonu(ntu)
double precision,parameter:: dlu=twoPi/dble(ntu),dlui=dble(ntu)/(twoPi+small)

 !Weight arrays for grid to ultra-fine grid interpolation:
double precision:: w00(mgu,mgu),w10(mgu,mgu)
double precision:: w01(mgu,mgu),w11(mgu,mgu)
integer:: ixfw(ntu),iyfw(ngu)
integer:: ix0w(ntu),iy0w(ngu)
integer:: ix1w(ntu),iy1w(ngu)

 !Grid -> Contour arrays:
double precision:: xgu(ntu+1),ygu(0:ngu)
double precision:: qlev(2*nlevm),qoff
integer(kind=dbleint):: ibx(ntu,0:1)

 !Half-grid -> Full-grid tri-diagonal arrays:
double precision:: etd(numLongGridPts),htd(numLongGridPts),ptd(numLongGridPts),xndeno

 !Coriolis frequency:
double precision:: fcor(numLatGridPts)

!Constants used in contour node placement (Fontane & Dritschel JCP 2009):
double precision,parameter:: amu=0.2d0, ell=6.25d0*dl
double precision,parameter:: elf=one/ell**2, densf=one/(amu*sqrt(ell))
double precision,parameter:: dm=amu**2*ell/four, dm2=dm**2, dmsq=four*dm2
double precision,parameter:: dmi=two/dm, d4small=small*dl**4


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains 

!=======================================================================

subroutine init_contours
! Initialises all quantities needed for contour advection.

implicit none

double precision:: fac,pxu,pyu,pxc,pyc,rlonu,rlonf
integer(kind=dbleint):: ngulong
integer:: ixu,iyu,ix,iy,ixx,iyy,popQIdx,i,j

!--------------------------------------------------------------------
 !Initialise area weights for interpolation of residual PV (qResid) 
 !onto the ultra-fine horizontal grid:
fac=one/dble(mgu)
do ixu=1,mgu
  pxu=fac*dble(ixu-1)
  pxc=one-pxu
  do iyu=1,mgu
    pyu=fac*dble(iyu-1)
    pyc=one-pyu
    w00(iyu,ixu)=pyc*pxc
    w10(iyu,ixu)=pyu*pxc
    w01(iyu,ixu)=pyc*pxu
    w11(iyu,ixu)=pyu*pxu
  enddo
enddo

 !modulo values to access above weights:
do ix=1,ntu
  ixx=(ix-1)/mgu
  ix0w(ix)=1+ixx
  ix1w(ix)=2+ixx-numLongGridPts*(ix0w(ix)/numLongGridPts)
  ixfw(ix)=ix-mgu*ixx
enddo

do iy=1,ngu-1
  iyy=iy/mgu
  iy0w(iy)=iyy
  iy1w(iy)=1+iyy
  iyfw(iy)=1+iy-mgu*iyy
enddo

!--------------------------------------------------------------
 !Quantities needed in generating new contours (ugrid2con):
qoff=dq*dble(nlevm)
! qoff: should be a large integer multiple of the 
!       contour interval, dq.  The multiple should exceed 
!       the maximum expected number of contour levels.

do popQIdx=1,2*nlevm
  qlev(popQIdx)=(dble(popQIdx)-f12)*dq-qoff
enddo

 !Fixed arrays:
do i=1,ntu
  rlonu=dlu*dble(i-1)-pi
  clonu(i)=cos(rlonu)
  slonu(i)=sin(rlonu)
enddo

 !Coordinates of grid lines (longitudes and latitudes):
do ix=1,ntu+1
  xgu(ix)=glxu*dble(ix-1)-pi
enddo
do iy=0,ngu
  ygu(iy)=glyu*dble(iy)-hpi
enddo

 !Grid box reference indices:
ngulong=ngu
do ix=1,ntu
  ibx(ix,1)=ngulong*(ix-1)
enddo
do ix=2,ntu
  ibx(ix,0)=ibx(ix-1,1)
enddo
ibx(1,0)=ibx(ntu,1)

!-----------------------------------------------
 !Fixed arrays:
do i=1,ntf
  rlonf=dlf*dble(i-1)-pi
  clonf(i)=cos(rlonf)
  slonf(i)=sin(rlonf)
enddo

 !Coriolis frequency (f):
do j=1,numLatGridPts
  fcor(j)=fpole*sin((dble(j)-f12)*dl-hpi)
enddo

!----------------------------------------------------------------------
 !Tri-diagonal arrays for half grid -> full grid interpolation of qResid
 !(used in congen.f90)

 !Initialise periodic tridiagonal problem:
htd(1)=one
ptd(1)=-f16*htd(1)
etd(1)=ptd(1)

do j=2,numLongGridPts
  htd(j)=one/(one+f16*etd(j-1))
  ptd(j)=-f16*ptd(j-1)*htd(j)
  etd(j)=-f16*htd(j)
enddo

ptd(ntm1)=etd(ntm1)+ptd(ntm1)
do j=ntm2,1,-1
  ptd(j)=etd(j)*ptd(j+1)+ptd(j)
enddo

xndeno=one/(one-etd(numLongGridPts)*ptd(1)-ptd(numLongGridPts))

return
end subroutine init_contours

!==========================================================================

subroutine velint(u1,u2,v1,v2,u,v,w)
! Interpolates the gridded velocity field (uu,vv) in each layer to 
! the PV contours as (u,v,w) in Cartesian components.

implicit none

 !Passed variables:
double precision:: u1(numLatGridPts,numLongGridPts),u2(numLatGridPts,numLongGridPts)
double precision:: v1(numLatGridPts,numLongGridPts),v2(numLatGridPts,numLongGridPts)
double precision:: u(numNodesTot),v(numNodesTot),w(numNodesTot)

 !Local variables:
double precision:: uu(0:numLatGridPts+1,numLongGridPts,numLayers),vv(0:numLatGridPts+1,numLongGridPts,numLayers)
double precision:: ri,rj,abl,bbl,cbl,dbl
double precision:: ulon,ulat,rr,rri,zula
integer:: i,ic,ip1,j,jp1,k,iz

 !Extend the velocity field adjacent to poles (j = 1 and numLatGridPts)
 !with a pi shift in longitude to simplify interpolation below:
uu(1:numLatGridPts,:,1)=u1
uu(1:numLatGridPts,:,2)=u2
vv(1:numLatGridPts,:,1)=v1
vv(1:numLatGridPts,:,2)=v2
do i=1,numLatGridPts
  ic=i+numLatGridPts
  uu(0,i,1)=-u1(1,ic)
  uu(0,i,2)=-u2(1,ic)
  vv(0,i,1)=-v1(1,ic)
  vv(0,i,2)=-v2(1,ic)

  uu(0,ic,1)=-u1(1,i)
  uu(0,ic,2)=-u2(1,i)
  vv(0,ic,1)=-v1(1,i)
  vv(0,ic,2)=-v2(1,i)

  uu(ngp1,i,1)=-u1(numLatGridPts,ic)
  uu(ngp1,i,2)=-u2(numLatGridPts,ic)
  vv(ngp1,i,1)=-v1(numLatGridPts,ic)
  vv(ngp1,i,2)=-v2(numLatGridPts,ic)

  uu(ngp1,ic,1)=-u1(numLatGridPts,i)
  uu(ngp1,ic,2)=-u2(numLatGridPts,i)
  vv(ngp1,ic,1)=-v1(numLatGridPts,i)
  vv(ngp1,ic,2)=-v2(numLatGridPts,i)
enddo

 !Loop over layers:
do iz=1,numLayers
   !Only process if the layer contains contours:
  if (layerLastContIdx(iz) > 0) then
     !Bi-linearly zonal and meridional velocity at the contour nodes (x,y,z):
    do k=layerFirstNodeIdx(iz),layerLastNodeIdx(iz)
      ri=dli*(pi+atan2(y(k),x(k)))
      i=1+int(ri)
      ip1=1+mod(i,numLongGridPts)
      bbl=dble(i)-ri
      abl=one-bbl

      rj=dli*(hpidl+asin(z(k)))
      j=int(rj)
      jp1=j+1
      cbl=rj-dble(j)
      dbl=one-cbl

      ulon=bbl*(dbl*uu(j,i,iz)  +cbl*uu(jp1,i,iz))+  &
           abl*(dbl*uu(j,ip1,iz)+cbl*uu(jp1,ip1,iz))
      ulat=bbl*(dbl*vv(j,i,iz)  +cbl*vv(jp1,i,iz))+  &
           abl*(dbl*vv(j,ip1,iz)+cbl*vv(jp1,ip1,iz))

      rr=sqrt(x(k)**2+y(k)**2+small)
      rri=one/rr
      zula=-z(k)*ulat
       !Convert to Cartesians at the contour nodes:
      u(k)=(x(k)*zula-y(k)*ulon)*rri
      v(k)=(y(k)*zula+x(k)*ulon)*rri
      w(k)=rr*ulat
    enddo
  endif
enddo

return
end subroutine velint

!==========================================================================

subroutine con2grid(qc)
! Calculates the PV anomaly field (qc) from the PV contours
! (x,y,z) in each layer.  Takes away Coriolis frequency (fcor).

implicit none

 !Passed array:
double precision:: qc(numLatGridPts,numLongGridPts,numLayers)

 !Local variables:
double precision:: qa(0:ngf+1,ntf),qaend(ngf/2)
double precision:: cx(numNodesTot),cy(numNodesTot),cz(numNodesTot),sq(numNodesTot)
double precision:: sig,rlatc,p
integer:: ilm1(numNodesTot),ntc(numNodesTot)
integer:: iz,k,ka,jump,ncr,i,j,ioff,ip1
integer:: ic,io,ie,je,ngh,nth,nghh,nthh,nghp1

!----------------------------------------------------------------
 !Begin a major loop over layers
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do iz=1,numLayers
  if (layerLastContIdx(iz) > 0) then
     !There are contours in this layer.  Initialise crossing information:
    do k=layerFirstNodeIdx(iz),layerLastNodeIdx(iz)
      ilm1(k)=int(dlfi*(pi+atan2(y(k),x(k))))
    enddo

    do k=layerFirstNodeIdx(iz),layerLastNodeIdx(iz)
      ka=next(k)
      cx(k)=z(k)*y(ka)-y(k)*z(ka)
      cy(k)=x(k)*z(ka)-z(k)*x(ka)
      cz(k)=x(k)*y(ka)-y(k)*x(ka)
      ntc(k)=ilm1(ka)-ilm1(k)
    enddo

    do k=layerFirstNodeIdx(iz),layerLastNodeIdx(iz)
      sig=sign(one,cz(k))
      sq(k)=dq*sig
      ntc(k)=ntc(k)-ntf*((2*ntc(k))/ntf)
      if (sig*dble(ntc(k)) < zero) ntc(k)=-ntc(k)
      if (abs(cz(k)) > zero) then
        cx(k)=cx(k)/cz(k)
        cy(k)=cy(k)/cz(k)
      endif
    enddo

     !Initialise PV jump array:
    qa=zero

     !Determine crossing indices:
    do k=layerFirstNodeIdx(iz),layerLastNodeIdx(iz)
      if (ntc(k) .ne. 0) then
        jump=sign(1,ntc(k))
        ioff=ntf+ilm1(k)+(1+jump)/2
        ncr=0
        do while (ncr .ne. ntc(k))
          i=1+mod(ioff+ncr,ntf)
          rlatc=dlfi*(hpi+atan(cx(k)*clonf(i)+cy(k)*slonf(i)))
          j=int(rlatc)+1
          p=rlatc-dble(j-1)
          qa(j,i)=  qa(j,i)+(one-p)*sq(k)
          qa(j+1,i)=qa(j+1,i)+    p*sq(k)
          ncr=ncr+jump
        enddo
      endif
    enddo

     !Get PV values, at half latitudes, by sweeping through latitudes:
    do i=1,ntf
      do j=2,ngf
        qa(j,i)=qa(j,i)+qa(j-1,i)
      enddo
    enddo
     !Here, qa(j,i) stands for the PV at latitude j-1/2,
     !from j = 1, ..., ngf.

     !----------------------------------------------------------------------
     !Average PV values on the fine grid to get corresponding 
     !values on the inversion grid (numLatGridPts,numLongGridPts):
    ngh=ngf
    nth=ntf

    do while (ngh > numLatGridPts)
       !Pre-store PV adjacent to poles at complementary longitudes (+pi):
      nthh=nth/2
      nghp1=ngh+1
      do i=1,nthh
        ic=i+nthh
        qa(0,i)=qa(1,ic)
        qa(0,ic)=qa(1,i)
        qa(nghp1,i)=qa(ngh,ic)
        qa(nghp1,ic)=qa(ngh,i)
      enddo

       !Work from SP to NP to define PV at full latitudes from averages
       !at adjacent half latitudes:
      do i=1,nth
        do j=0,ngh
          qa(j,i)=f12*(qa(j+1,i)+qa(j,i))
        enddo
      enddo

       !Now qa(j,i) is the PV at latitude j*(pi/ngh)-pi/2

       !Next 1-2-1 average these values to define PV at half latitudes
       !on a grid twice as coarse:
      nghh=ngh/2
      do i=1,nth
        do j=1,nghh
          je=2*j
          qa(j,i)=f12*qa(je-1,i)+f14*(qa(je-2,i)+qa(je,i))
        enddo
      enddo

       !Now perform analogous longitudinal 1-2-1 average:
      do j=1,nghh
        qaend(j)=f12*(qa(j,nth)+qa(j,1))
      enddo
      do i=1,nth-1
        ip1=i+1
        do j=1,nghh
          qa(j,i)=f12*(qa(j,i)+qa(j,ip1))
        enddo
      enddo
      do j=1,nghh
        qa(j,nth)=qaend(j)
      enddo
       !Now qa(j,i) gives the PV at the half-longitudes i + 1/2.

       !Average these on the twice coarser grid:
      do j=1,nghh
        qa(j,1)=f12*(qa(j,nth)+qa(j,1))
      enddo
      do i=2,nthh
        io=2*i-1
        ie=io-1
        do j=1,nghh
          qa(j,i)=f12*(qa(j,ie)+qa(j,io))
        enddo
      enddo

      ngh=nghh
      nth=nthh

    enddo

     !Finalise and take away f to define PV anomaly:
    do i=1,numLongGridPts
      do j=1,numLatGridPts
        qc(j,i,iz)=qa(j,i)-fcor(j)
      enddo
    enddo

  else
     !There are no contours in this layer; assign qc = 0:
    qc(:,:,iz)=zero
  endif
enddo
 !Ends loop over layers
    
return
end subroutine con2grid

!========================================================================
subroutine surgery
! Performs surgery and, afterwards, redistributes nodes

! Major revision 01/01/2001 by D. G. Dritschel to accelerate surgery
! using local boxes to minimise search costs.  All surgery is now
! done in a single way.

! Here we use rectangular, equal area z-longitude boxes, except at
! each pole where we use a circular patch of radius bwz, where bwz 
! is the z spacing of the grid lines.  This special treatment of
! the poles was introduced 07/07/2012 by dgd @ Wildwood Crest, NJ.

! Extended 7 November 2026 to multiple layers by DGD @ Kyoto
  
implicit none

 !Local arrays and parameters:
integer,parameter:: ngbs=numLongGridPts*numLatGridPts+1
! nlevm: max number of distinct vorticity levels
! ngbs:  (+1) max number of boxes used in surgery below
integer,parameter:: nsegm=npm
! nsegm: maximum number of segments belonging to contours of a 
!        given q level crossing through all the boxes.
!        The factor multiplying npm should be proportional to 
!        C*((dbar+amu*ell/2)/bw)^2, where C is equal to the 
!        ratio of the maximum number of nodes of a given q
!        level to npm, dbar is the average segment length,
!        and bw is the box width.

double precision:: dx(numNodesTot),dy(numNodesTot),dz(numNodesTot)
double precision:: xd(nprm),yd(nprm),zd(nprm)
double precision:: lonc(numNodesTot),dlonc(numNodesTot),dsq(numNodesTot)
double precision:: xa(npm),ya(npm),za(npm)
double precision:: xx,numNodesOnQLvl,approxNumLongBoxes,approxNumZBoxes,bw,longBoxWidth,recipLongBoxWidth,fac,zbm,bwz,bwzi
double precision:: c0,adl,adz,delx,dely,delz,aa,cc
double precision:: dxa,dya,dza
integer:: popQsFirstContIdx(nlevm),popQsLastContIdx(nlevm),popQsFirstNodeIdx(nlevm),popQsLastNodeIdx(nlevm),orderedPopQsIdxs(nlevm)
integer:: nspb(0:ngbs),kb1(0:ngbs),kb2(0:ngbs)
integer:: loc(nsegm),list(nsegm),node(nsegm)
integer:: i1a(nm),i2a(nm),nexta(npm)
integer:: i,ia,iz,firstContIdx,lastContIdx,oldQIdx,numPopulatedQLvls,currContIdx,newQIdx,qIdx,firstNodeIdx,lastNodeIdx
integer:: numLongBoxes,numZBoxes,nblp1,numBoxesTot,numSegs,ib
integer:: mbl,mbl1,mbl2,mbz1,mbz2,mbz,mmbz,mb,k,ks,isb,is
integer:: jbeg,npd
 !Logicals:
logical:: avail(numNodesTot)

!------------------------------------------------------------------------
 !Get work arrays for efficient surgery:
do i=1,numNodesTot
  xx=atan2(y(i),x(i))
  lonc(i)=xx-twoPi*dble(int(xx*hlxi))
enddo

do i=1,numNodesTot
  ia=next(i)
  xx=lonc(ia)-lonc(i)
  dlonc(i)=xx-twoPi*dble(int(xx*hlxi))
   !dlonc: the change in longitude from node i to node ia
  dx(i)=x(ia)-x(i)
  dy(i)=y(ia)-y(i)
  dz(i)=z(ia)-z(i)
enddo

do i=1,numNodesTot
  dsq(i)=dx(i)**2+dy(i)**2+dz(i)**2
enddo

 !Zero out counters for contours and nodes for building new contours:
numContoursTot=0
numNodesTot=0

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !              Begin a major loop over layers
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do iz=1,numLayers
  if (layerLastContIdx(iz) .gt. 0) then
    ! Then there are contours in this layer

    ! for this layer
    firstContIdx=layerFirstContIdx(iz) ! get its first contour index
    lastContIdx=layerLastContIdx(iz)   ! get its final contour index
    numPopulatedQLvls=1                           ! init number of populated q levels

    ! The following calculates the first and final contours of each populated 
    ! (i.e. contour containing) discrete q level
    ! Note: the following code assumes (thus requires) that 
    !       the contours are ordered by increasing q level 
    !       We will loop through contours which, under the ordering assumption,
    !       loops through the q levels from lowest to highest index
    !       As we go we store the first and final contour indices for each populated q level
    !       (Any discrete q level which has no contours is not stored)

    ! Start with the first contour in this layer
    oldQIdx=contsQIdx(firstContIdx)  ! store its q level index so we can tell when we have moved to a new q level
    popQsFirstContIdx(1)=firstContIdx    ! get its first contour index
    orderedPopQsIdxs(1)=oldQIdx           ! Since this q level has contours, store its index

    ! Loop through the remaining contours
    do currContIdx=firstContIdx+1,lastContIdx 
      newQIdx=contsQIdx(currContIdx)                   ! find what q level this contour is in
      if (newQIdx .gt. oldQIdx) then                            ! if the q level has increased, we are on a new level
        popQsLastContIdx(numPopulatedQLvls)=currContIdx-1  ! so we know to store the previous contour as the final of the previous q level
        numPopulatedQLvls=numPopulatedQLvls+1                           ! and increment the number of populated q levels
        popQsFirstContIdx(numPopulatedQLvls)=currContIdx   ! Now store the first contour of the new q level
        orderedPopQsIdxs(numPopulatedQLvls)=newQIdx         ! Since this q level has contours, store its index
        oldQIdx=newQIdx                                         ! store the current q level index so we know when we have moved to a new one
      endif
    enddo

    ! store the last contour as the final of the last q level
    popQsLastContIdx(numPopulatedQLvls)=lastContIdx

     !Note: orderedPopQsIdxs(newQIdx) gives the q level (contsQIdx) of contours.

    ! for populated q levels, create arrays for mapping Qidx -> first/last node index 
    do popQIdx=1,numPopulatedQLvls
      popQsFirstNodeIdx(popQIdx)=firstNodeIdxList(popQsFirstContIdx(popQIdx))
      popQsLastNodeIdx(popQIdx)=lastNodeIdxList(popQsLastContIdx(popQIdx))
    enddo

     !Reset initial layer counters for contours and nodes:
    layerFirstContIdx(iz)=numContoursTot+1
    layerFirstNodeIdx(iz)=numNodesTot+1

     !-----------------------------------------------------------------
     !              Begin a major loop over q levels
     !-----------------------------------------------------------------

    do popQIdx=1,numPopulatedQLvls
      ! for the current populated q level get: 
      
      ! the q level index
      qIdx=orderedPopQsIdxs(popQIdx)

      ! its first/last node indices
      firstNodeIdx=popQsFirstNodeIdx(popQIdx)
      lastNodeIdx=popQsLastNodeIdx(popQIdx)

      ! the number of nodes
      numNodesOnQLvl=dble(lastNodeIdx-firstNodeIdx+1)

      ! Contour segments are binned into boxes to reduce the search costs in surgery.  
      ! The optimal number of boxes is chosen by balancing the cost of boxing (creating boxes 
      ! and binning segments into them) with the cost of surgery search (finding which 
      ! segments to perform surgery on).

      ! any segment which crosses a box boundary will be binned into both boxes
      ! The total number of segments processed in surgery, numSegs, is 
      ! the true number of segments (equal to numNodesOnQLvl since there is one segment per node) 
      ! plus duplicates from segments which cross box boundaries.
      ! we estimate numSegs empirically, based on a realistically complex tests, as: 
      !   (1+3*(avgNodeSpacing/domainArea)*sqrt(numBoxesTot))*numNodesOnQLvl;

      ! boxing costs are estimated as:
      !   (3+15*numSegs/numBoxesTot)*numBoxesTot 

      ! surgery search costs are estimated as:
      !   24*numNodesOnQLvl*numSegs/numBoxesTot, to work out optimal box number;
      
      ! This balance results in a quartic equation for a = sqrt(numBoxesTot)/sqrt(numNodesOnQLvl)
      ! whose solution only depends on r = 3*(avgNodeSpacing/domainArea)*sqrt(numNodesOnQLvl)
      ! Note: (avgNodeSpacing/domainArea) was previously called eps (epsilon?) in this comment

      ! Fortunately, a only varies from 1.129 to 1.265 when r ranges from 0 to infinity.  
      ! Thus we take a = 1.2
      ! which gives 
      !    sqrt(numBoxesTot) = 1.2*sqrt(numNodesOnQLvl),
      ! => numBoxesTot = 1.44*numNodesOnQLvl

      ! For boxing (i.e. gridding the sphere), 
      ! we construct approximately square boxes in the plane of longitude and cartesian z (sine(latitude))
      ! We use circular patches at the poles to circumvent the pole problem.

      ! Since longitude spans a range of 2*pi, and z spans a range of 2, to achieve square boxes, 
      ! there should be pi times as many boxes in longitude as in z
      ! this is satisfied by:
      !   numLongBoxes = 1.2*sqrt(numNodesOnQLvl*pi)
      !   numZBoxes    = 1.2*sqrt(numNodesOnQLvl/pi) = numLongBoxes/pi
      approxNumLongBoxes=1.2d0*sqrt(numNodesOnQLvl*pi)
      approxNumZBoxes=approxNumLongBoxes/pi

      ! Round to the nearest ninteger
      numLongBoxes=nint(approxNumLongBoxes)
      numZBoxes=nint(approxNumZBoxes)

      ! we bound the number of boxes to be at least 1 and at most the number of grid points in that direction
      numLongBoxes=max(min(numLongBoxes,numLongGridPts),1) !enforce 1 <= numLongBoxes <= numLongGridPts
      numZBoxes=max(min(numZBoxes,numLatGridPts),1)        !enforce 1 <= numZBoxes <= numLatGridPts

      ! The total number of boxes (including polar patches) is then:
      numBoxesTot=numLongBoxes*numZBoxes+1 
      !   where the +1 accounts for one polar patch and the other comes from 0 indexing

      ! Calculate the width of the boxes in longitude
      ! and precompute its reciprocal to minimise division (compiler likely does this automatically)
      longBoxWidth=twoPi/numLongBoxes
      recipLongBoxWidth=one/longBoxWidth 

      !----------
      ! I don't understand the following...

      ! fac = tan(latitude)
      ! zbm = cos(latitude)
      ! bwz = sin(latitude)

      ! let zbm be the distance of the polar patch (boundary) from the origin in z
      ! the z thickness of the polar patch is then 1-zbm

      ! let bwz be the radius of the polar patch (boundary) x-y plane
      ! then, since we are working on a unit sphere, 
      ! zbm**2 + bwz**2 = 1

      ! Ensure radius of polar "boxes" equals the z box division:
      ! i.e. set bwz ~= fac
      ! for small fac this is satisfied by
      ! bwz = fac/sqrt(1+fac**2)

      ! get the nominal spacing of the boxes in the cartesian z direction
      ! z=sin(latitude) spans a range of 2 (i.e. [-1,1]) 
      ! which we divide by the number of boxes in z
      fac=two/numZBoxes 

      ! let the distance in z, from the origin, 
      ! at which the polar patch (boundary) is located be zbm
      zbm=one/sqrt(one+fac**2)

      bwz=zbm*fac
      bwzi=one/bwz !precompute the reciprocal to minimise division (compiler likely does this automatically)
      !----------

       !Box all segments [i,next(i)] to reduce search costs in surgery:
      do mb=0,numBoxesTot
        nspb(mb)=0
      enddo
      numSegs=0

      c0=pi+longBoxWidth*dble(numLongBoxes+1)
      !above, pi is the domain half-width (in longitude).
      do ib=firstNodeIdx,lastNodeIdx
        i=next(ib)
         !find range of boxes spanned by segment (i,next(i)):
        adl=c0+lonc(i)
        if (dlonc(i) > zero) then
          mbl1=int(recipLongBoxWidth*adl)-numLongBoxes
          mbl2=int(recipLongBoxWidth*(adl+dlonc(i)))-numLongBoxes
        else
          mbl1=int(recipLongBoxWidth*(adl+dlonc(i)))-numLongBoxes
          mbl2=int(recipLongBoxWidth*adl)-numLongBoxes
        endif
         !mbl1,mbl2 is the longitude range spanned by the segment

        adz=zbm+z(i)
        if (dz(i) > zero) then
          mbz1=int(bwzi*adz+one)
          mbz2=int(bwzi*(adz+dz(i))+one)
        else
          mbz1=int(bwzi*(adz+dz(i))+one)
          mbz2=int(bwzi*adz+one)
        endif
         !mbl1,mbl2 is the z range spanned by the segment

        do mbz=mbz1,mbz2
          if (mbz == 0) then
             !This is the "box" at the south pole:
            nspb(0)=nspb(0)+1
            numSegs=numSegs+1
            loc(numSegs)=0
            list(numSegs)=ib
          else if (mbz == (numZBoxes+1)) then
             !This is the "box" at the north pole:
            nspb(numBoxesTot)=nspb(numBoxesTot)+1
            numSegs=numSegs+1
            loc(numSegs)=numBoxesTot
            list(numSegs)=ib
          else
             !All other boxes:
            mmbz=mbz-1
            do mbl=mbl1,mbl2
              mb=numLongBoxes*mmbz+mod(mbl+(numLongBoxes-1),numLongBoxes)+1

              nspb(mb)=nspb(mb)+1
              numSegs=numSegs+1
              loc(numSegs)=mb
              list(numSegs)=ib
            enddo
          endif
        enddo
      enddo
       !nspb counts number of segments that cross through box mb.
       !numSegs counts the total number of segments crossing all boxes.
       !loc gives the box location of this segment
       !list gives the node before the one at the segment origin (i).

      kb1(0)=1
      do mb=0,numBoxesTot-1
        kb1(mb+1)=kb1(mb)+nspb(mb)
      enddo
      do mb=0,numBoxesTot
        kb2(mb)=kb1(mb)-1
      enddo

      do k=1,numSegs
        mb=loc(k)
        ks=kb2(mb)+1
        node(ks)=list(k)
        kb2(mb)=ks
      enddo

       !With the above, segments [node(kb1(mb)),next(node(kb1(mb))],
       ![node(kb1(mb)+1),next(node(kb1(mb)+1)], ..., 
       ![node(kb2(mb)),next(node(kb2(mb))] cross through box mb.

       !Note: roughly (3+15*numSegs/numBoxesTot)*numBoxesTot operations are required
       !      up to this point (in the loop over q levels) to do
       !      all the segment boxing.

       !Now search for surgery with segments crossing through the box 
       !containing i:
      do ib=firstNodeIdx,lastNodeIdx
        i=next(ib)

        mbz=int(bwzi*(zbm+z(i))+one)
        if (mbz == 0) then
           !This is the "box" at the south pole:
          mb=0
        else if (mbz == (numZBoxes+1)) then
           !This is the "box" at the north pole:
          mb=numBoxesTot
        else
           !All other boxes:
          mb=numLongBoxes*(mbz-1)+int(recipLongBoxWidth*(pi+lonc(i)))+1
           !mb: box containing the node i.
        endif

        do k=kb1(mb),kb2(mb)
          isb=node(k)
          is=next(isb)
           !Segment [is,next(is)] lies in box mb.  
           !Exclude the segments having node i at either endpoint:
          if ((is-ib)*(is-i) == 0) cycle
           !We will see if the node i gets within a distance dm 
           !to the line segment between is and next(is):
          delx=x(is)-x(i)
          dely=y(is)-y(i)
          delz=z(is)-z(i)
          aa=delx*dx(is)+dely*dy(is)+delz*dz(is)

           !Note: roughly 24*numNodesOnQLvl*(numSegs/numBoxesTot) operations are required 
           !up to the following statement (which is rarely satisfied); 
           !this is assumed to be the dominant cost of surgery.  
           !Even if each node i surgically reconnects on average 
           !once, the above estimate holds.

          if (aa*(aa+dsq(is))+d4small < zero) then
             !Passing this condition, a perpendicular can be dropped 
             !onto the line segment from is to next(is).  
             ![Tests show that this is satisfied only 6.4% of the time.]
             !Check distance to line segment:
            cc=(delx*dy(is)-dely*dx(is))**2+ &
               (dely*dz(is)-delz*dy(is))**2+ &
               (delz*dx(is)-delx*dz(is))**2
              !NB: sqrt(cc/dsq(is)) is the distance to the segment and 
              !    dm2=dm**2 below.  This distance must be strictly 
              !    positive, yet less than dm, to permit surgery:
            if (cc*(cc-dm2*dsq(is)) < zero) then
               !Surgery is now possible between node i and the segment 
               ![is,next(is)].
               ![Tests show that this is satisfied only 0.45% of the time.]
               !Move node i to a point midway between i and either is or 
               !next(is), whichever is closer:
              if (dsq(is)+two*aa < zero) then
                dxa=f12*(delx+dx(is))
                dya=f12*(dely+dy(is))
                dza=f12*(delz+dz(is))
                isb=is
                is=next(is)
              else
                dxa=f12*delx
                dya=f12*dely
                dza=f12*delz
              endif
               !Move node i & is to a common node; first deal with node i:
              x(i)=x(i)+dxa
              y(i)=y(i)+dya
              z(i)=z(i)+dza
              dx(i)=dx(i)-dxa
              dy(i)=dy(i)-dya
              dz(i)=dz(i)-dza
              dsq(i)=dx(i)**2+dy(i)**2+dz(i)**2
               !isb becomes the node before i:
              next(isb)=i
              dx(isb)=dx(isb)-dxa
              dy(isb)=dy(isb)-dya
              dz(isb)=dz(isb)-dza
              dsq(isb)=dx(isb)**2+dy(isb)**2+dz(isb)**2

               !Now deal with node is:
              x(is)=x(i)
              y(is)=y(i)
              z(is)=z(i)
              dx(is)=dx(is)+dxa
              dy(is)=dy(is)+dya
              dz(is)=dz(is)+dza
              dsq(is)=dx(is)**2+dy(is)**2+dz(is)**2
               !ib becomes the node before is:
              next(ib)=is
              dx(ib)=dx(ib)+dxa
              dy(ib)=dy(ib)+dya
              dz(ib)=dz(ib)+dza
              dsq(ib)=dx(ib)**2+dy(ib)**2+dz(ib)**2
                !Now update i and look for further possible surgery:
              i=next(ib)
              if (i == ib) exit
               !i = ib when a contour consists of a single node
            endif
          endif
        enddo
      enddo 

       !---------------------------------------------------------------
       !It remains to rebuild the contours using the next() information
       !and to eliminate small contours (having numNodesOnContourList < 5):
      do i=firstNodeIdx,lastNodeIdx
        avail(i)=.true.
      enddo

      jbeg=numContoursTot+1
      do i=firstNodeIdx,lastNodeIdx
        if (avail(i)) then
           !avail(i) is true if node i has not yet been associated with
           !a contour. Start a new contour:
          numContoursTot=numContoursTot+1

           !First point on the contour:
          npd=1
          fac=one/sqrt(x(i)**2+y(i)**2+z(i)**2)
          xd(1)=fac*x(i)
          yd(1)=fac*y(i)
          zd(1)=fac*z(i)

          avail(i)=.false.
          is=next(i)
           !Generate the contour using next() and finish when the
           !original node (i) is reached:
          do while (is .ne. i) 
            npd=npd+1
            fac=one/sqrt(x(is)**2+y(is)**2+z(is)**2)
            xd(npd)=fac*x(is)
            yd(npd)=fac*y(is)
            zd(npd)=fac*z(is)
            avail(is)=.false.
            is=next(is)
          enddo

           !Re-distribute nodes on this contour:
          numNodesOnContourList(numContoursTot)=npd
          if (npd > 3) &
            call renode(xd,yd,zd,npd,xa(numNodesTot+1),ya(numNodesTot+1),za(numNodesTot+1),numNodesOnContourList(numContoursTot))
          if (numNodesOnContourList(numContoursTot) > 3) then
            i1a(numContoursTot)=numNodesTot+1
            numNodesTot=numNodesTot+numNodesOnContourList(numContoursTot)
            i2a(numContoursTot)=numNodesTot
            contsQIdx(numContoursTot)=qIdx
            contourToLayerIdx(numContoursTot)=iz
            do is=i1a(numContoursTot),numNodesTot-1
              nexta(is)=is+1
            enddo
            nexta(numNodesTot)=i1a(numContoursTot)
          else
             !Delete contour if deemed too small (see renode):
            numContoursTot=numContoursTot-1
          endif
        endif
      enddo
       !Now go back and consider the next q level:
    enddo
     !--------------------------------------------------------------------
     !This ends the loop over q levels; surgery is complete for this layer
     !--------------------------------------------------------------------

     !Reset final layer counters for contours and nodes:
    layerLastContIdx(iz)=numContoursTot
    layerLastNodeIdx(iz)=numNodesTot
  endif
enddo
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !          The above ends the major loop over layers
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 !Copy nodes back to those in the argument of the subroutine:
do i=1,numNodesTot
  x(i)=xa(i)
  y(i)=ya(i)
  z(i)=za(i)
  next(i)=nexta(i)
enddo
 !next(i) is not altered until surgery.

do j=1,numContoursTot
  firstNodeIdxList(j)=i1a(j)
  lastNodeIdxList(j)=i2a(j)
enddo

return
end subroutine surgery

!=======================================================================

subroutine renode(xd,yd,zd,npd,xr,yr,zr,npr)
! Re-nodes a single closed contour (xd(i),yd(i),zd(i)), i = 1,...,npd 
!   and returns the new contour in (xr(i),yr(i),zr(i)), i = 1,...,npr

! Note: a closed contour closes on itself

! If npr = 0, the contour is too small and should be removed

implicit none

 !Passed variables:
double precision:: xd(nprm),yd(nprm),zd(nprm),xr(nprm),yr(nprm),zr(nprm)
integer:: npd,npr

 !Local variables:
double precision:: dx(npd),dy(npd),dz(npd),a(npd),b(npd),c(npd),d(npd),e(npd)
double precision:: dsa(npd),dsb(npd)
double precision:: bdif,ww,sum,fac,acc,p,eta,del,xtmp,ytmp,ztmp
double precision:: afac,ax,ay,az,sx,sy,sz
integer:: node(npd+1),nextr(0:npd)
integer:: i,ia,ib,ncorn,im,k,firstNodeIdx,lastNodeIdx,npseg
logical:: corner(npd)

!------------------------------------------------------------------
 !Define the node following a node:
do i=0,npd-1
  nextr(i)=i+1
enddo
nextr(npd)=1

 !Compute the contour increments:
do i=1,npd
  ia=nextr(i)
  dx(i)=xd(ia)-xd(i)
  dy(i)=yd(ia)-yd(i)
  dz(i)=zd(ia)-zd(i)
enddo

do i=1,npd
  dsa(i)=dx(i)**2+dy(i)**2+dz(i)**2
  e(i)=sqrt(dsa(i))
enddo

do ib=1,npd
  i=nextr(ib)
  dsb(i)=dsa(ib)
  a(i)=-dx(ib)
  b(i)=-dy(ib)
  c(i)=-dz(ib)
enddo

ncorn=0
do i=1,npd
  corner(i)=dx(i)*a(i)+dy(i)*b(i)+dz(i)*c(i) > zero
  if (corner(i)) then 
     !Keep track of corner locations for use in renoding below:
    ncorn=ncorn+1
    node(ncorn)=i
     !Set curvature to zero at corners:
    d(i)=zero
  else
    d(i)=(xd(i)*(dy(i)*c(i)-b(i)*dz(i))+ &
          yd(i)*(dz(i)*a(i)-c(i)*dx(i))+ &
          zd(i)*(dx(i)*b(i)-a(i)*dy(i)))/ &
          sqrt((a(i)*dsa(i)-dx(i)*dsb(i))**2+ &
               (b(i)*dsa(i)-dy(i)*dsb(i))**2+ &
               (c(i)*dsa(i)-dz(i)*dsb(i))**2+small3)
  endif
enddo

 !Calculate the cubic interpolation coefficients:
do i=1,npd
  ia=nextr(i)
  dsb(i)=e(i)*(d(ia)+d(i))
  bdif=e(i)*(d(ia)-d(i))
  a(i)=f16*bdif-f12*dsb(i)
  b(i)=f12*(dsb(i)-bdif)
  c(i)=f13*bdif
enddo

!------------------------------------------------------------------------
 !Use the spherical curvature expression (radius of the sphere = ell)
 !to ensure an adequate node density in low curvature regions.
do i=1,npd
  ww=one/(dsa(i)+dmsq)
  dsb(i)=ww*sqrt(elf*dsa(i)+dsb(i)**2)
  dsa(i)=ww*e(i)
enddo
 !NB: elf = 1/ell**2; v(i) = |xx_{i+1}-xx_{i}|**2; e(i)=sqrt{v(i)};
 !    dsb(i)/e(i) = (kappa_{i}+kappa_{i+1})/2; dmsq = (2*dm)**2

 !Re-assign curvature at a node from weighted average on either side
 !(v above is the weight):
do ib=1,npd
  i=nextr(ib)
  d(i)=(dsb(ib)+dsb(i))/(dsa(ib)+dsa(i))
enddo

 !Re-average to get interval value (effectively, four curvature
 !values go into getting the final interval value, dsb(i)):
do i=1,npd
  ia=nextr(i)
  dsb(i)=f12*(d(i)+d(ia))
enddo

 !Compute fractional number of nodes to be placed between old
 !nodes i and i+1:
do i=1,npd
  e(i)=e(i)*min(dmi,densf*sqrt(dsb(i))+dsb(i))
enddo
 !NB: dmi = 2/delta; densf = 1/(amu*sqrt{ell})

sum=zero
do i=1,npd
  sum=sum+e(i)
enddo

 !Number of points on renoded contour:
npr=nint(sum)+1
if (npr < 3) then
   !Contour too small - remove it:
  npr=0
  return
endif

!-------------------------------------------------------------
 !Redistribute nodes making sure to preserve corner locations:
if (ncorn == 0) then
   !No corners - simplest case

   !Make the sum of e(i) equal to npr:
  fac=dble(npr)/sum
  do i=1,npd
    e(i)=fac*e(i)
  enddo

   !The first node is fixed; find the remaining node positions:
  xr(1)=xd(1)
  yr(1)=yd(1)
  zr(1)=zd(1)
  acc=zero
  i=0
  do im=2,npr
    do while (acc < one)
      i=i+1
      acc=acc+e(i)
    enddo
    acc=acc-one
    p=one-acc/e(i)
    eta=p*(a(i)+p*(b(i)+p*c(i)))
    del=f12*p*(one-p)
    xtmp=yd(i)*dz(i)-zd(i)*dy(i)
    ytmp=zd(i)*dx(i)-xd(i)*dz(i)
    ztmp=xd(i)*dy(i)-yd(i)*dx(i)
    afac=sqrt((dx(i)**2+dy(i)**2+dz(i)**2)/(xtmp**2+ytmp**2+ztmp**2))
    ax=afac*xtmp
    ay=afac*ytmp
    az=afac*ztmp
    sx=dy(i)*az-dz(i)*ay
    sy=dz(i)*ax-dx(i)*az
    sz=dx(i)*ay-dy(i)*ax
    xr(im)=xd(i)+p*dx(i)+eta*ax+del*sx
    yr(im)=yd(i)+p*dy(i)+eta*ay+del*sy
    zr(im)=zd(i)+p*dz(i)+eta*az+del*sz
  enddo

else if (ncorn == 1) then
   !A single corner - start new contour at the corner:

   !Make the sum of e(i) equal to npr:
  fac=dble(npr)/sum
  do i=1,npd
    e(i)=fac*e(i)
  enddo

   !The first node (the corner) is fixed; find the remaining nodes:
  i=node(1)
  xr(1)=xd(i)
  yr(1)=yd(i)
  zr(1)=zd(i)
  acc=zero
  i=i-1
  do im=2,npr
    do while (acc < one)
      i=nextr(i)
      acc=acc+e(i)
    enddo
    acc=acc-one
    p=one-acc/e(i)
    eta=p*(a(i)+p*(b(i)+p*c(i)))
    del=f12*p*(one-p)
    xtmp=yd(i)*dz(i)-zd(i)*dy(i)
    ytmp=zd(i)*dx(i)-xd(i)*dz(i)
    ztmp=xd(i)*dy(i)-yd(i)*dx(i)
    afac=sqrt((dx(i)**2+dy(i)**2+dz(i)**2)/(xtmp**2+ytmp**2+ztmp**2))
    ax=afac*xtmp
    ay=afac*ytmp
    az=afac*ztmp
    sx=dy(i)*az-dz(i)*ay
    sy=dz(i)*ax-dx(i)*az
    sz=dx(i)*ay-dy(i)*ax
    xr(im)=xd(i)+p*dx(i)+eta*ax+del*sx
    yr(im)=yd(i)+p*dy(i)+eta*ay+del*sy
    zr(im)=zd(i)+p*dz(i)+eta*az+del*sz
  enddo

else
   !Multiple corners - start new contour at the first corner and
   !preserve locations of all corners:
  npr=0
  node(ncorn+1)=node(1)
  firstNodeIdx=node(1)

  do k=1,ncorn
    lastNodeIdx=node(k+1)
    i=firstNodeIdx
    sum=zero
    do while (i .ne. lastNodeIdx)
      sum=sum+e(i)
      i=nextr(i)
    enddo

    if (sum > small) then
       !Adequate spacing exists between corners to renode this segment

       !Number of points on renoded contour segment:
      npseg=nint(sum)+1

       !Make the sum of e(i) equal to npseg
      fac=dble(npseg)/sum
      i=firstNodeIdx
      do while (i .ne. lastNodeIdx)
        e(i)=fac*e(i)
        i=nextr(i)
      enddo

       !Fix the first node along the segment as the corner:
      xr(npr+1)=xd(firstNodeIdx)
      yr(npr+1)=yd(firstNodeIdx)
      zr(npr+1)=zd(firstNodeIdx)

      if (npseg > 1) then
         !Find the remaining points along this segment:
        acc=zero
        i=firstNodeIdx-1
        do im=npr+2,npr+npseg
          do while (acc < one)
            i=nextr(i)
            acc=acc+e(i)
          enddo
          acc=acc-one
          p=one-acc/e(i)
          eta=p*(a(i)+p*(b(i)+p*c(i)))
          del=f12*p*(one-p)
          xtmp=yd(i)*dz(i)-zd(i)*dy(i)
          ytmp=zd(i)*dx(i)-xd(i)*dz(i)
          ztmp=xd(i)*dy(i)-yd(i)*dx(i)
          afac=sqrt((dx(i)**2+dy(i)**2+dz(i)**2)/(xtmp**2+ytmp**2+ztmp**2))
          ax=afac*xtmp
          ay=afac*ytmp
          az=afac*ztmp
          sx=dy(i)*az-dz(i)*ay
          sy=dz(i)*ax-dx(i)*az
          sz=dx(i)*ay-dy(i)*ax
          xr(im)=xd(i)+p*dx(i)+eta*ax+del*sx
          yr(im)=yd(i)+p*dy(i)+eta*ay+del*sy
          zr(im)=zd(i)+p*dz(i)+eta*az+del*sz
        enddo
      endif

      npr=npr+npseg

    endif
    firstNodeIdx=lastNodeIdx

     !go on and consider the next segment (if any) between corners:
  enddo
endif

 !Normalise coordinates to unit sphere:
do i=1,npr
  fac=one/sqrt(xr(i)**2+yr(i)**2+zr(i)**2)
  xr(i)=fac*xr(i)
  yr(i)=fac*yr(i)
  zr(i)=fac*zr(i)
enddo

return
end subroutine renode

!==========================================================================

 !Main end module
end module contours
