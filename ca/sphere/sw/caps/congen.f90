module congen

!---------------------------------------------------------------------------
!  Converts PV contours to gridded values on an ultra-fine grid of
!  dimensions mgu*nt x mgu*ng, adds the residual q (interpolated to the
!  ultra-fine grid), then creates new contours.  

!  Adapted from ~dgd/cs/spe/sources/mcongen.F on 26/2/13
!  by Stuart King & DG Dritschel @ St Andrews

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

use common

implicit none

 !Grid -> Contour arrays:
double precision:: qa(0:nLatUFGridPts+1,nLongUFGridPts+1)

contains
 
!==========================================================================

subroutine recontour
!a-h = a,b,c,d,e,f,g,h
!o-z = o,p,q,r,s,t,u,v,w,x,y,z
implicit double precision(a-h,o-z)
implicit integer(i-n)

!declaring types for renamed var that fall outside implicit rules

 !Local parameters and variables:
double precision:: qqFullGrid(0:nLatGridPts,nLongGridPts)
double precision:: xTriDiag(nLongGridPts),uTriDiag(nLongGridPts)

!==========================================================
 !What does this section do?? Does it create the fine grid?
 !half grid -> full grid interpolation of qAnomResid
 !Loop over great circles in latitude (only half the longitudes c.f. nLongGridPts=2*nLatGridPts):
do ix=1,nLatGridPts
  ic=ix+nLatGridPts

   !Source vector:
  uTriDiag(1)=f23*(qAnomResid(1,ix)+qAnomResid(1,ic))
  do j=2,nLatGridPts
    uTriDiag(j)=f23*(qAnomResid(j,ix)+qAnomResid(j-1,ix))
  enddo
  uTriDiag(nLatGridPtsPlus1)=f23*(qAnomResid(nLatGridPts,ic)+qAnomResid(nLatGridPts,ix))
  do j=numLatGridPtsPlus2,nLongGridPts
    uTriDiag(j)=f23*(qAnomResid(nLongGridPtsPlusTwo-j,ic)+qAnomResid(nLongGridPtsPlusOne-j,ic))
  enddo

   !Interpolate qr by 4th-order method (periodic):
  xTriDiag(1)=uTriDiag(1)*hTriDiag(1)
  do j=2,nLongGridPts
    xTriDiag(j)=(uTriDiag(j)-f16*xTriDiag(j-1))*hTriDiag(j)
  enddo
  do j=nLongGridPtsMin2,1,-1
    xTriDiag(j)=etd(j)*xTriDiag(j+1)+xTriDiag(j)
  enddo
  xTriDiag(nLongGridPts)=(etd(nLongGridPts)*xTriDiag(1)+xTriDiag(nLongGridPts))*xndeno
  xend=xTriDiag(nLongGridPts)

  do j=1,nLongGridPtsMin1
    xTriDiag(j)=ptd(j)*xend+xTriDiag(j)
  enddo

   !Copy back into full grid array (qq0):
  do j=0,nLatGridPts
    qqFullGrid(j,ix)=xTriDiag(j+1)
  enddo
  qqFullGrid(0,ic)=xTriDiag(1)
  do j=1,nLatGridPts
    qqFullGrid(j,ic)=xTriDiag(nLongGridPtsPlusOne-j)
  enddo

enddo
 !Ends loops over great circles.  Interpolation complete.

 !Obtain unique polar values of qr for use below:
 !qAnomResidSPole is the south pole qr value
qAnomResidSPole=zero
qAnomResidNPole=zero

do ix=1,nLongGridPts
  qAnomResidSPole=qAnomResidSPole+qqFullGrid(0 ,ix)
  qAnomResidNPole=qAnomResidNPole+qqFullGrid(nLatGridPts,ix)
enddo

qAnomResidSPole=qAnomResidSPole/dble(nLongGridPts)
qAnomResidNPole=qAnomResidNPole/dble(nLongGridPts)

do ix=1,nLongGridPts
  qqFullGrid(0 ,ix)=qAnomResidSPole
  qqFullGrid(nLatGridPts,ix)=qAnomResidNPole
enddo

!------------------------------------------------------------
 !Obtain gridded PV from contours (if present):
if (n .gt. 0) then
   !Determine the PV value at the south pole (qSPole):
  qSPole=zero
   !Form great circles to carry out half grid -> full grid
   !interpolation of qs:
  do ix=1,nLatGridPts
    ic=ix+nLatGridPts

     !Source vector:
    uTriDiag(1)=f23*(qSpec(1,ix)+qSpec(1,ic))
    do j=2,nLatGridPts
      uTriDiag(j)=f23*(qSpec(j,ix)+qSpec(j-1,ix))
    enddo
    uTriDiag(nLatGridPtsPlus1)=f23*(qSpec(nLatGridPts,ic)+qSpec(nLatGridPts,ix))
    do j=numLatGridPtsPlus2,nLongGridPts
      uTriDiag(j)=f23*(qSpec(nLongGridPtsPlusTwo-j,ic)+qSpec(nLongGridPtsPlusOne-j,ic))
    enddo

     !Interpolate qs by 4th-order method (periodic):
    xTriDiag(1)=uTriDiag(1)*hTriDiag(1)
    do j=2,nLongGridPts
      xTriDiag(j)=(uTriDiag(j)-f16*xTriDiag(j-1))*hTriDiag(j)
    enddo
    do j=nLongGridPtsMin2,1,-1
      xTriDiag(j)=etd(j)*xTriDiag(j+1)+xTriDiag(j)
    enddo
    xTriDiag(nLongGridPts)=(etd(nLongGridPts)*xTriDiag(1)+xTriDiag(nLongGridPts))*xndeno
  
     !Increment south pole PV value (averaged below):
    qSPole=qSPole+ptd(1)*xTriDiag(nLongGridPts)+xTriDiag(1)
  enddo

   !Obtain average qSPole:
  qSPole=qSPole/dble(nLatGridPts)
  
   !Convert contours to gridded values:
  call con2ufgrid ! Converts PV contours (x,y,z) to gridded values (qa).

   !Bi-linear interpolate qr to the fine grid and add to qa:
  do ix=1,nLongUFGridPts
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)
  
    qa(0,ix)=qa(0,ix)+qAnomResidSPole
    do iy=1,nLatUFGridPts-1
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)

      qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qqFullGrid(iy0,ix0)+w10(iyf,ixf)*qqFullGrid(iy1,ix0) &
                       & +w01(iyf,ixf)*qqFullGrid(iy0,ix1)+w11(iyf,ixf)*qqFullGrid(iy1,ix1)
    enddo
    qa(nLatUFGridPts,ix)=qa(nLatUFGridPts,ix)+qAnomResidNPole
  enddo
   !qAnomResidSPole & qAnomResidNPole are the polar qr values (necessarily uniform).

   !Next adjust qa by a constant so that the PV at the south 
   !pole is the same as that in pvgrid.dat:
  qinc=qSPole-qa(0,1)
   !note: qa does not vary with ix at either pole
  do ix=1,nLongUFGridPts
    do iy=0,nLatUFGridPts
      qa(iy,ix)=qa(iy,ix)+qinc
    enddo
  enddo

else
   !Recall: qa(0:nLatUFGridPts+1,nLongUFGridPts+1)
   !No contours: Bi-linear interpolate qr (which here contains the full q)
   !to the fine grid as qa:
  do ix=1,nLongUFGridPts
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    qa(0,ix)=qAnomResidSPole
    do iy=1,nLatUFGridPts-1
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)
      qa(iy,ix)=w00(iyf,ixf)*qqFullGrid(iy0,ix0)+w10(iyf,ixf)*qqFullGrid(iy1,ix0) &
                                     &   +w01(iyf,ixf)*qqFullGrid(iy0,ix1) &
                                     &   +w11(iyf,ixf)*qqFullGrid(iy1,ix1)
    enddo
    qa(nLatUFGridPts,ix)=qAnomResidNPole
  enddo
   !qAnomResidSPole & qAnomResidNPole are the polar q values (necessarily uniform).

endif

 !Add a periodic column at ix = ntu+1:
ix=nLongUFGridPts+1
do iy=0,nLatUFGridPts
  qa(iy,ix)=qa(iy,1)
enddo

 !Counters for total number of nodes and contours:                                           
npt=0
n=0

 !Generate new contours:
call ufgrid2con !Generates new contours (xd,yd) from the gridded data (qa).

return
end subroutine

!=======================================================================

subroutine con2ufgrid
! Converts PV contours (x,y,z) to gridded values (qa).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local arrays:
double precision:: cx(npt),cy(npt),cz(npt)
double precision:: sq(npt)
integer:: ntc(npt),ilm1(npt)

!----------------------------------------------------------------
 !Initialise crossing information:
do k=1,npt
  ilm1(k)=int(dLongUFInv*(pi+atan2(y(k),x(k))))
enddo

do k=1,npt
  ka=next(k)
  cx(k)=z(k)*y(ka)-y(k)*z(ka)
  cy(k)=x(k)*z(ka)-z(k)*x(ka)
  cz(k)=x(k)*y(ka)-y(k)*x(ka)
  ntc(k)=ilm1(ka)-ilm1(k)
enddo

do k=1,npt
  sig=sign(one,cz(k))
  sq(k)=dq*sig
  ntc(k)=ntc(k)-nLongUFGridPts*((2*ntc(k))/nLongUFGridPts)
  if (sig*dble(ntc(k)) .lt. zero) ntc(k)=-ntc(k)
  if (abs(cz(k)) .gt. zero) then
    cx(k)=cx(k)/cz(k)
    cy(k)=cy(k)/cz(k)
  endif
enddo

!----------------------------------------------------------------------
 !Initialise PV jump array:
do i=1,nLongUFGridPts
  do j=0,nLatUFGridPts+1
    qa(j,i)=zero
  enddo
enddo

 !Determine crossing indices:
do k=1,npt
  if (ntc(k) .ne. 0) then
    jump=sign(1,ntc(k))
    ioff=nLongUFGridPts+ilm1(k)+(1+jump)/2
    ncr=0
    do while (ncr .ne. ntc(k))
      i=1+mod(ioff+ncr,nLongUFGridPts)
      rlatc=dLongUFInv*(hpi+atan(cx(k)*cosLongUF(i)+cy(k)*sinLongUF(i)))
      j=int(rlatc)+1
      p=rlatc-dble(j-1)
      qa(j,i)=  qa(j,i)+(one-p)*sq(k)
      qa(j+1,i)=qa(j+1,i)+    p*sq(k)
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get PV values, at half latitudes, by sweeping through latitudes:
do i=1,nLongUFGridPts
  do j=2,nLatUFGridPts
    qa(j,i)=qa(j,i)+qa(j-1,i)
  enddo
enddo
 !Here, qa(j,i) stands for the PV at latitude j-1/2,
 !from j = 1, ..., ngu.

 !Determine unique polar values:
qasp=zero
qanp=zero
do i=1,nLongUFGridPts
  qasp=qasp+qa(1  ,i)
  qanp=qanp+qa(nLatUFGridPts,i)
enddo
qasp=qasp/dble(nLongUFGridPts)
qanp=qanp/dble(nLongUFGridPts)

 !Average half-grid PV to full grid:
do i=1,nLongUFGridPts
  qa(0,i)=qasp
  do j=1,nLatUFGridPts-1
    qa(j,i)=f12*(qa(j,i)+qa(j+1,i))
  enddo
  qa(nLatUFGridPts,i)=qanp
enddo

return
end subroutine

!=======================================================================
subroutine ufgrid2con
! Generates new contours (xd,yd) from the gridded data (qa).

implicit double precision(a-h,o-z)
implicit integer(i-n)

integer,parameter:: ncrm=3*(npm/2)
! ncrm: max number of contour crossings of a single field
!       level on the finest grid

 !Local Grid -> Contour arrays:
integer(kind=dbleint),parameter:: totalUFGridPts=int(nLongUFGridPts,kind=dbleint)*int(nLatUFGridPts,kind=dbleint)
integer(kind=dbleint):: kob,kib(ncrm)
double precision:: ycr(ncrm),xcr(ncrm)
double precision:: qrx(nLongUFGridPts+1),qry(0:nLatUFGridPts)
double precision:: xd(nprm),yd(nprm),zd(nprm)
integer:: isx(nLongUFGridPts+1),isy(0:nLatUFGridPts)
integer:: icrtab(totalUFGridPts,2)
integer(kind=halfint):: noctab(totalUFGridPts)
logical:: free(ncrm),keep

!--------------------------------------------------------
 !Recall: qa is the residual PV interpolated to the ultra-fine grid
 !and qa(0:nLatUFGridPts+1,nLongUFGridPts+1)

 !First get the beginning and ending contour (q) levels:
qamax=max(qa(0,1),qa(nLatUFGridPts,1)) !max of poles
qamin=min(qa(0,1),qa(nLatUFGridPts,1)) !min of poles

 !find the min and max of residual PV on the ultra-fine grid
 !(qa is uniform at the extended edges iy = 0 and ngu)
 !since column nLongUFGridPts+1 is periodic copy of column 1
 !so the loop can run to nLongUFGridPts only and not nLongUFGridPts+1
do ix=1,nLongUFGridPts
  do iy=1,nLatUFGridPts-1 !does not include poles as already considered
    qamax=max(qamax,qa(iy,ix))
    qamin=min(qamin,qa(iy,ix))
  enddo
enddo

!I think these lines center the contour levels about zero
!the +1/2 is so that int gives closest int not floor
levbeg=int((qoff+qamin)*dqi+f12)+1 
levend=int((qoff+qamax)*dqi+f12)

 !Return if no levels to process:
if (levbeg .gt. levend) return 

 !Loop over contour levels and process:
do lev=levbeg,levend
   !Integer index giving contour level (index):                                                        
  indq=lev-nlevm+(lev-1)/nlevm-1
   !as lev could be from 1 to 2nlevm
   !indq runs from -nlevm to int(nlevm +1-1/nlevm)=nlevm

   !Counter for total number of grid line crossings:
  ncr=0

   !Contour level (q) being sought: 
  qtmp=qlev(lev)

   !Initialise number of crossings per box:
  do kob=1,totalUFGridPts
    noctab(kob)=0
  enddo

   !Find x grid line crossings first:
  do ix=1,nLongUFGridPts

    do iy=0,nLatUFGridPts
      qry(iy)=qa(iy,ix)-qtmp    !difference of grid point value and contour value
      isy(iy)=sign(one,qry(iy)) !sign of the difference: +1 if qa >= qtmp; -1 if qa < qtmp
    enddo

    do iy=0,nLatUFGridPts-1
      if (isy(iy) .ne. isy(iy+1)) then ! criterion for crossing 
        ncr=ncr+1
        inc=(1-isy(iy))/2           !inc = 1 if qa < qtmp; inc = 0 if qa >= qtmp
        kib(ncr)=iy+1+ibx(ix,inc)   !if inc=0, kib = iy+1+ibx(ix,0); if inc=1, kib = iy+1+ibx(ix,1)
        kob=iy+1+ibx(ix,1-inc)      !if inc=0, kob = iy+1+ibx(ix,1); if inc=1, kob = iy+1+ibx(ix,0)
        noctab(kob)=noctab(kob)+1   !increment number of crossings in box kob
        icrtab(kob,noctab(kob))=ncr ! note where crossing number ncr occurs
        
        xcr(ncr)=xgu(ix)                                   !coordinate of crossed x grid line
        ycr(ncr)=ygu(iy)-glyu*qry(iy)/(qry(iy+1)-qry(iy))  !coordinate of crossed y grid line
      endif
    enddo

  enddo

!   Above, kib = grid box into which the contour (containing icr) is going
!          kob =   "   "  out of "    "     "         "       "    " coming
!     [kob -> icr -> kib:  icr lies at the boundary between kob & kib]

 !Find y grid line crossings next (no crossings can occur at iy=0,ngu):
  do iy=1,nLatUFGridPts-1
    ygt=ygu(iy)

    do ix=1,nLongUFGridPts+1
      qrx(ix)=qa(iy,ix)-qtmp
      isx(ix)=sign(one,qrx(ix))
    enddo

    do ix=1,nLongUFGridPts
      if (isx(ix) .ne. isx(ix+1)) then
        ncr=ncr+1
        inc=(1-isx(ix))/2
        kib(ncr)=ibx(ix,1)+iy+1-inc   !in x crossings: kib(ncr)=ibx(ix,inc)+iy+1
        kob=ibx(ix,1)+iy+inc          !in x crossings: kob=ibx(ix,1-inc)+iy+1
        noctab(kob)=noctab(kob)+1
        icrtab(kob,noctab(kob))=ncr
        ycr(ncr)=ygt
        xcr(ncr)=xgu(ix)-glxu*qrx(ix)/(qrx(ix+1)-qrx(ix))
      endif
    enddo

  enddo

!------------------------------------------------------------------------
   !Now re-build contours - converting to spherical geometry:
  do icr=1,ncr
    free(icr)=.true. !mark all crossings as not associated with a contour
  enddo

  do icr=1,ncr !for each crossing
    if (free(icr)) then 
       !A new contour (indexed n) starts here:
      n=n+1 !increment global contour counter
      ind(n)=indq !store the contour level index
      ibeg=npt+1 !the node index of the first point on this contour
      i1(n)=ibeg !store the first node index for this contour

       !First point on the contour:
       !in cartesian coords
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
       !Why last crossing?

      icrn=icrtab(kob,noc)
       !icrn gives the next point after icr (icrn is leaving box kob)
       !next point = the last crossing in box kob

      do while (icrn .ne. icr) !keep going until we return to starting point icr
        noctab(kob)=noc-1 ! decrement number of crossings in box kob
         !noctab is usually zero now except for boxes with a
         !maximum possible 2 crossings
        npd=npd+1

        coslat=cos(ycr(icrn))
        xd(npd)=coslat*cos(xcr(icrn))
        yd(npd)=coslat*sin(xcr(icrn))
        zd(npd)=sin(ycr(icrn))

        free(icrn)=.false. ! mark icrn as associated with a contour

        kob=kib(icrn)       !box the contour is entering next
        noc=noctab(kob)     !number of crossings in this box
        !Use last crossing (noc) in this box (kob) as the next node:
        icrn=icrtab(kob,noc)
      enddo

       !Re-distribute nodes on this contour 3 times to reduce complexity:
      keep=.false.
      do
        call renode(xd,yd,zd,npd,x(ibeg),y(ibeg),z(ibeg),np(n))
         !Delete contour if deemed too small (see renode):
        if (np(n) .eq. 0) exit
        call renode(x(ibeg),y(ibeg),z(ibeg),np(n),xd,yd,zd,npd)
         !Delete contour if deemed too small (see renode):
        if (npd .eq. 0) exit
        call renode(xd,yd,zd,npd,x(ibeg),y(ibeg),z(ibeg),np(n))
         !Delete contour if deemed too small (see renode):
        if (np(n) .eq. 0) exit
         !Contour is big enough to keep:
        keep=.true.
        exit
      enddo

      if (keep) then 
        npt=npt+np(n) !increment total node counter
        iend=ibeg+np(n)-1 !last node index on this contour
        i2(n)=iend
        do i=ibeg,iend-1
          next(i)=i+1 !link nodes on this contour
        enddo
        next(iend)=ibeg !link last node to first node
      else
        n=n-1 !undo the increment of global contour counter
      endif

      free(icr)=.false. ! mark icr as associated with a contour
    endif
  enddo

enddo
!End of loop over contour levels
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

return
end subroutine

!=======================================================================

 !Main end module
end module
