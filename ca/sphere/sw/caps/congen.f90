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
double precision:: qAdiab(0:nLatUFGridPts+1,nLongUFGridPts+1)

contains
 
!==========================================================================

subroutine recontour
!a-h = a,b,c,d,e,f,g,h
!o-z = o,p,q,r,s,t,u,v,w,x,y,z
implicit double precision(a-h,o-z)
implicit integer(i-n)

!declaring types for renamed var that fall outside implicit rules

 !Local parameters and variables:
double precision:: qAnomFull0idx(0:nLatGridPts,nLongGridPts)
double precision:: xTriDiag(nLongGridPts),uTriDiag(nLongGridPts)

!==========================================================
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
    qAnomFull0idx(j,ix)=xTriDiag(j+1)
  enddo
  qAnomFull0idx(0,ic)=xTriDiag(1)
  do j=1,nLatGridPts
    qAnomFull0idx(j,ic)=xTriDiag(nLongGridPtsPlusOne-j)
  enddo

enddo
 !Ends loops over great circles.  Interpolation complete.

 !Obtain unique polar values of qr for use below:
 !qAnomResidSPole is the south pole qr value
qAnomResidSPole=zero
qAnomResidNPole=zero

do ix=1,nLongGridPts
  qAnomResidSPole=qAnomResidSPole+qAnomFull0idx(0 ,ix)
  qAnomResidNPole=qAnomResidNPole+qAnomFull0idx(nLatGridPts,ix)
enddo

qAnomResidSPole=qAnomResidSPole/dble(nLongGridPts)
qAnomResidNPole=qAnomResidNPole/dble(nLongGridPts)

do ix=1,nLongGridPts
  qAnomFull0idx(0 ,ix)=qAnomResidSPole
  qAnomFull0idx(nLatGridPts,ix)=qAnomResidNPole
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
  call con2ufgrid

   !Bi-linear interpolate qr to the fine grid and add to qa:
  do ix=1,nLongUFGridPts
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)
  
    qAdiab(0,ix)=qAdiab(0,ix)+qAnomResidSPole
    do iy=1,nLatUFGridPts-1
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)

      qAdiab(iy,ix)=qAdiab(iy,ix)+w00(iyf,ixf)*qAnomFull0idx(iy0,ix0)+w10(iyf,ixf)*qAnomFull0idx(iy1,ix0) &
                       & +w01(iyf,ixf)*qAnomFull0idx(iy0,ix1)+w11(iyf,ixf)*qAnomFull0idx(iy1,ix1)
    enddo
    qAdiab(nLatUFGridPts,ix)=qAdiab(nLatUFGridPts,ix)+qAnomResidNPole
  enddo
   !qAnomResidSPole & qAnomResidNPole are the polar qr values (necessarily uniform).

   !Next adjust qa by a constant so that the PV at the south 
   !pole is the same as that in pvgrid.dat:
  qinc=qSPole-qAdiab(0,1)
   !note: qa does not vary with ix at either pole
  do ix=1,nLongUFGridPts
    do iy=0,nLatUFGridPts
      qAdiab(iy,ix)=qAdiab(iy,ix)+qinc
    enddo
  enddo

else

   !No contours: interpolate qr (which here contains the full q)
   !to the fine grid as qa:
  do ix=1,nLongUFGridPts
    ixf=ixfw(ix)
    ix0=ix0w(ix)
    ix1=ix1w(ix)

    qAdiab(0,ix)=qAnomResidSPole
    do iy=1,nLatUFGridPts-1
      iyf=iyfw(iy)
      iy0=iy0w(iy)
      iy1=iy1w(iy)
      qAdiab(iy,ix)=w00(iyf,ixf)*qAnomFull0idx(iy0,ix0)+w10(iyf,ixf)*qAnomFull0idx(iy1,ix0) &
                                     &   +w01(iyf,ixf)*qAnomFull0idx(iy0,ix1) &
                                     &   +w11(iyf,ixf)*qAnomFull0idx(iy1,ix1)
    enddo
    qAdiab(nLatUFGridPts,ix)=qAnomResidNPole
  enddo
   !qAnomResidSPole & qAnomResidNPole are the polar q values (necessarily uniform).

endif

 !Add a periodic column at ix = ntu+1:
ix=nLongUFGridPts+1
do iy=0,nLatUFGridPts
  qAdiab(iy,ix)=qAdiab(iy,1)
enddo

 !Counters for total number of nodes and contours:                                           
npt=0
n=0

 !Generate new contours:
call ufgrid2con

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
  ilm1(k)=int(dlui*(pi+atan2(y(k),x(k))))
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
    qAdiab(j,i)=zero
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
      rlatc=dlui*(hpi+atan(cx(k)*clonu(i)+cy(k)*slonu(i)))
      j=int(rlatc)+1
      p=rlatc-dble(j-1)
      qAdiab(j,i)=  qAdiab(j,i)+(one-p)*sq(k)
      qAdiab(j+1,i)=qAdiab(j+1,i)+    p*sq(k)
      ncr=ncr+jump
    enddo
  endif
enddo

 !Get PV values, at half latitudes, by sweeping through latitudes:
do i=1,nLongUFGridPts
  do j=2,nLatUFGridPts
    qAdiab(j,i)=qAdiab(j,i)+qAdiab(j-1,i)
  enddo
enddo
 !Here, qa(j,i) stands for the PV at latitude j-1/2,
 !from j = 1, ..., ngu.

 !Determine unique polar values:
qasp=zero
qanp=zero
do i=1,nLongUFGridPts
  qasp=qasp+qAdiab(1  ,i)
  qanp=qanp+qAdiab(nLatUFGridPts,i)
enddo
qasp=qasp/dble(nLongUFGridPts)
qanp=qanp/dble(nLongUFGridPts)

 !Average half-grid PV to full grid:
do i=1,nLongUFGridPts
  qAdiab(0,i)=qasp
  do j=1,nLatUFGridPts-1
    qAdiab(j,i)=f12*(qAdiab(j,i)+qAdiab(j+1,i))
  enddo
  qAdiab(nLatUFGridPts,i)=qanp
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
integer(kind=dbleint),parameter:: ntng=int(nLongUFGridPts,kind=dbleint)*int(nLatUFGridPts,kind=dbleint)
integer(kind=dbleint):: kob,kib(ncrm)
double precision:: ycr(ncrm),xcr(ncrm)
double precision:: qrx(nLongUFGridPts+1),qry(0:nLatUFGridPts)
double precision:: xd(nprm),yd(nprm),zd(nprm)
integer:: isx(nLongUFGridPts+1),isy(0:nLatUFGridPts)
integer:: icrtab(ntng,2)
integer(kind=halfint):: noctab(ntng)
logical:: free(ncrm),keep

!--------------------------------------------------------
 !First get the beginning and ending contour levels:
qamax=max(qAdiab(0,1),qAdiab(nLatUFGridPts,1))
qamin=min(qAdiab(0,1),qAdiab(nLatUFGridPts,1))
 !(qa is uniform at the extended edges iy = 0 and ngu)
do ix=1,nLongUFGridPts
  do iy=1,nLatUFGridPts-1
    qamax=max(qamax,qAdiab(iy,ix))
    qamin=min(qamin,qAdiab(iy,ix))
  enddo
enddo

levbeg=int((qoff+qamin)*dqi+f12)+1
levend=int((qoff+qamax)*dqi+f12)

 !Return if no levels to process:
if (levbeg .gt. levend) return 

 !Loop over contour levels and process:
do lev=levbeg,levend
   !Integer index giving contour level:                                                        
  indq=lev-nlevm+(lev-1)/nlevm-1
   !Counter for total number of grid line crossings:
  ncr=0

   !Contour level being sought:
  qtmp=qlev(lev)

   !Initialise number of crossings per box:
  do kob=1,ntng
    noctab(kob)=0
  enddo

   !Find x grid line crossings first:
  do ix=1,nLongUFGridPts
    xgt=xgu(ix)

    do iy=0,nLatUFGridPts
      qry(iy)=qAdiab(iy,ix)-qtmp
      isy(iy)=sign(one,qry(iy))
    enddo

    do iy=0,nLatUFGridPts-1
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
  do iy=1,nLatUFGridPts-1
    ygt=ygu(iy)

    do ix=1,nLongUFGridPts+1
      qrx(ix)=qAdiab(iy,ix)-qtmp
      isx(ix)=sign(one,qrx(ix))
    enddo

    do ix=1,nLongUFGridPts
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
  do icr=1,ncr
    free(icr)=.true.
  enddo

  do icr=1,ncr
    if (free(icr)) then
       !A new contour (indexed n) starts here:
      n=n+1
      ind(n)=indq
      ibeg=npt+1
      i1(n)=ibeg

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
        npt=npt+np(n)
        iend=ibeg+np(n)-1
        i2(n)=iend
        do i=ibeg,iend-1
          next(i)=i+1
        enddo
        next(iend)=ibeg
      else
        n=n-1
      endif

      free(icr)=.false.
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
