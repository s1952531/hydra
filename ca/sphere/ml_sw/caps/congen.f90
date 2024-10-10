module congen

!---------------------------------------------------------------------------
!         Converts PV contours to gridded values on an
!         ultra-fine grid of dimensions mgu*nt x mgu*ng (periodic 
!         in x but free slip walls in y), adds the residual q 
!         (interpolated to the ultra-fine grid), then creates new contours.  

!         This routine processes one contour level at a time.

!         This avoids storing large arrays in memory.  After
!         running this programme.

!         Adapted from ~dgd/cs/spe/sources/mcongen.F on 15/2/13
!         by Stuart King @ St Andrews

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

use common

implicit none

 !Common parameters for the congen module:
integer,parameter:: nplm=npm,nlm=nplm/20+nplm/200
integer,parameter:: nlevm=2000
!      nplm:     maximum number of  nodes  
!      nlm:      maximum number of contours

integer:: i1d(nlm),i2d(nlm),npd(nlm)

 !Grid -> Contour arrays:
double precision:: qa(0:ngu+1,ntu+1)

 !Basic parameters:
integer:: na,nd,npta,nptd

!==========================================================================

subroutine recontour

implicit double precision(a-h,o-z)
implicit integer(i-n)

integer,parameter:: ntm1=nt-1,ntm2=nt-2,ntp1=nt+1,ntp2=nt+2
integer,parameter:: ngp1=ng+1,ngp2=ng+2

 !Local parameters and variables:
double precision:: qq0(0:ng,nt)
double precision:: etd(nt),htd(nt),ptd(nt)
double precision:: xtd(nt),utd(nt)
integer:: i1(nlm),i2(nlm),np(nlm),lay(nlm)
integer:: jl1(nlay),jl2(nlay),il1(nlay),il2(nlay)
logical:: contours


!--------------------------------------------------------------------
 !Allow for the possibility of no contours:
contours=(n .gt. 0)
if (contours) then

   !Find beginning and ending nodes and contours in each layer:
  if (nlay .eq. 1) then
     !Only one layer present:
    jl1(1)=1
    jl2(1)=n
    il1(1)=1
    il2(1)=npt
  else
    !Multiple layers present:
    do lz=1,nlay
      jl2(lz)=0
    enddo
    do j=1,n
      il2(lay(j))=i2(j)
      jl2(lay(j))=j
    enddo
    ibeg=1
    jbeg=1
    do lz=1,nlay
      il1(lz)=ibeg
      ibeg=il2(lz)+1
      jl1(lz)=jbeg
      jbeg=jl2(lz)+1
    enddo
  endif

endif

 !Counters for total number of nodes and contours:
npta=0
na=0

!==========================================================
 !Process data layer by layer:
do lz=1,nlay
   !Form great circles to carry out half grid -> full grid
   !interpolation of qd -> qq0:

   !Initialise periodic tridiagonal problem:
  htd(1)=one
  ptd(1)=-f16*htd(1)
  etd(1)=ptd(1)

  do j=2,nt
    htd(j)=one/(one+f16*etd(j-1))
    ptd(j)=-f16*ptd(j-1)*htd(j)
    etd(j)=-f16*htd(j)
  enddo

  ptd(ntm1)=etd(ntm1)+ptd(ntm1)
  do j=ntm2,1,-1
    ptd(j)=etd(j)*ptd(j+1)+ptd(j)
  enddo

  xndeno=one/(one-etd(nt)*ptd(1)-ptd(nt))

   !Loop over great circles in latitude (only half the longitudes):
  do ix=1,ng
    ic=ix+ng

     !Source vector:
    utd(1)=f23*(qd(1,ix)+qd(1,ic))
    do j=2,ng
      utd(j)=f23*(qd(j,ix)+qd(j-1,ix))
    enddo
    utd(ngp1)=f23*(qd(ng,ic)+qd(ng,ix))
    do j=ngp2,nt
      utd(j)=f23*(qd(ntp2-j,ic)+qd(ntp1-j,ic))
    enddo

     !Interpolate qd by 4th-order method (periodic):
    xtd(1)=utd(1)*htd(1)
    do j=2,nt
      xtd(j)=(utd(j)-f16*xtd(j-1))*htd(j)
    enddo
    do j=ntm2,1,-1
      xtd(j)=etd(j)*xtd(j+1)+xtd(j)
    enddo
    xtd(nt)=(etd(nt)*xtd(1)+xtd(nt))*xndeno
    xend=xtd(nt)

    do j=1,ntm1
      xtd(j)=ptd(j)*xend+xtd(j)
    enddo

     !Copy back into full grid array (qq0):
    do j=0,ng
      qq0(j,ix)=xtd(j+1)
    enddo
    qq0(0,ic)=xtd(1)
    do j=1,ng
      qq0(j,ic)=xtd(ntp1-j)
    enddo

  enddo
   !Ends loops over great circles.  Interpolation complete.

   !Obtain unique polar values of qd for use below:
  qdsp=zero
  qdnp=zero
  do ix=1,nt
    qdsp=qdsp+qq0(0 ,ix)
    qdnp=qdnp+qq0(ng,ix)
  enddo
  qdsp=qdsp/dble(nt)
  qdnp=qdnp/dble(nt)
  do ix=1,nt
    qq0(0 ,ix)=qdsp
    qq0(ng,ix)=qdnp
  enddo

  !------------------------------------------------------------
   !Process contours layer by layer (if multiple layers occur):
  if (contours) then
    ioff=il1(lz)-1
    nptd=il2(lz)-ioff

     !Get relative contour indices:
    joff=jl1(lz)-1
    nd  =jl2(lz)-joff
  
    do j=1,nd
      jj=joff+j
      i1d(j)=i1(jj)-ioff
      i2d(j)=i2(jj)-ioff
      npd(j)=np(jj)
    enddo

     !Define next(i), which  gives the node following node i:
    do i=1,nptd-1
      next(i)=i+1
    enddo
    do j=1,nd
      next(i2d(j))=i1d(j)
    enddo

     !Determine the PV value at the south pole (qsp):
    qsp=zero
     !Form great circles to carry out half grid -> full grid
     !interpolation of qs:
    do ix=1,ng
      ic=ix+ng

       !Source vector:
      utd(1)=f23*(qs(1,ix)+qs(1,ic))
      do j=2,ng
        utd(j)=f23*(qs(j,ix)+qs(j-1,ix))
      enddo
      utd(ngp1)=f23*(qs(ng,ic)+qs(ng,ix))
      do j=ngp2,nt
        utd(j)=f23*(qs(ntp2-j,ic)+qs(ntp1-j,ic))
      enddo

       !Interpolate qs by 4th-order method (periodic):
      xtd(1)=utd(1)*htd(1)
      do j=2,nt
        xtd(j)=(utd(j)-f16*xtd(j-1))*htd(j)
      enddo
      do j=ntm2,1,-1
        xtd(j)=etd(j)*xtd(j+1)+xtd(j)
      enddo
      xtd(nt)=(etd(nt)*xtd(1)+xtd(nt))*xndeno
  
       !Increment south pole PV value (averaged below):
      qsp=qsp+ptd(1)*xtd(nt)+xtd(1)
    enddo

     !Obtain average qsp:
    qsp=qsp/dble(ng)
  
     !Convert contours to gridded values:
    call con2ugrid

     !Bi-linear interpolate qs to the fine grid and add to qa:
    do ix=1,ntu
      ixf=ixfw(ix)
      ix0=ix0w(ix)
      ix1=ix1w(ix)
  
      qa(0,ix)=qa(0,ix)+qdsp
      do iy=1,ngu-1
        iyf=iyfw(iy)
        iy0=iy0w(iy)
        iy1=iy1w(iy)
  
        qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qq0(iy0,ix0)+w10(iyf,ixf)*qq0(iy1,ix0) &
                         & +w01(iyf,ixf)*qq0(iy0,ix1)+w11(iyf,ixf)*qq0(iy1,ix1)
      enddo
      qa(ngu,ix)=qa(ngu,ix)+qdnp
    enddo
     !qdsp & qdnp are the polar qd values (necessarily uniform).

     !Next adjust qa by a constant so that the PV at the south 
     !pole is the same as that in pvgrid.dat:
    qinc=qsp-qa(0,1)
     !note: qa does not vary with ix at either pole
    do ix=1,ntu
      do iy=0,ngu
        qa(iy,ix)=qa(iy,ix)+qinc
      enddo
    enddo

  else

     !No contours: interpolate qd (which here contains the full q) to the fine grid as qa:
    do ix=1,ntu
      ixf=ixfw(ix)
      ix0=ix0w(ix)
      ix1=ix1w(ix)

      qa(0,ix)=qdsp
      do iy=1,ngu-1
        iyf=iyfw(iy)
        iy0=iy0w(iy)
        iy1=iy1w(iy)

        qa(iy,ix)=w00(iyf,ixf)*qq0(iy0,ix0)+w10(iyf,ixf)*qq0(iy1,ix0) &
                                       &   +w01(iyf,ixf)*qq0(iy0,ix1) &
                                       &   +w11(iyf,ixf)*qq0(iy1,ix1)
      enddo
      qa(ngu,ix)=qdnp
    enddo
     !qdsp & qdnp are the polar q values (necessarily uniform).

  endif

   !Add a periodic column at ix = ntu+1:
  ix=ntu+1
  do iy=0,ngu
    qa(iy,ix)=qa(iy,1)
  enddo

   !Generate new contours in layer lz:
  call ugrid2con(lz)

enddo
!============================================================================
 !Ends loop over layers.

 !Close files:
close(10)
close(11)
close(12)
close(21)
close(22)

open(13,file='head_congen.dat',status='unknown')
 !Write one line header for congen.dat:
write(13,'(i7,1x,i8,1x,f12.5,1x,f16.12)') na,npta,t,dq
close(13)

write(*,*) 
write(*,*) ' All done.  To create pvcont.dat, type:'
write(*,*) ' cat head_congen.dat all_contours.dat all_nodes.dat > pvcont.dat'
write(*,*) ' rm  head_congen.dat all_contours.dat all_nodes.dat'

return
end subroutine

!=======================================================================

subroutine con2ugrid
! Converts PV contours (x,y,z) to gridded values (qa).

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Local arrays:
integer:: ntc(nplm),ilm1(nplm)

!----------------------------------------------------------------
 !Initialise crossing information:
do k=1,nptd
  ilm1(k)=int(dlui*(pi+atan2(yd(k),xd(k))))
enddo

do k=1,nptd
  ka=next(k)
  v(k)=xd(k)*yd(ka)-yd(k)*xd(ka)
  c(k)=zd(k)*yd(ka)-yd(k)*zd(ka)
  d(k)=xd(k)*zd(ka)-zd(k)*xd(ka)
  ntc(k)=ilm1(ka)-ilm1(k)
enddo

do k=1,nptd
  sig=sign(one,v(k))
  a(k)=dq*sig
  ntc(k)=ntc(k)-ntu*((2*ntc(k))/ntu)
  if (sig*dble(ntc(k)) .lt. zero) ntc(k)=-ntc(k)
  if (abs(v(k)) .gt. zero) then
    c(k)=c(k)/v(k)
    d(k)=d(k)/v(k)
  endif
enddo

!----------------------------------------------------------------------
 !Initialise PV jump array:
do i=1,ntu
  do j=0,ngu+1
    qa(j,i)=zero
  enddo
enddo

 !Determine crossing indices:
do k=1,nptd
  if (ntc(k) .ne. 0) then
    jump=sign(1,ntc(k))
    ioff=ntu+ilm1(k)+(1+jump)/2
    ncr=0
    do while (ncr .ne. ntc(k))
      i=1+mod(ioff+ncr,ntu)
      rlatc=dlui*(hpi+atan(c(k)*clonu(i)+d(k)*slonu(i)))
      j=int(rlatc)+1
      p=rlatc-dble(j-1)
      qa(j,i)=  qa(j,i)+(one-p)*a(k)
      qa(j+1,i)=qa(j+1,i)+    p*a(k)
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
qasp=zero
qanp=zero
do i=1,ntu
  qasp=qasp+qa(1  ,i)
  qanp=qanp+qa(ngu,i)
enddo
qasp=qasp/dble(ntu)
qanp=qanp/dble(ntu)

 !Average half-grid PV to full grid:
do i=1,ntu
  qa(0,i)=qasp
  do j=1,ngu-1
    qa(j,i)=f12*(qa(j,i)+qa(j+1,i))
  enddo
  qa(ngu,i)=qanp
enddo

return
end subroutine

!=======================================================================
subroutine ugrid2con(lz)
! Generates new contours (xd,yd) from the gridded data (qa) 
! in layer lz.

implicit double precision(a-h,o-z)
implicit integer(i-n)

integer,parameter:: nreno=2
! nreno: number of times renode is called to reduce point
!        density on contours.

integer,parameter:: ncrm=3*nplm/2
! ncrm: max number of contour crossings of a single field
!       level on the finest grid

 !Local Grid -> Contour arrays:
integer,parameter:: ntng=ntu*ngu
double precision:: ycr(ncrm),xcr(ncrm)
double precision:: qdx(ntu+1),qdy(0:ngu)
integer:: kib(ncrm),kob(ncrm),ipo(ncrm)
integer:: isx(ntu+1),isy(0:ngu)
integer:: icrtab(ntng,2)
integer:: noctab(ntng)
logical:: free(ncrm)

!--------------------------------------------------------
 !First get the beginning and ending contour levels:
qamax=max(qa(0,1),qa(ngu,1))
qamin=min(qa(0,1),qa(ngu,1))
 !(qa is uniform at the extended edges iy = 0 and ngu)
do ix=1,ntu
  do iy=1,ngu-1
    qamax=max(qamax,qa(iy,ix))
    qamin=min(qamin,qa(iy,ix))
  enddo
enddo

levbeg=int((qoff+qamin)*dqi+f12)+1
levend=int((qoff+qamax)*dqi+f12)

inda=0

if (levbeg .le. levend) then
 !Loop over contour levels and process:
  do lev=levbeg,levend
     !Counter for total number of grid line crossings:
    ncr=0

     !Contour level being sought:
    qtmp=qlev(lev)

     !Find x grid line crossings first:
    do ix=1,ntu
      xgt=xgu(ix)

      do iy=0,ngu
        qdy(iy)=qa(iy,ix)-qtmp
        isy(iy)=sign(one,qdy(iy))
      enddo

      do iy=0,ngu-1
        if (isy(iy) .ne. isy(iy+1)) then
          ncr=ncr+1
          inc=(1-isy(iy))/2
          kib(ncr)=iy+1+ibx(ix,inc)
          kob(ncr)=iy+1+ibx(ix,1-inc)
          xcr(ncr)=xgt
          ycr(ncr)=ygu(iy)-dlu*qdy(iy)/(qdy(iy+1)-qdy(iy))
        endif
      enddo

    enddo

!     Above, kib = grid box into which the contour (containing icr) is going
!            kob =   "   "  out of "    "     "         "       "    " coming
!       [kob -> icr -> kib:  icr lies at the boundary between kob & kib]

   !Find y grid line crossings next (no crossings can occur at iy=0,ngu):
    do iy=1,ngu-1
      ygt=ygu(iy)

      do ix=1,ntu+1
        qdx(ix)=qa(iy,ix)-qtmp
        isx(ix)=sign(one,qdx(ix))
      enddo

      do ix=1,ntu
        if (isx(ix) .ne. isx(ix+1)) then
          ncr=ncr+1
          inc=(1-isx(ix))/2
          kib(ncr)=ibx(ix,1)+iy+1-inc
          kob(ncr)=ibx(ix,1)+iy+inc
          ycr(ncr)=ygt
          xcr(ncr)=xgu(ix)-dlu*qdx(ix)/(qdx(ix+1)-qdx(ix))
        endif
      enddo

    enddo

     !Initialise number of crossings per box:
    do i=1,ncr
      noctab(kob(i))=0
    enddo

    do icr=1,ncr
        !icr is the index of the current crossing at level lev.
       k=kob(icr)
        !Accumulate number of crossings in this box:
       noctab(k)=noctab(k)+1
        !Assign crossing to box, permitting 2 crossings:
       icrtab(k,noctab(k))=icr
    enddo

    do icr=1,ncr
      k=kib(icr)
      noc=noctab(k)
       !Use last crossing in this box as the next node:
      kob(icr)=icrtab(k,noc)
       !kob(icr) now gives the next point after icr
      noctab(k)=noc-1
       !This will normally be zero, except for boxes with 2 crossings;
       !this allows a second use of this box.
    enddo

!-----------------
     !Now re-build contours:
    j=0
    i=0
  
    do icr=1,ncr
      free(icr)=.true.
    enddo

    do icr=1,ncr
      if (free(icr)) then
         !A new contour (j) starts here:
        i=i+1
        j=j+1
        i1d(j)=i
        ipo(i)=icr
        icrn=kob(icr)

        do while (icrn .ne. icr)
           !Find remaining points on contour j:
          i=i+1
          ipo(i)=icrn
          free(icrn)=.false.
          icrn=kob(icrn)
        enddo
        i2d(j)=i
        npd(j)=i2d(j)-i1d(j)+1
      endif
    enddo

    nptd=0
    nrem=0
    ndt=j

    do j=1,ndt
      if (npd(j) .lt. 5) then
        nrem=nrem+1
      else
        ibeg=nptd+1
        ioff=ibeg-i1d(j)
        do i=i1d(j),i2d(j)
          icr=ipo(i)
          xd(ioff+i)=xcr(icr)
          yd(ioff+i)=ycr(icr)
        enddo
        nd=j-nrem
        npd(nd)=npd(j)
        i1d(nd)=ibeg
        nptd=nptd+npd(j)
        i2d(nd)=nptd
      endif
    enddo
    nd=ndt-nrem
     !Done rebuilding contours.

    if (nd .gt. 0) then
!-----------------------------------------------------------
       !Convert to spherical coordinates: 
       !Note that xd above is longitude while yd is latitude:
      do i=1,nptd
        rlon=xd(i)
        rlat=yd(i)
        clat=cos(rlat)
        xd(i)=clat*cos(rlon)
        yd(i)=clat*sin(rlon)
        zd(i)=sin(rlat)
      enddo
!-------------------------------------------------------
       !Remove points that are extremely close together:
      do i=1,nptd-1
        v(i)=(xd(i+1)-xd(i))**2+(yd(i+1)-yd(i))**2+(zd(i+1)-zd(i))**2
      enddo
  
      do j=1,nd
        i=i2d(j)
        ip1=i1d(j)
        v(i)=(xd(ip1)-xd(i))**2+(yd(ip1)-yd(i))**2+(zd(ip1)-zd(i))**2
      enddo

      ndt=nd
      nd=0
      nptd=0

      do j=1,ndt
        nptbeg=nptd
        do i=i1d(j),i2d(j)
          if (v(i) .gt. small) then
            nptd=nptd+1
            xd(nptd)=xd(i)
            yd(nptd)=yd(i)
            zd(nptd)=zd(i)
          endif
        enddo
        npdiff=nptd-nptbeg
        if (npdiff .lt. 5) then
          nptd=nptbeg
        else
          nd=nd+1
          i1d(nd)=nptbeg+1
          i2d(nd)=nptd
          npd(nd)=npdiff
        endif
      enddo
  
      if (nd .gt. 0) then    !!!!Really is this needed already in a block ? 
!------------------------------------------
         !Redistribute points on the contours:
        do ireno=1,nreno
          call renode
        enddo
!-------------------------------------------
         !Write contour indices for this level:
        inda=lev-nlevm+(lev-1)/nlevm-1
        do j=1,nd
          i1a=npta+i1d(j)
          write(11,'(i6,1x,i8,2(1x,i4))') npd(j),i1a,inda,lz
        enddo
         !Write nodes:
        do i=1,nptd
          write(12,'(f12.9,1x,f12.9,1x,f12.9)') xd(i),yd(i),zd(i)
        enddo

         !Augment number of contours and nodes:
        na=na+nd
        npta=npta+nptd

      endif
    endif
  enddo
   !End of loop over contour levels
endif

return
end subroutine

!=======================================================================

 !Main end module
end module
