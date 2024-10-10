!------------------------------------------------------------------------
!   General purpose contouring routine based on the method developed
!   originally by Dritschel & Ambaum, Mon. Wea. Rev. 134(9), 
!   2503-2514 (2006).
!------------------------------------------------------------------------

program g2c

 !Import contants and parameters:
use constants

implicit none

 !Define common space:
integer,parameter:: nx=ng, ny=ng
integer,parameter:: nlevm=2000, npm=400*ng*ng, nplm=npm, nm=npm/20

integer,parameter:: nxp1=nx+1,nyp1=ny+1
integer,parameter:: nxp2=nx+2 nyp2=ny+2

double precision,parameter:: small=1.d-12, small3=small*small*small

 !Contour arrays:
double precision:: xa(npm), ya(npm)
double precision:: xd(nplm),yd(nplm),dx(nplm),dy(nplm)
double precision:: a(nplm), b(nplm), c(nplm), d(nplm)
double precision:: u(nplm), v(nplm), q(nplm)
integer:: next(nplm)
integer:: i1d(nm),i2d(nm),npd(nm),indd(nm)
integer:: i1a(nm),i2a(nm),npa(nm),inda(nm),laya(nm)
integer:: il1a(0:nz),il2a(0:nz),jl1a(0:nz),jl2a(0:nz)

 !Grid -> Contour arrays:
double precision:: qa(0:nyp1,0:nxp1)
double precision:: xg(0:nx),yg(0:ny)
double precision:: qlev(2*nlevm)
integer:: ibx(0:nx,0:1),iby(0:ny,0:1)

 !Basic parameters:
double precision:: amu,ell,elf,densf,dm,dmsq,dmi
double precision:: glx,gly,dq,dqi,qoff
double precision:: xbeg,xend,ybeg,yend
double precision:: xwid,hlxi,xavg,ywid,hlyi,yavg
integer:: ixbeg,iybeg
integer:: na,nd,npta,nptd
integer:: levbeg,levend

 !Control variables:
logical:: corn(npm),last,perix,periy

 !Local variables:
integer:: ix,iy,iz

!===================================================================
 !Start of executable statements:


      character infile*80

c--------------------------------------------------------------------
c      Open input data file:
      write(*,*) ' File containing the gridded data? '
      read(*,'(a)') infile
      open(20,file=infile,status='old')
      rewind 20

      write(*,*) ' (1) x inner dimension, y outer, or (2) vice-versa?'
      read(*,*) iord

      write(*,*) ' Column containing the data? '
      read(*,*) kol

      write(*,*) ' Is the data periodic in x (1 = yes, 0 = no)?'
      read(*,*) ixper
      perix=ixper .eq. 1

      write(*,*) ' Is the data periodic in y (1 = yes, 0 = no)?'
      read(*,*) iyper
      periy=iyper .eq. 1

c      Define fixed arrays and constants:
      call init

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c      Process each layer in turn:
      do lz=0,nz

c---------------------------------
c        Read data for this layer:
        if (iord .eq. 1) then
          do iy=1,ny
            do ix=1,nx
              read(20,*) (qvar(k),k=1,kol)
              qa(iy,ix)=qvar(kol)
            enddo
          enddo
        else
          do ix=1,nx
            do iy=1,ny
              read(20,*) (qvar(k),k=1,kol)
              qa(iy,ix)=qvar(kol)
            enddo
          enddo
        endif

c-----------------------------------------------------------
c        Add a row or rows to deal with boundary conditions:
        if (perix) then
c          The data are periodic in x:
          do iy=1,ny
            qa(iy,nxp1)=qa(iy,1)
          enddo

          if (periy) then
c            The data are periodic in y also:
            do ix=1,nxp1
              qa(nyp1,ix)=qa(1,ix)
            enddo
          else
c            The data are not periodic in y:
            qbot=zero
            qtop=zero
            do ix=1,nx
              qbot=qbot+qa(1,ix)
              qtop=qtop+qa(ny,ix)
            enddo
            qbot=qbot/dble(nx)
            qtop=qtop/dble(nx)
            do ix=1,nxp1
              qa(0,ix)   =qbot
              qa(nyp1,ix)=qtop
            enddo
          endif

        else
c          The data are not periodic in x:

          if (periy) then
c            The data are periodic in y:
            qlft=zero
            qrgt=zero
            do iy=1,ny
              qlft=qlft+qa(iy,1)
              qrgt=qrgt+qa(iy,nx)
            enddo
            qlft=qlft/dble(ny)
            qrgt=qrgt/dble(ny)
            do iy=1,ny
              qa(iy,0)   =qlft
              qa(iy,nxp1)=qrgt
            enddo

            do ix=0,nxp1
              qa(nyp1,ix)=qa(1,ix)
            enddo
          else
c            The data are not periodic in y either:
            qavg=zero
            do ix=1,nx
              qavg=qavg+qa(1,ix)+qa(ny,ix)
            enddo
            do iy=1,ny
              qavg=qavg+qa(iy,1)+qa(iy,nx)
            enddo
            qavg=qavg/dble(2*(nx+ny))

            do iy=1,ny
              qa(iy,0)   =qavg
              qa(iy,nxp1)=qavg
            enddo

            do ix=0,nxp1
              qa(0,ix)   =qavg
              qa(nyp1,ix)=qavg
            enddo
          endif

        endif

c----------------------------------
c        Generate contours (xd,yd):
        call grid2con

c------------------------------------------------
c        Store the contours (xd,yd) into (xa,ya):
        call storecon(lz)

      enddo
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c      Ends loop over layers
      close(20)

c      Write contours in congen.dat:
      call writecon

      write(*,*) 
      write(*,*) ' The contours are now available in congen.dat'
      write(*,*) '   and the levels are available in levels.dat'

      end

c=======================================================================
      subroutine init
c      Defines various fixed arrays and constants.

      include 'g2c_com.i'

c------------------------------------------------------
      write(*,*) ' Domain limits, x_min and x_max?'
      read(*,*) xbeg,xend
      xwid=xend-xbeg
      if (perix) then
        ixbeg=1
        hlxi=one/(xwid/two+small)
        xavg=f12*(xbeg+xend)
      else
        ixbeg=0
      endif

      write(*,*) ' Domain limits, y_min and y_max?'
      read(*,*) ybeg,yend
      ywid=yend-ybeg
      if (periy) then
        iybeg=1
        hlyi=one/(ywid/two+small)
        yavg=f12*(ybeg+yend)
      else
        iybeg=0
      endif

c      Save domain limits for use in pcon.F (prints the contours):
      open(22,file='domain_limits',status='unknown')
      write(22,'(2(1x,f16.9))') xbeg,xend
      write(22,'(2(1x,f16.9))') ybeg,yend
      close(22)

      write(*,*) ' Contour interval, dq?'
      read(*,*) dq

      write(*,*) ' Contours are found for the field values '
      write(*,*) ' qbar +/- dq/2, qbar +/- 3*dq/2, etc.'
      write(*,*) ' Enter qbar:'
      read(*,*) qbar

      dqi=one/dq
      qoff=dq*dble(nlevm)-qbar
c      qoff: should be a large integer multiple of the 
c            contour interval, dq.  The multiple should exceed 
c            the maximum expected number of contour levels.

      do lev=1,2*nlevm
        qlev(lev)=(dble(lev)-f12)*dq-qoff
      enddo

      glx=xwid/dble(nx-1+ixbeg)
      gly=ywid/dble(ny-1+iybeg)
c      glx, gly: the x & y grid lengths

c--------------------------------------------------
c      Node redistribution parameters (see renode):

c      Maximum node separation / dx:
      sepmax=1.25d0
c      Node spacing parameter:
      amu=0.2d0
c      Large-scale length:
      ell=four*sepmax**2*sqrt(glx*gly)
c      Surgical scale:
      dm=ell*amu**2/four
      write(*,*)
      write(*,'(2(a,f9.7))') ' Using mu = ',amu,' and L = ',ell

      dmsq=four*dm**2
      dmi=two/dm
      elf=one/ell**2
      densf=one/(amu*sqrt(ell))

c---------------------------------------------
c      Define x & y grid lines (see grid2con):
      do ix=ixbeg,nx
        xg(ix)=glx*dble(ix-1)+xbeg
      enddo

      do iy=iybeg,ny
        yg(iy)=gly*dble(iy-1)+ybeg
      enddo
c      Note: contours are squeezed back into the domain after.

c----------------------------------------------------------------
c      Define grid-box reference indices taking into account BCs:
      if (perix) then
c        The data are periodic in x:

        if (periy) then
c          The data are periodic in y also:
          do ix=1,nx
            ibx(ix,1)=ny*(ix-1)
          enddo
          do ix=2,nx
            ibx(ix,0)=ibx(ix-1,1)
          enddo
          ibx(1,0)=ibx(nx,1)

          do iy=1,ny
            iby(iy,1)=iy
          enddo
          do iy=2,ny
            iby(iy,0)=iby(iy-1,1)
          enddo
          iby(1,0)=iby(ny,1)

        else
c          The data are not periodic in y:
          do ix=1,nx
            ibx(ix,1)=nyp2*(ix-1)
          enddo
          do ix=2,nx
            ibx(ix,0)=ibx(ix-1,1)
          enddo
          ibx(1,0)=ibx(nx,1)

          do iy=0,ny
            iby(iy,1)=iy+2
            iby(iy,0)=iy+1
          enddo

        endif

      else
c        The data are not periodic in x:

        if (periy) then
c          The data are periodic in y:
          do ix=0,nx
            ibx(ix,1)=ny*(ix+1)
            ibx(ix,0)=ny*ix
          enddo

          do iy=1,ny
            iby(iy,1)=iy
          enddo
          do iy=2,ny
            iby(iy,0)=iby(iy-1,1)
          enddo
          iby(1,0)=iby(ny,1)

        else
c          The data are not periodic in y either:
          do ix=0,nx
            ibx(ix,1)=nyp2*ix
          enddo
          do ix=1,nx
            ibx(ix,0)=ibx(ix-1,1)
          enddo

          do iy=0,ny
            iby(iy,1)=iy+2
            iby(iy,0)=iy+1
          enddo

        endif

      endif

      return
      end

c=======================================================================
      subroutine grid2con
c      Generates new contours (xd,yd) from the gridded data (qa) 

      include 'g2c_com.i'

c      Local Grid -> Contour arrays:
      parameter (ncrm=2*npm/nz)
c      ncrm: max number of contour crossings on the finest grid

      dimension ycr(ncrm),xcr(ncrm)
      dimension lcr(ncrm),kib(ncrm),kob(ncrm),ipo(ncrm)
      dimension lgrx(0:nxp1),isix(0:nx)
      dimension lgry(0:nyp1),isiy(0:ny)
      dimension icrtab(nyp2*nxp2,2)
      dimension ncrpl(2*nlevm),ibpl(2*nlevm)
      integer*1 noctab(nyp2*nxp2)
      logical crox(0:nx),croy(0:ny),free(ncrm)
      logical keep(nplm)

c------------------------------------------------------
c      Counter for total number of grid line crossings:
      ncr=0
c      Counters for each level:
      do lev=1,2*nlevm
        ncrpl(lev)=0
      enddo

c------------------------------------------------------
c      Find x grid line crossings first (x = constant):
      do ix=1,nx
        xgt=xg(ix)

        do iy=iybeg,nyp1
          lgry(iy)=nint((qoff+qa(iy,ix))*dqi)
        enddo

        do iy=iybeg,ny
          lgrd=lgry(iy+1)-lgry(iy)
          croy(iy)=(lgrd .ne. 0)
          isiy(iy)=sign(1,lgrd)
        enddo

        do iy=iybeg,ny
          if (croy(iy)) then
            lgrt=lgry(iy)
            jump=isiy(iy)
            inc=(1+jump)/2
            kibt=iby(iy,1)+ibx(ix,inc)
            kobt=iby(iy,1)+ibx(ix,1-inc)
            fac=gly/(qa(iy+1,ix)-qa(iy,ix))
            do while (lgrt .ne. lgry(iy+1))
              ncr=ncr+1
              xcr(ncr)=xgt
              lev=lgrt+inc
              lcr(ncr)=lev
              kib(ncr)=kibt
              kob(ncr)=kobt
              ycr(ncr)=yg(iy)+fac*(qlev(lev)-qa(iy,ix))
              ncrpl(lev)=ncrpl(lev)+1
              lgrt=lgrt+jump
            enddo
          endif
        enddo

      enddo

c      Above, kib = grid box into which the contour (containing icr) is going
c             kob =   "   "  out of "    "     "         "       "    " coming
c        [kob -> icr -> kib:  icr lies at the boundary between kob & kib]

c-----------------------------------------------------
c      Find y grid line crossings next (y = constant):
      do iy=1,ny
        ygt=yg(iy)

        do ix=ixbeg,nxp1
          lgrx(ix)=nint((qoff+qa(iy,ix))*dqi)
        enddo

        do ix=ixbeg,nx
          lgrd=lgrx(ix+1)-lgrx(ix)
          crox(ix)=(lgrd .ne. 0)
          isix(ix)=sign(1,lgrd)
        enddo

        do ix=ixbeg,nx
          if (crox(ix)) then
            lgrt=lgrx(ix)
            jump=isix(ix)
            inc=(1+jump)/2
            kibt=ibx(ix,1)+iby(iy,1-inc)
            kobt=ibx(ix,1)+iby(iy,inc)
            fac=glx/(qa(iy,ix+1)-qa(iy,ix))
            do while (lgrt .ne. lgrx(ix+1))
              ncr=ncr+1
              ycr(ncr)=ygt
              lev=lgrt+inc
              lcr(ncr)=lev
              kib(ncr)=kibt
              kob(ncr)=kobt
              xcr(ncr)=xg(ix)+fac*(qlev(lev)-qa(iy,ix))
              ncrpl(lev)=ncrpl(lev)+1
              lgrt=lgrt+jump
            enddo
          endif
        enddo

      enddo

c      Find the first and last level having non-zero number of crossings:
      lev=1
      do while (ncrpl(lev) .eq. 0) 
        lev=lev+1
      enddo
      levbeg=lev
      do while (ncrpl(lev) .ne. 0) 
        lev=lev+1
      enddo
      levend=lev-1

      ibpl(levbeg)=1
      do lev=levbeg+1,levend
        ibpl(lev)=ibpl(lev-1)+ncrpl(lev-1)
      enddo

      if (ncr .gt. ncrm) then
        write(*,'(a,i8)') ' ncr  = ',ncr
        write(*,'(a,i8)') ' ncrm = ',ncrm
        write(*,*) ' ncr > ncrm!  Stopping!'
        stop
      endif

      do icr=1,ncr
        lev=lcr(icr)
        ibt=ibpl(lev)
        ipo(ibt)=icr
        ibpl(lev)=ibt+1
      enddo

      do lev=levbeg,levend
        ibeg=ibpl(lev)-ncrpl(lev)
        iend=ibpl(lev)-1

c        Initialise number of crossings per box:
        do i=ibeg,iend
          noctab(kob(ipo(i)))=0
        enddo

        write(*,*) ' # of crossings = ',iend-ibeg+1
        do i=ibeg,iend
          icr=ipo(i)
c          icr is the index of the current crossing at level lev.
          k=kob(icr)
c          accumulate number of crossings in this box:
          noctab(k)=noctab(k)+1
c          assign crossing to box, permitting 2 crossings:
          icrtab(k,noctab(k))=icr
        enddo

        do i=ibeg,iend
          icr=ipo(i)
          k=kib(icr)
          noc=noctab(k)
c          Use last crossing in this box as the next node:
          kob(icr)=icrtab(k,noc)
c          kob(icr) now gives the next point after icr
          noctab(k)=noc-1
c          This will normally be zero, except for boxes with 2 crossings;
c          this allows a second use of this box.
        enddo

      enddo

c-----------------
c      Now re-build contours:
      j=0
      i=0

      do icr=1,ncr
        free(icr)=.true.
      enddo

      do icr=1,ncr
        if (free(icr)) then
c          A new contour (j) starts here:
          i=i+1
          j=j+1
          i1d(j)=i
          indd(j)=lcr(icr)
          ipo(i)=icr
          icrn=kob(icr)

          do while (icrn .ne. icr)
c            Find remaining points on contour j:
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
          indd(nd)=indd(j)
          i1d(nd)=ibeg
          nptd=nptd+npd(j)
          i2d(nd)=nptd
        endif
      enddo
c      Done rebuilding contours.

c-------------------------------------------------------------------
      if (perix .and. periy) then
c        Redistribute points on the contours to reduce point density:
        call renode

      else
c        For non-periodic data, remove points along domain boundaries:

c        Define the point following a node i by next(i):
        do i=1,nptd-1
          next(i)=i+1
        enddo

        do j=1,nd
          next(i2d(j))=i1d(j)
        enddo

        if (perix) then
          do i=1,nptd
            xx=xd(next(i))-xd(i)
            dx(i)=xx-xwid*int(xx*hlxi)
          enddo
        else
          do i=1,nptd
            xd(i)=min(xend,max(xbeg,xd(i)))
          enddo
          do i=1,nptd
            dx(i)=xd(next(i))-xd(i)
            corn(i)=abs(dx(i)) .gt. small
c            corn = false when the segment runs along the left or right edge
          enddo
          do i=1,nptd
            ia=next(i)
            keep(ia)=corn(i) .or. corn(ia)
c            keep = false means node ia will be eliminated below
          enddo
        endif

        if (periy) then
          do i=1,nptd
            yy=yd(next(i))-yd(i)
            dy(i)=yy-ywid*int(yy*hlyi)
          enddo
        else
          do i=1,nptd
            yd(i)=min(yend,max(ybeg,yd(i)))
          enddo
          do i=1,nptd
            dy(i)=yd(next(i))-yd(i)
            corn(i)=abs(dy(i)) .gt. small
c            corn = false when the segment runs along the top or bottom edge
          enddo
          do i=1,nptd
            ia=next(i)
            keep(ia)=corn(i) .or. corn(ia)
c            keep = false means node ia will be eliminated below
          enddo
        endif

        ndt=nd
        nd=0
        nptd=0

        do j=1,ndt
          nptbeg=nptd
          do i=i1d(j),i2d(j)
            if (keep(i)) then
              nptd=nptd+1
              xd(nptd)=xd(i)
              yd(nptd)=yd(i)
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
      endif

      return
      end
c
c=======================================================================
      subroutine storecon(lz)
c      Stores the shifted contours xd, yd, i1d, i2d ... 
c                  into the arrays xa, ya, i1a, i2a ....  
c      Also updates na & npta and the node & contours layer indices 
c      il1a, il2a, jl1a & jl2a.

      include 'g2c_com.i'

c-----------------------------------------------------
c      Store starting layer indices:
      jl1a(lz)=na+1
      il1a(lz)=npta+1

c      Store nodes:
      do i=1,nptd
        xa(npta+i)=xd(i)
        ya(npta+i)=yd(i)
      enddo

c      Store contour indices:
      do j=1,nd
        jn=na+j
        i1a(jn) =npta+i1d(j)
        i2a(jn) =npta+i2d(j)
        npa(jn) =npd(j)
        inda(jn)=indd(j)
        laya(jn)=lz
      enddo

c      Augment number of contours and nodes:
      na=na+nd
      npta=npta+nptd

c      Store ending layer indices:
      jl2a(lz)=na
      il2a(lz)=npta

      return
      end

c=======================================================================
      subroutine renode
c      Re-nodes each contour while preserving corner locations.
c      Uses square-root dependence on a weighted sum of nearby 
c      curvature values.

      include 'g2c_com.i'

c------------------------------------------------------------------------
c      Define the point following a node i by next(i):
      do i=1,nptd-1
        next(i)=i+1
      enddo

      do j=1,nd
        next(i2d(j))=i1d(j)
      enddo

c------------------------------------------------------------------------
c      Get the updated cubic interpolation coefficients:
      call cubic

c------------------------------------------------------------------------
c      Use the spherical curvature expression (radius of the sphere = ell)
c      to ensure an adequate node density in low curvature regions.
      do i=1,nptd
        ww=one/(v(i)+dmsq)
        u(i)=ww*sqrt(elf*v(i)+u(i)**2)
        v(i)=ww*d(i)
      enddo
c      NB: elf = 1/ell**2; v(i) = |xx_{i+1}-xx_{i}|**2; d(i)=sqrt{v(i)};
c          u(i)/d(i) = (kappa_{i}+kappa_{i+1})/2; dmsq = (2*dm)**2

c      Re-assign curvature at a node from weighted average on either side
c      (v above is the weight):
      do ib=1,nptd
        i=next(ib)
        q(i)=(u(ib)+u(i))/(v(ib)+v(i))
      enddo

c      Re-average to get interval value (effectively, four curvature
c      values go into getting the final interval value, u(i)):
      do i=1,nptd
        ia=next(i)
        u(i)=f12*(q(i)+q(ia))
      enddo

c      Compute fractional number of nodes to be placed between old
c      nodes i and i+1:
      do i=1,nptd
        d(i)=d(i)*min(dmi,densf*sqrt(u(i))+u(i))
      enddo
c      NB: dmi = 2/delta; densf = 1/(amu*sqrt{ell})

c------------------------------------------------------------------------
c      Now begin the redistribution of nodes contour by contour,
c      making sure to preserve corner locations:
      nptd=0
      do j=1,nd
        inew=1
        i1t=i1d(j)
        i1d(j)=nptd+1
300       u(nptd+inew)=xd(i1t)
          v(nptd+inew)=yd(i1t)
          sum=zero
          i=i1t
310         sum=sum+d(i)
            i=i+1
            last=i .gt. i2d(j)
            if (last) goto 330
            if (corn(i)) goto 320
            goto 310
320       if (sum .lt. small) then
            i1t=i
            goto 300
          else
            i2t=i-1
            goto 340
          endif
330       if (sum .lt. small) then
            inew=inew-1
            goto 390
          else
            i2t=i-1
          endif
340       npseg=nint(sum)+1
c          npseg-1 is the number of nodes to be placed on the contour segment.
          fac=dble(npseg)/sum
          do i=i1t,i2t
            d(i)=fac*d(i)
          enddo
c          Now, the sum of d(i) is equal to npseg.
c          The first node along a contour (segment) is fixed;
c          find the new node positions:
          acc=zero
          i=i1t-1
          do im=nptd+inew+1,nptd+inew+npseg-1
            if (acc .ge. one) goto 370
360           acc=acc+d(i+1)
              i=i+1
              if (acc .lt. one) goto 360
370         acc=acc-one
            p=one-acc/d(i)
            eta=p*(a(i)+p*(b(i)+p*c(i)))
            u(im)=xd(i)+p*dx(i)-eta*dy(i)
            v(im)=yd(i)+p*dy(i)+eta*dx(i)
          enddo
          if (last) then
            inew=inew+npseg-1
            goto 390
          else
            inew=inew+npseg
            i1t=i2t+1
            goto 300
          endif
390     npd(j)=inew
        nptd=nptd+inew
      enddo

c------------------------------------------------------------------------
c      Switch arrays around again:
      if (perix) then
        do i=1,nptd
          xx=u(i)-xavg
          xd(i)=xavg+xx-xwid*int(xx*hlxi)
        enddo
      else
        do i=1,nptd
          xd(i)=min(xend,max(xbeg,u(i)))
        enddo
      endif

      if (periy) then
        do i=1,nptd
          yy=v(i)-yavg
          yd(i)=yavg+yy-ywid*int(yy*hlyi)
        enddo
      else
        do i=1,nptd
          yd(i)=min(yend,max(ybeg,v(i)))
        enddo
      endif

c      Reset ending contour indices:
      do j=1,nd
        i2d(j)=i1d(j)+npd(j)-1
      enddo

      return
      end

c=======================================================================
      subroutine cubic
c      Calculates the interpolation coefficients between the nodes 
c      [xd(i),yd(i)] and [xd(next(i)),yd(next(i))], i = 1, ..., nptd.

c      The interpolation approximately enforces continuity of curvature 
c      (except at corners which have effectively infinite curvature).

      include 'g2c_com.i'

c----------------------------------------------------------------------
      if (perix) then
        do i=1,nptd
          ia=next(i)
          xx=xd(ia)-xd(i)
          dx(i)=xx-xwid*int(xx*hlxi)
        enddo
      else
        do i=1,nptd
          ia=next(i)
          dx(i)=xd(ia)-xd(i)
        enddo
      endif

      if (periy) then
        do i=1,nptd
          ia=next(i)
          yy=yd(ia)-yd(i)
          dy(i)=yy-ywid*int(yy*hlyi)
        enddo
      else
        do i=1,nptd
          ia=next(i)
          dy(i)=yd(ia)-yd(i)
        enddo
      endif

      do i=1,nptd
        v(i)=dx(i)*dx(i)+dy(i)*dy(i)+small
        d(i)=sqrt(v(i))
      enddo

      do ib=1,nptd
        i=next(ib)
        u(i)=v(ib)
        c(i)=-dx(ib)
        q(i)=-dy(ib)
      enddo

      do i=1,nptd
        corn(i)=dx(i)*c(i)+dy(i)*q(i) .gt. zero
        if (corn(i)) then
c          Set curvature to zero at corners:
          b(i)=zero
        else
          b(i)=(dx(i)*q(i)-c(i)*dy(i))/
     .     sqrt((c(i)*v(i)-dx(i)*u(i))**2+
     .          (q(i)*v(i)-dy(i)*u(i))**2+small3)
        endif
      enddo

      do i=1,nptd
        ia=next(i)
        u(i)=d(i)*(b(ia)+b(i))
        c(i)=d(i)*(b(ia)-b(i))
      enddo

c      Calculate the cubic interpolation coefficients:
      do i=1,nptd
        a(i)=f16*c(i)-f12*u(i)
        b(i)=f12*(u(i)-c(i))
        c(i)=f13*c(i)
      enddo

      return
      end

c=======================================================================
      subroutine writecon
c      Writes contours to congen.dat and contour levels to levels.dat:

      include 'g2c_com.i'

c----------------------------------------------------------------------
c      Adjust the levels so that - values have negative indices
c                            and + values have positive indices:
      do j=1,na
        lev=inda(j)
        inda(j)=lev-nlevm+(lev-1)/nlevm-1
      enddo

      open(12,file='levels.dat',status='unknown')

      write(*,*) 
      write(*,*) ' Level      Value '
      write(*,*) ' -----    ---------- '
      do lev=levbeg,levend
        levnew=lev-nlevm+(lev-1)/nlevm-1
        write(12,'(2x,i5,4x,f13.9)') levnew,qlev(lev)
        write( *,'(2x,i5,4x,f13.9)') levnew,qlev(lev)
      enddo

      close(12)

c----------------------------------------------------------------------
      open(11,file='congen.dat',status='unknown')

      write(*,*)
      write(*, '(a,i6,a,i7)') '   # contours = ',na,
     .                           '   # nodes = ',npta

      write(11,'(i6,1x,i7,1x,f7.2)') na,npta,zero
      do j=1,na
        write(11,'(i5,1x,i7,2(1x,i5),f14.10)') 
     .    npa(j),i1a(j),inda(j),laya(j),dq
      enddo
      do i=1,npta
        write(11,'(f12.9,1x,f12.9)') xa(i),ya(i)
      enddo

      close(11)

      return
      end
