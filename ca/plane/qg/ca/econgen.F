c---------------------------------------------------------------------------
c         Converts contours in pvcont.bin to gridded values on an
c         ultra-fine grid of dimensions mgu*nh x mgu*nh, adds the
c         the residual PV in pvresi.bin (interpolated to the ultra-
c         fine grid), then creates new contours.

c         ==> This version is for a non-zero planetary vorticity
c             gradient (beta);  use mcongen.F otherwise.

c         ==> It is ESSENTIAL that the PV jump be of the form
c             2*pi*beta/n, where n is a positive integer.

c         ==> Dynamic memory allocation is used; the alternative
c             is to use mergepv.F, ylimits.F & bcongen.F separately
c             (see e.g. the script oldf2d).

c         This routine processes one contour level at a time, and
c         writes data to disk, in three files:

c            head_congen.bin   : one-line header containing total
c                                # of contours, nodes and t = 0
c            all_contours.bin  : all contour index information
c            all_nodes.bin     : x,y coordinates of all nodes

c         This avoids storing large arrays in memory.  After
c         running this programme, the conventional pvcont.bin
c         can be obtained by issuing the command:

c       cat head_congen.bin all_contours.bin all_nodes.bin > pvcont.bin

c         This is normally run via a script such as batch_clam2d and
c         batch_caps2d. When the algorithm used is CASL, this routine
c         also computes min, max and rms values of the PV field for
c         further adaptation of the sepctral damping rate.

c         Adapted mcongen.F & bcongen.F on 16/9/2008 by D.G. Dritschel
c         @ St Andrews
c         Modified on 24/9/2008 by J. Fontane @ St Andrews

c         Replaced bi-linear interpolation of the residual PV by a
c         much more accurate spectral interpolation on 22/9/2009
c         (D.G. Dritschel @ St Andrews)

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      program congen

      include 'econgen_com.i'

c      Local parameters and variables:

c      Used for copying qd into finer spectral array qdf:
      parameter (kxoff=nhf-nh,kyoff=kxoff)

      dimension qdf(nhf,nhf),qd(nh,nh),bety(nhf)
      dimension qqmin(nhf+1),qqmax(nhf+1)
      dimension i1(nm),i2(nm),np(nm),lay(nm)
      dimension jl1(nz),jl2(nz),il1(nz),il2(nz)
      logical conts

c--------------------------------------------------------------------
c      Read PV jump (dq), node redistribution parameters (amu,ell)
c      and planetary vorticity gradient (beta) from parameter file:
      open(20,file='params.dat',status='old')
      rewind 20
      read(20,*)
      read(20,*) tsave,tsim,dtmax,dq
      read(20,*) amu,ell
      read(20,*) beta
      close(20)

c      Contrast in beta*y across the domain in y:
      ppvdif=twopi*beta

c      Check that ppvdif is an integer multiple of dq:
      brat=ppvdif/dq
      ibrat=nint(brat)
      if (abs(brat-dble(ibrat)) .gt. 1.d-7) then
        write(*,'(a,f12.7)') ' 2*pi*beta/dq = ',brat
        write(*,*) ' This should be an integer!'
        stop
      endif

c-----------------------------------------------------------------
c      Quantities needed in contour -> grid conversion (con2grid):
      fnwhf=dble(nhf/2)

c      Grid lengths and their inverses:
      gl=twopi/dble(nhf)
      gli=one/gl

c--------------------------------------------------------------
c      Quantities needed in generating new contours (grid2con):
      dqi=one/dq
      qoff=dq*dble(nlevm)
c      qoff: should be a large integer multiple of the
c            contour interval, dq.  The multiple should exceed
c            the maximum expected number of contour levels.

      do lev=1,2*nlevm
        qlev(lev)=(dble(lev)-half)*dq-qoff
      enddo

c      Used for extending array in y:
      bfac=beta*gl

c      x coordinates of grid lines:
      do ix=1,nhf+1
        xg(ix)=gl*dble(ix-1)-pi
      enddo

c-----------------------------------------------
c      Used in node redistribution (see renode):
      dm=amu**2*ell/four
      dmsq=four*dm**2
      dmi=two/dm
      elf=one/ell**2
      densf=one/(amu*sqrt(ell))

c--------------------------------------------------------------------
c      Open contour data file:
      open(10,file='pvcont.bin',status='old',form='unformatted')
      rewind 10

      read(10) n,npt,t
      write(*, '(a,i8,a,i10,a,f12.5)') '   # contours = ',n,
     .                                 '   # nodes = ',npt,
     .                                 '   time = ',t

c      Allow for the possibility of no contours:
      conts=(n .gt. 0)
      if (conts) then

c        Read only contour indices for now, to save memory:
        do j=1,n
          read(10) np(j),i1(j),ind,lay(j)
        enddo

c        Get starting and ending layer indices, as well as i2(j):
        do lz=1,nz
          jl2(lz)=0
        enddo

        do j=1,n
          i2(j)=i1(j)+np(j)-1
          il2(lay(j))=i2(j)
          jl2(lay(j))=j
        enddo

        ibeg=1
        jbeg=1

        do lz=1,nz
          if (jl2(lz) .gt. 0) then
            il1(lz)=ibeg
            ibeg=il2(lz)+1
            jl1(lz)=jbeg
            jbeg=jl2(lz)+1
          endif
        enddo

      endif

c--------------------------------------------------------------------
c      Open residual PV file containing qd on the inversion grid:
      open(21,file='pvresi.bin',status='old',form='unformatted')
      rewind 21

c      Used below for spectral interpolation:
      fmgu=dble(mgu)

c--------------------------------------------------------------------
c      Open output files:
      open(11,file='all_contours.bin',status='unknown',
     .        form='unformatted')
      open(12,file='all_nodes.bin',status='unknown',
     .        form='unformatted')

c      Counters for total number of nodes and contours:
      npta=0
      na=0

c+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%
c      Process data layer by layer:
      do lz=1,nz

c        Read gridded residual PV from pvresi.bin for this layer:
        read(21) qd

c        Convert qd to spectral space:
        call ptospc(qd)

c        Pad qd array with zeros and write into qdf:
        do kx=1,nhf
          do ky=1,nhf
            qdf(kx,ky)=zero
          enddo
        enddo

        do ky=1,nwh+1
          do kx=1,nwh+1
            qdf(kx,ky)=fmgu*qd(kx,ky)
          enddo
          do kx=nh,nwh+2,-1
            qdf(kxoff+kx,ky)=fmgu*qd(kx,ky)
          enddo
        enddo

        do ky=nh,nwh+2,-1
          do kx=1,nwh+1
            qdf(kx,kyoff+ky)=fmgu*qd(kx,ky)
          enddo
          do kx=nh,nwh+2,-1
            qdf(kxoff+kx,kyoff+ky)=fmgu*qd(kx,ky)
          enddo
        enddo

c        Return qdf to physical space:
        call spctopf(qdf)

        if (conts) then
c          Convert contours in this layer to a fine gridded field, qa:

          ioff=il1(lz)-1
          nptd=il2(lz)-ioff

c          Read nodes in this layer:
          do i=1,nptd
            read(10) xd(i),yd(i)
          enddo

c          Get relative contour indices:
          joff=jl1(lz)-1
          nd  =jl2(lz)-joff

          do j=1,nd
            jj=joff+j
            i1d(j)=i1(jj)-ioff
            i2d(j)=i2(jj)-ioff
            npd(j)=np(jj)
          enddo

c          Convert contours to gridded values:
          call con2grid

c          Add fine-grid qd to qa:
          do ix=1,nhf
            do iy=1,nhf
              qa(iy,ix)=qa(iy,ix)+qdf(iy,ix)
            enddo
          enddo

        else

c          No contours: interpolate qd to the fine grid as qa
c                       and add beta*y:
          do iy=1,nhf
            bety(iy)=beta*(gl*dble(iy-1)-pi)
          enddo

          do ix=1,nhf
            do iy=1,nhf
              qa(iy,ix)=bety(iy)+qdf(iy,ix)
            enddo
          enddo

        endif

c        Find min & max PV in every grid line:
        do iy=1,nhf
          qqmin(iy)=qa(iy,1)
          qqmax(iy)=qa(iy,1)
        enddo

        do ix=2,nhf
          do iy=1,nhf
            qqmin(iy)=min(qqmin(iy),qa(iy,ix))
            qqmax(iy)=max(qqmax(iy),qa(iy,ix))
          enddo
        enddo

c        Define values at y = +pi:
        qqmin(nhf+1)=qqmin(1)+ppvdif
        qqmax(nhf+1)=qqmax(1)+ppvdif
c        ppvdif=twopi*beta is the contrast of beta*y across the domain.

c        Find overall maximum:
        qamin=qqmin(1)
        qamax=qqmax(1)
        do iy=2,nhf+1
          qamin=min(qamin,qqmin(iy))
          qamax=max(qamax,qqmax(iy))
        enddo

c        Range of PV levels to contour:
        levbeg=int((qoff+qamin)/dq+half)+1
        levend=int((qoff+qamax)/dq+half)

c        Find minimum grid line iy = iymin containing all contours
c        within basic periodic box:
        levmin=int(qamin/dq-half)
        qc=(dble(levmin)+half)*dq
        iy=1
1       iy=iy-1
        iy0=1+mod(100*nhf+iy-1,nhf)
c        Adequate for contours extending over < 100 periodic boxes!
        qadd=bfac*dble(iy-iy0)
        if (qadd+qqmax(iy0) .gt. qc) goto 1
        iymin=iy

c        Find maximum grid line iy = iymax containing all contours
c        within basic periodic box:
        levmax=int(qamax/dq+half)
        qc=(dble(levmax)-half)*dq
        iy=nhf+1
2       iy=iy+1
        iy0=1+mod(iy-1,nhf)
        qadd=bfac*dble(iy-iy0)
        if (qadd+qqmin(iy0) .lt. qc) goto 2
        iymax=iy

c        For dimensioning allocatable arrays in grid2con:
        ngps=nhf*(iymax-iymin+1)

c        Generate new contours in layer lz:
        call grid2con(lz,iymin,iymax,ngps)

      enddo
c+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%+%
c      Ends loop over layers.

      close(10)
      close(11)
      close(12)
      close(21)

c      Write one line header for congen.bin:
      open(13,file='head_congen.bin',status='unknown',
     .        form='unformatted')
      write(13) na,npta,t
      close(13)

c      Write one line header needed in batch_clam2d:
      open(13,file='head_congen.dat',status='unknown')
      write(13,'(i8,1x,i10,1x,f12.5,1x,f16.12)') na,npta,t,dq
      close(13)

      write(*,*)
      write(*,*) ' All done.  To create pvcont.bin, type:'
      write(*,*) ' cat head_congen.bin all_contours.bin all_nodes.bin'
     .          ,' > pvcont.bin'
      write(*,*) ' rm  head_congen.bin all_contours.bin all_nodes.bin'

c--------------------------------------------------------------------
c      Possibly computes min, max and rms values of the fine grid PV
      write(*,*)
      write(*,*) ' Do you want to compute min, max and rms values?'
      write(*,*) ' yes (0) or no (1). Choice? '
      read(*,*) iopt

      if (iopt .eq. 0) then
c        Get the rms value (min and max are already computed above):
        qarms=0.d0

        do ix=1,nhf
          do iy=1,nhf
            qarms=qarms+qa(iy,ix)**2
          enddo
        enddo

        qarms=sqrt(qarms/dble(nhf*nhf))

        write(*,*) ' The rms field value is'
        write(*,'(1x,f11.5)') qarms
        write(*,*)
        write(*,*) ' The min and max field values are'
        write(*,'(2(1x,f11.5))') qamin,qamax
        write(*,*)
        write(*,*) ' The max absolute field value is'
        write(*,'(1x,f11.5)') max(-qamin,qamax)

      endif

      end

c=======================================================================
      subroutine con2grid
c      Converts PV contours (x,y) to gridded values (qa)
c      on the horizontal grid (nhf,nhf).

      include 'econgen_com.i'

c      Local parameters and arrays:
      parameter (thrpi=three*pi,pinv=one/pi)

      dimension qjy(nhf),qright(nhf)
      dimension nyc(nplm),iyc(nplm)

c-----------------------------------------------------
c      Ensure points all lie inside the domain:
      do i=1,nptd
        xd(i)=xd(i)-twopi*int(xd(i)*hlxi)
        yd(i)=yd(i)-twopi*int(yd(i)*hlyi)
      enddo

c-----------------------------------------------------
c      Get coordinate differences and form logical array telling
c      if (i,i+1) crosses y = -pi:
      do i=1,nptd
        u(i)=gli*(yd(i)+pi)
      enddo
      do i=1,nptd-1
        ia=i+1
        b(i)=u(ia)
        xx=xd(ia)-xd(i)
        dx(i)=xx-twopi*int(xx*hlxi)
        yy=yd(ia)-yd(i)
        dy(i)=yy-twopi*int(yy*hlyi)
        corn(i)=(abs(xx) .gt. pi)
      enddo
      do j=1,nd
        i =i2d(j)
        ia=i1d(j)
        b(i)=u(ia)
        xx=xd(ia)-xd(i)
        dx(i)=xx-twopi*int(xx*hlxi)
        yy=yd(ia)-yd(i)
        dy(i)=yy-twopi*int(yy*hlyi)
        corn(i)=(abs(xx) .gt. pi)
      enddo

c----------------------------------------------------------------
c      Fill PV jump array at the left edge:
c      initialise array:
      do iy=1,nhf
        qjy(iy)=zero
      enddo

      do i=1,nptd
        if (corn(i)) then
          cc=sign(one,dx(i))
          p=-(xd(i)-pi*cc)/(dx(i)+small)
          iyc(i)=1+int(gli*(mod(thrpi+yd(i)+p*dy(i),twopi)))
          c(i)=dq*cc
        endif
      enddo

c      Fill jumps:
      do i=1,nptd
        if (corn(i)) then
          iy=iyc(i)
          qjy(iy)=qjy(iy)+c(i)
        endif
      enddo

c----------------------------------------------------------------
c      Construct grid PV:
c      initialise crossing information:
      do i=1,nptd
        iym1=int(u(i))
        idyt=int(b(i))-iym1
        corn(i)=(idyt .ne. 0)
        if (corn(i)) then
          nyc(i)=idyt-nhf*((2*idyt)/nhf)
          jump=sign(1,nyc(i))
          iyc(i)=iym1+(1+jump)/2
          fac=gl/dy(i)
          p=fac*(iyc(i)-u(i))
          u(i)=pinv*(xd(i)+p*dx(i))
          b(i)=pinv*abs(fac)*dx(i)
          c(i)=dq*jump
        endif
      enddo

c----------------------------------------------------------------------
c      do a vertical sweep to initialize the PV along this "boundary"
      qright(1)=zero
      do iy=1,nhf-1
        qright(iy+1)=qright(iy)+qjy(iy)
      enddo

c      initialise PV jump array:
      do ix=1,nhf
        do iy=1,nhf
          qa(iy,ix)=zero
        enddo
      enddo

      do i=1,nptd
        if (corn(i)) then
c        the contour crosses through at least one horizontal grid line
          jump=sign(1,nyc(i))
          ncr=0
240         iy=1+mod(nhf+iyc(i)+ncr,nhf)
            ix=1+int(fnwhf*(u(i)+one-two*int(u(i))))
            qa(iy,ix)=qa(iy,ix)-c(i)
            u(i)=u(i)+b(i)
            ncr=ncr+jump
          if (ncr .ne. nyc(i)) goto 240
        endif
      enddo

c      get PV values by horizontal sweeping:
      do iy=1,nhf
        qa(iy,nhf)=qright(iy)-qa(iy,nhf)
      enddo
      do ix=nhf-1,1,-1
        do iy=1,nhf
          qa(iy,ix)=qa(iy,ix+1)-qa(iy,ix)
        enddo
      enddo

c      Remove domain average:
      qavg=zero
      do ix=1,nhf
        do iy=1,nhf
          qavg=qavg+qa(iy,ix)
        enddo
      enddo
      qavg=qavg/dble(nhf**2)
      do ix=1,nhf
        do iy=1,nhf
          qa(iy,ix)=qa(iy,ix)-qavg
        enddo
      enddo

      return
      end

c=======================================================================
      subroutine grid2con(lz,iymin,iymax,ngps)

c      Generates new contours (xd,yd) from the gridded data (qa)
c      on the horizontal grid (nhf,nhf).

      include 'econgen_com.i'

      parameter (nreno=2)
c      nreno: number of times renode is called to reduce point
c             density on contours.

      parameter (ncrm=3*nplm/4)
c      ncrm: max number of contour crossings of a single field
c            level on the finest grid

c      Local (fixed) Grid -> Contour arrays:
      dimension kib(ncrm),kob(ncrm),ipo(ncrm)
      logical free(ncrm),below,above

c      *Allocatable* arrays:
      real*8, dimension(:), allocatable :: qd,yg
      integer, dimension(:), allocatable :: isi
      integer, dimension(:,:), allocatable :: iby
      integer, dimension(:,:), allocatable :: icrtab
      integer*1, dimension(:), allocatable :: noctab

      allocate(qd(iymin:iymax),yg(iymin:iymax))
      allocate(isi(iymin:iymax),iby(iymin:iymax,0:1),icrtab(ngps,2))
      allocate(noctab(ngps))

c      For diagnostic purposes, keep track of max number of crossings,
c      nodes and contours across all PV levels:
      maxncr=0
      maxnpt=0
      maxnc =0

c      Initialisation:
      nby=1-iymin/nhf
      yoff=pi+twopi*dble(nby)

c      Grid box reference indices:
      nyf=iymax-iymin+1
      do ix=1,nhf
        ibx(ix,1)=nyf*(ix-1)
      enddo
      do ix=2,nhf
        ibx(ix,0)=ibx(ix-1,1)
      enddo
      ibx(1,0)=ibx(nhf,1)

      do iy=iymin,iymax
        iby(iy,0)=iy-iymin
        iby(iy,1)=iy-iymin+1
        yg(iy)=gl*dble(iy-1)-pi
      enddo

c---------------------------------------------------------------
      if (levbeg .le. levend) then
c      Loop over contour levels and process:
      do lev=levbeg,levend

c      Counter for total number of grid line crossings:
      ncr=0

c      Contour level being sought:
      qtmp=(dble(lev)-half)*dq-qoff

c      Find x grid line crossings first:
      do ix=1,nhf
        xgt=xg(ix)

        if (iymin .lt. 1) then
        do iy=iymin,0
          iy0=1+mod(nby*nhf+iy-1,nhf)
          qd(iy)=qa(iy0,ix)+bfac*dble(iy-iy0)-qtmp
        enddo
        endif

        do iy=1,nhf
          qd(iy)=qa(iy,ix)-qtmp
        enddo

        if (iymax .gt. nhf) then
        do iy=nhf+1,iymax
          iy0=1+mod(iy-1,nhf)
          qd(iy)=qa(iy0,ix)+bfac*dble(iy-iy0)-qtmp
        enddo
        endif

        do iy=iymin,iymax
          isi(iy)=sign(one,qd(iy))
        enddo

        do iy=iymin,iymax-1
          if (isi(iy) .ne. isi(iy+1)) then
            ncr=ncr+1
            inc=(1-isi(iy))/2
            kib(ncr)=iby(iy,1)+ibx(ix,inc)
            kob(ncr)=iby(iy,1)+ibx(ix,1-inc)
            u(ncr)=xgt
            v(ncr)=yg(iy)-gl*qd(iy)/(qd(iy+1)-qd(iy))
          endif
        enddo

      enddo

c      Above, kib = grid box into which the contour (containing icr) is going
c             kob =   "   "  out of "    "     "         "       "    " coming
c        [kob -> icr -> kib:  icr lies at the boundary between kob & kib]

c      Find y grid line crossings next:
      if (iymin .lt. 0) then
c      This is the region below the central domain:
      do iy=iymin+1,0
        ygt=yg(iy)
        iy0=1+mod(nby*nhf+iy-1,nhf)
        qadd=bfac*dble(iy-iy0)-qtmp

        do ix=1,nhf
          qd(ix)=qa(iy0,ix)+qadd
          isi(ix)=sign(one,qd(ix))
        enddo
        qd(nhf+1)=qa(iy0,1)+qadd
        isi(nhf+1)=sign(one,qd(nhf+1))

        do ix=1,nhf
          if (isi(ix) .ne. isi(ix+1)) then
            ncr=ncr+1
            inc=(1-isi(ix))/2
            kib(ncr)=ibx(ix,1)+iby(iy,1-inc)
            kob(ncr)=ibx(ix,1)+iby(iy,inc)
            v(ncr)=ygt
            u(ncr)=xg(ix)-gl*qd(ix)/(qd(ix+1)-qd(ix))
          endif
        enddo

      enddo
      endif

c      This is the central domain:
      do iy=1,nhf
        ygt=yg(iy)

        do ix=1,nhf
          qd(ix)=qa(iy,ix)-qtmp
          isi(ix)=sign(one,qd(ix))
        enddo
        qd(nhf+1)=qa(iy,1)-qtmp
        isi(nhf+1)=sign(one,qd(nhf+1))

        do ix=1,nhf
          if (isi(ix) .ne. isi(ix+1)) then
            ncr=ncr+1
            inc=(1-isi(ix))/2
            kib(ncr)=ibx(ix,1)+iby(iy,1-inc)
            kob(ncr)=ibx(ix,1)+iby(iy,inc)
            v(ncr)=ygt
            u(ncr)=xg(ix)-gl*qd(ix)/(qd(ix+1)-qd(ix))
          endif
        enddo

      enddo

      if (iymax .gt. nhf+1) then
c      This is the region above the central domain:
      do iy=nhf+1,iymax-1
        ygt=yg(iy)
        iy0=1+mod(iy-1,nhf)
        qadd=bfac*dble(iy-iy0)-qtmp

        do ix=1,nhf
          qd(ix)=qa(iy0,ix)+qadd
          isi(ix)=sign(one,qd(ix))
        enddo
        qd(nhf+1)=qa(iy0,1)+qadd
        isi(nhf+1)=sign(one,qd(nhf+1))

        do ix=1,nhf
          if (isi(ix) .ne. isi(ix+1)) then
            ncr=ncr+1
            inc=(1-isi(ix))/2
            kib(ncr)=ibx(ix,1)+iby(iy,1-inc)
            kob(ncr)=ibx(ix,1)+iby(iy,inc)
            v(ncr)=ygt
            u(ncr)=xg(ix)-gl*qd(ix)/(qd(ix+1)-qd(ix))
          endif
        enddo

      enddo
      endif

c      Initialise number of crossings per box:
      do i=1,ncr
        noctab(kob(i))=0
      enddo

      do icr=1,ncr
c        icr is the index of the current crossing at level lev.
        k=kob(icr)
c        accumulate number of crossings in this box:
        noctab(k)=noctab(k)+1
c        assign crossing to box, permitting 2 crossings:
        icrtab(k,noctab(k))=icr
      enddo

      do icr=1,ncr
        k=kib(icr)
        noc=noctab(k)
c        Use last crossing in this box as the next node:
        kob(icr)=icrtab(k,noc)
c        kob(icr) now gives the next point after icr
        noctab(k)=noc-1
c        This will normally be zero, except for boxes with 2 crossings;
c        this allows a second use of this box.
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
        ibeg=nptd+1
        ioff=ibeg-i1d(j)

c        See if contour at least partly lies in original domain:
        below=.false.
        above=.true.
        do i=i1d(j),i2d(j)
          icr=ipo(i)
          xd(ioff+i)=u(icr)
          yd(ioff+i)=v(icr)
          if (v(icr) .lt. -pi) below=.true.
          if (v(icr) .lt.  pi) above=.false.
        enddo

c        All contours having any portion with y < -pi are eliminated,
c        as are contours lying entirely above y = +pi:
        if ((npd(j) .lt. 5) .or. below .or. above) then
          nrem=nrem+1
        else
          nd=j-nrem
          npd(nd)=npd(j)
          i1d(nd)=ibeg
          nptd=nptd+npd(j)
          i2d(nd)=nptd
        endif

      enddo
      nd=ndt-nrem

c      Done rebuilding contours.
      if (nrem .lt. ndt) then
c      Remove contours which have periodic copies:

c-------------------------------------------------------
c      Diagnose max # of crossings, nodes and contours:
      maxncr=max(maxncr,ncr)
      maxnpt=max(maxnpt,nptd)
      maxnc =max(maxnc, nd)

c-------------------------------------------------------
c      Put all y coordinates between -pi and pi:
      do i=1,nptd
        yd(i)=mod(yd(i)+yoff,twopi)-pi
      enddo
c      yoff is an odd multiple of pi - see top of routine.

c-------------------------------------------------------
c      Remove points that are extremely close together:
      do i=1,nptd-1
        xx=xd(i+1)-xd(i)
        dx(i)=xx-twopi*int(xx*hlxi)
        yy=yd(i+1)-yd(i)
        dy(i)=yy-twopi*int(yy*hlyi)
      enddo

      do j=1,nd
        ibeg=i1d(j)
        iend=i2d(j)
        xx=xd(ibeg)-xd(iend)
        dx(iend)=xx-twopi*int(xx*hlxi)
        yy=yd(ibeg)-yd(iend)
        dy(iend)=yy-twopi*int(yy*hlyi)
      enddo

      do i=1,nptd
        v(i)=dx(i)**2+dy(i)**2
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

      if (nd .gt. 0) then
c------------------------------------------------------------------------
c      Redistribute points on the contours:
      do ireno=1,nreno
        call renode
      enddo

c      Write contour indices for this level:
      inda=lev-nlevm+(lev-1)/nlevm-1
      do j=1,nd
        i1a=npta+i1d(j)
        write(11) npd(j),i1a,inda,lz
      enddo
c      Write nodes:
      do i=1,nptd
        write(12) xd(i),yd(i)
      enddo

c      Augment number of contours and nodes:
      na=na+nd
      npta=npta+nptd

      endif
      endif
      enddo
c      End of loop over contour levels
      endif

c      For diagnostic purposes only:
      write(*,*) ' Max # crossings, nodes and contours of any level:'
      write(*,'(3(1x,i10))') maxncr,maxnpt,maxnc

      deallocate(qd,yg)
      deallocate(isi,iby,icrtab)
      deallocate(noctab)

      return
      end

c=======================================================================
      subroutine renode
c      Re-nodes each contour while preserving corner locations.
c      Uses square-root dependence on a weighted sum of nearby
c      curvature values.

      include 'econgen_com.i'

c------------------------------------------------------------------------
c      Set up next array:

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
        u(i)=half*(q(i)+q(ia))
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
            eta=p*(p*(b(i)+p*c(i))-b(i)-c(i))
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
      do i=1,nptd
        xd(i)=u(i)-twopi*int(u(i)*hlxi)
        yd(i)=v(i)-twopi*int(v(i)*hlyi)
      enddo

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

      include 'econgen_com.i'

c----------------------------------------------------------------------
      do i=1,nptd
        ia=next(i)
        xx=xd(ia)-xd(i)
        dx(i)=xx-twopi*int(xx*hlxi)
        yy=yd(ia)-yd(i)
        dy(i)=yy-twopi*int(yy*hlyi)
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
c      Set curvature to zero at corners:
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
        b(i)=half*(u(i)-c(i))
        c(i)=third*c(i)
      enddo

      return
      end

c==========================================================================
      subroutine ptospc(var)
c      Physical -> Spectral FFT

      implicit double precision(a-h,o-z)
      include 'econgen_dim.i'

      parameter (nhsq=nh*nh)
      dimension trig(2*nh),work(nhsq)

      dimension var(nh,nh)

      ifail=0
      call c06fpf(nh,nh,var,'i',trig,work,ifail)
      do ky=2,nh
        do ix=1,ky-1
          dum=var(ky,ix)
          var(ky,ix)=var(ix,ky)
          var(ix,ky)=dum
        enddo
      enddo
      call c06fpf(nh,nh,var,'s',trig,work,ifail)

      return
      end

c==========================================================================
      subroutine spctopf(var)
c      Spectral -> Physical FFT

      include 'econgen_com.i'

      dimension trigf(2*nhf),var(nhf,nhf)

      ifail=0
      call c06gqf(nhf,nhf,var,ifail)
      call c06fqf(nhf,nhf,var,'i',trigf,qa,ifail)
      do kx=2,nhf
        do iy=1,kx-1
          dum=var(iy,kx)
          var(iy,kx)=var(kx,iy)
          var(kx,iy)=dum
        enddo
      enddo
      call c06gqf(nhf,nhf,var,ifail)
      call c06fqf(nhf,nhf,var,'s',trigf,qa,ifail)

      return
      end
