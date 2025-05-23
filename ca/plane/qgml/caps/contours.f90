module contours

! Module containing all subroutines related to contour advection.

! Revised by D G Dritschel on 31 July 2020 to use any number of layers.
  
use constants

implicit none

 !PV contours:
double precision:: xq(npm),yq(npm),qavg(nz),qjump(nz)
integer:: indq(nm),npq(nm),i1q(nm),i2q(nm),layq(nm)
integer:: il1q(nz),il2q(nz),jl1q(nz),jl2q(nz)
integer:: nextq(npm),nptq,nq

 !Contour to grid conversion quantities (fine grid used in this module):
double precision:: xgf(ngf)
double precision,parameter:: glf=twopi/dble(ngf),glfi=dble(ngf)/twopi

 !Contour to grid conversion quantities (ultra-fine grid used in congen):
double precision:: xgu(ngu),ygu(ngu)
double precision,parameter:: glu=twopi/dble(ngu),glui=dble(ngu)/twopi

 !Next grid points (accounting for periodicity) used in interpolation:
integer:: igp1(ng)

 !Area weights used for interpolation (used only in module congen):
double precision:: w00(mgu,mgu),w10(mgu,mgu),w01(mgu,mgu),w11(mgu,mgu)
integer:: igfw(ngu),ig0w(ngu),ig1w(ngu)

 !Grid box reference indices used in congen:
integer:: ibg(ngu+1)

 !Surgery & node redistribution quantities:
double precision,parameter:: amu=0.2d0,ell=6.25d0*gl
double precision,parameter:: dm=amu**2*ell/four,dm2=dm**2,d4small=small*gl**4
double precision,parameter:: dmsq=four*dm2,dmi=two/dm
double precision,parameter:: elf=one/ell**2,densf=one/(amu*sqrt(ell))


!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains 

!=======================================================================

subroutine init_contours
! Initialises all quantities needed for contour advection.

implicit none

double precision:: fac,pxf,pxc,pyf,pyc
integer:: ig,igg,ix,ixf,iyf

!------------------------------------------------------------------------
 !Next grid points used in velocity interpolation (velint) and elsewhere; 
 !enforces periodicity in x and y:
do ig=1,ngm1
  igp1(ig)=ig+1
enddo
igp1(ng)=1

 !Fine grid lines needed for contour-to-grid conversion (con2grid):
do ix=1,ngf
  xgf(ix)=glf*dble(ix-1)-pi
enddo

!==========================================================================
 !Initialise ultra-fine grid lines needed in contour regeneration (congen):
do ix=1,ngu
  xgu(ix)=glu*dble(ix-1)-pi
enddo
ygu=xgu

 !Grid box reference indices used in congen:
ibg(1)=ngu-1
do ig=2,ngu+1
  ibg(ig)=ig-2
enddo

!----------------------------------------------------------------------
 !Area weights for interpolation of a gridded field onto the ultra-fine 
 !horizontal grid (congen):
fac=one/dble(mgu)

do ixf=1,mgu
  pxf=fac*dble(ixf-1)
  pxc=one-pxf

  do iyf=1,mgu
    pyf=fac*dble(iyf-1)
    pyc=one-pyf

    w00(iyf,ixf)=pyc*pxc
    w10(iyf,ixf)=pyf*pxc
    w01(iyf,ixf)=pyc*pxf
    w11(iyf,ixf)=pyf*pxf
  enddo
enddo

do ig=1,ngu-mgu
  igg=(ig-1)/mgu
  ig0w(ig)=1+igg
  ig1w(ig)=2+igg
  igfw(ig)=ig-mgu*igg
enddo
do ig=ngu-mgu+1,ngu
  ig0w(ig)=ng
  ig1w(ig)=1
  igfw(ig)=ig-mgu*ngm1
enddo

return
end subroutine init_contours

!=======================================================================
      
subroutine velint(uu,vv,uq,vq)

! Bi-linearly interpolates the velocity (uu,vv) in each layer j to the
! contour nodes (xq,yq) in that layer and stores the result in (uq,vq).

implicit none

double precision:: xx,px,pxc,yy,py,pyc
integer:: iz,i,ix0,ix1,iy0,iy1

 !Passed arrays:
double precision:: uu(ng,ng,nz),vv(ng,ng,nz)
double precision:: uq(nptq),vq(nptq)

do iz=1,nz
  if (jl2q(iz) > 0) then
    do i=il1q(iz),il2q(iz)
      xx=gli*(xq(i)+pi)
      ix0=1+int(xx)
      pxc=dble(ix0)-xx
      px=one-pxc
      ix1=igp1(ix0)

      yy=gli*(yq(i)+pi)
      iy0=1+int(yy)
      pyc=dble(iy0)-yy
      py=one-pyc
      iy1=igp1(iy0)

      uq(i)=pyc*(pxc*uu(iy0,ix0,iz)+px*uu(iy0,ix1,iz)) &
            +py*(pxc*uu(iy1,ix0,iz)+px*uu(iy1,ix1,iz))

      vq(i)=pyc*(pxc*vv(iy0,ix0,iz)+px*vv(iy0,ix1,iz)) &
            +py*(pxc*vv(iy1,ix0,iz)+px*vv(iy1,ix1,iz))
    enddo
  endif
enddo

return
end subroutine velint

!=======================================================================

subroutine con2grid(qc)
! Contour -> grid conversion in each layer.
! The gridded field is returned in the array qc.

implicit none

 !Passed array:
double precision:: qc(ng,ng,nz)

 !Local arrays and variables:
double precision:: qod0(ngf/2),qod1(ngf/2),qod2(ngf/2)
double precision:: qev0(0:ngf/2),qev1(0:ngf/2),qev2(0:ngf/2)
double precision:: qa(ngf+1,ngf)
double precision:: qjx(ngf),qbot(ngf)
double precision:: dx(nptq),dy(nptq)
double precision:: xx,yy,cc,p,sqjump,px0,qavg0,qadd
integer:: ixc(nptq),ngc(nptq)
integer:: iz,i,ia,ix,iy,ixdif,ixbeg,jump,ncr,ngh,ngff,mix,mixm1,miy
logical:: crossx(nptq)

if (nptq == 0) then 
   !No contours to convert: return qc = 0:
  qc=zero
  return
endif

 !Begin a major loop over layers
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do iz=1,nz
  if (jl2q(iz) > 0) then
     !There are contours in this layer.  Initialise interior x grid line 
     !crossing information and fill the PV jump array along lower boundary:
    do i=il1q(iz),il2q(iz)
      ixc(i)=1+int(glfi*(xq(i)+pi))
    enddo

    qjx=zero
    do i=il1q(iz),il2q(iz)
      ia=nextq(i)
      xx=xq(ia)-xq(i)
      dx(i)=xx-twopi*dble(int(xx*pinv))
      yy=yq(ia)-yq(i)
      dy(i)=yy-twopi*dble(int(yy*pinv))
      ixdif=ixc(ia)-ixc(i)
      ngc(i)=ixdif-ngf*((2*ixdif)/ngf)
      crossx(i)=(ngc(i) /= 0)
      if (abs(yy) > pi) then
         !The contour segment (i,ia) crosses y = -pi; find x location:
        cc=sign(one,dy(i))
        p=-(yq(i)-pi*cc)/(dy(i)+small)
        ix=1+int(glfi*(mod(thrpi+xq(i)+p*dx(i),twopi)))
        qjx(ix)=qjx(ix)-qjump(iz)*cc
         !Note: qjx gives the jump going from ix to ix+1
      endif
    enddo

     !Sum q jumps to obtain the gridded q along lower boundary:
    qbot(1)=zero
     !Corner value cannot be determined a priori; qavg is used for this below
    do ix=1,ngf-1
      qbot(ix+1)=qbot(ix)+qjx(ix)
    enddo

     !Initialise interior q jump array:
    qa(2:ngf+1,:)=zero

     !Determine x grid line crossings and accumulate q jumps:
    do i=il1q(iz),il2q(iz)
      if (crossx(i)) then
        jump=sign(1,ngc(i))
        ixbeg=ixc(i)+(jump-1)/2+ngf
        sqjump=qjump(iz)*sign(one,dx(i))
        ncr=0
        do while (ncr /= ngc(i)) 
          ix=1+mod(ixbeg+ncr,ngf)
          xx=xgf(ix)-xq(i)
          px0=(xx-twopi*dble(int(xx*pinv)))/dx(i)
           !The contour crossed the fine grid line ix at the point
           !   x = xq(i) + px0*dx(i) and y = yq(i) + px0*dy(i):
          yy=yq(i)+px0*dy(i)
          iy=2+int(glfi*(yy-twopi*dble(int(yy*pinv))+pi))
           !Increment q jump between the grid lines iy-1 & iy:
          qa(iy,ix)=qa(iy,ix)+sqjump
           !Go on to consider next x grid line (if there is one):
          ncr=ncr+jump
        enddo
      endif
    enddo

     !Get q values by sweeping through y:
    do ix=1,ngf 
      qa(1,ix)=qbot(ix)
      do iy=2,ngf
        qa(iy,ix)=qa(iy,ix)+qa(iy-1,ix)
      enddo
    enddo

     !Average the PV field in qa to the coarser grid (ng,ng):
    ngh=ngf
    do while (ngh > ng)
      ngff=ngh
      ngh=ngh/2
       !Perform nine-point averaging:
      do iy=1,ngh
        miy=2*iy
        qod2(iy)=qa(miy-1,ngff)
        qev2(iy)=qa(miy,ngff)
      enddo
      qev2(0)=qa(ngff,ngff)
      do ix=1,ngh
        mix=2*ix
        mixm1=mix-1
        do iy=1,ngh
          miy=2*iy
          qod1(iy)=qa(miy-1,mixm1)
          qod0(iy)=qa(miy-1,mix)
          qev1(iy)=qa(miy,mixm1)
          qev0(iy)=qa(miy,mix)
        enddo
        qev1(0)=qev1(ngh)
        qev0(0)=qev0(ngh)
        do iy=1,ngh
          qa(iy,ix)=0.0625d0*(qev0(iy)+qev0(iy-1)+qev2(iy)+qev2(iy-1)) &
                    +0.125d0*(qev1(iy)+qev1(iy-1)+qod0(iy)+qod2(iy)) &
                      +0.25d0*qod1(iy)
        enddo
        do iy=1,ngh
          qod2(iy)=qod0(iy)
          qev2(iy)=qev0(iy)
        enddo
        qev2(0)=qev0(0)
      enddo
    enddo

     !Restore average (qavg(iz)):
    qavg0=danorm*sum(qa(1:ng,1:ng))
    qadd=qavg(iz)-qavg0
    qc(:,:,iz)=qa(1:ng,1:ng)+qadd
     !Now qc has the correct average in this layer

  else
     !There are no contours in this layer; assign qc = qavg(iz):
    qc(:,:,iz)=qavg(iz)
  endif
enddo
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Ends the major loop over layers

return
end subroutine con2grid

!=======================================================================

subroutine surgery
! Contour surgery (D. G. Dritschel, J. Comput. Phys. 77, 477--483, 1988).

! Contours are represented by nodes (xq(i),yq(i)), i = 1, ..., nptq, 
! and nextq(i) gives the index of the node following i.
! indq(j) is an integer index distinguishing different contour levels.
! layq(j) gives the layer containing contour j, with contours j = jl1q(iz)
! to jl2q(iz) occupying layer iz and nodes i = il1q(iz) to il2q(iz) also 
! occupying this layer.
! npq(j) is the number of points on a contour, i1q(j) is the starting index 
! of a contour and i2q(j) is the ending index.
! nq is the number of contours and nptq is the total number of nodes.

! After surgery is complete for a given level indq and layer layq, 
! this routine redistributes nodes (renode).

! Major revision 01/01/2001 by D. G. Dritschel to accelerate surgery
! using local boxes to minimise search costs.  All surgery is now
! done in a single way.

! Revised from the channel geometry stratified case 30 Sep 2013. 

implicit none

 !Local variables:
integer,parameter:: ngbs=ng*ng
double precision:: xa(npm),ya(npm)
double precision:: xd(nprm),yd(nprm)
double precision:: dx(nptq),dy(nptq),dsq(nptq)
double precision:: xx,yy,fnq,fnb,bwi,dnb
double precision:: delx,dely,aa,cc,dxa,dya
integer:: i1a(nm),i2a(nm),inda(nm)
integer:: jq1(nlevm),jq2(nlevm),iq1(nlevm),iq2(nlevm),levq(nlevm)
integer:: nspb(ngbs),kb1(ngbs),kb2(ngbs)
integer:: loc(nplm),list(nplm),node(nplm)
integer:: nexta(npm),icre(nm)
integer:: lev,levp,levt,nlev,iz,j,j1,j2
integer:: i,ib,ie,is,isb,ibeg,iend,nptlev,nseg
integer:: nb,nbox,mbx1,mbx2,mby1,mby2,mbx,mmbx,mb,mby,k,ks,npd
logical:: avail(nptq)

 !If there are no contours, there is nothing to do:
if (nq == 0) return

 !Otherwise, first get work arrays for efficient surgery:
do i=1,nptq-1
  xx=xq(i+1)-xq(i)
  dx(i)=xx-twopi*dble(int(xx*pinv))
  yy=yq(i+1)-yq(i)
  dy(i)=yy-twopi*dble(int(yy*pinv))
  dsq(i)=dx(i)**2+dy(i)**2
  avail(i)=.true.
enddo
avail(nptq)=.true.

do j=1,nq
  ie=i2q(j)
  is=i1q(j)
  xx=xq(is)-xq(ie)
  dx(ie)=xx-twopi*dble(int(xx*pinv))
  yy=yq(is)-yq(ie)
  dy(ie)=yy-twopi*dble(int(yy*pinv))
  dsq(ie)=dx(ie)**2+dy(ie)**2
enddo

 !Zero out counters for contours and nodes for building new contours:
nq=0
nptq=0

 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !              Begin a major loop over layers
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do iz=1,nz
  if (jl2q(iz) > 0) then
     !There are contours in this layer.  Calculate beginning and 
     !ending contours (jq1,jq2) for each distinct value of ind:
    j1=jl1q(iz)
    j2=jl2q(iz)

    levp=indq(j1)
    nlev=1
    jq1(1)=j1
    levq(1)=levp
    do j=j1+1,j2
      lev=indq(j)
      if (lev > levp) then
        jq2(nlev)=j-1
        nlev=nlev+1
        jq1(nlev)=j
        levq(nlev)=lev
        levp=lev
      endif
    enddo
    jq2(nlev)=j2
     !Note: levq(lev) gives the q level (indq) of contours.

    do lev=1,nlev
      iq1(lev)=i1q(jq1(lev))
      iq2(lev)=i2q(jq2(lev))
    enddo

     !Reset initial layer counters for contours and nodes:
    jl1q(iz)=nq+1
    il1q(iz)=nptq+1

     !-----------------------------------------------------------------
     !              Begin a major loop over q levels
     !-----------------------------------------------------------------
    do lev=1,nlev
      ibeg=iq1(lev)
      iend=iq2(lev)

       !Work out optimal box number, nbox, for fast surgery:
      nptlev=iend-ibeg+1
       !nptlev: number of nodes having this q level.
      fnq=dble(nptlev)

       !Balance boxing costs, (3+15*nseg/nbox)*nbox, with surgery 
       !search costs, 24*nptlev*nseg/nbox, to work out optimal box number;
       !here we estimate nseg, the total number of segments counted
       !in all boxes, as (1+3*eps*nb)*nptlev, where nb=sqrt(nbox);
       !this has been verified in realistically complex tests).
       !This balance results in a quartic equation for a = nb/sqrt(nptlev)
       !whose solution only depends on r = 3*eps*sqrt(nptlev), where
       !eps=(average distance between adjacent nodes)/(domain area).
       !Fortunately, a only varies from 1.129 to 1.265 over the entire
       !range of r (from 0 to infinity).  Here, therefore, we simply
       !take a = 1.2, i.e. nb = 1.2*sqrt(nptlev), as an approximation:
      fnb=1.2d0*sqrt(fnq)
      nb=max(min(nint(fnb),ng),1)
       !nb: number of boxes in the x & y directions; nb <= ng.
      nbox=nb*nb
       !nbox: total number of boxes
      bwi=oms*dble(nb)/twopi

       !Box all segments [i,nextq(i)] to reduce search costs in surgery:
      do mb=1,nbox
        nspb(mb)=0
      enddo
      nseg=0

      dnb=dble(nb)
      do ib=ibeg,iend
        i=nextq(ib)
         !Find range of boxes spanned by the segment (i,nextq(i)):
        if (dx(i) > zero) then
          mbx1=int(dnb+bwi*(xq(i)+pi))
          mbx2=int(dnb+bwi*(xq(i)+dx(i)+pi))
        else
          mbx1=int(dnb+bwi*(xq(i)+dx(i)+pi))
          mbx2=int(dnb+bwi*(xq(i)+pi))
        endif
        if (dy(i) > zero) then
          mby1=int(dnb+bwi*(yq(i)+pi))
          mby2=int(dnb+bwi*(yq(i)+dy(i)+pi))
        else
          mby1=int(dnb+bwi*(yq(i)+dy(i)+pi))
          mby2=int(dnb+bwi*(yq(i)+pi))
        endif
         !mbx1+1,mbx2+1 is the x range of boxes spanned by the segment while
         !mby1+1,mby2+1 is the y range.  

        do mbx=mbx1,mbx2
          mmbx=mod(mbx,nb)
          do mby=mby1,mby2
            mb=nb*mmbx+mod(mby,nb)+1
            nspb(mb)=nspb(mb)+1
             !nspb counts number of segments that cross through box mb.
            nseg=nseg+1
             !nseg counts the total number of segments crossing all boxes.
            loc(nseg)=mb
             !loc gives the box location of the segment
            list(nseg)=ib
             !list gives the node before the one at the segment origin (i).
          enddo
        enddo
      enddo
       !nseg: the total number of segment crossings found for all boxes.
       !**NB: nseg/nbox = average number of segments crossing a box.

      kb1(1)=1
      do mb=1,nbox-1
        kb1(mb+1)=kb1(mb)+nspb(mb)
      enddo
      do mb=1,nbox
        kb2(mb)=kb1(mb)-1
      enddo

      do k=1,nseg
        mb=loc(k)
        ks=kb2(mb)+1
        node(ks)=list(k)
        kb2(mb)=ks
      enddo

       !With the above, segments [node(kb1(mb)),nextq(node(kb1(mb))],
       ![node(kb1(mb)+1),nextq(node(kb1(mb)+1)], ..., 
       ![node(kb2(mb)),nextq(node(kb2(mb))] cross through box mb.
    
       !Note: roughly (3+15*nseg/nbox)*nbox operations are required
       !      up to this point (in the loop over q levels) to do
       !      all the segment boxing.

       !-------------------------------------------------------------
       !Now search for surgery with segments crossing through the box 
       !containing i:
      do ib=ibeg,iend
        i=nextq(ib)
        mb=nb*int(bwi*(xq(i)+pi))+int(bwi*(yq(i)+pi))+1
         !mb: box containing the node i.
        do k=kb1(mb),kb2(mb)
          isb=node(k)
          is=nextq(isb)
           !Segment [is,nextq(is)] lies in box mb.  Exclude segments 
           !ending in a edge or having node i at either endpoint:
          if ((is-ib)*(is-i) == 0) cycle
           !nextq(is) could lie at an edge as a result of previous surgery.

           !Next see if the node i lies within a distance dm to the line 
           !segment between is and nextq(is):
          xx=xq(is)-xq(i)
          delx=xx-twopi*dble(int(xx*pinv))
          yy=yq(is)-yq(i)
          dely=yy-twopi*dble(int(yy*pinv))
          aa=delx*dx(is)+dely*dy(is)

           !Note: roughly 24*nptlev*(nseg/nbox) operations are required 
           !up to the following statement (which is rarely satisfied); 
           !this is assumed to be the dominant cost of surgery.  
           !Even if each node i surgically reconnects on average 
           !once, the above estimate holds.

          if (aa*(aa+dsq(is))+d4small < zero) then
             !Passing this condition, a perpendicular can be dropped 
             !onto the line segment from is to nextq(is).  
             ![Tests show that this is satisfied only 6.4% of the time.]
             !Check distance to line segment:
            cc=(delx*dy(is)-dely*dx(is))**2
             !NB: sqrt(cc/dsq(is)) is the distance to the segment and 
             !    dm2=dm**2 below.  This distance must be strictly 
             !    positive, yet less than dm, to permit surgery:
            if (cc*(cc-dm2*dsq(is)) < zero) then
               !Surgery is now possible between node i and the segment 
               ![is,nextq(is)].
               ![Tests show that this is satisfied only 0.45% of the time.]
               !Move node i to a point midway between i and either is or 
               !nextq(is), whichever is closer:
              if (dsq(is)+two*aa < zero) then
                 !node i is closest to node nextq(is):
                dxa=f12*(delx+dx(is))
                dya=f12*(dely+dy(is))
                isb=is
                is=nextq(is)
              else
                 !node i is closest to node is:
                dxa=f12*delx
                dya=f12*dely
              endif
               !Move node i & is to a common node; first deal with node i:
              xx=xq(i)+dxa
              xq(i)=oms*(xx-twopi*dble(int(xx*pinv)))
              yy=yq(i)+dya
              yq(i)=oms*(yy-twopi*dble(int(yy*pinv)))
              dx(i)=dx(i)-dxa
              dy(i)=dy(i)-dya
              dsq(i)=dx(i)**2+dy(i)**2
               !isb becomes the node before i:
              nextq(isb)=i
              dx(isb)=dx(isb)-dxa
              dy(isb)=dy(isb)-dya
              dsq(isb)=dx(isb)**2+dy(isb)**2

               !Now deal with node is:
              xq(is)=xq(i)
              yq(is)=yq(i)
              dx(is)=dx(is)+dxa
              dy(is)=dy(is)+dya
              dsq(is)=dx(is)**2+dy(is)**2
               !ib becomes the node before is:
              nextq(ib)=is
              dx(ib)=dx(ib)+dxa
              dy(ib)=dy(ib)+dya
              dsq(ib)=dx(ib)**2+dy(ib)**2

               !Now update i and look for further possible surgery:
              i=is
              if (i == ib) exit
               !i = ib when a contour consists of a single node
            endif
          endif
        enddo
      enddo

       !----------------------------------------------------------------
       !It remains to rebuild the contours using the nextq() information

       !Current contour level:
      levt=levq(lev)

      do i=ibeg,iend
        if (avail(i)) then
         !avail(i) is true if node i has not yet been associated with
         !a contour.  Start a new contour:
          nq=nq+1

         !First point on the contour:
          npd=1
          xd(1)=xq(i)
          yd(1)=yq(i)

          avail(i)=.false.
          is=nextq(i)
           !Generate the contour using nextq() and finish when the
           !original node (i) is reached:
          do while (is /= i) 
            npd=npd+1
            xd(npd)=xq(is)
            yd(npd)=yq(is)
            avail(is)=.false.
            is=nextq(is)
          enddo

           !Re-distribute nodes on this contour:
          npq(nq)=npd
          if (npd > 3) &
               call renode(xd,yd,npd,xa(nptq+1),ya(nptq+1),npq(nq))
          if (npq(nq) > 3) then
            i1a(nq)=nptq+1
            nptq=nptq+npq(nq)
            i2a(nq)=nptq
            inda(nq)=levt
            layq(nq)=iz
            do is=i1a(nq),nptq-1
              nexta(is)=is+1
            enddo
            nexta(nptq)=i1a(nq)
          else
             !Delete contour if deemed too small (see renode):
            nq=nq-1
          endif

        endif
      enddo

       !Now go back and consider the next q level:
    enddo
     !--------------------------------------------------------------------
     !This ends the loop over q levels; surgery is complete for this layer
     !--------------------------------------------------------------------

     !Reset final layer counters for contours and nodes:
    jl2q(iz)=nq
    il2q(iz)=nptq
  endif
enddo
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !          The above ends the major loop over layers
 !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 !Copy nodes back to those in the argument of the subroutine:
do i=1,nptq
  xq(i)=xa(i)
  yq(i)=ya(i)
  nextq(i)=nexta(i)
enddo
 !nextq(i) is not altered until surgery.

do j=1,nq
  i1q(j)=i1a(j)
  i2q(j)=i2a(j)
  indq(j)=inda(j)
enddo

return
end subroutine surgery

!==========================================================================

subroutine renode(xd,yd,npd,xr,yr,npr)
! Re-nodes a single closed contour (xd(i),yd(i)), i = 1,...,npd 
!   and returns the new contour in (xr(i),yr(i)), i = 1,...,npr

! Note: a closed contour closes on itself

! If npr = 0, the contour is too small and should be removed

! -----------------------------------------------------------
! Intended for a rectangular domain with free-slip boundaries
! -----------------------------------------------------------

implicit none

 !Passed arrays and indices:
double precision:: xd(nprm),yd(nprm),xr(nprm),yr(nprm)
integer:: npd,npr

 !Local variables:
double precision:: dx(npd),dy(npd),a(npd),b(npd),c(npd),d(npd),e(npd)
double precision:: u(npd),v(npd)
double precision:: xx,yy,ww,sumden,fac,acc,p,eta
integer:: node(npd+1),next(0:npd)
integer:: i,ia,ib,ncorn,ibeg,iend,k,npseg,im
logical:: corner(npd)

!------------------------------------------------------------------
 !Define the node following a node:
do i=0,npd-1
  next(i)=i+1
enddo
next(npd)=1

 !Compute the contour increments:
do i=1,npd
  ia=next(i)
  xx=xd(ia)-xd(i)
  dx(i)=xx-twopi*dble(int(xx*pinv))
  yy=yd(ia)-yd(i)
  dy(i)=yy-twopi*dble(int(yy*pinv))
enddo

 !Calculate the cubic interpolation coefficients:
do i=1,npd
  v(i)=dx(i)*dx(i)+dy(i)*dy(i)
  e(i)=sqrt(v(i))
enddo

do ib=1,npd
  i=next(ib)
  u(i)=v(ib)
  a(i)=-dx(ib)
  c(i)=-dy(ib)
enddo

ncorn=0
do i=1,npd
  corner(i)=dx(i)*a(i)+dy(i)*c(i) > zero
  if (corner(i)) then 
     !Keep track of corner locations for use in renoding below:
    ncorn=ncorn+1
    node(ncorn)=i
     !Set curvature to zero at corners:
    b(i)=zero
  else
    b(i)=(dx(i)*c(i)-a(i)*dy(i))/ &
     & sqrt((a(i)*v(i)-dx(i)*u(i))**2+(c(i)*v(i)-dy(i)*u(i))**2+small3)
  endif
enddo

do i=1,npd
  ia=next(i)
  u(i)=e(i)*(b(ia)+b(i))
  c(i)=e(i)*(b(ia)-b(i))
enddo

 !Calculate the cubic interpolation coefficients:
do i=1,npd
  a(i)=f16*c(i)-f12*u(i)
  b(i)=f12*(u(i)-c(i))
  c(i)=f13*c(i)
enddo

!------------------------------------------------------------------------
 !Use the spherical curvature expression (radius of the sphere = ell)
 !to ensure an adequate node density in low curvature regions.
do i=1,npd
  ww=one/(v(i)+dmsq)
  u(i)=ww*sqrt(elf*v(i)+u(i)**2)
  v(i)=ww*e(i)
enddo
 !NB: elf = 1/ell**2; v(i) = |xx_{i+1}-xx_{i}|**2; e(i)=sqrt{v(i)};
 !    u(i)/e(i) = (kappa_{i}+kappa_{i+1})/2; dmsq = (2*dm)**2

 !Re-assign curvature at a node from weighted average on either side
 !(v above is the weight):
do ib=1,npd
  i=next(ib)
  d(i)=(u(ib)+u(i))/(v(ib)+v(i)+small)
enddo

 !Re-average to get interval value (effectively, four curvature
 !values go into getting the final interval value, u(i)):
do i=1,npd
  ia=next(i)
  u(i)=f12*(d(i)+d(ia))
enddo

 !Compute fractional number of nodes to be placed between old
 !nodes i and i+1:
do i=1,npd
  e(i)=e(i)*min(dmi,densf*sqrt(u(i))+u(i))
enddo
 !NB: dmi = 2/delta; densf = 1/(amu*sqrt{ell})

sumden=zero
do i=1,npd
  sumden=sumden+e(i)
enddo
 !Number of points on renoded contour:
npr=nint(sumden)+1
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
  fac=dble(npr)/sumden
  do i=1,npd
    e(i)=fac*e(i)
  enddo

   !The first node is fixed; find the remaining node positions:
  xr(1)=xd(1)
  yr(1)=yd(1)
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
    xx=xd(i)+p*dx(i)-eta*dy(i)
    xr(im)=oms*(xx-twopi*dble(int(xx*pinv)))
    yy=yd(i)+p*dy(i)+eta*dx(i)
    yr(im)=oms*(yy-twopi*dble(int(yy*pinv)))
  enddo

else if (ncorn == 1) then
   !A single corner - start new contour at the corner:

   !Make the sum of e(i) equal to npr:
  fac=dble(npr)/sumden
  do i=1,npd
    e(i)=fac*e(i)
  enddo

   !The first node (the corner) is fixed; find the remaining nodes:
  i=node(1)
  xr(1)=xd(i)
  yr(1)=yd(i)
  acc=zero
  i=i-1
  do im=2,npr
    do while (acc < one)
      i=next(i)
      acc=acc+e(i)
    enddo
    acc=acc-one
    p=one-acc/e(i)
    eta=p*(a(i)+p*(b(i)+p*c(i)))
    xx=xd(i)+p*dx(i)-eta*dy(i)
    xr(im)=oms*(xx-twopi*dble(int(xx*pinv)))
    yy=yd(i)+p*dy(i)+eta*dx(i)
    yr(im)=oms*(yy-twopi*dble(int(yy*pinv)))
  enddo

else
   !Multiple corners - start new contour at the first corner and
   !preserve locations of all corners:
  npr=0
  node(ncorn+1)=node(1)
  ibeg=node(1)

  do k=1,ncorn
    iend=node(k+1)
    i=ibeg
    sumden=zero
    do while (i /= iend)
      sumden=sumden+e(i)
      i=next(i)
    enddo

    if (sumden > small) then
       !Adequate spacing exists between corners to renode this segment

       !Number of points on renoded contour segment:
      npseg=nint(sumden)+1

       !Make the sum of e(i) equal to npseg
      fac=dble(npseg)/sumden
      i=ibeg
      do while (i /= iend)
        e(i)=fac*e(i)
        i=next(i)
      enddo

       !Fix the first node along the segment as the corner:
      xr(npr+1)=xd(ibeg)
      yr(npr+1)=yd(ibeg)

      if (npseg > 1) then
         !Find the remaining points along this segment:
        acc=zero
        i=ibeg-1
        do im=npr+2,npr+npseg
          do while (acc < one)
            i=next(i)
            acc=acc+e(i)
          enddo
          acc=acc-one
          p=one-acc/e(i)
          eta=p*(a(i)+p*(b(i)+p*c(i)))
          xx=xd(i)+p*dx(i)-eta*dy(i)
          xr(im)=oms*(xx-twopi*dble(int(xx*pinv)))
          yy=yd(i)+p*dy(i)+eta*dx(i)
          yr(im)=oms*(yy-twopi*dble(int(yy*pinv)))
        enddo
      endif

      npr=npr+npseg

    endif
    ibeg=iend

     !go on and consider the next segment (if any) between corners:
  enddo
endif

return
end subroutine renode

!=======================================================================

 !Main end module
end module contours
