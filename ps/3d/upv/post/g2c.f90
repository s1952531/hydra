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
integer,parameter:: npm=100*ng*ng, nplm=npm, nm=npm/20
integer,parameter:: nlevm=2000

integer,parameter:: nxp1=nx+1, nyp1=ny+1
integer,parameter:: nxp2=nx+2, nyp2=ny+2

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
double precision:: qa(nyp1,nxp1)
double precision:: xg(0:nx),yg(0:ny)
double precision:: qlev(2*nlevm),qamin,qamax
integer:: ibx(0:nx,0:1),iby(0:ny,0:1)

 !Basic parameters:
double precision:: amu,ell,elf,densf,dm,dmsq,dmi
double precision:: glx,gly,dq,dqi,qoff
double precision:: xbeg,xend,ybeg,yend
double precision:: xwid,hlxi,xavg,ywid,hlyi,yavg
real:: tr4,q3dr4(ng,ng,0:nz)
integer:: na,nd,npta,nptd
integer:: levbegall,levendall
integer:: levbeg,levend
integer:: loop,lz,iopt

 !Control variables:
logical:: corn(npm),last

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

write(*,*) ' Choose one of the following 3D fields to contour:'
write(*,*) ' (1) ql'
write(*,*) ' (2) d'
write(*,*) ' (3) g'
write(*,*) ' (4) r'
write(*,*) ' (5) w'
write(*,*) ' (6) pn'
write(*,*)
write(*,*) ' Option?'
read(*,*) iopt

write(*,*) ' Time to image?'
read(*,*) tr4
loop=nint(tr4/tgsave)+1

 !Open data file and read:
if (iopt .eq. 1) then
  open(31,file='3d/ql.r4',form='unformatted',access='direct', &
                        status='old',recl=ntbytes)
else if (iopt .eq. 2) then
  open(31,file= '3d/d.r4',form='unformatted',access='direct', &
                        status='old',recl=ntbytes)
else if (iopt .eq. 3) then
  open(31,file= '3d/g.r4',form='unformatted',access='direct', &
                        status='old',recl=ntbytes)
else if (iopt .eq. 4) then
  open(31,file= '3d/r.r4',form='unformatted',access='direct', &
                        status='old',recl=ntbytes)
else if (iopt .eq. 5) then
  open(31,file= '3d/w.r4',form='unformatted',access='direct', &
                        status='old',recl=ntbytes)
else if (iopt .eq. 6) then
  open(31,file='3d/pn.r4',form='unformatted',access='direct', &
                        status='old',recl=ntbytes)
else
  write(*,*) ' Not a valid option!  Exiting!'
  stop
endif

read(31,rec=loop) tr4,q3dr4
close(31)

qamin=minval(q3dr4)
qamax=maxval(q3dr4)

 !Define fixed arrays and constants:
call init(qamin,qamax)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Process each layer in turn:
do lz=0,nz
   !Extract field at this level:
  qa(1:ny,1:nx)=dble(q3dr4(:,:,lz))

   !Add an extra periodic column and row:
  qa(1:ny,nxp1)=qa(1:ny,1)
  qa(nyp1,1:nxp1)=qa(1,1:nxp1)

   !Generate contours (xd,yd):
  call grid2con

   !Store the contours (xd,yd) into (xa,ya):
  call storecon(lz)

enddo
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !Ends loop over layers

 !Write contours to congen.asc:
call writecon

write(*,*) 
write(*,*) ' To image these data, type'
write(*,*) 
write(*,*) ' python pc3d.py'
write(*,*) 

contains
 !Subroutine definitions follow:

!=======================================================================

subroutine init(fmin,fmax)
 !Defines various fixed arrays and constants.

double precision:: fmin,fmax,qbar,sepmax
integer:: ix,iy,lev

xbeg=-pi
xend=pi
xwid=xend-xbeg

hlxi=one/(f12*xwid+small)
xavg=f12*(xbeg+xend)

ybeg=-pi
yend=pi
ywid=yend-ybeg

hlyi=one/(f12*ywid+small)
yavg=f12*(ybeg+yend)

write(*,*)
write(*,'(a,f15.9)') ' Minimum field value = ',fmin
write(*,'(a,f15.9)') ' Maximum field value = ',fmax

write(*,*)
write(*,*) ' Contours are found for the field values qbar +/- dq/2,'
write(*,*) ' qbar +/- 3*dq/2, qbar +/- 5*dq/2, etc.'
write(*,*)
write(*,*) ' Enter qbar and dq:'
read(*,*) qbar,dq

dqi=one/dq
qoff=dq*dble(nlevm)-qbar
 !qoff: should be a large integer multiple of the 
 !      contour interval, dq.  The multiple should exceed 
 !      the maximum expected number of contour levels.
do lev=1,2*nlevm
  qlev(lev)=(dble(lev)-f12)*dq-qoff
enddo

glx=xwid/dble(nx)
gly=ywid/dble(ny)

 !Node redistribution parameters (see renode):
sepmax=1.25d0
amu=0.2d0
ell=four*sepmax**2*sqrt(glx*gly)
dm=ell*amu**2/four
dmsq=four*dm**2
dmi=two/dm
elf=one/ell**2
densf=one/(amu*sqrt(ell))

 !Define x & y grid lines (see grid2con):
do ix=1,nx
  xg(ix)=glx*dble(ix-1)+xbeg
enddo

do iy=1,ny
  yg(iy)=gly*dble(iy-1)+ybeg
enddo

 !Define grid-box reference indices:
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

 !Initialise counters:
na=0
npta=0
levbegall=2*nlevm
levendall=-2*nlevm

return
end subroutine init

!=======================================================================

subroutine grid2con
 !Generates new contours (xd,yd) from the 2D gridded data (qa) 

integer,parameter:: ncrm=2*npm/nz
 !ncrm: max number of contour crossings with grid lines

double precision:: ycr(ncrm),xcr(ncrm)
double precision:: xgt,ygt,fac
integer:: lcr(ncrm),kib(ncrm),kob(ncrm),ipo(ncrm)
integer:: lgrx(0:nxp1),isix(0:nx)
integer:: lgry(0:nyp1),isiy(0:ny)
integer:: icrtab(nyp2*nxp2,2)
integer:: ncrpl(2*nlevm),ibpl(2*nlevm)
integer:: ncr,lgrd,lgrt,jump,inc,kibt,kobt
integer:: lev,icr,ibt,ibeg,iend,k,noc
integer:: ix,iy,i,j,icrn,nrem,ioff,ndt
integer,parameter:: short = selected_int_kind(1)
integer(kind=short):: noctab(nyp2*nxp2)
logical:: crox(0:nx),croy(0:ny),free(ncrm)
logical:: keep(nplm)

 !Counter for total number of grid line crossings:
ncr=0
 !Counters for each level:
ncrpl=0

 !Find x grid line crossings first (x = constant):
do ix=1,nx
  xgt=xg(ix)

  do iy=1,nyp1
    lgry(iy)=nint((qoff+qa(iy,ix))*dqi)
  enddo

  do iy=1,ny
    lgrd=lgry(iy+1)-lgry(iy)
    croy(iy)=(lgrd .ne. 0)
    isiy(iy)=sign(1,lgrd)
  enddo

  do iy=1,ny
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

 !Above, kib = grid box into which the contour (containing icr) is going
 !       kob =   "   "  out of "    "     "         "       "    " coming
 !       [kob -> icr -> kib:  icr lies at the boundary between kob & kib]

 !Find y grid line crossings next (y = constant):
do iy=1,ny
  ygt=yg(iy)

  do ix=1,nxp1
    lgrx(ix)=nint((qoff+qa(iy,ix))*dqi)
  enddo

  do ix=1,nx
    lgrd=lgrx(ix+1)-lgrx(ix)
    crox(ix)=(lgrd .ne. 0)
    isix(ix)=sign(1,lgrd)
  enddo

  do ix=1,nx
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

if (ncr .eq. 0) then
  nd=0
  nptd=0
  return
endif

 !Find the first and last level having non-zero number of crossings:
lev=1
do while (ncrpl(lev) .eq. 0) 
  lev=lev+1
enddo
levbeg=lev
do while (ncrpl(lev) .ne. 0) 
  lev=lev+1
enddo
levend=lev-1

levbegall=min(levbegall,levbeg)
levendall=max(levendall,levend)

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

   !Initialise number of crossings per box:
  do i=ibeg,iend
    noctab(kob(ipo(i)))=0
  enddo

  do i=ibeg,iend
    icr=ipo(i)
     !icr is the index of the current crossing at level lev.
    k=kob(icr)
     !accumulate number of crossings in this box:
    noctab(k)=noctab(k)+1
     !assign crossing to box, permitting 2 crossings:
    icrtab(k,noctab(k))=icr
  enddo

  do i=ibeg,iend
    icr=ipo(i)
    k=kib(icr)
    noc=noctab(k)
     !Use last crossing in this box as the next node:
    kob(icr)=icrtab(k,noc)
     !kob(icr) now gives the next point after icr
    noctab(k)=noc-1
     !This will normally be zero, except for boxes with 2 crossings;
     !this allows a second use of this box.
  enddo
enddo

 !Now re-build contours:
j=0
i=0

do icr=1,ncr
  free(icr)=.true.
enddo

do icr=1,ncr
  if (free(icr)) then
     !new contour (j) starts here:
    i=i+1
    j=j+1
    i1d(j)=i
    indd(j)=lcr(icr)
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
    indd(nd)=indd(j)
    i1d(nd)=ibeg
    nptd=nptd+npd(j)
    i2d(nd)=nptd
  endif
enddo
 !Done rebuilding contours.

 !Redistribute points on the contours to reduce point density:
call renode

return
end subroutine grid2con

!=======================================================================

subroutine storecon(iz)
 !Stores the shifted contours xd, yd, i1d, i2d ... 
 !into the arrays xa, ya, i1a, i2a ....  
 !Also updates na & npta and the node & contours layer indices 
 !il1a, il2a, jl1a & jl2a.

integer:: iz,i,j,jn

 !Store starting layer indices:
jl1a(iz)=na+1
il1a(iz)=npta+1

 !Store nodes:
do i=1,nptd
  xa(npta+i)=xd(i)
  ya(npta+i)=yd(i)
enddo

 !Store contour indices:
do j=1,nd
  jn=na+j
  i1a(jn) =npta+i1d(j)
  i2a(jn) =npta+i2d(j)
  npa(jn) =npd(j)
  inda(jn)=indd(j)
  laya(jn)=iz
enddo

 !Augment number of contours and nodes:
na=na+nd
npta=npta+nptd

 !Store ending layer indices:
jl2a(iz)=na
il2a(iz)=npta

return
end subroutine storecon

!=======================================================================

subroutine renode
 !Re-nodes each contour while preserving corner locations.
 !Uses square-root dependence on a weighted sum of nearby 
 !curvature values.

double precision:: ww,sum,xx,yy,acc,eta,fac,p
integer:: iz,i,j,jn,ia,ib,inew,i1t,i2t,im,npseg

 !Define the point following a node i by next(i):
do i=1,nptd-1
  next(i)=i+1
enddo

do j=1,nd
  next(i2d(j))=i1d(j)
enddo

 !Get the updated cubic interpolation coefficients:
call cubic

 !Use the spherical curvature expression (radius of the sphere = ell)
 !to ensure an adequate node density in low curvature regions.
do i=1,nptd
  ww=one/(v(i)+dmsq)
  u(i)=ww*sqrt(elf*v(i)+u(i)**2)
  v(i)=ww*d(i)
enddo
 !NB: elf = 1/ell**2; v(i) = |xx_{i+1}-xx_{i}|**2; d(i)=sqrt{v(i)};
 !         u(i)/d(i) = (kappa_{i}+kappa_{i+1})/2; dmsq = (2*dm)**2

 !Re-assign curvature at a node from weighted average on either side
 !(v above is the weight):
do ib=1,nptd
  i=next(ib)
  q(i)=(u(ib)+u(i))/(v(ib)+v(i))
enddo

 !Re-average to get interval value (effectively, four curvature
 !values go into getting the final interval value, u(i)):
do i=1,nptd
  ia=next(i)
  u(i)=f12*(q(i)+q(ia))
enddo

 !Compute fractional number of nodes to be placed between old
 !nodes i and i+1:
do i=1,nptd
  d(i)=d(i)*min(dmi,densf*sqrt(u(i))+u(i))
enddo
 !NB: dmi = 2/delta; densf = 1/(amu*sqrt{ell})

 !Now begin the redistribution of nodes contour by contour,
 !making sure to preserve corner locations:
nptd=0
do j=1,nd
  inew=1
  i1t=i1d(j)
  i1d(j)=nptd+1
300 u(nptd+inew)=xd(i1t)
    v(nptd+inew)=yd(i1t)
    sum=zero
    i=i1t
310   sum=sum+d(i)
      i=i+1
      last=i .gt. i2d(j)
      if (last) goto 330
      if (corn(i)) goto 320
      goto 310
320 if (sum .lt. small) then
      i1t=i
      goto 300
    else
      i2t=i-1
      goto 340
    endif
330 if (sum .lt. small) then
      inew=inew-1
      goto 390
    else
      i2t=i-1
    endif
340 npseg=nint(sum)+1
     !npseg-1 is the number of nodes to be placed on the contour segment.
    fac=dble(npseg)/sum
    do i=i1t,i2t
      d(i)=fac*d(i)
    enddo
     !Now, the sum of d(i) is equal to npseg.
     !The first node along a contour (segment) is fixed;
     !find the new node positions:
    acc=zero
    i=i1t-1
    do im=nptd+inew+1,nptd+inew+npseg-1
      if (acc .ge. one) goto 370
360     acc=acc+d(i+1)
        i=i+1
        if (acc .lt. one) goto 360
370   acc=acc-one
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
390 npd(j)=inew
  nptd=nptd+inew
enddo

 !Switch arrays around again:
do i=1,nptd
  xx=u(i)-xavg
  xd(i)=xavg+xx-xwid*int(xx*hlxi)
enddo

do i=1,nptd
  yy=v(i)-yavg
  yd(i)=yavg+yy-ywid*int(yy*hlyi)
enddo

 !Reset ending contour indices:
do j=1,nd
  i2d(j)=i1d(j)+npd(j)-1
enddo

return
end subroutine renode

!=======================================================================

subroutine cubic
 !Calculates the interpolation coefficients between the nodes 
 ![xd(i),yd(i)] and [xd(next(i)),yd(next(i))], i = 1, ..., nptd.

 !The interpolation approximately enforces continuity of curvature 
 !(except at corners which have effectively infinite curvature).

double precision:: xx,yy
integer:: i,ia,ib

do i=1,nptd
  ia=next(i)
  xx=xd(ia)-xd(i)
  dx(i)=xx-xwid*int(xx*hlxi)
enddo

do i=1,nptd
  ia=next(i)
  yy=yd(ia)-yd(i)
  dy(i)=yy-ywid*int(yy*hlyi)
enddo

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
   !Set curvature to zero at corners:
    b(i)=zero
  else
    b(i)=(dx(i)*q(i)-c(i)*dy(i))/ &
         sqrt((c(i)*v(i)-dx(i)*u(i))**2+(q(i)*v(i)-dy(i)*u(i))**2+small3)
  endif
enddo

do i=1,nptd
  ia=next(i)
  u(i)=d(i)*(b(ia)+b(i))
  c(i)=d(i)*(b(ia)-b(i))
enddo

 !Calculate the cubic interpolation coefficients:
do i=1,nptd
  a(i)=f16*c(i)-f12*u(i)
  b(i)=f12*(u(i)-c(i))
  c(i)=f13*c(i)
enddo

return
end subroutine cubic

!=======================================================================

subroutine writecon
 !Writes contours to congen.asc and contour levels to levels.asc:

integer:: i,j,lev,levnew
character(len=20):: levname,conname
character(len=4):: pind

 !Save number of vertical grid points and scaled mean depth (so the 
 !total rescaled depth is L_D = c/f) to a file for use by pc3d.py:
open(9,file='3d/vstruct.asc',status='replace')
write(9,'(i4,1x,f12.8)') nz+1,cgw/cof
close(9)

 !Create filenames for output:
write(pind,'(i4.4)') loop
if (iopt .eq. 1) then
  levname='ql'//pind//'levels.asc'
  conname='ql'//pind//'contours.asc'
else if (iopt .eq. 2) then
  levname='d'//pind//'levels.asc'
  conname='d'//pind//'contours.asc'
else if (iopt .eq. 3) then
  levname='g'//pind//'levels.asc'
  conname='g'//pind//'contours.asc'
else if (iopt .eq. 4) then
  levname='r'//pind//'levels.asc'
  conname='r'//pind//'contours.asc'
else if (iopt .eq. 5) then
  levname='w'//pind//'levels.asc'
  conname='w'//pind//'contours.asc'
else if (iopt .eq. 6) then
  levname='pn'//pind//'levels.asc'
  conname='pn'//pind//'contours.asc'
endif

 !Adjust the levels so that - values have negative indices
 !                      and + values have positive indices:
do j=1,na
  lev=inda(j)
  inda(j)=lev-nlevm+(lev-1)/nlevm-1
enddo

open(12,file='3d/'//levname,status='replace')

write(*,*) 
write(*,*) ' Level      Value '
write(*,*) ' -----   ----------- '
do lev=levbegall,levendall
  levnew=lev-nlevm+(lev-1)/nlevm-1
  write(12,'(2x,i5,3x,f14.9)') levnew,qlev(lev)
  write( *,'(2x,i5,3x,f14.9)') levnew,qlev(lev)
enddo

close(12)

open(11,file='3d/'//conname,status='replace')

write(*,*)
write(*, '(a,i6,a,i7)') '   # contours = ',na,'   # nodes = ',npta

write(11,'(i6,1x,i7,1x,f7.2)') na,npta,zero
do j=1,na
  write(11,'(i5,1x,i7,2(1x,i5),f14.10)') npa(j),i1a(j),inda(j),laya(j),dq
enddo
do i=1,npta
  write(11,'(f12.9,1x,f12.9)') xa(i),ya(i)
enddo

close(11)

return
end subroutine writecon

!================================================================

end program g2c
