program r4toc2
!  ---------------------------------------------------------------------
!  |   Creates character maps for producing images of various fields   |
!  |   either at a fixed time (using data in the grid subdirectory or  |
!  |   in the fine subdirectory (only for PV)), or over a selected     |
!  |   range of times (using data in the main job directory).          |
!  |                                                                   |
!  |   Images may be created either on a longitude-latitude grid or    |
!  |   in orthographic projection rotated about a specified axis.      |
!  |                                                                   |
!  |   All input data files are assumed to be unfortmatted and direct  |
!  |   access.  Furthermore, they should have a ".r4" extension to     |
!  |   to indicate single-precision real.                              |
!  |                                                                   |
!  |   Output files (ending in .c2) are also unformatted and direct    |
!  |   access.  Here, they contain a pair of characters (hence "c2")   |
!  |   per pixel to enable practically continuous true colour;         |
!  |   images can be converted to various formats using the script     |
!  |   c2image, and they can be made into movies using c2movie.        |
!  |                                                                   |
!  |       Completed 27 June 2013 by D G Dritschel @ St Andrews        |
!  ---------------------------------------------------------------------

use constants

implicit none

integer,parameter:: nbpc=1
! npbc: number of bytes used to represent a single character;
!       normally 1 but can be 4 (as when using Intel Fortran)

real,parameter:: cfrac=1024./(256.*256)
 !This parameter is used to access the fraction of the colourmap 
 !extending from cfrac to 1-cfrac (to avoid oversaturation)

real:: crange,coff,cfac,qqmin,qqmax
integer:: fopt,iopt,jopt,sopt,kbeg,kend,kint,period
integer(kind=dbleint):: nx,ny,nbread,nbprec

 !Input variable names (in data files read):
character(len=2),dimension(7),parameter:: ivar= &
               & ['hh','uu','vv','zz','qq','qq','dd']
 !Output variable names (in data files created):
character(len=2),dimension(7),parameter:: ovar= &
               & ['hh','uu','vv','zz','qq','qa','dd']
character(len=2),dimension(2),parameter:: gformat=['ll','op']
character(len=3):: suffix

!---------------------------------------------------------
 !Quantities required for converting data into characters:
crange=float(256**2)
coff=cfrac*crange

!------------------------------------------------------
 !Select time evolution or fixed time:
write(*,*) ' Image sequence type?'
write(*,*)
write(*,*) '    (1) time evolution (all fields)'
write(*,*) '    (2) fixed time (all fields)'
write(*,*) '    (3) ultra-fine resolution PV or PV anomaly'
read(*,*) sopt

if (sopt .lt. 1 .or. sopt .gt. 3) then
  write(*,*) ' Option not recognised.  *** Stopping ***'
  stop
endif

!------------------------------------------------------
 !Select image format:
write(*,*) ' Format?'
write(*,*)
write(*,*) '    (1) longitude-latitude grid'
write(*,*) '    (2) orthographic perspective'
read(*,*) fopt

if (fopt .lt. 1 .or. fopt .gt. 2) then
  write(*,*) ' Option not recognised.  *** Stopping ***'
  stop
endif

!------------------------------------------------------
 !Select field to image:
write(*,*) ' Field to display?'
write(*,*)
if (sopt .lt. 3) then
   !Data are read from hh, uu, vv, qq or dd.r4 either in the current
   !directory (sopt = 1) or in the grid subdirectory (sopt = 2)
  nx=nt
  ny=ng
   !Here the data dimensions are nt longitudes and ng latitudes
  write(*,*) '    (1)       depth anomaly, h'
  write(*,*) '    (2)      zonal velocity, u'
  write(*,*) '    (3) meridional velocity, v'
  write(*,*) '    (4)  relative vorticity, zeta'
  write(*,*) '    (5)             full PV, q'
  write(*,*) '    (6)          PV anomaly, q - f'
  write(*,*) '    (7)          divergence, delta'
  read(*,*) iopt

  if (iopt .lt. 1 .or. iopt .gt. 7) then
    write(*,*) ' Option not recognised.  *** Stopping ***'
    stop
  endif

else
   !The PV field is read from fine/pnnn where nnn is the selected period
  nx=ntu
  ny=ngu
   !Here the data dimensions are ntu longitudes and ngu latitudes
  write(*,*) '    (1)    full PV, q'
  write(*,*) '    (2) PV anomaly, q - f'
  read(*,*) iopt

  if (iopt .lt. 1 .or. iopt .gt. 2) then
    write(*,*) ' Option not recognised.  *** Stopping ***'
    stop
  endif
   !For use in the subroutines below:
  iopt=iopt+4
endif

!------------------------------------------------------
 !Define input and output record lengths (in bytes):
nbread=4*(nx*ny+1)
if (fopt .eq. 1) then
  nbprec=2*nx*ny/nbpc
else
  nbprec=2*ny*ny/nbpc
endif

!------------------------------------------------------
 !Select method of imaging:
write(*,*)
write(*,*) ' Let "f" stand for the field shown.'
write(*,*)

write(*,*) ' Choose on of the following options:'
write(*,*) ' (1) image between -|f|_max and |f|_max,'
write(*,*) ' (2) image between   f_min  and  f_max,  or'
write(*,*) ' (3) specify min and max value?'
read(*,*) jopt

if (jopt .eq. 3) then
  write(*,'(a,i1,a)') '  Min & max f to image? '
  read(*,*) qqmin,qqmax
   !For creating colours in subroutines below:
  cfac=(crange-2.*coff)/(qqmax-qqmin)
endif

!-----------------------------------------------------------------------
 !Determine input file to image, open it and the associated output file:
if (sopt .eq. 3) then
   !This is for imaging the ultra-fine grid PV field:
  write(*,*) ' Time period (choose 0 for the initial one)?'
  read(*,*) period
  write(suffix(1:3),'(i3.3)') period
  kbeg=1
  kend=1
  kint=1

   !Open selected input data file from the fine subdirectory:
  open(11,file='fine/qq'//suffix//'.r4',form='unformatted', &
      & access='direct',status='old',recl=nbread)

   !Open output data file:
  open(22,file='fine/'//ovar(iopt)//gformat(fopt)//suffix//'.c2', &
        & form='unformatted',access='direct',status='unknown',recl=nbprec)
else
   !This is for imaging a field on the regular inversion grid:
  if (sopt .eq. 1) then
     !Image a time sequence from data in the main job directory:
    write(*,*) & 
     & ' Enter the beginning & ending frames, and the interval (e.g. 1):'
    read(*,*) kbeg,kend,kint

     !Open input data file:
    open(11,file=ivar(iopt)//'.r4',form='unformatted', &
        & access='direct',status='old',recl=nbread)

     !Open output data file:
    open(22,file=ovar(iopt)//gformat(fopt)//'.c2',form='unformatted', & 
        & access='direct',status='unknown',recl=nbprec)
  else
     !Image a fixed period from data in the grid subdirectory:
    write(*,*) ' Time period (choose 0 for the initial one)?'
    read(*,*) period
    write(suffix(1:3),'(i3.3)') period
    kbeg=period+1
    kend=kbeg
    kint=1

     !Open input data file:
    open(11,file='grid/'//ivar(iopt)//'.r4',form='unformatted', & 
        & access='direct',status='old',recl=nbread)

     !Open output data file:
    open(22,file='grid/'//ovar(iopt)//gformat(fopt)//suffix//'.c2', &
          & form='unformatted',access='direct',status='unknown',recl=nbprec)
  endif
endif

!-------------------------------------------------------------
 !Convert data to character maps:
if (fopt .eq. 1) then
  call llgrid(ny,nx,kbeg,kend,kint,iopt)
else
  call opgrid(ny,nx,kbeg,kend,kint,iopt,sopt)
endif

close(11)
close(22)

!============================================================================

 ! Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine llgrid(ny,nx,kbeg,kend,kint,iopt)
! Converts data to a character map on a longitude-latitude grid

implicit none

integer(kind=dbleint):: nx,ny,ii

integer:: kbeg,kend,kint,iopt,i,j,k,loop,ich,ich1,ich2

real:: qq(0:ny+1,nx),cof(ny),var(5)
real:: t,aqqmin,aqqmax,bqqmin,bqqmax

character(len=2):: bqq(ny,nx)

!---------------------------------------------------------------
!Create character image(s):
bqqmin=1.e20
bqqmax=-bqqmin
loop=0
do k=kbeg,kend,kint
  loop=loop+1
   !Read each frame of the data:
  read(11,rec=k) t,qq(1:ny,1:nx)

   !If PV anomaly, subtract f from PV:
  if (iopt .eq. 6) then
    do j=1,ny
      cof(j)=fpole*sin((float(j)-f12)*dl-hpi)
    enddo
    do i=1,nx
      do j=1,ny
        qq(j,i)=qq(j,i)-cof(j)
      enddo
    enddo
  endif

   !Work out min & max values for information:
  aqqmin=qq(1,1)
  aqqmax=qq(1,1)
  do i=1,nx
    do j=1,ny
      aqqmin=min(aqqmin,qq(j,i))
      aqqmax=max(aqqmax,qq(j,i))
    enddo
  enddo
  write(*,'(3(a,f12.6))') ' t = ',t,' min value = ',aqqmin, & 
                                  & ' max value = ',aqqmax

   !Keep track of min & max over all time:
  bqqmin=min(bqqmin,aqqmin)
  bqqmax=max(bqqmax,aqqmax)

   !Possibly create dynamic colourmap:
  if (jopt .eq. 1) then
    qqmax=max(abs(aqqmin),abs(aqqmax))
    qqmin=-qqmax
    cfac=(crange-2.*coff)/(qqmax-qqmin)
  else if (jopt .eq. 2) then
    qqmax=aqqmax
    qqmin=aqqmin
    cfac=(crange-2.*coff)/(qqmax-qqmin)
  endif

   !Convert each data value to a pair of characters:
  do i=1,nx
    do j=1,ny
      ich=nint(coff+cfac*(min(qqmax,max(qqmin,qq(j,i)))-qqmin))
      ich1=ich/256+1
      ich2=ich-(ich1-1)*256+1
      bqq(j,i)=char(ich1)//char(ich2)
    enddo
  enddo

   !Write character image:
  write(22,rec=loop) bqq

enddo

if (kend .gt. kbeg) then
  write(*,*)
  write(*,*) ' Over all time, the min & max values of the field are'
  write(*,'(f12.6,a,f12.6)') bqqmin,' & ',bqqmax
endif

write(*,*)
write(*,'(2(a,i5))') ' Note, the image has dimensions ',nx,' by ',ny

end subroutine

!=======================================================================

subroutine opgrid(ny,nx,kbeg,kend,kint,iopt,sopt)
! Converts data to a character map in an orthographic projection

implicit none

integer(kind=dbleint):: nx,ny,ii

integer:: kbeg,kend,kint,iopt,sopt,i,j,k,loop,ich,ich1,ich2
integer:: iy,iz,ic,ip1,jp1,nalp

real:: qq(0:ny+1,nx),cof(ny),var(5)
real:: xg(ny,ny),yg(ny),zg(ny)
real:: t,aqqmin,aqqmax,bqqmin,bqqmax
real:: dy,dz,rlatc,clatc,slatc,rlonc,clonc,slonc
real:: xp,yp,zp,xm,xt,yt,zt,ri,rj,aa,bb,cc,dd,qqt
real:: t11,t12,t21,t22,t31,alp,dalp,calp,salp
real:: a11,a12,a13,a21,a22,a23,a31,a32,a33

character(len=2):: bqq(ny,ny),white

logical:: inside(ny,ny)

!------------------------------------------------------------------
 !For creating a white surround of each image:
white=char(0)//char(0)

 !Work out points inside disk of view:
dy=two/float(ny)
dz=dy
do iy=1,ny
  yg(iy)=dy*(float(iy)-f12)-one
  do iz=1,ny
    zg(iz)=dz*(float(iz)-f12)-one
    inside(iz,iy)=(yg(iy)**2+zg(iz)**2 .lt. one)
    if (inside(iz,iy)) xg(iz,iy)=sqrt(one-yg(iy)**2-zg(iz)**2)
  enddo
enddo

!------------------------------------------------------------------
 !Either show a time sequence from a fixed perspective or a fixed 
 !time rotated about a specified axis:
if (sopt .eq. 1) then
   !Show a time sequence from a fixed perspective:

  write(*,*) ' Latitude & longitude of the direction of view (degrees)?'
  read(*,*) rlatc,rlonc
  rlatc=rlatc*pi/180.
  rlonc=rlonc*pi/180.

  clonc=cos(rlonc)
  slonc=sin(rlonc)
  clatc=cos(rlatc)
  slatc=sin(rlatc)

  bqqmin=1.e20
  bqqmax=-bqqmin
  loop=0
  do k=kbeg,kend,kint
    loop=loop+1
     !Read each frame of the data:
    read(11,rec=k) t,qq(1:ny,1:nx)

     !If PV anomaly, subtract f from PV:
    if (iopt .eq. 6) then
      do j=1,ny
        cof(j)=fpole*sin((float(j)-f12)*dl-hpi)
      enddo
      do i=1,nx
        do j=1,ny
          qq(j,i)=qq(j,i)-cof(j)
        enddo
      enddo
    endif

     !Copy latitudes adjacent to poles (j = 1 and ny) with a pi
     !shift in longitude to simplify interpolation below:
    if (iopt .eq. 2 .or. iopt .eq. 3) then
       !Here, we are imaging u or v; flip sign crossing the pole:
      do i=1,ny
        ic=i+ny
        qq(0,i)=-qq(1,ic)
        qq(0,ic)=-qq(1,i)
        qq(ny+1,i)=-qq(ny,ic)
        qq(ny+1,ic)=-qq(ny,i)
      enddo
    else
       !Here, we are imaging a scalar field; preserve sign crossing the pole:
      do i=1,ny
        ic=i+ny
        qq(0,i)=qq(1,ic)
        qq(0,ic)=qq(1,i)
        qq(ny+1,i)=qq(ny,ic)
        qq(ny+1,ic)=qq(ny,i)
      enddo
    endif

     !Work out min & max values for information:
    aqqmin=qq(1,1)
    aqqmax=qq(1,1)
    do i=1,nx
      do j=1,ny
        aqqmin=min(aqqmin,qq(j,i))
        aqqmax=max(aqqmax,qq(j,i))
      enddo
    enddo
    write(*,'(3(a,f12.6))') ' t = ',t,' min value = ',aqqmin, &
                                    & ' max value = ',aqqmax

     !Keep track of min & max over all time:
    bqqmin=min(bqqmin,aqqmin)
    bqqmax=max(bqqmax,aqqmax)

     !Possibly create dynamic colourmap:
    if (jopt .eq. 1) then
      qqmax=max(abs(aqqmin),abs(aqqmax))
      qqmin=-qqmax
      cfac=(crange-2.*coff)/(qqmax-qqmin)
    else if (jopt .eq. 2) then
      qqmax=aqqmax
      qqmin=aqqmin
      cfac=(crange-2.*coff)/(qqmax-qqmin)
    endif

     !Do interpolation in chosen perspective and convert each data value
     !to a pair of characters:
    do iy=1,ny
      yp=yg(iy)
      do iz=1,ny
        if (inside(iz,iy)) then
          xp=xg(iz,iy)
          zp=zg(iz)

          xm=xp*clatc-zp*slatc
          zt=zp*clatc+xp*slatc
          yt=yp*clonc+xm*slonc
          xt=xm*clonc-yp*slonc

           !Find lat & lon then bi-linearly interpolate qq:
          ri=dli*(pi+atan2(yt,xt))
          i=1+int(ri)
          ip1=1+mod(i,nx)
          bb=float(i)-ri
          aa=one-bb

          rj=dli*(hpidl+asin(zt))
          j=int(rj)
          jp1=j+1
          cc=rj-float(j)
          dd=one-cc

          qqt=bb*(dd*qq(j,i)+cc*qq(jp1,i))+aa*(dd*qq(j,ip1)+cc*qq(jp1,ip1))
          ich=nint(coff+cfac*(min(qqmax,max(qqmin,qqt))-qqmin))
          ich1=ich/256+1
          ich2=ich-(ich1-1)*256+1
          bqq(iz,iy)=char(ich1)//char(ich2)
        else
          bqq(iz,iy)=white
        endif
      enddo
    enddo

     !Write character image:
    write(22,rec=loop) bqq
  enddo

else
   !Show a fixed time rotated about a specified axis:
  write(*,*) ' Latitude & longitude of the axis of rotation (degrees)?'
  read(*,*) rlatc,rlonc
  rlatc=rlatc*pi/180.
  rlonc=rlonc*pi/180.

   !Prepare rotation matrix:
  clonc=cos(rlonc)
  slonc=sin(rlonc)
  clatc=cos(rlatc)
  slatc=sin(rlatc)
  t11=clonc*slatc
  t12=-slonc
  a13=clonc*clatc
  t21=slonc*slatc
  t22=clonc
  a23=slonc*clatc
  t31=-clatc
  a33=slatc

  write(*,*) ' Number of frames to display in one full rotation?'
  read(*,*) nalp
  dalp=twopi/float(nalp)

   !Read a single frame of data:
  read(11,rec=kbeg) t,qq(1:ny,1:nx)

   !Copy latitudes adjacent to poles (j = 1 and ny) with a pi
   !shift in longitude to simplify interpolation below:
  if (iopt .eq. 2 .or. iopt .eq. 3) then
     !Here, we are imaging u or v; flip sign crossing the pole:
    do i=1,ny
      ic=i+ny
      qq(0,i)=-qq(1,ic)
      qq(0,ic)=-qq(1,i)
      qq(ny+1,i)=-qq(ny,ic)
      qq(ny+1,ic)=-qq(ny,i)
    enddo
  else
     !Here, we are imaging a scalar field; preserve sign crossing the pole:
    do i=1,ny
      ic=i+ny
      qq(0,i)=qq(1,ic)
      qq(0,ic)=qq(1,i)
      qq(ny+1,i)=qq(ny,ic)
      qq(ny+1,ic)=qq(ny,i)
    enddo
  endif

   !Work out min & max values for information:
  aqqmin=qq(1,1)
  aqqmax=qq(1,1)
  do i=1,nx
    do j=1,ny
      aqqmin=min(aqqmin,qq(j,i))
      aqqmax=max(aqqmax,qq(j,i))
    enddo
  enddo
  write(*,'(2(a,f12.6))') ' min field value = ',aqqmin, &
                        & ' max field value = ',aqqmax

   !Possibly create dynamic colourmap:
  if (jopt .eq. 1) then
    qqmax=max(abs(aqqmin),abs(aqqmax))
    qqmin=-qqmax
    cfac=(crange-2.*coff)/(qqmax-qqmin)
  else if (jopt .eq. 2) then
    qqmax=aqqmax
    qqmin=aqqmin
    cfac=(crange-2.*coff)/(qqmax-qqmin)
  endif

   !Loop over rotation angle alpha
  do loop=1,nalp
    alp=dalp*float(loop-1)
    calp=cos(alp)
    salp=sin(alp)

    a11=t11*calp-t12*salp
    a12=t12*calp+t11*salp
    a21=t21*calp-t22*salp
    a22=t22*calp+t21*salp
    a31=t31*calp
    a32=t31*salp

     !Do interpolation in chosen perspective and convert each data value
     !to a pair of characters:
    do iy=1,ny
      yp=yg(iy)
      do iz=1,ny
        if (inside(iz,iy)) then
          xp=xg(iz,iy)
          zp=zg(iz)

          xt=a11*xp+a12*yp+a13*zp
          yt=a21*xp+a22*yp+a23*zp
          zt=a31*xp+a32*yp+a33*zp

           !Find lat & lon then bi-linearly interpolate qq:
          ri=dli*(pi+atan2(yt,xt))
          i=1+int(ri)
          ip1=1+mod(i,nx)
          bb=float(i)-ri
          aa=one-bb

          rj=dli*(hpidl+asin(zt))
          j=int(rj)
          jp1=j+1
          cc=rj-float(j)
          dd=one-cc

          qqt=bb*(dd*qq(j,i)+cc*qq(jp1,i))+aa*(dd*qq(j,ip1)+cc*qq(jp1,ip1))
          ich=nint(coff+cfac*(min(qqmax,max(qqmin,qqt))-qqmin))
          ich1=ich/256+1
          ich2=ich-(ich1-1)*256+1
          bqq(iz,iy)=char(ich1)//char(ich2)
        else
          bqq(iz,iy)=white
        endif
      enddo
    enddo

     !Write character image:
    write(22,rec=loop) bqq
  enddo
endif

if (kend .gt. kbeg) then
  write(*,*)
  write(*,*) ' Over all time, the min & max values of the field are'
  write(*,'(f12.6,a,f12.6)') bqqmin,' & ',bqqmax
endif

write(*,*)
write(*,'(2(a,i5))') ' Note, the image has dimensions ',ny,' by ',ny

end subroutine

!=========================================================

end program
