program ortho
!  ---------------------------------------------------------------------
!  |   Creates an orthographic projection of the data in ff.r4 where   |
!  |   ff is any selected field in the evolution subdirectory.         |
!  |   Writes ffop.r4 which can be imaged using odv.                   |
!  |                                                                   |
!  |   Adapted from bt code on => 12.02.20.21 <= by dgd @ St Andrews   |
!  ---------------------------------------------------------------------

use constants

implicit none

real:: rlatc,rlonc
integer:: ny,nx,nbprec
integer:: iopt,jopt
character(len=8):: vfile1,vfile2

!---------------------------------------------------------------------
write(*,*) ' Which field do you wish to project?'
write(*,*) '   (1)  PV anomaly (q),'
write(*,*) '   (2)  dimensionless height anomaly (h)?'
write(*,*) '   (3)  relative vorticity (zeta),'
write(*,*) '   (4)  velocity divergence (delta), or'
write(*,*) '   (5)  acceleration divergence (gamma), or'
write(*,*) ' Option?'
read(*,*) iopt

write(*,*) ' (1) Full or (2) balanced fields?'
read(*,*) jopt
if (jopt*(3-jopt) <= 0) then
  write(*,*) ' *** Not a valid option; exiting!'
  stop
endif

!write(*,*) ' Width (in pixels) of output image?'
!read(*,*) ny
ny=numLatGridPts !need this to be compatible with the odv viewing script [could use ny=numLongGridPts]
nx=ny

 !Define output record length (in bytes):
nbprec=4*(ny*nx+1)

!---------------------------------------------------------------------
if (iopt == 1) then
  if (jopt == 1) then
     !Open input files for q:
    open(11,file='evolution/qAnom1.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
    open(12,file='evolution/qAnom2.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
     !Open output files:      
    open(21,file='evolution/q1op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    open(22,file='evolution/q2op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    vfile1='q1op'
    vfile2='q2op'
  else
     !Open input files (balanced):
    open(11,file='evolution/bq1.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
    open(12,file='evolution/bq2.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
     !Open output files:      
    open(21,file='evolution/bq1op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    open(22,file='evolution/bq2op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    vfile1='bq1op'
    vfile2='bq2op'
  endif
 
else if (iopt == 2) then
  if (jopt == 1) then
     !Open input files for h:
    open(11,file='evolution/h1.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
    open(12,file='evolution/h2.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
     !Open output files:      
    open(21,file='evolution/h1op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    open(22,file='evolution/h2op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    vfile1='h1op'
    vfile2='h2op'
  else
     !Open input files (balanced):
    open(11,file='evolution/bh1.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
    open(12,file='evolution/bh2.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
     !Open output files:      
    open(21,file='evolution/bh1op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    open(22,file='evolution/bh2op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    vfile1='bh1op'
    vfile2='bh2op'
  endif

else if (iopt == 3) then
  if (jopt == 1) then
     !Open input files for z:
    open(11,file='evolution/z1.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
    open(12,file='evolution/z2.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
     !Open output files:      
    open(21,file='evolution/z1op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    open(22,file='evolution/z2op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    vfile1='z1op'
    vfile2='z2op'
  else
     !Open input files (balanced):
    open(11,file='evolution/bz1.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
    open(12,file='evolution/bz2.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
     !Open output files:      
    open(21,file='evolution/bz1op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    open(22,file='evolution/bz2op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    vfile1='bz1op'
    vfile2='bz2op'
  endif

else if (iopt == 4) then
  if (jopt == 1) then
     !Open input files for d:
    open(11,file='evolution/d1.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
    open(12,file='evolution/d2.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
     !Open output files:      
    open(21,file='evolution/d1op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    open(22,file='evolution/d2op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    vfile1='d1op'
    vfile2='d2op'
  else
     !Open input files (balanced):
    open(11,file='evolution/bd1.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
    open(12,file='evolution/bd2.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
     !Open output files:      
    open(21,file='evolution/bd1op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    open(22,file='evolution/bd2op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    vfile1='bd1op'
    vfile2='bd2op'
  endif

else if (iopt == 5) then
  if (jopt == 1) then
     !Open input files for g:
    open(11,file='evolution/g1.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
    open(12,file='evolution/g2.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
     !Open output files:      
    open(21,file='evolution/g1op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    open(22,file='evolution/g2op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    vfile1='g1op'
    vfile2='g2op'
  else
     !Open input files (balanced):
    open(11,file='evolution/bg1.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
    open(12,file='evolution/bg2.r4', form='unformatted', &
          access='direct',status='old',recl=nbytes)
     !Open output files:      
    open(21,file='evolution/bg1op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    open(22,file='evolution/bg2op.r4',form='unformatted', &
          access='direct',status='replace',recl=nbprec)
    vfile1='bg1op'
    vfile2='bg2op'
  endif
else
  write(*,*) ' *** Not a valid option; exiting!'
  stop
endif

!---------------------------------------------------------------------
write(*,*) ' Latitude & longitude of the direction of view (degrees)?'
read(*,*) rlatc,rlonc

 !Project data orthographically in each layer:
call opgrid(rlatc,rlonc,ny,nx,1)
call opgrid(rlatc,rlonc,ny,nx,2)

 !Close files:
close(11)
close(12)
close(21)
close(22)

write(*,*)
write(*,*) ' *** View the results using'
write(*,*)
write(*,*) ' odv '//vfile1
write(*,*) ' odv '//vfile2
write(*,*)

!============================================================================

 ! Internal subroutine definitions (inherit global variables):

contains 

!============================================================================

subroutine opgrid(rlatc,rlonc,ny,nx,contourIdxToLayerIdx)
! Projects lat-lon data in an r4 file orthographically

implicit none

integer:: ny,nx,contourIdxToLayerIdx

integer:: i,j,loop,iread
integer:: iy,iz,ic,ip1,jp1

real:: qq(0:ngp1,numLongGridPts),qqop(nx,ny)
real:: yg(ny),zg(nx)
real:: qqbg,t,dy,dz,det
real:: rlatc,clatc,slatc,rlonc,clonc,slonc
real:: xp,yp,zp,xm,xt,yt,zt,ri,rj,aa,botTopo,cc,dd,qqt

!------------------------------------------------------------------
!write(*,*) ' Background field level (for points off the surface)?'
!read(*,*) qqbg
qqbg=zero

!------------------------------------------------------------------
 !Define points (yg,zg) in cross-sectional plane of view:
dy=2.0/float(ny)
do iy=1,ny
  yg(iy)=dy*(float(iy)-0.5)-1.0
enddo
dz=dy
do iz=1,nx
  zg(iz)=dz*(float(iz-nx/2)-0.5)
enddo

clatc=cos(rlatc*pi/180.)
slatc=sin(rlatc*pi/180.)

clonc=cos(rlonc*pi/180.)
slonc=sin(rlonc*pi/180.)

!------------------------------------------------------------------
 !Process selected frames:
loop=0
do
  loop=loop+1
   !Read each frame of the data:
  read(10+contourIdxToLayerIdx,rec=loop,iostat=iread) t,qq(1:numLatGridPts,1:numLongGridPts)
  if (iread .ne. 0) exit 

   !Copy latitudes adjacent to poles (j = 1 and numLatGridPts) with a pi
   !shift in longitude to simplify interpolation below:
  do i=1,numLatGridPts
    ic=i+numLatGridPts
    qq(0,i)=qq(1,ic)
    qq(0,ic)=qq(1,i)
    qq(ngp1,i)=qq(numLatGridPts,ic)
    qq(ngp1,ic)=qq(numLatGridPts,i)
  enddo

   !Do interpolation in chosen perspective:
  do iy=1,ny
    yp=yg(iy)
    do iz=1,nx
      zp=zg(iz)
      det=1.-yp**2-zp**2
      if (det .gt. zero) then
        xp=sqrt(det)
        xm=xp*clatc-zp*slatc
        zt=zp*clatc+xp*slatc
        xm=xp*clatc-zp*slatc
        yt=yp*clonc+xm*slonc
        xt=xm*clonc-yp*slonc

         !Find longitude & latitude then bi-linearly interpolate qq:
        ri=dli*(pi+atan2(yt,xt))
        i=1+int(ri)
        ip1=1+mod(i,numLongGridPts)
        botTopo=float(i)-ri
        aa=one-botTopo

        rj=dli*(hpidl+asin(zt))
        j=int(rj)
        jp1=j+1
        cc=rj-float(j)
        dd=one-cc

        qqop(iz,iy)=botTopo*(dd*qq(j,i)  +cc*qq(jp1,i))   &
                   +aa*(dd*qq(j,ip1)+cc*qq(jp1,ip1))
      else
        qqop(iz,iy)=qqbg
      endif
    enddo
  enddo

   !Write character image:
  write(20+contourIdxToLayerIdx,rec=loop) qqop
enddo

end subroutine opgrid

!=========================================================

end program ortho
