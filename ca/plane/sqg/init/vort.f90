program ellipse
! Initialises an elliptical vortex.

use constants

implicit double precision(a-h,o-z)

!integer,parameter:: ngridf=nxf*nyf,nbytesf=4*(ngridf+1)
double precision:: qod0(nyu/2),qod1(nyu/2),qod2(nyu/2)
double precision:: qev0(0:nyu/2),qev1(0:nyu/2),qev2(0:nyu/2)
double precision:: qa(nyu,nxu),qq(ny,nx)

integer :: iopt

glxf=ellx/dble(nxu)
glyf=elly/dble(nyu)

open(11,file='data.dat',status='replace',action='write')

do ix=1,nxu
   do iy=1,nyu
      qa(iy,ix) = 0.d0
   enddo   
enddo

write(*,*) ' Enter option: '
write(*,*) '   - 1: parabolic profile '
write(*,*) '   - 2: nearly uniform with a smooth edge (exp(-r**(2n))) '
write(*,*) '   - 3: explicitly uniform '
write(*,*) '   - 4: edge-intensified (x0+ar^2n)e^(-r^2n)'
read(*,*) iopt

if (iopt == 2) then
   write(*,*) ' Enter n'
   read(*,*) nsq
endif

if (iopt == 4) then
   write(*,*) ' Enter x0,a,n'
   read(*,*) x0,aa,nsq
endif

write(*,*) ' Enter the number of vortices '
read(*,*) nv

write(11,*) 'x0,aa,nn',x0,aa,nsq
write(11,*) 'nv',nv
write(11,*) ' xv,yv,rv,qv'




do iv=1,nv
   write(*,*) ' Enter x,y,r,q for vortex ',iv
   read(*,*) xv,yv,rv,qv
   write(11,*) xv,yv,rv,qv
   do ix=1,nxu
      x = xmin+glxf*dble(ix-1) - xv
      do iy=1,nyu
         y = ymin+glxf*dble(iy-1) - yv
         rsq = (x**2 + y**2)/rv**2
         rl = 1.d0 - rsq
!         OPTION 1: parabolic
         if (iopt == 1) then
            if ( rl > 0.d0 ) qa(iy,ix) = qa(iy,ix) + qv*sqrt(rl)
!         OPTION 2 sharp smooth edge
         else if (iopt == 2) then
            qa(iy,ix) = qa(iy,ix) + qv*exp(-rsq**nsq)
!         OPTION 3 explicitly uniform
         else if (iopt == 3) then
            if ( rl > 0.d0 ) qa(iy,ix) = qa(iy,ix) + qv
         else if (iopt == 4) then
            rrr = rsq**nsq
            qa(iy,ix) = qa(iy,ix) + (x0+aa*rrr)*exp(-rrr)
         endif
     enddo   
   enddo   
enddo

close(11)

!------------------------------------------------------------------------
 !Average the buoyancy field in qa to the coarser grid (ny,nx):
nxh=nxu
nyh=nyu
do while (nxh .gt. nx)
  nxuf=nxh
  nxh=nxh/2
  nyuf=nyh
  nyh=nyh/2
   !Perform nine-point averaging:
  do iy=1,nyh
    miy=2*iy
    qod2(iy)=qa(miy-1,nyuf)
    qev2(iy)=qa(miy,nyuf)
  enddo
  qev2(0)=qa(nyuf,nyuf)
  do ix=1,nxh
    mix=2*ix
    mixm1=mix-1
    do iy=1,nyh
      miy=2*iy
      qod1(iy)=qa(miy-1,mixm1)
      qod0(iy)=qa(miy-1,mix)
      qev1(iy)=qa(miy,mixm1)
      qev0(iy)=qa(miy,mix)
    enddo
    qev1(0)=qev1(nyh)
    qev0(0)=qev0(nyh)
    do iy=1,nyh
      qa(iy,ix)=0.0625d0*(qev0(iy)+qev0(iy-1)+qev2(iy)+qev2(iy-1)) &
              & +0.125d0*(qev1(iy)+qev1(iy-1)+qod0(iy)+qod2(iy)) &
              &   +0.25d0*qod1(iy)
    enddo
    do iy=1,nyh
      qod2(iy)=qod0(iy)
      qev2(iy)=qev0(iy)
    enddo
    qev2(0)=qev0(0)
  enddo
enddo

 !Calculate and remove average qq:
qavg=zero
do ix=1,nx
  do iy=1,ny
    qavg=qavg+qa(iy,ix)
  enddo
enddo
qavg=qavg/dble(nx*ny)

do ix=1,nx
  do iy=1,ny
    qq(iy,ix)=qa(iy,ix)-qavg
  enddo
enddo

open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program
