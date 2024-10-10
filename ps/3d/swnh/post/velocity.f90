!#########################################################################
!  Computes the vertically-averaged velocity (u2d,v2d) where the vertical
!  average is 1/h*int_0^h (.)dz.

!  Writes (u2d,v2d) to 2d/u.r4 & v.r4.

!  Written 5/8/2019 by D G Dritschel @ New York
!#########################################################################

program velocity

 !Import spectral module:
use spectral

implicit none

 !Various arrays needed below:
double precision:: u(ng,ng,0:nz),v(ng,ng,0:nz),r(ng,ng,0:nz)
double precision:: d(ng,ng,0:nz),g(ng,ng,0:nz),zeta(ng,ng,0:nz)

double precision:: qs(ng,ng,0:nz),ds(ng,ng,0:nz),gs(ng,ng,0:nz)

double precision:: h2d(ng,ng),u2d(ng,ng),v2d(ng,ng)

!Other local variables:
real:: tr4,q2dr4(ng,ng),q3dr4(ng,ng,0:nz)
integer:: loop,iread,iz

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open input data files:
open(31,file='3d/ql.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
open(32,file= '3d/d.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
open(33,file= '3d/g.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
open(44,file= '2d/h.r4',form='unformatted',access='direct', &
                      status='old',recl=nhbytes)

 !Open output files:
open(51,file='2d/u.r4',form='unformatted',access='direct', &
                     status='replace',recl=nhbytes)
open(52,file='2d/v.r4',form='unformatted',access='direct', &
                     status='replace',recl=nhbytes)

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0

  read(31,rec=loop,iostat=iread) tr4,q3dr4
  if (iread .ne. 0) exit 
  zeta=dble(q3dr4)
  call ptospc3d(zeta,qs,0,nz)

  read(32,rec=loop,iostat=iread) tr4,q3dr4
  if (iread .ne. 0) exit 
  d=dble(q3dr4)
  call ptospc3d(d,ds,0,nz)

  read(33,rec=loop,iostat=iread) tr4,q3dr4
  if (iread .ne. 0) exit 
  g=dble(q3dr4)
  call ptospc3d(g,gs,0,nz)

  read(44,rec=loop,iostat=iread) tr4,q2dr4
  if (iread .ne. 0) exit 
  h2d=dble(q2dr4)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

  !Obtain velocity field by inversion:
  call main_invert(qs,ds,gs,r,u,v,zeta)
  !Note: qs, ds & gs are in spectral space while 
  !      r, u, v and zeta are in physical space.

  !Compute vertically-averaged velocity:
  u2d=zero
  v2d=zero
  do iz=0,nz
    u2d=u2d+weight(iz)*(one+r(:,:,iz))*u(:,:,iz)
    v2d=v2d+weight(iz)*(one+r(:,:,iz))*v(:,:,iz)
  enddo
  u2d=u2d/(one+h2d)
  v2d=v2d/(one+h2d)

  !Write field data to the 2d subdirectory:
  write(51,rec=loop) tr4,real(u2d)
  write(52,rec=loop) tr4,real(v2d)
enddo

 !Close files:
close(31)
close(32)
close(33)
close(44)
close(51)
close(52)

write(*,*)
write(*,*) ' Vertically averaged horizontal velocity written to 2d/u.r4 & v.r4'

 !End main program
end program velocity
!=======================================================================
