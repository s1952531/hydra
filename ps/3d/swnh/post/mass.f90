!########################################################################
!  Computes the vertically-averaged local rate of change of h_tilde
!  due to vertical variations and due to the vertically-averaged flow.
!  Writes these quantities as 2d/m3.r4 and 2d/m2.r4 respectively.

!  Also, the rms values are written to m2_rms.asc and m3_rms.asc.

!          Written 7/8/2019 by D G Dritschel @ New York
!########################################################################

program mass

 !Import spectral module:
use spectral

implicit none

 !Various arrays needed below:
double precision:: u(ng,ng,0:nz),v(ng,ng,0:nz),r(ng,ng,0:nz)
double precision:: d(ng,ng,0:nz),g(ng,ng,0:nz),zeta(ng,ng,0:nz)

double precision:: qs(ng,ng,0:nz),ds(ng,ng,0:nz),gs(ng,ng,0:nz)

double precision:: h2d(ng,ng),u2d(ng,ng),v2d(ng,ng)
double precision:: dm2(ng,ng),dm3(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wke(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng)
double precision:: htoti(ng,ng),wfac(ng,ng)

 !Other local variables:
double precision:: anrms,a3rms,agrms,ahrms
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
open(34,file= '3d/r.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
open(44,file= '2d/h.r4',form='unformatted',access='direct', &
                      status='old',recl=nhbytes)

 !Open output files:
open(51,file='2d/m2.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(52,file='2d/m3.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(81,file='m2_rms.asc',status='replace')
open(82,file='m3_rms.asc',status='replace')

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

  read(34,rec=loop,iostat=iread) tr4,q3dr4
  if (iread .ne. 0) exit 
  r=dble(q3dr4)

  read(44,rec=loop,iostat=iread) tr4,q2dr4
  if (iread .ne. 0) exit 
  h2d=dble(q2dr4)
  htoti=one/(one+h2d)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

  !Obtain velocity field by inversion:
  call main_invert(qs,ds,gs,r,u,v,zeta)
  !Note: qs, ds & gs are in spectral space while 
  !      r, u, v and zeta are in physical space.

  !Compute vertically-averaged velocity:
  u2d=zero
  v2d=zero
  do iz=0,nz
    wfac=weight(iz)*(one+r(:,:,iz))
    u2d=u2d+wfac*u(:,:,iz)
    v2d=v2d+wfac*v(:,:,iz)
  enddo
  u2d=htoti*u2d
  v2d=htoti*v2d

  !Compute mass changes due to vertically varying flow:
  dm3=zero
  do iz=0,nz
    wkp=(r(:,:,iz)-h2d)*(u(:,:,iz)-u2d)
    wkq=(r(:,:,iz)-h2d)*(v(:,:,iz)-v2d)
    call divs(wkp,wkq,wka)
    call spctop(ng,ng,wka,wkp,xfactors,yfactors,xtrig,ytrig)
    dm3=dm3-weight(iz)*wkp
  enddo

  !Compute mass changes due to vertically-averaged flow:
  wkp=(one+h2d)*u2d
  wkq=(one+h2d)*v2d
  call divs(wkp,wkq,wka)
  call spctop(ng,ng,wka,dm2,xfactors,yfactors,xtrig,ytrig)
  dm2=-dm2

  !Write data:
  write(81,'(1x,f12.5,1x,e14.7)') tr4,sqrt(dsumi*sum(dm2**2))
  write(82,'(1x,f12.5,1x,e14.7)') tr4,sqrt(dsumi*sum(dm3**2))

  !Write field data to the 2d subdirectory:
  write(51,rec=loop) tr4,real(dm2)
  write(52,rec=loop) tr4,real(dm3)
enddo

 !Close files:
close(31)
close(32)
close(33)
close(34)
close(44)
close(51)
close(52)
close(81)
close(82)

write(*,*)
write(*,*) ' Fields of 2d & 3d mass changes written to 2d/m2.r4 & 2d/m3.r4'
write(*,*)
write(*,*) ' Rms norms of m2 & m2 are written to m2_rms.asc & m3_rms.asc'

!End main program
end program mass
!=======================================================================
