!#########################################################################
!  Computes the vertically-averaged ageostrophic acceleration (agx,agy),
!  the non-hydrostatic acceleration (anx,any), and the acceleration due
!  to 3D variations (a3x,a3y).   Writes (agx,agy) to 2d/agx.r4 & agy.r4,
!  (anx,any) to 2d/anx.r4 & any.r4, and (a3x,a3y) to 2d/a3x.r4 & a3y.r4.

!  Computes rms norms of the hydrostatic acceleration (ahx,ahy), as well
!  as of (agx,agy), (anx,any) and (a3x,a3y).  Calling these A_h, A_g, A_n
!  and A_3 respectively, then A_g/A_h is written to ag_rms.asc, A_n/A_h
!  is written to an_rms.asc, and A_3/A_h is written to a3_rms.asc.  The
!  unscaled norms A_g, A_n and A_3 are written to ag-an-a3_rms.asc.

!       Originally written 30/4/2019 by D G Dritschel @ St Andrews
!                  Revised  5/8/2019 by D G Dritschel @ New York
!#########################################################################

program accelnh

 !Import spectral module:
use spectral

implicit none

 !Various arrays needed below:
double precision:: u(ng,ng,0:nz),v(ng,ng,0:nz),r(ng,ng,0:nz)
double precision:: d(ng,ng,0:nz),g(ng,ng,0:nz),zeta(ng,ng,0:nz)
double precision:: pn(ng,ng,0:nz)

double precision:: qs(ng,ng,0:nz),ds(ng,ng,0:nz),gs(ng,ng,0:nz)

double precision:: h2d(ng,ng),u2d(ng,ng),v2d(ng,ng),pn2d(ng,ng)
double precision:: ahx(ng,ng),ahy(ng,ng)
double precision:: anx(ng,ng),any(ng,ng)
double precision:: a2x(ng,ng),a2y(ng,ng)
double precision:: a3x(ng,ng),a3y(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wke(ng,ng),wkp(ng,ng)
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
open(36,file='3d/pn.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
open(44,file= '2d/h.r4',form='unformatted',access='direct', &
                      status='old',recl=nhbytes)

 !Open output files:
open(51,file= '2d/an.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(52,file='2d/anx.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(53,file='2d/any.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(61,file= '2d/a3.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(62,file='2d/a3x.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(63,file='2d/a3y.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(71,file= '2d/ag.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(72,file='2d/agx.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(73,file='2d/agy.r4',form='unformatted',access='direct', &
                       status='replace',recl=nhbytes)
open(81,file='an_rms.asc',status='replace')
open(82,file='a3_rms.asc',status='replace')
open(83,file='ag_rms.asc',status='replace')
open(84,file='ag-an-a3_rms.asc',status='replace')

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

  read(36,rec=loop,iostat=iread) tr4,q3dr4
  if (iread .ne. 0) exit 
  pn=dble(q3dr4)

  read(44,rec=loop,iostat=iread) tr4,q2dr4
  if (iread .ne. 0) exit 
  h2d=dble(q2dr4)
  htoti=one/(one+h2d)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

  !Obtain velocity field by inversion:
  call main_invert(qs,ds,gs,r,u,v,zeta)
  !Note: qs, ds & gs are in spectral space while 
  !      r, u, v and zeta are in physical space.

  !Compute vertically-averaged velocity and non-hydrostatic pressure:
  u2d=zero
  v2d=zero
  pn2d=zero
  do iz=0,nz
    wfac=weight(iz)*(one+r(:,:,iz))
    u2d=u2d+wfac*u(:,:,iz)
    v2d=v2d+wfac*v(:,:,iz)
    pn2d=pn2d+wfac*pn(:,:,iz)
  enddo
  u2d=htoti*u2d
  v2d=htoti*v2d
  !Note: pn is not divided by 1+h2d; only the acceleration is below

  !Compute acceleration due to vertically varying velocity:
  a3x=zero
  a3y=zero
  do iz=0,nz
    wfac=weight(iz)*(one+r(:,:,iz))
    !Accumulate the vertical average of u*du/dx+v*du/dy in a3x:
    wkp=u(:,:,iz)
    call ptospc(ng,ng,wkp,wke,xfactors,yfactors,xtrig,ytrig)
    call xderiv(ng,ng,hrkx,wke,wka)
    call yderiv(ng,ng,hrky,wke,wkb)
    call spctop(ng,ng,wka,anx,xfactors,yfactors,xtrig,ytrig) !anx = du/dx
    call spctop(ng,ng,wkb,any,xfactors,yfactors,xtrig,ytrig) !any = du/dy
    a3x=a3x+wfac*(u(:,:,iz)*anx+v(:,:,iz)*any)
    !Accumulate the vertical average of u*dv/dx+v*dv/dy in a3y:
    wkp=v(:,:,iz)
    call ptospc(ng,ng,wkp,wke,xfactors,yfactors,xtrig,ytrig)
    call xderiv(ng,ng,hrkx,wke,wka)
    call yderiv(ng,ng,hrky,wke,wkb)
    call spctop(ng,ng,wka,anx,xfactors,yfactors,xtrig,ytrig) !anx = dv/dx
    call spctop(ng,ng,wkb,any,xfactors,yfactors,xtrig,ytrig) !any = dv/dy
    a3y=a3y+wfac*(u(:,:,iz)*anx+v(:,:,iz)*any)
  enddo
  a3x=htoti*a3x
  a3y=htoti*a3y

  !Compute acceleration due to vertically-averaged velocity:
  wkp=u2d
  call ptospc(ng,ng,wkp,wke,xfactors,yfactors,xtrig,ytrig)
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)
  call spctop(ng,ng,wka,anx,xfactors,yfactors,xtrig,ytrig) !anx = du_bar/dx
  call spctop(ng,ng,wkb,any,xfactors,yfactors,xtrig,ytrig) !any = du_bar/dy
  a2x=u2d*anx+v2d*any
  wkp=v2d
  call ptospc(ng,ng,wkp,wke,xfactors,yfactors,xtrig,ytrig)
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)
  call spctop(ng,ng,wka,anx,xfactors,yfactors,xtrig,ytrig) !anx = dv_bar/dx
  call spctop(ng,ng,wkb,any,xfactors,yfactors,xtrig,ytrig) !any = dv_bar/dy
  a2y=u2d*anx+v2d*any

  !Difference in the above accelerations is the eddy term:
  a3x=a2x-a3x
  a3y=a2y-a3y

  !Compute non-hydrostatic acceleration, -(1+h_tilde)^{-1}*grad(pn2d):
  call ptospc(ng,ng,pn2d,wke,xfactors,yfactors,xtrig,ytrig)
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)
  call spctop(ng,ng,wka,anx,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkb,any,xfactors,yfactors,xtrig,ytrig)
  !Normalise to complete definition of vertically-averaged quantities:
  anx=-htoti*anx
  any=-htoti*any

  !Compute hydrostatic acceleration, -c^2*grad(h_tilde):
  call ptospc(ng,ng,h2d,wke,xfactors,yfactors,xtrig,ytrig)
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)
  call spctop(ng,ng,wka,ahx,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkb,ahy,xfactors,yfactors,xtrig,ytrig)
  ahx=-csq*ahx
  ahy=-csq*ahy

  !Compute various rms norms:
  ahrms=sum(ahx**2+ahy**2)
  anrms=sum(anx**2+any**2)
  a3rms=sum(a3x**2+a3y**2)

  !Store vertically-averaged ageostrophic acceleration (reuse ahx & ahy):
  ahx=ahx+cof*v2d
  ahy=ahy-cof*u2d
  agrms=sum(ahx**2+ahy**2)
  !Write data:
  write(81,'(1x,f12.5,1x,e14.7)') tr4,sqrt(anrms/ahrms)
  write(82,'(1x,f12.5,1x,e14.7)') tr4,sqrt(a3rms/ahrms)
  write(83,'(1x,f12.5,1x,e14.7)') tr4,sqrt(agrms/ahrms)
  write(84,'(1x,f12.5,3(1x,e14.7))') tr4,sqrt(dsumi*agrms), &
                       sqrt(dsumi*anrms),sqrt(dsumi*a3rms)

  !Write field data to the 2d subdirectory:
  wkb=sqrt(anx**2+any**2)
  write(51,rec=loop) tr4,real(wkb)
  write(52,rec=loop) tr4,real(anx)
  write(53,rec=loop) tr4,real(any)

  wkb=sqrt(a3x**2+a3y**2)
  write(61,rec=loop) tr4,real(wkb)
  write(62,rec=loop) tr4,real(a3x)
  write(63,rec=loop) tr4,real(a3y)

  wkb=sqrt(ahx**2+ahy**2)
  write(71,rec=loop) tr4,real(wkb)
  write(72,rec=loop) tr4,real(ahx)
  write(73,rec=loop) tr4,real(ahy)
enddo

 !Close files:
close(31)
close(32)
close(33)
close(34)
close(36)
close(44)
close(51)
close(52)
close(53)
close(61)
close(62)
close(63)
close(71)
close(72)
close(73)
close(81)
close(82)
close(83)
close(84)

write(*,*)
write(*,*) ' Ageostrophic acceleration and magnitude written to'
write(*,*) ' 2d/agx.r4, 2d/agy.r4 and 2d/ag.r4'
write(*,*)
write(*,*) ' Non-hydrostatic acceleration and magnitude written to'
write(*,*) ' 2d/anx.r4, 2d/any.r4 and 2d/an.r4'
write(*,*)
write(*,*) ' 3D variable acceleration and magnitude written to'
write(*,*) ' 2d/a3x.r4, 2d/a3y.r4 and 2d/a3.r4'
write(*,*)
write(*,*) ' Rms norms (relative to hydrostatic acceleration) are written to'
write(*,*) ' ag_rms.asc, an_rms.asc and a3_rms.asc'
write(*,*)
write(*,*) ' Raw rms norms are written to'
write(*,*) ' ag-an-a3_rms.asc'

 !End main program
end program accelnh
!=======================================================================
