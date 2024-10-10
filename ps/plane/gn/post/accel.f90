!#########################################################################
!  Computes the magnitude of the non-hydrostatic acceleration and writes
!  aa.r4.  Also writes the components to ax.r4 and ay.r4.

!  Computes also the rms non-hydrostatic acceleration relative to the 
!  rms total acceleration and writes an_rms.asc.

!           Written 30/4/2019 by D G Dritschel @ St Andrews
!#########################################################################

program accelnh

 !Import spectral module:
use spectral

implicit none

 !Physical arrays:
double precision:: uu(ng,ng),vv(ng,ng),hh(ng,ng),htot(ng,ng)
double precision:: dd(ng,ng),zz(ng,ng),gg(ng,ng),aa(ng,ng)

 !Spectral arrays:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng)

 !Other local variables:
double precision:: uio,vio,anrms
real:: tr4,qqr4(ng,ng)
integer:: loop,iread

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Open input data files:
open(31,file='hh.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(32,file='dd.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(33,file='zz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(34,file='gg.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)

 !Open output data files:
open(41,file='aa.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(42,file='ax.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(43,file='ay.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(61,file='an_rms.asc',status='replace')

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit 
  hh=dble(qqr4)

  read(32,rec=loop,iostat=iread) tr4,qqr4
  dd=dble(qqr4)

  read(33,rec=loop,iostat=iread) tr4,qqr4
  zz=dble(qqr4)

  read(34,rec=loop,iostat=iread) tr4,qqr4
  gg=dble(qqr4)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

  !Obtain velocity field from dd & zz:
  aa=dd
  call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
  wke=rlap*wke
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)

  aa=zz
  call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
  wke=rlap*wke
  call xderiv(ng,ng,hrkx,wke,wkc)
  call yderiv(ng,ng,hrky,wke,wkd)

  wka=wka-wkd
  wkb=wkb+wkc
  call spctop(ng,ng,wka,uu,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkb,vv,xfactors,yfactors,xtrig,ytrig)

  !Add mean flow (uio,vio):
  uio=-sum(hh*uu)*dsumi
  vio=-sum(hh*vv)*dsumi
  uu=uu+uio
  vv=vv+vio

  !Compute J(u,v):
  call jacob(uu,vv,aa)

  !Complete definition of gamma-tilde:
  aa=gg+two*(aa-dd**2)

  !Filter:
  call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
  wke=filt*wke
  call spctop(ng,ng,wke,aa,xfactors,yfactors,xtrig,ytrig)

  !Define 1+hh:
  htot=one+hh

  !Multiply by htot and filter:
  aa=aa*htot
  call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
  wke=filt*wke
  call spctop(ng,ng,wke,aa,xfactors,yfactors,xtrig,ytrig)

  !Square htot and filter:
  wkb=htot**2
  call ptospc(ng,ng,wkb,wke,xfactors,yfactors,xtrig,ytrig)
  wke=filt*wke
  call spctop(ng,ng,wke,wkb,xfactors,yfactors,xtrig,ytrig)

  !Multiply htot^2 and htot*gamma-tilde, filter and obtain gradient:
  aa=aa*wkb
  call ptospc(ng,ng,aa,wke,xfactors,yfactors,xtrig,ytrig)
  wke=filt*wke
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)
  call spctop(ng,ng,wka,uu,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkb,vv,xfactors,yfactors,xtrig,ytrig)

  !Normalise to find non-hydrostatic acceleration (use uu & vv):
  uu=hbsq3*uu/htot
  vv=hbsq3*vv/htot

  !Next obtain hydrostatic acceleration (use dd & zz):
  call ptospc(ng,ng,hh,wke,xfactors,yfactors,xtrig,ytrig)
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkb)
  call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkb,zz,xfactors,yfactors,xtrig,ytrig)
  dd=-csq*dd
  zz=-csq*zz

  !Form relative rms acceleration and write data:
  anrms=sqrt(sum(uu**2+vv**2)/sum((dd+uu)**2+(zz+vv)**2))
  write(61,'(1x,f12.5,1x,e14.7)') tr4,anrms

  !Obtain non-hydrostatic acceleration magnitude and write data:
  aa=sqrt(uu**2+vv**2)
  write(41,rec=loop) tr4,real(aa)

  !Write components of acceleration also:
  write(42,rec=loop) tr4,real(uu)
  write(43,rec=loop) tr4,real(vv)
enddo

 !Close files:
close(31)
close(32)
close(33)
close(34)
close(41)
close(42)
close(43)
close(61)

write(*,*)
write(*,*) ' Non-hydrostatic acceleration magnitude written to aa.r4'
write(*,*) ' while the components are written to ax.r4 and ay.r4'

 !End main program
end program accelnh
!=======================================================================
