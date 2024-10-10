program energy
!  -------------------------------------------------------------------------
!  |   Computes the various components making up the total energy from     |
!  |   data in zz.r4, dd.r4 and hh.r4.                                     |
!  -------------------------------------------------------------------------

 !Import modules:
use common
use spectral

implicit none

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

 !Physical fields:
double precision:: aa(ng,ng),dd(ng,ng),htot(ng,ng)

 !Spectral arrays:
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng)

 !Other local variables:
real:: tr4,qqr4(ng,ng)

 !Diagnostic quantities:
double precision:: ekin,epot,ediv,etot,uio,vio,eu
integer:: loop,iread

!---------------------------------------------------------------
 !Open files containing height anomaly and divergence fields:
open(31,file='hh.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(32,file='dd.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(35,file='zz.r4',form='unformatted', & 
      access='direct', status='old',recl=nbytes)

 !Open output file:
open(22,file='e.asc',status='replace')

!---------------------------------------------------------------
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
  read(35,rec=loop,iostat=iread) tr4,qqr4
  zz=dble(qqr4)

  write(*,'(a,f12.5)') ' Processing t = ',tr4

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

  htot=one+hh
  ekin=f12*garea*sum(htot*(uu**2+vv**2))
  epot=f12*garea*csq*sum(hh**2)
  dd=htot*dd
  ediv=f12*garea*hbsq3*sum(htot*dd**2)
  etot=ekin+epot+ediv

  write(22,'(f9.2,5(1x,f16.9))') tr4,ediv,ekin,ekin+ediv,epot,etot
enddo

 !Close files:
close(22)
close(31)
close(32)
close(35)

write(*,*)
write(*,*) ' t vs E_div, E_kin, E_kin+E_div, E_pot & E_tot are ready in e.asc.'

return
end subroutine

 !End main program
end program
!=======================================================================
