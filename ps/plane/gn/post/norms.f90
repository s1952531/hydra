!#########################################################################
!    Computes the rms values of height, vorticity, divergence and 
!    acceleration divergence - for the full, balanced and imbalanced
!    fields - and writes the data to hnorm.asc, znorm.asc, dnorm.asc
!    and gnorm.asc.

!    The output can be viewed using plotcol or compare_norms.py 
!    (when comparing different simulation results).

!                   *** Must run dgbal first ***

!           Written 3/5/2018 by D G Dritschel @ St Andrews
!#########################################################################

program ispectra

 !Import constants:
use constants

implicit none

 !Physical arrays:
double precision::  hh(ng,ng), zz(ng,ng), dd(ng,ng), gg(ng,ng)
double precision:: bhh(ng,ng),bzz(ng,ng),bdd(ng,ng),bgg(ng,ng)

 !Other local variables:
double precision:: rms,rmsb,rmsi
real:: tr4,qqr4(ng,ng)
integer:: loop,iread

!----------------------------------------------------------------------
 !Open input data files:
open(32,file='dd.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(33,file='gg.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(34,file='hh.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(35,file='zz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(42,file='bdd.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(43,file='bgg.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(44,file='bhh.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(45,file='bzz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)

 !Open output data files:
open(51,file='hnorms.asc',status='replace')
open(52,file='znorms.asc',status='replace')
open(53,file='dnorms.asc',status='replace')
open(54,file='gnorms.asc',status='replace')

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  !Read original fields:
  read(32,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit 
  dd=dble(qqr4)
  read(33,rec=loop,iostat=iread) tr4,qqr4
  gg=dble(qqr4)
  read(34,rec=loop,iostat=iread) tr4,qqr4
  hh=dble(qqr4)
  read(35,rec=loop,iostat=iread) tr4,qqr4
  zz=dble(qqr4)

  !Read balanced fields:
  read(42,rec=loop,iostat=iread) tr4,qqr4
  bdd=dble(qqr4)
  read(43,rec=loop,iostat=iread) tr4,qqr4
  bgg=dble(qqr4)
  read(44,rec=loop,iostat=iread) tr4,qqr4
  bhh=dble(qqr4)
  read(45,rec=loop,iostat=iread) tr4,qqr4
  bzz=dble(qqr4)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

  !Compute and write rms h values:
  rms=sqrt(dsumi*sum(hh**2))
  rmsb=sqrt(dsumi*sum(bhh**2))
  rmsi=sqrt(dsumi*sum((hh-bhh)**2))
  write(51,'(f9.2,3(1x,e14.7))') tr4,rms,rmsb,rmsi

  !Compute and write rms zeta values:
  rms=sqrt(dsumi*sum(zz**2))
  rmsb=sqrt(dsumi*sum(bzz**2))
  rmsi=sqrt(dsumi*sum((zz-bzz)**2))
  write(52,'(f9.2,3(1x,e14.7))') tr4,rms,rmsb,rmsi

  !Compute and write rms delta values:
  rms=sqrt(dsumi*sum(dd**2))
  rmsb=sqrt(dsumi*sum(bdd**2))
  rmsi=sqrt(dsumi*sum((dd-bdd)**2))
  write(53,'(f9.2,3(1x,e14.7))') tr4,rms,rmsb,rmsi

  !Compute and write rms gamma values:
  rms=sqrt(dsumi*sum(gg**2))
  rmsb=sqrt(dsumi*sum(bgg**2))
  rmsi=sqrt(dsumi*sum((gg-bgg)**2))
  write(54,'(f9.2,3(1x,e14.7))') tr4,rms,rmsb,rmsi

enddo

 !Close files:
close(32)
close(33)
close(34)
close(35)
close(42)
close(43)
close(44)
close(45)
close(51)
close(52)
close(53)
close(54)

write(*,*)
write(*,*) ' The full, balanced and imbalanced rms field values are now'
write(*,*) ' available in hnorm.asc, znorm.asc, dnorm.asc and gnorm.asc'
write(*,*) ' for h, zeta, delta and gamma, respectively.'

 !End main program
end program
!=======================================================================
