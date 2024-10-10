!#########################################################################
!  Computes the imbalanced spectra of height, vorticity, divergence and 
!  acceleration divergence and writes the data to ispectra.asc.
!  The output can be viewed using ispec_view.

!  Also saves the balanced spectra, which can be viewed with bspec_view.

!  *** Must run dgbal first ***

!           Written 6/4/2018 by D G Dritschel @ St Andrews
!#########################################################################

program ispectra

 !Import contants, parameters and common arrays:
use constants
use spectral

implicit none

 !Physical arrays:
double precision::  hh(ng,ng), zz(ng,ng), dd(ng,ng), gg(ng,ng)
double precision:: bhh(ng,ng),bzz(ng,ng),bdd(ng,ng),bgg(ng,ng)
double precision::  aa(ng,ng)

 !Spectral arrays:
double precision:: wka(ng,ng)

 !Other local variables:
double precision:: hspec(0:ng),zspec(0:ng),dspec(0:ng),gspec(0:ng)
real:: tr4,qqr4(ng,ng)
integer:: loop,iread,k

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

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
open(50,file='bspectra.asc',status='replace')
open(60,file='alt-bspectra.asc',status='replace')
open(51,file='ispectra.asc',status='replace')
open(61,file='alt-ispectra.asc',status='replace')

!----------------------------------------------------------------------
 !Read data and process:
loop=0
do
  loop=loop+1
  iread=0
  read(32,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit 
  dd=dble(qqr4)
  read(33,rec=loop,iostat=iread) tr4,qqr4
  gg=dble(qqr4)
  read(34,rec=loop,iostat=iread) tr4,qqr4
  hh=dble(qqr4)
  read(35,rec=loop,iostat=iread) tr4,qqr4
  zz=dble(qqr4)

  read(42,rec=loop,iostat=iread) tr4,qqr4
  bdd=dble(qqr4)
  read(43,rec=loop,iostat=iread) tr4,qqr4
  bgg=dble(qqr4)
  read(44,rec=loop,iostat=iread) tr4,qqr4
  bhh=dble(qqr4)
  read(45,rec=loop,iostat=iread) tr4,qqr4
  bzz=dble(qqr4)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

  !Compute 1d spectra of balanced fields:
  aa=bhh
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,hspec)
  aa=bdd
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,dspec)
  aa=bgg
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,gspec)
  aa=bzz
  call ptospc(ng,ng,aa,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,zspec)
  hspec=spmf*hspec
  dspec=spmf*dspec
  gspec=spmf*gspec
  zspec=spmf*zspec

  write(50,'(f9.2,1x,i5)') tr4,kmaxred
  write(60,'(f9.2,1x,i5)') tr4,kmaxred
  do k=1,kmaxred
    write(50,'(4(1x,f12.8))') alk(k),log10(zspec(k)+1.d-32), &
                                     log10(dspec(k)+1.d-32), &
                                     log10(gspec(k)+1.d-32)
    write(60,'(4(1x,f12.8))') alk(k),log10(hspec(k)+1.d-32), &
                                     log10(dspec(k)+1.d-32), &
                                     log10(gspec(k)+1.d-32)
  enddo

  !Obtain imbalanced fields:
  hh=hh-bhh
  dd=dd-bdd
  gg=gg-bgg
  zz=zz-bzz

  !Compute 1d spectra:
  call ptospc(ng,ng,hh,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,hspec)
  call ptospc(ng,ng,dd,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,dspec)
  call ptospc(ng,ng,gg,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,gspec)
  call ptospc(ng,ng,zz,wka,xfactors,yfactors,xtrig,ytrig)
  call spec1d(wka,zspec)
  hspec=spmf*hspec
  dspec=spmf*dspec
  gspec=spmf*gspec
  zspec=spmf*zspec

  write(51,'(f9.2,1x,i5)') tr4,kmaxred
  write(61,'(f9.2,1x,i5)') tr4,kmaxred
  do k=1,kmaxred
    write(51,'(4(1x,f12.8))') alk(k),log10(zspec(k)+1.d-32), &
                                     log10(dspec(k)+1.d-32), &
                                     log10(gspec(k)+1.d-32)
    write(61,'(4(1x,f12.8))') alk(k),log10(hspec(k)+1.d-32), &
                                     log10(dspec(k)+1.d-32), &
                                     log10(gspec(k)+1.d-32)
  enddo
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
close(50)
close(60)
close(51)
close(61)

write(*,*)
write(*,*) ' Spectra of imbalanced zeta, delta and gamma now in ispectra.asc.'
write(*,*) ' Spectra of *balanced* zeta, delta and gamma now in bspectra.asc.'

 !End main program
end program
!=======================================================================
