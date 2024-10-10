program flowprop
!  ------------------------------------------------------------------
!  |   Computes various flow properties (max Froude number, min &   |
!  |   max polar Rossby number, etc, from data in grid/pnnn         |
!  |   (with nnn = 000, 001, etc).                                  |
!  ------------------------------------------------------------------

use constants

implicit double precision(a-h,o-z)

real:: tin,uu(ng,nt),vv(ng,nt),hh(ng,nt),zz(ng,nt)
double precision:: wka(ng,nt),clat(ng)

!--------------------------------------------------------------
 !Define various constants and useful arrays:
cgwi=one/cgw
fpolei=one/fpole
do j=1,ng
  rlat=(dble(j)-f12)*dl-hpi
  clat(j)=cos(rlat)
enddo

rsum=f1112*(clat(1)+clat(ng))
do j=2,ngm1
  rsum=rsum+clat(j)
enddo
 !Note, rsum = 1/sin(dl/2) - sin(dl/2)/6 exactly.
rsumi=one/rsum
dsumi=rsumi/dble(nt)

!---------------------------------------------------------------
write(*,*) '  Enter the time frame to process (0, 1, ...):'
read(*,*) kfr

t=tsim*dble(kfr)
write(*,'(a,f9.2)') ' Processing t = ',t

irec=kfr+1
 !Open & read gridded data files: 
open(40,file='grid/hh.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
open(41,file='grid/uu.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
open(42,file='grid/vv.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
open(43,file='grid/zz.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
read(40,rec=irec) tin,hh
read(41,rec=irec) tin,uu
read(42,rec=irec) tin,vv
read(43,rec=irec) tin,zz
close(40)
close(41)
close(42)
close(43)

!-----------------------------------------------------------
 !Compute diagnostics:
hhmin=zero
hhmax=zero
hhrms=zero

do i=1,nt
  do j=1,ng
    hhmin=min(hhmin,hh(j,i))
    hhmax=max(hhmax,hh(j,i))
    wka(j,i)=clat(j)*hh(j,i)**2
  enddo
enddo

do i=1,nt
  hhrms=hhrms+f1112*(wka(1,i)+wka(ng,i))
  do j=2,ngm1
    hhrms=hhrms+wka(j,i)
  enddo
enddo
hhrms=sqrt(hhrms*dsumi)

!---------------
romin=zero
romax=zero
rorms=zero

do i=1,nt
  do j=1,ng
    romin=min(romin,zz(j,i))
    romax=max(romax,zz(j,i))
    wka(j,i)=clat(j)*zz(j,i)**2
  enddo
enddo
romin=fpolei*romin
romax=fpolei*romax

do i=1,nt
  rorms=rorms+f1112*(wka(1,i)+wka(ng,i))
  do j=2,ngm1
    rorms=rorms+wka(j,i)
  enddo
enddo
rorms=fpolei*sqrt(rorms*dsumi)

!---------------
frmax=zero
frrms=zero

do i=1,nt
  do j=1,ng
    frsq=(uu(j,i)**2+vv(j,i)**2)/(one+hh(j,i))
    frmax=max(frmax,sqrt(frsq))
    wka(j,i)=clat(j)*frsq
  enddo
enddo
frmax=cgwi*frmax

do i=1,nt
  frrms=frrms+f1112*(wka(1,i)+wka(ng,i))
  do j=2,ngm1
    frrms=frrms+wka(j,i)
  enddo
enddo
frrms=cgwi*sqrt(frrms*dsumi)

!------------------------
 !Write out results:
write(*,*)
write(*,*) ' Min, max and rms depth anomaly h:'
write(*,'(3(3x,f9.5))') hhmin,hhmax,hhrms
write(*,*)
write(*,*) ' Min, max and rms polar Rossby number zeta/f_pole:'
write(*,'(3(3x,f9.5))') romin,romax,rorms

write(*,*)
write(*,*) ' Max and rms Froude number |u|/[c(1+h)^{1/2}]:'
write(*,'(2(3x,f9.5))') frmax,frrms

end program
