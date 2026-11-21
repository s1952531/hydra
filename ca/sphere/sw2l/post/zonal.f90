program zonal
!  -------------------------------------------------------------------------
!  |   Computes the zonal mean zonal velocity u, <u>, height anomaly, <h>, |
!  |   and zonal mean PV, <q>, as a function of latitude from data         |
!  |   generated previously by caps.                                       |
!  |                                                                       |
!  |   View the results using the script zv                                |
!  |                                                                       |
!  |   Written => 12.02.20.21 <= by D G Dritschel @ St Andrews             |
!  |   Extended to 2 layers on 9 Nov 2026 by D G Dritschel @ Kyoto         |
!  -------------------------------------------------------------------------

 ! Import constants and parameters:
use constants
 ! Import spectral module:
use spectral

 ! Declarations:
implicit none

double precision:: qAnom1(numLatGridPts,numLongGridPts),qAnom2(numLatGridPts,numLongGridPts)
double precision:: h1(numLatGridPts,numLongGridPts),h2(numLatGridPts,numLongGridPts)
double precision:: z1(numLatGridPts,numLongGridPts),z2(numLatGridPts,numLongGridPts)
double precision:: u1(numLatGridPts,numLongGridPts),u2(numLatGridPts,numLongGridPts)

double precision:: wka(numLatGridPts,numLongGridPts),wkb(numLatGridPts,numLongGridPts),wkc(numLatGridPts,numLongGridPts)
double precision:: zu1(numLatGridPts),zu2(numLatGridPts),zh1(numLatGridPts),zh2(numLatGridPts),zq1(numLatGridPts),zq2(numLatGridPts)
double precision:: phi(numLatGridPts)
double precision:: zfac,t

real:: qqr4(numLatGridPts,numLongGridPts),tr4
integer:: i,j,loop,iread

!---------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Define latitude and scaling constant for averaging:
zfac=one/dble(numLongGridPts)
do j=1,numLatGridPts
  phi(j)=(dble(j)-f12)*dl-hpi
enddo

!---------------------------------------------------------------
 !Open all input data files:
open(31,file='evolution/h1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(41,file='evolution/h2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

open(32,file='evolution/z1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(42,file='evolution/z2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open all output data files:
open(51,file='evolution/zu1.asc',status='replace')
open(61,file='evolution/zu2.asc',status='replace')
open(52,file='evolution/zh1.asc',status='replace')
open(62,file='evolution/zh2.asc',status='replace')
open(53,file='evolution/zq1.asc',status='replace')
open(63,file='evolution/zq2.asc',status='replace')

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,qqr4
  if (iread .ne. 0) exit
  t=dble(tr4)
  h1=dble(qqr4)
  read(41,rec=loop,iostat=iread) tr4,qqr4
  h2=dble(qqr4)

  read(32,rec=loop,iostat=iread) tr4,qqr4
  z1=dble(qqr4)
  read(42,rec=loop,iostat=iread) tr4,qqr4
  z2=dble(qqr4)

  write(*,'(a,f9.2)') ' Processing t = ',t

  !-----------------------------------------------------------------
   !Invert Laplace's operator on relative vorticity:

   !::::::::::::::::::::::::::::::::::::::::::::::::
   !Layer 1:
  call laplinv(z1,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Switch sign of wkb to define the zonal velocity:
  u1=-wkb

   !Return to physical space:
  call revfft(numLatGridPts,numLongGridPts,u1,trig,factors) 

   !::::::::::::::::::::::::::::::::::::::::::::::::
   !Layer 2:

   !Invert Laplace's operator on relative vorticity:
  call laplinv(z2,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Switch sign of wkb to define the zonal velocity:
  u2=-wkb

   !Return to physical space:
  call revfft(numLatGridPts,numLongGridPts,u2,trig,factors) 

  !-----------------------------------------------------------------
   !Define PV:
  do i=1,numLongGridPts
    qAnom1(:,i)=(cof+z1(:,i))/(one+h1(:,i))
    qAnom2(:,i)=(cof+z2(:,i))/(one+h2(:,i))
  enddo
  
   !Compute zonal averages:
  zu1=zero
  zu2=zero
  zh1=zero
  zh2=zero
  zq1=zero
  zq2=zero

  do j=1,numLatGridPts
    zu1(j)=zu1(j)+zfac*sum(u1(j,:))
    zu2(j)=zu2(j)+zfac*sum(u2(j,:))
    zh1(j)=zh1(j)+zfac*sum(h1(j,:))
    zh2(j)=zh2(j)+zfac*sum(h2(j,:))
    zq1(j)=zq1(j)+zfac*sum(qAnom1(j,:))
    zq2(j)=zq2(j)+zfac*sum(qAnom2(j,:))
  enddo

   !Write diagnostic data for this time:
  write(51,'(f13.6,1x,i5)') t,numLatGridPts
  write(61,'(f13.6,1x,i5)') t,numLatGridPts
  write(52,'(f13.6,1x,i5)') t,numLatGridPts
  write(62,'(f13.6,1x,i5)') t,numLatGridPts
  write(53,'(f13.6,1x,i5)') t,numLatGridPts
  write(63,'(f13.6,1x,i5)') t,numLatGridPts
  do j=1,numLatGridPts
    write(51,'(2(1x,f12.8))') zu1(j),phi(j)
    write(61,'(2(1x,f12.8))') zu2(j),phi(j)
    write(52,'(2(1x,f12.8))') zh1(j),phi(j)
    write(62,'(2(1x,f12.8))') zh2(j),phi(j)
    write(53,'(2(1x,f12.8))') zq1(j),phi(j)
    write(63,'(2(1x,f12.8))') zq2(j),phi(j)
  enddo

enddo

 !Close all files:
close(31)
close(41)
close(32)
close(42)

close(51)
close(61)
close(52)
close(62)
close(53)
close(63)

write(*,*)
write(*,*) ' To view the results type'
write(*,*)
write(*,*) ' zv'
write(*,*)

end program zonal
