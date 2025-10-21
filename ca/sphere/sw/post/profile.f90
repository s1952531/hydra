program profile
!  -------------------------------------------------------------------------
!  |   Computes the zonal mean and time mean dimensionless height anomaly  |
!  |   h, zonal velocity u, and potential vorticity q.                     |
!  |                                                                       |
!  |   The quantities are output vs latitude in the file havg.asc,         |
!  |   uavg.asc and qavg.asc.                                              |
!  -------------------------------------------------------------------------

use spectral

implicit none

real:: t,qqr4(nLatGridPts,nLongGridPts)
double precision:: uu(nLatGridPts,nLongGridPts),vv(nLatGridPts,nLongGridPts),hh(nLatGridPts,nLongGridPts),zz(nLatGridPts,nLongGridPts)
double precision:: qq(nLatGridPts,nLongGridPts),dd(nLatGridPts,nLongGridPts),gg(nLatGridPts,nLongGridPts)
double precision:: qs(nLatGridPts,nLongGridPts),ds(nLatGridPts,nLongGridPts),gs(nLatGridPts,nLongGridPts)
double precision:: zuu(nLatGridPts),zhh(nLatGridPts),zqq(nLatGridPts)
double precision:: auu(nLatGridPts),ahh(nLatGridPts),aqq(nLatGridPts)
double precision:: t1,t2,tmax,zfac,phi
integer:: i,j,loop,navg,iread

!---------------------------------------------------------------
 !Open all input data files:
open(31,file='evolution/qq.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(32,file='evolution/dd.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(33,file='evolution/gg.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)
open(34,file='evolution/hh.r4',form='unformatted',access='direct', &
                             status='old',recl=nbytes)

!---------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral
zfac=one/dble(nLongGridPts)

!---------------------------------------------------------------
 !Read final time in data:
open(15,file='evolution/ecomp.asc',status='old',access='append')
backspace(15)
read(15,*) tmax
close(15)

write(*,'(a,f9.2)') ' The final time in the simulation is ',tmax
write(*,*) ' Time averaging is done over t1 <= t <= t2.  Enter t1 & t2:'
read(*,*) t1,t2

if (t2 .gt. tmax) t2=tmax

 !Compute and save initial profiles:
read(31,rec=1) t,qqr4
qq=dble(qqr4)
read(34,rec=1) t,qqr4
hh=dble(qqr4)
zhh=zero
zqq=zero
do i=1,nLongGridPts
  zhh=zhh+hh(:,i)
  zqq=zqq+qq(:,i)
enddo
zhh=zhh*zfac
zuu=zero
zqq=zqq*zfac-corFreq

open(21,file='evolution/hini.asc',status='replace')
open(22,file='evolution/uini.asc',status='replace')
open(23,file='evolution/qini.asc',status='replace')

do j=1,nLatGridPts
  phi=(dble(j)-f12)*dl-hpi
  write(21,'(2(1x,f14.10))') zhh(j),phi
  write(22,'(2(1x,f14.10))') zuu(j),phi
  write(23,'(2(1x,f14.10))') zqq(j),phi
enddo
close(21)
close(22)
close(23)

 !Initialise time and zonal averages:
ahh=zero
auu=zero
aqq=zero
navg=0

 !Read data and process:
loop=0
do
  loop=loop+1
  read(31,rec=loop,iostat=iread) t,qqr4
  if (iread .ne. 0) exit 
  if ((t+1.d-6 .gt. t1) .and. (t-1.d-6 .lt. t2)) then
    qq=dble(qqr4)
    read(32,rec=loop) t,qqr4
    dd=dble(qqr4)
    read(33,rec=loop) t,qqr4
    gg=dble(qqr4)
    read(34,rec=loop) t,qqr4
    hh=dble(qqr4)
    write(*,'(a,f9.2)') ' Processing t = ',t

    navg=navg+1

     !FFT some fields:
    do i=1,nLongGridPts
      qs(:,i)=qq(:,i)-corFreq
    enddo
    call forfft(nLatGridPts,nLongGridPts,qs,trig,factors) 
    ds=dd
    call forfft(nLatGridPts,nLongGridPts,ds,trig,factors) 
    gs=gg
    call forfft(nLatGridPts,nLongGridPts,gs,trig,factors) 
  
     !Find velocity field (uu,vv):  
    call main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
     !Note: qs, ds & gs are in semi-spectral space while 
     !      hh, uu, vv, qq and zz are in physical space.
    
     !Initialise averages:
    zhh=zero
    zuu=zero
    zqq=zero

     !Average zonally and accumulate in time:
    do i=1,nLongGridPts
      zhh=zhh+hh(:,i)
      zuu=zuu+uu(:,i)
      zqq=zqq+qq(:,i)
    enddo

    zhh=zhh*zfac
    zuu=zuu*zfac
    zqq=zqq*zfac
    ahh=ahh+zhh
    auu=auu+zuu
    aqq=aqq+zqq
  endif
enddo
close(31)
close(32)
close(33)
close(34)

 !Finish time averaging:
zfac=one/dble(navg)
ahh=ahh*zfac
auu=auu*zfac
aqq=aqq*zfac-corFreq

 !Open and write output files:
open(21,file='evolution/havg.asc',status='replace')
open(22,file='evolution/uavg.asc',status='replace')
open(23,file='evolution/qavg.asc',status='replace')

do j=1,nLatGridPts
  phi=(dble(j)-f12)*dl-hpi
  write(21,'(2(1x,f14.10))') ahh(j),phi
  write(22,'(2(1x,f14.10))') auu(j),phi
  write(23,'(2(1x,f14.10))') aqq(j),phi
enddo
close(21)
close(22)
close(23)

write(*,*)
write(*,*) &
     ' The results are ready in evolution/havg.asc, uavg.asc and qavg.asc'
write(*,*)
write(*,*) ' View using e.g. plotcol evolution/havg.asc evolution/hini.asc &'

end program profile
