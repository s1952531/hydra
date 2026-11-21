program hov

! Creates hovmuller diagrams of the data created by running zonal.f90  

! Written 4 September 2021 by D G Dritschel @ St Andrews
! Extended to two layers on 9 November 2026 by D G Dritschel @ Kyoto

use spectral

implicit none

double precision:: zuu(numLatGridPts),zhh(numLatGridPts),zek(numLatGridPts),zqq(numLatGridPts),zvq(numLatGridPts),zqg(0:numLatGridPts)
double precision:: t,ombar
integer:: iread, loop, j

!---------------------------------------------------------------------
 !Open all input data files:
open(51,file='evolution/zu1.asc',status='old')
open(61,file='evolution/zu2.asc',status='old')
open(52,file='evolution/zh1.asc',status='old')
open(62,file='evolution/zh2.asc',status='old')
open(53,file='evolution/zq1.asc',status='old')
open(63,file='evolution/zq2.asc',status='old')

open(71,file='hovu1.r8',form='unformatted',access='stream',status='replace')
open(81,file='hovu2.r8',form='unformatted',access='stream',status='replace')
open(72,file='hovh1.r8',form='unformatted',access='stream',status='replace')
open(82,file='hovh2.r8',form='unformatted',access='stream',status='replace')
open(73,file='hovq1.r8',form='unformatted',access='stream',status='replace')
open(83,file='hovq2.r8',form='unformatted',access='stream',status='replace')

write(71) 0.d0
write(81) 0.d0
write(72) 0.d0
write(82) 0.d0
write(73) 0.d0
write(83) 0.d0

loop=0
do
  iread=0
  read(51,*,iostat=iread) t
  if (iread .ne. 0) exit 
  read(61,*) t
  read(52,*) t
  read(62,*) t
  read(53,*) t
  read(63,*) t
  do j=1,numLatGridPts
    read(51,*) zu1(j)
    read(61,*) zu2(j)
    read(52,*) zh1(j)
    read(62,*) zh2(j)
    read(53,*) zq1(j)
    read(63,*) zq2(j)
  enddo

   !Remove mean angular velocity from zuu:
  ombar=(f1112*(zu1(1)+zu1(numLatGridPts))+sum(zu1(2:ngm1)))*rsumi
  zu1=zu1-ombar*cosLat
  ombar=(f1112*(zu2(1)+zu2(numLatGridPts))+sum(zu2(2:ngm1)))*rsumi
  zu2=zu2-ombar*cosLat

  write(71) zu1
  write(81) zu2
  write(72) zh1
  write(82) zh2
  write(73) zq1
  write(83) zq2

  loop=loop+1
   !End of loop over time:
enddo

 !Close files:
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


open(88,file='hov-plotting',status='replace')
write(88,*) ' To plot Hovmuller diagrams of u_bar(t,phi), etc, type:'
write(* ,*) ' To plot Hovmuller diagrams of u_bar(t,phi), etc, type:'
write(* ,*)
write(88,'(a,i5,1x,i4)') ' dataview hovu1.r8 -ndim ',loop,numLatGridPts
write(* ,'(a,i5,1x,i4)') ' dataview hovu1.r8 -ndim ',loop,numLatGridPts
write(88,'(a,i5,1x,i4)') ' dataview hovu2.r8 -ndim ',loop,numLatGridPts
write(* ,'(a,i5,1x,i4)') ' dataview hovu2.r8 -ndim ',loop,numLatGridPts
write(88,'(a,i5,1x,i4)') ' dataview hovh1.r8 -ndim ',loop,numLatGridPts
write(* ,'(a,i5,1x,i4)') ' dataview hovh1.r8 -ndim ',loop,numLatGridPts
write(88,'(a,i5,1x,i4)') ' dataview hovh2.r8 -ndim ',loop,numLatGridPts
write(* ,'(a,i5,1x,i4)') ' dataview hovh2.r8 -ndim ',loop,numLatGridPts
write(88,'(a,i5,1x,i4)') ' dataview hovq1.r8 -ndim ',loop,numLatGridPts
write(* ,'(a,i5,1x,i4)') ' dataview hovq1.r8 -ndim ',loop,numLatGridPts
write(88,'(a,i5,1x,i4)') ' dataview hovq2.r8 -ndim ',loop,numLatGridPts+1
write(* ,'(a,i5,1x,i4)') ' dataview hovq2.r8 -ndim ',loop,numLatGridPts+1
write(88,*)
write(* ,*)
close(88)

 !End main program
end program hov
