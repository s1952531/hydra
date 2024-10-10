program pvsep
!  -------------------------------------------------------------------------
!  |  Divides the PV anomaly, q, into even and odd parts, writes these     |
!  |  for imaging as qe & qo.r4, respectively, and computes the associated |
!  |  the associated enstrophy, writing this (and the enstrophy in q) to   |
!  |  enstrophy.asc.                                                       |
!  -------------------------------------------------------------------------

 !Import constants and parameters:
use constants

implicit none

integer,parameter:: nyh=ny/2,nyc=nyh+1,nhbytes=4*(nx*nyh+1)

!---------------------------------------------------------
 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

 !Physical fields:
double precision:: qq(ny,nx),qo(nyh,nx),qe(nyh,nx)

 !Diagnostic quantities:
double precision:: ensa,enso,ense
real:: qqr4(ny,nx)
real:: tr4

integer:: loop,iread,iy

!---------------------------------------------------------------
 !Open file containing PV anomaly field:
open(31,file='qq.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output files:
open(22,file='enstrophy.asc',status='replace')

open(41,file='qo.r4',form='unformatted', &
        access='direct',status='replace',recl=nhbytes)

open(42,file='qe.r4',form='unformatted', &
        access='direct',status='replace',recl=nhbytes)

!---------------------------------------------------------------
 !Read data and process:
loop=0
do
   loop=loop+1
   iread=0
   read(31,rec=loop,iostat=iread) tr4,qqr4
   if (iread .ne. 0) exit 
   write(*,'(a,f12.5)') ' Processing t = ',tr4

   !Convert PV to double precision as qq:
   qq=dble(qqr4)

   !Define qo & qe:
   qo(1,:)=zero
   qe(1,:)=qq(nyc,:)
   do iy=1,nyh-2
      qo(iy+1,:)=f12*(qq(nyc+iy,:)-qq(nyc-iy,:))
      qe(iy+1,:)=f12*(qq(nyc+iy,:)+qq(nyc-iy,:))
   enddo
   qo(nyh,:)=zero
   qe(nyh,:)=qq(1,:)
   !nyh = ny/2 and nyc = ny/2 + 1 (i.e. y = 0) above.

   !Calculate various enstrophies:
   ensa=f12*garea*sum(qq**2)
   enso=garea*sum(qo**2)
   ense=garea*sum(qe**2)

   !Write data:
   write(22,'(1x,f13.5,3(1x,e14.7))') tr4,ensa,enso,ense

   write(41,rec=loop) tr4,real(qo)
   write(42,rec=loop) tr4,real(qe)

enddo

 !Close output files:
close(22)
close(41)
close(42)

write(*,*)
write(*,*) ' The results are ready in qo.r4, qe.r4 and enstrophy.asc'
write(*,*)

return
end subroutine diagnose

 !End main program
end program pvsep
!=======================================================================
