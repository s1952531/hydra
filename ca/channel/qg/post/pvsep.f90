program pvsep
!  -------------------------------------------------------------------------
!  |   Divides the PV anomaly, q_a = q - beta*y, into even and odd parts,  |
!  |   writes these for imaging as qe & qo.r4, respectively, and computes  |
!  |   the associated enstrophy, writing this (and the enstrophy in qa)    |
!  |   to enstrophy.asc. Also writes q_a to qa.r4.                         |
!  -------------------------------------------------------------------------

 !Import constants and parameters:
use constants

implicit none

integer,parameter:: nyh=ny/2,nhbytes=4*(nx*(nyh+1)+1)

!---------------------------------------------------------
 !Read data and process:
call diagnose

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine diagnose

implicit none

 !Physical fields:
double precision:: qq(0:ny,0:nxm1),qo(0:nyh,0:nxm1),qe(0:nyh,0:nxm1)

 !Diagnostic quantities:
double precision:: bety(0:ny),ensa,enso,ense
real:: qqr4(0:ny,0:nxm1)
real:: tr4

integer:: loop,iread,ix,iy

!---------------------------------------------------------------
 !Open file containing PV field:
open(31,file='evolution/qq.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output files:
open(22,file='evolution/enstrophy.asc',status='replace')

open(40,file='evolution/qa.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

open(41,file='evolution/qo.r4',form='unformatted', &
        access='direct',status='replace',recl=nhbytes)

open(42,file='evolution/qe.r4',form='unformatted', &
        access='direct',status='replace',recl=nhbytes)

 !Define beta*y:
do iy=0,ny
  bety(iy)=beta*(ymin+gly*dble(iy))
enddo

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

   !Remove beta*y:
  do ix=0,nxm1
    do iy=0,ny
      qq(iy,ix)=qq(iy,ix)-bety(iy)
    enddo
  enddo

   !Define qo & qe:
  do ix=0,nxm1
    do iy=0,nyh
      qo(iy,ix)=f12*(qq(nyh+iy,ix)-qq(nyh-iy,ix))
      qe(iy,ix)=f12*(qq(nyh+iy,ix)+qq(nyh-iy,ix))
    enddo
  enddo

   !Calculate various enstrophies:
  ensa=f12*garea*sum(qq**2)
  enso=garea*sum(qo**2)
  ense=garea*sum(qe**2)

   !Write data:
  write(22,'(1x,f13.5,3(1x,e14.7))') tr4,ensa,enso,ense

  write(40,rec=loop) tr4,real(qq)
  write(41,rec=loop) tr4,real(qo)
  write(42,rec=loop) tr4,real(qe)

enddo

 !Close output files:
close(22)
close(40)
close(41)
close(42)

write(*,*)
write(*,*) ' The results are ready in evolution/qo.r4, evolution/qe.r4,'
write(*,*) ' evolution/qa.r4 and evolution/enstrophy.asc'

return
end subroutine diagnose

 !End main program
end program pvsep
!=======================================================================
