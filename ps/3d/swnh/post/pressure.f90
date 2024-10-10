!#########################################################################
!  Computes the vertically-averaged non-hydrostatic pressure and writes
!  2d/pn.r4.

!        Written 4 August 2019 by D G Dritschel @ Wildwood Crest
!#########################################################################

program pressnh

 !Import constants and parameters:
use constants

implicit none

 !Various arrays needed below:
double precision:: r(ng,ng,0:nz),pn(ng,ng,0:nz)
double precision:: pn2d(ng,ng)
double precision:: weight(0:nz)

 !Other local variables:
real:: tr4,q3dr4(ng,ng,0:nz)
integer:: loop,iread,iz

!---------------------------------------------------------
 !Initialise weights:
weight(0)=f12
do iz=1,nz-1
  weight(iz)=one
enddo
weight(nz)=f12
weight=weight/dble(nz)

!---------------------------------------------------------------
 !Open input data files:
open(31,file='3d/r.r4' ,form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
open(36,file='3d/pn.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)

 !Open output file:
open(51,file='2d/pn.r4',form='unformatted',access='direct', &
                      status='replace',recl=nhbytes)

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,q3dr4
  if (iread .ne. 0) exit 
  r=dble(q3dr4)

  read(36,rec=loop,iostat=iread) tr4,q3dr4
  if (iread .ne. 0) exit 
  pn=dble(q3dr4)

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

  !Compute vertically-averaged non-hydrostatic pressure:
  pn2d=zero
  do iz=0,nz-1
    !(pn = 0 when iz = nz)
    pn2d=pn2d+weight(iz)*pn(:,:,iz)*(one+r(:,:,iz))
  enddo

  !Write data:
  write(51,rec=loop) tr4,real(pn2d)
enddo

 !Close files:
close(31)
close(36)
close(51)

write(*,*)
write(*,*) ' Vertically-averaged non-hydrostatic pressure written to 2d/pn.r4'

 !End main program
end program pressnh
!=======================================================================
