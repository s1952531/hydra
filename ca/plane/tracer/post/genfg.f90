program genfg

! Program for post-processing the contour and residual files created by a run with 
! the casl suite of f90 code.

use contours
use congen

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Declarations:
real:: tr4,qdr4(ny,nx)
character(len=3):: pind
double precision:: qq(ny,nx)

!-----------------------------------------------------------------
write(*,*) ' Preparing...'
call init_contours

write(*,*) ' What is the numbered suffix of the file?'
read(*,*) iind
write(pind(1:3),'(i3.3)') iind

 !Open residual needed to build ultra-fine-grid vorticity with congen:
open(40,file='cont/qqresi.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
read(40,rec=iind) tr4,qdr4
close(40) 
do ix=1,nx
   do iy=1,ny
    qq(iy,ix)=dble(qdr4(iy,ix))
  enddo
enddo

 !Open buoyancy contours:
open(40,file='cont/qqsynopsis.asc',status='old')
do i=1,iind
  read(40,*) nq,nptq,tr8,qjump
enddo
close(40)

open(40,file='cont/qqindex'//pind,form='unformatted',status='old',access='stream')
read(40) npq(1:nq),i1q(1:nq),indq(1:nq)
close(40)

open(40,file='cont/qqnodes'//pind,form='unformatted',status='old',access='stream')
read(40) xq(1:nptq),yq(1:nptq)
close(40)

 !Reconstruct nextz array:
do j=1,nq
  ibeg=i1q(j)
  iend=ibeg+npq(j)-1
  do i=ibeg,iend-1
    nextq(i)=i+1
  enddo
  nextq(iend)=ibeg
enddo 

call ugridsave(qq,xq,yq,qjump,nextq,indq,npq,i1q,nq,nptq,tr8)


contains 

!-----------------------------------------------------------------
subroutine ugridsave(qq,xq,yq,dq,nextq,indq,npq,i1q,nq,nptq,t)

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
real:: qar4(nyu,nxu)
double precision:: qq(ny,nx)
double precision:: xq(npm),yq(npm)
integer:: nextq(npm),indq(nm),npq(nm),i1q(nm)

call con2ugrid


!Bi-linear interpolate the residual qq to the fine grid and add to qa:
do ix=1,nxu
  ixf=ixfw(ix)
  ix0=ix0w(ix)
  ix1=ix1w(ix)

  do iy=1,nyu
    iyf=iyfw(iy)
    iy0=iy0w(iy)
    iy1=iy1w(iy)

    qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*qq(iy0,ix0) &
                     & +w10(iyf,ixf)*qq(iy1,ix0) &
                     & +w01(iyf,ixf)*qq(iy0,ix1) &
                     & +w11(iyf,ixf)*qq(iy1,ix1)
  enddo
enddo

do ix=1,nxu
  do iy=1,nyu
    qar4(iy,ix)=real(qa(iy,ix))
  enddo
enddo

write(*,*) ' Writing output file...'

open(44,file='fine/qq'//pind//'.r4',form='unformatted',status='replace')
write(44) real(t),qar4
close(44)

write(*,'(a,i8,a,i8)') ' Fine grid dimensions are: ',nxu,', by ',nyu

return
end subroutine

!==============================================

end program
