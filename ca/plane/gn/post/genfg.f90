program genfg

! Program to create an image of the PV field on the ultra-fine grid

use common

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Declarations:
real:: tr4,qdr4(ng,ng)
character(len=3):: pind
double precision:: q0(ng,ng)

!-----------------------------------------------------------------
write(*,*) ' Initialising ...'
call init_contours

write(*,*) ' What is the numbered suffix of the file to process?'
read(*,*) iind
write(pind(1:3),'(i3.3)') iind

 !Open and read PV contours:
open(40,file='cont/qqsynopsis.asc',status='old')
i=0
do while (i .lt. iind)
  read(40,*) nq,nptq,tr8,qjump
  i=int(tr8/tcsave)
enddo
close(40)

 !See if file exists:
if (i .ne. iind) then
  write(*,*) ' *** File not found!  Exiting.'
  stop
endif

write(*,*) ' nq = ',nq,' nptq = ',nptq,' t = ',tr8

 !Open and read residual needed to build ultra-fine-grid PV:
open(40,file='cont/qqresi.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
read(40,rec=iind) tr4,qdr4
close(40) 
do ix=1,ng
   do iy=1,ng
    q0(iy,ix)=dble(qdr4(iy,ix))
  enddo
enddo

write(*,*) ' Read residual PV at t = ',tr4

open(40,file='cont/qqindex'//pind,form='unformatted',status='old')
read(40) npq(1:nq),i1q(1:nq),indq(1:nq)
close(40)

open(40,file='cont/qqnodes'//pind,form='unformatted',status='old')
read(40) xq(1:nptq),yq(1:nptq)
close(40)

 !Reconstruct nextq array:
do j=1,nq
  ibeg=i1q(j)
  iend=ibeg+npq(j)-1
  do i=ibeg,iend-1
    nextq(i)=i+1
  enddo
  nextq(iend)=ibeg
enddo 

 !Generate and save fine grid field:
call ugridsave(q0,tr4)

contains 

!======================================================================

subroutine ugridsave(q0,tr4)

use congen

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: q0(ng,ng)
real tr4

 !Local arrays:
real:: rqa(ngu,ngu),rrqa(ngu/2,ngu/2)

!----------------------------------------------------------------------
! Perform contour->grid conversion:
call con2ugrid

! Bi-linear interpolate the residual q0 to the fine grid and add to qa:
do ix=1,ngu
  ixf=ixfw(ix)
  ix0=ix0w(ix)
  ix1=ix1w(ix)

  do iy=1,ngu
    iyf=ixfw(iy)
    iy0=ix0w(iy)
    iy1=ix1w(iy)

    qa(iy,ix)=qa(iy,ix)+w00(iyf,ixf)*q0(iy0,ix0) &
                     & +w10(iyf,ixf)*q0(iy1,ix0) &
                     & +w01(iyf,ixf)*q0(iy0,ix1) &
                     & +w11(iyf,ixf)*q0(iy1,ix1)
  enddo
enddo

 !Convert to real*4:
do ix=1,ngu
  do iy=1,ngu
    rqa(iy,ix)=real(qa(iy,ix))
  enddo
enddo

 !Open output file:
open(44,file='fine/qq'//pind//'.r4',form='unformatted', &
      access='stream',status='replace')

write(*,*)
write(*,*) ' Choose the type of image:'
write(*,*) '   (1) Full domain view using all data points'
write(*,*) '   (2) As above, but using every nth point'
write(*,*) '   (3) Partial domain view'
read(*,*) iopt

write(*,*) ' Writing fine/qq'//pind//'.r4'
if (iopt .eq. 1) then
  ngv=ngu
  ngv=ngu
  write(44) tr4,rqa
else if (iopt .eq. 2) then
  write(*,*) ' Enter the stride, n:'
  read(*,*) inc
  ngv=ngu/inc
  ngv=ngu/inc
  do ix=1,ngv
    ixo=(ix-1)*inc+1
    do iy=1,ngv
      iyo=(iy-1)*inc+1
      rrqa(iy,ix)=rqa(iyo,ixo)
    enddo
  enddo
  write(44) tr4,rrqa(1:ngv,1:ngv)
else
  write(*,'(a,i5,a)') ' There are ',ngu,' x grid points.'
  write(*,*) ' Range of x grid points (i1,i2) on the inversion grid?'
  write(*,*) ' (Note: periodic wrapping is allowed using i2 < i1)'
  read(*,*) ix1,ix2
  ix1v=(ix1-1)*mgu+1
  ix2v=ix2*mgu
   !Allow periodic wrapping:
  ngv=mod(ix2v-ix1v+ngu,ngu)+1

  write(*,'(a,i5,a)') ' There are ',ngu,' y grid points.'
  write(*,*) ' Range of y grid points (i1,i2) on the inversion grid?'
  write(*,*) ' (Note: periodic wrapping is allowed using i2 < i1)'
  read(*,*) iy1,iy2
  iy1v=(iy1-1)*mgu+1
  iy2v=iy2*mgu
   !Allow periodic wrapping:
  ngv=mod(iy2v-iy1v+ngu,ngu)+1
  if (iy1v .lt. iy2v) then
    if (ix1v .lt. ix2v) then
      write(44) tr4,rqa(iy1v:iy2v,ix1v:ix2v)
    else
      write(44) tr4,rqa(iy1v:iy2v,ix1v:ngu),rqa(iy1v:iy2v,1:ix2v)
    endif
  else
    if (ix1v .lt. ix2v) then
      write(44) tr4,(rqa(iy1v:ngu,ix),rqa(1:iy2v,ix),ix=ix1v,ix2v)
    else
      write(44) tr4,(rqa(iy1v:ngu,ix),rqa(1:iy2v,ix),ix=ix1v,ngu), &
                    (rqa(iy1v:ngu,ix),rqa(1:iy2v,ix),ix=1,ix2v)
    endif
  endif
endif

close(44)

write(*,*)
write(*,*) ' Image by typing the command'
write(*,*)
write(*,'(a,i5,1x,i5)') ' dataview fine/qq'//pind//'.r4 -ndim ',ngv,ngv
write(*,*)

return
end subroutine

!==============================================

end program
