program check

! Reads coords.f90 to extract a selected grid curve; used to check that
! different resolutions have corresponding grid points.

use constants

implicit none
  
 !Arrays related to the conformal map:
double precision:: xo(0:ny,0:nx),yo(0:ny,0:nx)

integer:: iopt,ix,iy

!-------------------------------------------------------------
 !Read in conformal map:
open(11,file='coords.r8',form='unformatted', &
    & access='direct',status='old',recl=2*(nbytes-4))
read(11,rec=1) xo
read(11,rec=2) yo
close(11)

!-------------------------------------------------------------
write(*,*) ' Which kind of grid line do you wish to extract:'
write(*,*) ' (1) y = const (in conformal domain), or'
write(*,*) ' (2) x = const (in conformal domain)?'
read(*,*) iopt

if (iopt .eq. 1) then
  write(*,*) ' Enter the y index of the coordinate curve:'
  read(*,*) iy
  open(22,file='ycoord.asc',status='replace')
  do ix=0,nx
    write(22,*) xo(iy,ix),yo(iy,ix)
  enddo
  write(*,*) ' You can plot the coordinate curve using'
  write(*,*) ' plotcol ycoord.asc'
else
  write(*,*) ' Enter the x index of the coordinate curve:'
  read(*,*) ix
  open(22,file='xcoord.asc',status='replace')
  do iy=0,ny
    write(22,*) xo(iy,ix),yo(iy,ix)
  enddo
  write(*,*) ' You can plot the coordinate curve using'
  write(*,*) ' plotcol xcoord.asc'
endif


end program check
