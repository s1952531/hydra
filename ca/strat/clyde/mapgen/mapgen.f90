!#####################################################################
!                      Conformal Map Test Code
!#####################################################################

program mapgen

 !Import contants, parameters and common arrays needed for inversion etc:
use spectral

implicit none

 !Arrays needed for comparison with original map (X,Y):
double precision:: xx(0:ny,0:nx),yy(0:ny,0:nx)
double precision:: dyydx(0:ny,0:nx),dyydy(0:ny,0:nx)

 !Lower and upper edge values of Y needed to generate new map:
double precision:: ybot(0:nx),ytop(0:nx)

 !Work array:
double precision:: zz(0:ny,0:nx)

!-----------------------------------------------------------------
 !Initialise inversion constants and arrays and read original map:
call init_spectral

 !Copy original map:
xx=xori
yy=yori
dyydx=dyoridx
dyydy=dyoridy

 !Store edge values of Y:
ybot=yori(0,:)
ytop=yori(ny,:)

 !Generate new map:
call conform(ybot,ytop)

 !Write new map and derivatives:
open(41,file='newyori.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(41,rec=1) real(yori)
write(41,rec=2) real(dyoridx)
write(41,rec=3) real(dyoridy)
close(41)

open(41,file='newxori.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(41,rec=1) real(xori)
!write(41,rec=2) real(dxoridx)
!write(41,rec=3) real(dxoridy)
close(41)

 !Write difference fields:
yori=yori-yy
dyoridx=dyoridx-dyydx
dyoridy=dyoridy-dyydy
open(41,file='diffnewyori.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(41,rec=1) real(yori)
write(41,rec=2) real(dyoridx)
write(41,rec=3) real(dyoridy)
close(41)

 !Write difference fields:
xori=xori-xx
!dxoridx=dxoridx-dxxdx
!dxoridy=dxoridy-dxxdy
open(41,file='diffnewxori.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(41,rec=1) real(xori)
!write(41,rec=2) real(dxoridx)
!write(41,rec=3) real(dxoridy)
close(41)

 !Compute rms and max differences:
write(*,*) ' Max abs error in X = ',maxval(abs(xori))
zz=xori**2
write(*,*) ' R.m.s.  error in X = ', &
   (f14*(zz(0,0)+zz(0,nx)+zz(ny,0)+zz(ny,nx)) &
   +f12*(sum(zz(1:nym1,0)+zz(1:nym1,nx))+sum(zz(0,1:nxm1)+zz(ny,1:nxm1))) &
   +sum(zz(1:nym1,1:nxm1)))*dsumi

write(*,*) ' Max abs error in Y = ',maxval(abs(yori))
zz=yori**2
write(*,*) ' R.m.s.  error in Y = ', &
   (f14*(zz(0,0)+zz(0,nx)+zz(ny,0)+zz(ny,nx)) &
   +f12*(sum(zz(1:nym1,0)+zz(1:nym1,nx))+sum(zz(0,1:nxm1)+zz(ny,1:nxm1))) &
   +sum(zz(1:nym1,1:nxm1)))*dsumi

write(*,*) ' Max abs error in Y_x = ',maxval(abs(dyoridx))
zz=dyoridx**2
write(*,*) ' R.m.s.  error in Y_x = ', &
   (f14*(zz(0,0)+zz(0,nx)+zz(ny,0)+zz(ny,nx)) &
   +f12*(sum(zz(1:nym1,0)+zz(1:nym1,nx))+sum(zz(0,1:nxm1)+zz(ny,1:nxm1))) &
   +sum(zz(1:nym1,1:nxm1)))*dsumi

write(*,*) ' Max abs error in Y_y = ',maxval(abs(dyoridy))
zz=dyoridy**2
write(*,*) ' R.m.s.  error in Y_y = ', &
   (f14*(zz(0,0)+zz(0,nx)+zz(ny,0)+zz(ny,nx)) &
   +f12*(sum(zz(1:nym1,0)+zz(1:nym1,nx))+sum(zz(0,1:nxm1)+zz(ny,1:nxm1))) &
   +sum(zz(1:nym1,1:nxm1)))*dsumi

end program mapgen
