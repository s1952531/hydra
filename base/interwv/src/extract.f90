program extract

use constants
use sta2dfft

implicit double precision (a-h,o-z)

 !Parameters and array declarations:
integer,parameter:: nbytes=8+nx*(ny+1)*8
integer,parameter:: nxout=nx,nyout=ny
integer,parameter:: nbytesout=8+nxout*(nyout+1)*8
 !Gridded arrays:
double precision:: zz(0:ny,nx),bb(0:ny,nx)
double precision:: zzspc(nx,0:ny),bbspc(nx,0:ny)
double precision:: zzout(0:nyout,nxout),bbout(0:nyout,nxout)
double precision:: zzoutspc(nxout,0:nyout),bboutspc(nxout,0:nyout)
 !FFT arrays:
integer:: xfactors(5),yfactors(5)
double precision:: xtrig(2*nx),ytrig(2*ny)
double precision:: hrkx(nx),hrky(ny)
integer:: xfacout(5),yfacout(5)
double precision:: xtrigout(2*nxout),ytrigout(2*nyout)
double precision:: horkx(nxout),horky(nyout)

!------------------------------------------------

 !Initialise FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)
call init2dfft(nxout,nyout,ellx,elly,xfacout,yfacout,xtrigout,ytrigout,horkx,horky)

open(10,file='speed.asc',status='unknown')
open(24,file='bb.r8',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
open(25,file='zz.r8',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

 !Get state of interest:
write(*,*) ' Enter the saved state number to be extracted: '
read(*,*) nstate

do i=1,nstate
  read(10,*) uref,dum
enddo
close(10)

read(24,rec=nstate) amp,bb
read(25,rec=nstate) amp,zz
close(24)
close(25)

 !Use transforms to interpolate accurately to the output grid:
call ptospc_fc(nx,ny,bb,bbspc,xfactors,yfactors,xtrig,ytrig)
call ptospc_fc(nx,ny,zz,zzspc,xfactors,yfactors,xtrig,ytrig)

 !Pad bbspc and put into bboutspc:
do ky=0,nyout
  do kx=1,nxout
    bboutspc(kx,ky)=zero
    zzoutspc(kx,ky)=zero
  enddo
enddo

kxoff=nxout-nx
fmgu=sqrt(dble(nxout)/dble(nx))*sqrt(dble(nyout)/dble(ny))
do ky=0,ny
  do kx=1,nwx+1
    bboutspc(kx,ky)=fmgu*bbspc(kx,ky)
    zzoutspc(kx,ky)=fmgu*zzspc(kx,ky)
  enddo
  do kx=nx,nwx+2,-1
    bboutspc(kxoff+kx,ky)=fmgu*bbspc(kx,ky)
    zzoutspc(kxoff+kx,ky)=fmgu*zzspc(kx,ky)
  enddo
enddo

call spctop_fc(nxout,nyout,bboutspc,bbout,xfacout,yfacout,xtrigout,ytrigout)
call spctop_fc(nxout,nyout,zzoutspc,zzout,xfacout,yfacout,xtrigout,ytrigout)

open(24,file='bb_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytesout)
open(25,file='zz_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytesout)
write(24,rec=1) zero,bbout
write(25,rec=1) zero,zzout
close(24)
close(25)

open(10,file='extract_consts.asc',status='unknown')
write(10,'(i8,i8)') nxout,nyout
write(10,'(3(1x,f14.9))') ellx,zero,elly
write(10,'(1x,f14.9)') uref
close(10)


end program
