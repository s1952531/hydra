program setuptest

use wavesetup

implicit double precision (a-h,o-z)

 !Parameters and array declarations:
integer,parameter:: nbytes=8+nx*(ny+1)*8
double precision:: output(0:ny,nx),yg(0:ny)

!------------------------------------------------------------
call wave_init

 !Define y on the coarse grid:
do iy=0,ny
  yg(iy)=gly*dble(iy)
enddo

write(*,'(a,2(1x,f16.9))') 'cp, cp_long = ',cp,cplong

open(10,file='basic.asc',status='unknown')
do iy=0,nyf
  write(10,'(3(1x,f20.14))') glyf*dble(iy),bbarf(iy),bfsf(iy)
enddo
close(10)

psirms=zero
do ix=1,nx
  do iy=1,nym1
    psirms=psirms+psi(iy,ix)**2
  enddo
enddo
amp=sqrt(psirms/dble(ny*nx))/cp

do ix=1,nx
  output(0,ix)=zero
  do iy=1,nym1
    output(iy,ix)=psi(iy,ix)+cp*yg(iy)
  enddo
  output(ny,ix)=cp*elly
enddo

open(11,file='pp_ini.r8',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
write(11,rec=1) amp,output
close(11)

end program


