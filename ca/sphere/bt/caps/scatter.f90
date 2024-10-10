program scatter
!  -------------------------------------------------------------------------
!  |   Computes psi from zeta and lists this for all grid points in the    |
!  |   output file psi-vor.asc.                                            |
!  -------------------------------------------------------------------------

use spectral

implicit none

real:: tr4,zzr4(ng,nt)
double precision:: zz(ng,nt),qq(ng,nt),pp(ng,nt),wkb(ng,nt)
integer:: i,j,period

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open input data file:
open(44,file='grid/zz.r4',form='unformatted', & 
    & access='direct',status='old',recl=nbytes)

write(*,*) ' Time period (choose 0 for the initial one)?'
read(*,*) period

 !Read data:
read(44,rec=period+1) tr4,zzr4(1:ng,1:nt)
close(44)
write(*,'(a,f12.5)') ' t = ',tr4

 !Convert to double precision:
zz=dble(zzr4)
qq=zz

 !FFT in longitude (semi-spectral):
call forfft(ng,nt,qq,trig,factors) 

 !Invert Laplace's operator on the vorticity:
call laplinv(qq,pp,wkb)

 !Get psi (pp) in physical space:
call revfft(ng,nt,pp,trig,factors)

 !Write ascii data:
open(22,file='psi-vor.asc',status='replace')

do i=1,nt
  do j=1,ng
    write(22,'(2(1x,f14.10))') pp(j,i),zz(j,i)
  enddo
enddo

close(22)

write(*,*)
write(*,*) ' psi vs zeta is listed in psi-vor.asc'

end program
