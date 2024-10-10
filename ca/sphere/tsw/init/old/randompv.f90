program ranpv

!=====================================================
!         Sets up a random PV anomaly field
!=====================================================

 !Import contants, parameters and common arrays:
use constants
use variables
use spectral

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: qq(ng,nt)

!------------------------------------------------------------
 !Initialise spectral module:
call init_spectral

 !Initialize random # generator:
do i=1,iseed
  uni=rand(0)
enddo

 !Generate squared wavenumber arrays used in the forcing:
rksri=one/dble(ksr)
do m=1,nt
  do j=1,ng
    plon(j,m)=(rksri*wave(m)*clati(j))**2
    glon(j,m)=exp(-plon(j,m))
  enddo
enddo
do k=1,nt
  plat(k)=(rksri*wave(k))**2
  glat(k)=exp(-plat(k))
enddo

 !--------------------------------------------------------
write(*,*) ' Enter the rms PV anomaly relative to f_pole:'
read(*,*) eps

 !generate random PV anomaly field qq with this rms value:
call ranspec(qq,eps)

 !Add f:
do i=1,nt
  do j=1,ng
    qq(j,i)=qq(j,i)+cof(j)
  enddo
enddo

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program
