program adjust

!============================================================
!  Sets up a random height anomaly field in a state of rest
!============================================================

 !Import spectral module:
use spectral

 !Declarations:
implicit none

double precision:: wka(nLatGridPts,nLongGridPts),wkb(nLatGridPts,nLongGridPts),wkc(nLatGridPts,nLongGridPts),wkd(nLatGridPts,nLongGridPts)
double precision:: var(nLatGridPts,nLongGridPts)
double precision:: eps,rksri,fac,hmax,ahmax,r

integer, dimension(:), allocatable :: seed
integer:: i,j,k,ksr,jseed,ic,m

!------------------------------------------------------------
 !Initialise spectral module:
call init_spectral

write(*,*) ' The height anomaly is set to a random field with a variance'
write(*,*) ' spectrum proportional to k^5*exp(-2k^2/ksr^2), with zero'
write(*,*) ' global average and a maximum absolute value equal to h_max.'
write(*,*)

write(*,*) ' Enter h_max:'
read(*,*) hmax

write(*,*) ' Enter ksr:'
read(*,*) ksr

write(*,*) ' Enter a random seed (integer):'
read(*,*) jseed

 !Initialize random number generator:
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=jseed
do i=1,jseed
  call random_seed(put=seed)
enddo

 !Initialize squared wavenumber arrays used below:
rksri=one/dble(ksr)
do m=1,nLongGridPts
  do j=1,nLatGridPts
    plon(j,m)=(rksri*wave(m)*clati(j))**2
    glon(j,m)=exp(-plon(j,m))
  enddo
enddo
do k=1,nLongGridPts
  plat(k)=(rksri*wave(k))**2
  glat(k)=exp(-plat(k))
enddo

!------------------------------------------------------------------
 !Generate random values between -1 and +1 at all grid points:
call random_number(wkc)
wkc=two*wkc-one

 !Create great circles:
do i=1,nLatGridPts
  ic=i+nLatGridPts
  do j=1,nLatGridPts
    wka(i,j)     =wkc(j,i)
    wka(i,nLongGridPtsPlusOne-j)=wkc(j,ic)
  enddo
enddo

 !FFT in latitude:
call forfft(nLatGridPts,nLongGridPts,wka,trig,factors)

 !Apply latitudinal spectrum:
do k=1,nLongGridPts
  do i=1,nLatGridPts
    wka(i,k)=wka(i,k)*glat(k)
    wkb(i,k)=wka(i,k)*plat(k)
  enddo
enddo

 !Inverse FFT in latitude:
call revfft(nLatGridPts,nLongGridPts,wka,trig,factors)
call revfft(nLatGridPts,nLongGridPts,wkb,trig,factors)

 !Unpack arrays:
do i=1,nLatGridPts
  ic=i+nLatGridPts
  do j=1,nLatGridPts
    wkc(j,i) =wka(i,j)
    wkc(j,ic)=wka(i,nLongGridPtsPlusOne-j)
    wkd(j,i) =wkb(i,j)
    wkd(j,ic)=wkb(i,nLongGridPtsPlusOne-j)
  enddo
enddo

 !FFT in longitude:
call forfft(nLatGridPts,nLongGridPts,wkc,trig,factors)
call forfft(nLatGridPts,nLongGridPts,wkd,trig,factors)

 !Apply longitudinal spectrum and define var:
var=glon*(plon*wkc+wkd)

 !Return var to physical space:
call revfft(nLatGridPts,nLongGridPts,var,trig,factors)

 !Remove global mean:
call zeroavg(var)

 !Normalise:
ahmax=maxval(abs(var))
fac=hmax/ahmax
var=fac*var

 !Write initial height anomaly field:
open(20,file='hh_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,var
close(20)

 !Define PV for a rest state:
do i=1,nLongGridPts
  do j=1,nLatGridPts
    var(j,i)=fpole*slat(j)/(one+var(j,i))
  enddo
enddo

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,var
close(20)

 !Write initial divergence field:
var=zero
open(20,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,var
close(20)

end program adjust
