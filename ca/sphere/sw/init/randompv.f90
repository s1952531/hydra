program ranpv

!=====================================================
! Sets up a random PV anomaly field using the forcing
! function in module force.f90.
!=====================================================

 !Import necessary modules:
use constants
use force

 !Declarations:
implicit none

<<<<<<< Updated upstream
double precision:: qq(ng,nt),qp(ng)
double precision:: norm(nbeg:nend),kbi,wrat
integer:: order,i,j

!------------------------------------------------------------
!Initialise normalisation for spherical harmonics.  In general,
!norm(n) = sqrt(S(n)/(2*n+1)) where S(n) is the spectrum and
!n is the order.

if (nbeg > 1) then
   !Here we choose a normalisation corresponding to a flat spectrum
   !S(n) = 1 (usually narrow-band), where n = order:
   do order = nbeg, nend
      norm(order) = one/sqrt(dble(1+2*order))
   enddo
else
   !We consider the spectrum S(n) = wrat*exp(-wrat) where 
   !wrat = (n/kb)^3 and kb = 4*nend/9:
   kbi = 2.25d0/dble(nend)
   do order = nbeg, nend          !Note nbeg = 1 here
      wrat = (kbi*dble(order))**3
      norm(order) = wrat*exp(-wrat)/sqrt(dble(1+2*order))
   enddo
   !Note, nend is chosen so that norm < 10^{-4} when n = nend+1.
   !That is, 2.25 factor in the definition of kbi ensures this, or
   !equivalently the relation nend = 9*kb/4 ensures this.
endif

!Initialise forcing module (and spherical module):
call init_forcing(norm)

!Generate random PV anomaly field:
call generate_forcing(qq,brms)

!Add resting-state PV, 2*Omega*sin(latitude):
do j=1,ng
  qp(j)=fpole*sin((dble(j)-f12)*dl-hpi)
=======
double precision:: pslon(nLatGridPts,nLongGridPts),gslon(nLatGridPts,nLongGridPts)
double precision:: wka(nLatGridPts,nLongGridPts),wkb(nLatGridPts,nLongGridPts),wkc(nLatGridPts,nLongGridPts),wkd(nLatGridPts,nLongGridPts)
double precision:: var(nLatGridPts,nLongGridPts)
double precision:: pslat(nLongGridPts),gslat(nLongGridPts)
double precision:: qbar(nLatGridPts)
double precision:: eps,rk0i,fac,vrms,avrms,r

integer, dimension(:), allocatable :: seed
integer:: i,j,k,k0,jseed,ic,m

!------------------------------------------------------------
 !Initialise spectral module:
call init_spectral

write(*,*) ' The PV is set equal to f + a random field with variance'
write(*,*) ' spectrum proportional to k^5*exp(-2k^2/k_0^2), with zero'
write(*,*) ' global average and an rms value equal to vrms.'
write(*,*)

write(*,*) ' Enter the rms PV anomaly relative to f_pole:'
read(*,*) eps
vrms=eps*fpole

write(*,*) ' Enter k_0:'
read(*,*) k0

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
rk0i=one/dble(k0)
do m=1,nLongGridPts
  do j=1,nLatGridPts
    pslon(j,m)=(rk0i*wave(m)*clati(j))**2
    gslon(j,m)=exp(-pslon(j,m))
  enddo
enddo
do k=1,nLongGridPts
  pslat(k)=(rk0i*wave(k))**2
  gslat(k)=exp(-pslat(k))
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
    wka(i,k)=wka(i,k)*gslat(k)
    wkb(i,k)=wka(i,k)*pslat(k)
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
var=gslon*(pslon*wkc+wkd)

 !Return var to physical space:
call revfft(nLatGridPts,nLongGridPts,var,trig,factors)

 !Remove global mean:
call zeroavg(var)

 !Find rms value:
call getrms(var,avrms)

 !Normalise:
fac=vrms/avrms
var=fac*var

qbar=zero
do i=1,nLongGridPts
  do j=1,nLatGridPts
    qbar(j)=qbar(j)+var(j,i)**2
  enddo
enddo
qbar=sqrt(clati*qbar/dble(nLongGridPts))
open(77,file='qbar_init.asc',status='replace')
do j=1,nLatGridPts
  write(77,*) qbar(j),(dble(j)-f12)*dl-hpi
>>>>>>> Stashed changes
enddo

<<<<<<< Updated upstream
do i=1,nt
   qq(:,i)=qq(:,i)+qp
=======
 !Add f:
do i=1,nLongGridPts
  do j=1,nLatGridPts
    var(j,i)=var(j,i)+corFreq(j)
  enddo
>>>>>>> Stashed changes
enddo

!Write total initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
     access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program ranpv
