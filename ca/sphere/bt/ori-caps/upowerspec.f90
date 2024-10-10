program upowerspec
!  -------------------------------------------------------------------------
!  |  Computes the 1D power spectrum for vorticity from data in the fine   |
!  |  grid subdirectory, previously created using genfg.f90                |
!  |                                                                       |
!  |  Output is to the formatted "spectra<nnn>", where <nnn> = 000, 001,   |
!  |  etc is the period to be analysed.                                    |
!  |                                                                       |
!  |  Note, the maximum spherical harmonic degree is 2800.                 |
!  -------------------------------------------------------------------------

! Use Spherical Harmonic Package "SHTOOLS":
use SHTOOLS
! Use constants and parameters from local job directory:
use constants

implicit none

integer, parameter:: lmaxc=min(ngu/2-1,2800)

real:: tr4,qqr4(ngu,ntu)
double precision:: qq(ngu,ntu),qq0(ngu,ntu)

 !Half-grid -> Full-grid tri-diagonal arrays:
double precision:: utd(ntu),xtd(ntu),etd(ntu),htd(ntu),ptd(ntu),xndeno,xend

 !Spectral arrays:
double precision:: cilm(2,lmaxc+1,lmaxc+1),ss(lmaxc+1)

integer:: period,i,ic,j,l,lmax
character(len=3):: pind

!---------------------------------------------------------------
 !Select time frame:
write(*,*) ' Time period (choose 0 for the initial one)?'
read(*,*) period
write(pind(1:3),'(i3.3)') period

 !Open input data file:
open(44,file='fine/qq'//pind//'.r4',form='unformatted',status='old')

 !Read data:
read(44) tr4,qqr4
close(44)
write(*,'(a,f12.5)') ' t = ',tr4

 !Convert to double precision:
qq=dble(qqr4)

!------------------------------------------------------
 !Initialise periodic tridiagonal problem for half-grid to full-grid
 !interpolation:
htd(1)=one
ptd(1)=-f16*htd(1)
etd(1)=ptd(1)

do j=2,ntu
  htd(j)=one/(one+f16*etd(j-1))
  ptd(j)=-f16*ptd(j-1)*htd(j)
  etd(j)=-f16*htd(j)
enddo

ptd(ntu-1)=etd(ntu-1)+ptd(ntu-1)
do j=ntu-2,1,-1
  ptd(j)=etd(j)*ptd(j+1)+ptd(j)
enddo

xndeno=one/(one-etd(ntu)*ptd(1)-ptd(ntu))

 !-------------------------------------------------------------
 !Interpolate from the half grid to the full grid in latitude:
do i=1,ngu
  ic=i+ngu

   !Source vector:
  utd(1)=f23*(qq(1,i)+qq(1,ic))
  do j=2,ngu
    utd(j)=f23*(qq(j,i)+qq(j-1,i))
  enddo
  utd(ngu+1)=f23*(qq(ngu,ic)+qq(ngu,i))
  do j=ngu+2,ntu
    utd(j)=f23*(qq(ntu+2-j,ic)+qq(ntu+1-j,ic))
  enddo

   !Interpolate qq by 4th-order method (periodic):
  xtd(1)=utd(1)*htd(1)
  do j=2,ntu
    xtd(j)=(utd(j)-f16*xtd(j-1))*htd(j)
  enddo
  do j=ntu-2,1,-1
    xtd(j)=etd(j)*xtd(j+1)+xtd(j)
  enddo
  xtd(ntu)=(etd(ntu)*xtd(1)+xtd(ntu))*xndeno
  xend=xtd(ntu)
    
  do j=1,ntu-1
    xtd(j)=ptd(j)*xend+xtd(j)
  enddo

   !Copy back into full grid array (qq0), from North Pole to
   !South Pole, but ignoring South Pole (as required by SHTOOLS):
  do j=0,ngu-1
    qq0(ngu-j,i)=xtd(j+1)
  enddo
  qq0(ngu,ic)=xtd(1)
  do j=1,ngu-1
    qq0(ngu-j,ic)=xtd(ntu+1-j)
  enddo
enddo
 !Ends loops over great circles.  Interpolation complete.

 !-------------------------------------------------------------
 ! Expand qq0 in spherical harmonics:
call SHExpandDH(qq0, ngu, cilm, lmax, sampling = 2, csphase = 1, &
                lmax_calc = lmaxc)

 ! Compute the power spectrum of qq0:
call SHPowerSpectrum(cilm, lmaxc, ss)

 !Open and write output file:
open(20,file='fine/spectra'//pind,status='unknown')
do l=1,lmaxc
  write(20,'(2(1x,f12.8))') log10(dble(l)),log10(ss(l+1)+1.d-32)
enddo
close(20)

write(*,*)
write(*,*) ' The fine-grid vorticity spectra are in fine/spectra'//pind

end program
