program powerspec
!  -------------------------------------------------------------------------
!  |  Computes the 1D power spectrum of h, delta, zeta and q (the PV       |
!  |  anomaly) from data in the grid subdirectory, at all available times. |
!  |                                                                       |
!  |  Output is to the formatted "spectra.asc", which may be viewed using  |
!  |  spec_view (which takes this file as a default).                      |
!  -------------------------------------------------------------------------

! Use Spherical Harmonic Package "SHTOOLS":
use SHTOOLS
! Use constants and parameters from local job directory:
use constants

implicit none

integer, parameter:: nvar=4,lmaxc=ng/2-1

real:: tr4,qqr4(ng,nt)
double precision:: qq(ng,nt),qq0(ng,nt),cof(ng)

 !Half-grid -> Full-grid tri-diagonal arrays:
double precision:: utd(nt),xtd(nt),etd(nt),htd(nt),ptd(nt),xndeno,xend

 !Spectral arrays:
double precision:: cilm(2,ng+1,ng+1),ss(ng+1),spec(ng,nvar),lwv(ng)

integer:: loop,iopt,iread,i,ic,j,l,lmax

 !Input variable names (in data files read):
character(len=2),dimension(nvar),parameter:: ivar=['hh','dd','zz','qq']

!------------------------------------------------------
 !Initialise Coriolis frequency:
do j=1,ng
  cof(j)=fpole*sin((dble(j)-f12)*dl-hpi)
enddo

 !Initialise periodic tridiagonal problem for half-grid to full-grid
 !interpolation:
htd(1)=one
ptd(1)=-f16*htd(1)
etd(1)=ptd(1)

do j=2,nt
  htd(j)=one/(one+f16*etd(j-1))
  ptd(j)=-f16*ptd(j-1)*htd(j)
  etd(j)=-f16*htd(j)
enddo

ptd(ntm1)=etd(ntm1)+ptd(ntm1)
do j=ntm2,1,-1
  ptd(j)=etd(j)*ptd(j+1)+ptd(j)
enddo

xndeno=one/(one-etd(nt)*ptd(1)-ptd(nt))

 !Log-scaled wavenumbers:
do l=1,ng
  lwv(l)=log10(dble(l))
enddo

 !Open output file:
open(20,file='spectra.asc',status='unknown')

!------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1

   !Loop over input variables:
  do iopt=1,nvar
     !Open input data file & read data:
    open(11,file='grid/'//ivar(iopt)//'.r4',form='unformatted', & 
         & access='direct',status='old',recl=nbytes)
    iread=0
    read(11,rec=loop,iostat=iread) tr4,qqr4
    if (iread .ne. 0) exit 
    close(11)

    if (iopt .eq. 1) write(*,'(a,f9.2)') ' Processing t = ',tr4

     !Convert data to double precision:
    qq=dble(qqr4)

    if (iopt .eq. 4) then
       !Subtract f from q to define anomaly:
      do i=1,nt
        do j=1,ng
          qq(j,i)=qq(j,i)-cof(j)
        enddo
      enddo
    endif

     !-------------------------------------------------------------
     !Interpolate from the half grid to the full grid in latitude:
    do i=1,ng
      ic=i+ng

       !Source vector:
      utd(1)=f23*(qq(1,i)+qq(1,ic))
      do j=2,ng
        utd(j)=f23*(qq(j,i)+qq(j-1,i))
      enddo
      utd(ngp1)=f23*(qq(ng,ic)+qq(ng,i))
      do j=ngp2,nt
        utd(j)=f23*(qq(ntp2-j,ic)+qq(ntp1-j,ic))
      enddo

       !Interpolate qq by 4th-order method (periodic):
      xtd(1)=utd(1)*htd(1)
      do j=2,nt
        xtd(j)=(utd(j)-f16*xtd(j-1))*htd(j)
      enddo
      do j=ntm2,1,-1
        xtd(j)=etd(j)*xtd(j+1)+xtd(j)
      enddo
      xtd(nt)=(etd(nt)*xtd(1)+xtd(nt))*xndeno
      xend=xtd(nt)
    
      do j=1,ntm1
        xtd(j)=ptd(j)*xend+xtd(j)
      enddo

       !Copy back into full grid array (qq0), from North Pole to
       !South Pole, but ignoring South Pole (as required by SHTOOLS):
      do j=0,ngm1
        qq0(ng-j,i)=xtd(j+1)
      enddo
      qq0(ng,ic)=xtd(1)
      do j=1,ngm1
        qq0(ng-j,ic)=xtd(ntp1-j)
      enddo
    enddo
     !Ends loops over great circles.  Interpolation complete.

     !-------------------------------------------------------------
     ! Expand qq0 in spherical harmonics:
    call SHExpandDH(qq0, ng, cilm, lmax, sampling = 2, csphase = 1, &
                    lmax_calc = lmaxc)

     ! Compute the power spectrum of qq0:
    call SHPowerSpectrum(cilm, lmax, ss)

    do l=1,lmaxc
      spec(l,iopt)=log10(ss(l+1)+1.d-32)
    enddo

   !End loop over input variables:
  enddo

   ! Write log10(spectrum) to the output file:
  write(20,'(f12.5)') tr4
  do l=1,lmaxc
    write(20,'(5(1x,f12.8))') lwv(l),(spec(l,iopt),iopt=1,nvar)
  enddo

  if (iread .ne. 0) exit 
 !End loop over available time periods:
enddo

close(20)

write(*,*)
write(*,*) ' The h, delta, zeta & q spectra are in spectra.asc'
write(*,*) ' View them by typing'
write(*,*) '  spec_view'

end program
