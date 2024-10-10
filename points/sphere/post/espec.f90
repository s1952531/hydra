program spectrum
!==================================================================
! Computes the energy spectrum associated with a set of point
! vortices on a sphere

! Written 26 Aprile 2013 by D G Dritschel @ New York, USA

! See the "eprep" script for compilation
!==================================================================

use SHTOOLS
	
implicit none	
integer, parameter ::	nlat = 1024, nlon = 2*nlat, maxdeg = nlat
integer, parameter ::	nvm = 1000000
real*8, parameter ::    pi = 3.14159265358979d0, dl = pi/dble(nlat)
real*8 ::		cilm(2, maxdeg+1, maxdeg+1), psi(nlat, nlon)
real*8 ::		spectra(maxdeg+1)
real*8 ::		x(nvm), y(nvm), z(nvm), s(nvm)
real*8 ::		rlon, rlat, t, fac
real*8 ::		clon, slon, xg, yg, zg, rg
integer ::		nv, i, j, k, kfr, nfr
integer ::		lmax, lmaxc, l, m
character*4 ::		suffix

write(*,*) ' Number of point vortices?'
read(*,*) nv

! Open file containing vortex strengths and read them:
open(10,file='strengths.dat',status='old')
do i=1,nv
  read(10,*) s(i)
enddo
close(10)

! Open file containing vortex positions and read them:
open(11,file='points.dat',status='old')
write(*,*) ' Enter the frame to read (0 for initial):'
read(*,*) nfr
write(suffix(1:4),'(i4.4)') nfr

do kfr=0,nfr
  read(11,*) t
  do i=1,nv
    read(11,*) x(i),y(i),z(i)
  enddo
enddo
close(11)

! Ensure they are on the sphere:
do i=1,nv
  fac=1.d0/sqrt(x(i)**2+y(i)**2+z(i)**2)
  x(i)=fac*x(i)
  y(i)=fac*y(i)
  z(i)=fac*z(i)
enddo

write(*,'(a,f12.5)') ' Read vortices at t = ',t
write(*,*)

write(*,'(a,i4,a,i4,a)') ' Computing the streamfunction on a ',nlat, &
                         ' x ',nlon,' lat-lon grid...'

! Region between the poles:
do i=1,nlon
  rlon=dble(i-1)*dl-pi
  clon=cos(rlon)
  slon=sin(rlon)

  do j=1,nlat
     ! co-latitude starts from NP:
    rlat=dble(j-1)*dl
    zg=cos(rlat)
    rg=sin(rlat)
    xg=rg*clon
    yg=rg*slon
    psi(j,i)=0.d0
     ! Sum over point vortices and accumulate streamfunction:
    do k=1,nv
      psi(j,i)=psi(j,i)+s(k)*log(1.d0-x(k)*xg-y(k)*yg-z(k)*zg)
    enddo
  enddo
enddo

 ! Expand psi in spherical harmonics:
lmaxc=min(2800,nlat/2-1)
call SHExpandDH(psi, nlat, cilm, lmax, sampling = 2, csphase = 1, &
                lmax_calc = lmaxc)

!call SHExpandGLQ(cilm, lmax, gridglq, w, plx, zero, &
!                 norm, csphase, lmax_calc)


 ! Compute the power spectrum of psi:
call SHPowerSpectrum(cilm, lmax, spectra)	

 ! Write out data:
open(20,file='spec'//suffix,status='unknown')
do l=1, lmaxc
  write(20,'(1x,i4,1x,1p,e14.7)') l,spectra(l+1)*dble(l*(l+1))
enddo
close(20)

write(*,*) ' The energy spectrum, l vs E(l), is ready in spec'//suffix

end program spectrum
