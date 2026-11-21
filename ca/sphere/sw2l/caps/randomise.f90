module randomise

!===========================================================================
! Module which generates a random field from a selected range of spherical
! harmonics with degree (or order) n satisfying ordbeg <= n <= ordend.

! The random field ff is obtained on a latitude-longitude grid using:
!       call generate_random_field(ff, fmax, ordbeg, ordend)
! where fmax is the maximum value of |ff|, while ff is a double precision
! array dimensioned ff(numLatGridPts,numLongGridPts).

! Adapted from force.f90 by D. G. Dritschel, St Andrews UK, 24 October 2026
!===========================================================================

 !Import parameters and constants:
use constants

 !Import spherical module (stand alone):
use spherical

implicit none

contains

!=========================================================================
  
subroutine generate_random_field(ff, fmax, ordbeg, ordend)
 !Generates a field ff from a superposition of spherical harmonics.
 !The maximum value of |ff| is passed in fmax.

implicit none

 !Passed variables:
double precision:: ff(numLatGridPts,numLongGridPts), fmax
integer:: ordbeg, ordend

 !Local variable:
integer:: nsph

 !Total number of spherical harmonics used:
nsph=ordend+1-ordbeg+(ordend*(ordend+1)-ordbeg*(ordbeg-1))/2

 !Obtain the field ff:
call genff(ff, fmax, ordbeg, ordend, nsph)

return
end subroutine generate_random_field

!=========================================================================
  
subroutine genff(ff, fmax, ordbeg, ordend, nsph)
 !This does all the work apart from defining nsph

implicit none

 !Passed variables:
double precision:: ff(numLatGridPts,numLongGridPts), fmax
integer:: ordbeg, ordend, nsph

 !Local variables:
double precision:: norm(ordbeg:ordend)
double precision:: p(numLatGridPts), poly(numLatGridPts,nsph), wave(nsph)
double precision:: sinphi(numLatGridPts), cosphi(numLatGridPts), lambda(numLongGridPts)
double precision:: wka(numLatGridPts,numLongGridPts), cosfac(numLongGridPts)
double precision:: kbi, wrat, theta, sfac, pfac
integer, dimension(:), allocatable :: seed
integer:: i, j, k, m, n

!------------------------------------------------------------------
 !Initialise normalisation for spherical harmonics.  In general,
 !norm(n) = sqrt(S(n)/(2*n+1)) where S(n) is the spectrum and
 !n is the order.

if (ordbeg > 1) then
   !Here we choose a normalisation corresponding to a flat spectrum
   !S(n) = 1 (usually narrow-band), where n = order:
  do n = ordbeg, ordend
    norm(n) = one/sqrt(dble(1+2*n))
  enddo
else
   !We consider the spectrum S(n) = wrat*exp(-wrat) where 
   !wrat = (n/kb)^3 and kb = 4*ordend/9:
  kbi = 2.25d0/dble(ordend)
  do n = ordbeg, ordend          !Note ordbeg = 1 here
    wrat = (kbi*dble(n))**3
    norm(n) = wrat*exp(-wrat)/sqrt(dble(1+2*n))
  enddo
   !Note, ordend is chosen so that norm < 10^{-4} when n = ordend+1.
   !That is, 2.25 factor in the definition of kbi ensures this, or
   !equivalently the relation ordend = 9*kb/4 ensures this.
endif

!------------------------------------------------------------------
 !Define half-grid values of sinphi = sin(latitude):
do j = 1, numLatGridPts
  sinphi(j) = sin(dl*(dble(j)-0.5d0)-hpi)
enddo

 !Define cos(lat):
cosphi = sqrt(1.d0-sinphi**2)

 !Define longitudes:
do i = 1, numLongGridPts
  lambda(i) = dl*dble(i-1)-pi
enddo

 !Initialise random number generator:
call random_seed(size = k)
allocate(seed(1:k))
seed(:) = rngSeed
do i = 1, rngSeed
  call random_seed(put = seed)
enddo

!------------------------------------------------------------------
 !Generate normalised associated Legendre polynomials:
k = 0
do n = ordbeg, ordend
   !Normalisation factor:
  pfac = norm(n)

   !m = 0; include a factor of 1/2 (usual cosine series weight):
  k = k+1                         ! Index to store the polynomials
  wave(k) = 0.d0                  ! Stores azimuthal wavenumber
  call generate_spherical(sinphi,numLatGridPts,n,0,p)
   !p = normalised Y_n^0 without the azimuthal dependence
  poly(:,k) = 0.5d0*pfac*p        ! Stores re-scaled polynomial

   !m > 0:
  do m = 1, n
    k = k+1
    wave(k) = dble(m)
    call generate_spherical(sinphi,numLatGridPts,n,m,p)
     !p = normalised Y_n^m without the azimuthal dependence
    poly(:,k) = pfac*p
  enddo

enddo

!------------------------------------------------------------------
 !Generate random field:
ff = 0.d0
do k = 1, nsph
   !Random phase (theta):
  call random_number(theta)       !   0 < theta < 1 here
  theta = pi*(2.d0*theta-1.d0)    ! -pi < theta < pi now

   !cos(m*lambda + theta):
  cosfac = cos(wave(k)*lambda+theta)

   !Sum poly*cosfac into ff:
  do i = 1,numLongGridPts
     ff(:,i) = ff(:,i)+poly(:,k)*cosfac(i)
  enddo
enddo

 !Scale ff to have maximum absolute value:
sfac = fmax/maxval(abs(ff))
ff = sfac*ff

return
end subroutine genff

end module randomise
