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
enddo

do i=1,nt
   qq(:,i)=qq(:,i)+qp
enddo

!Write total initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
     access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program ranpv
