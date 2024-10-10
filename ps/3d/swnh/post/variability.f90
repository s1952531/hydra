!#########################################################################
!  Computes the "variability" in theta of the horizontal velocity spectra
!  by finding the vertical average Fourier coefficient for each wavevector
!  then the variance.  The variability is found by summing the variance
!  (mean square) over wavenumber shells and dividing this by the square 
!  of the mean also summed.  Finally, the square root of this ratio gives
!  the variability.

!           Written 4/7/2019 by D G Dritschel @ St Andrews
!#########################################################################

program variability

 !Import contants, parameters and common arrays from spectral module:
use spectral

implicit none

 !Physical arrays:
double precision:: u(ng,ng,0:nz),v(ng,ng,0:nz)
double precision:: r(ng,ng,0:nz),zeta(ng,ng,0:nz)
double precision:: wkp(ng,ng)

 !Spectral arrays:
double precision:: qs(ng,ng,0:nz),ds(ng,ng,0:nz),gs(ng,ng,0:nz)
double precision:: usbar(ng,ng),vsbar(ng,ng)
double precision:: usvar(ng,ng),vsvar(ng,ng)

 !Other local variables:
double precision:: sa(0:ng),sb(0:ng),sv(ng),t
real:: tr4,q3dr4(ng,ng,0:nz)
integer:: loop,iz,kx,ky,k
character(len=4):: pind

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

write(*,*) ' Enter the time you wish to analyse:'
read(*,*) t
loop=nint(t/tgsave)+1
write(pind,'(i4.4)') loop

 !Open input data files:
open(31,file='3d/ql.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
open(32,file= '3d/d.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
open(33,file= '3d/g.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)

 !Read data:
read(31,rec=loop) tr4,q3dr4
zeta=dble(q3dr4)
call ptospc3d(zeta,qs,0,nz)

read(32,rec=loop) tr4,q3dr4
zeta=dble(q3dr4)
call ptospc3d(zeta,ds,0,nz)

read(33,rec=loop) tr4,q3dr4
zeta=dble(q3dr4)
call ptospc3d(zeta,gs,0,nz)

close(31)
close(32)
close(33)

 !Obtain horizontal velocity field (u,v) by inversion:
call main_invert(qs,ds,gs,r,u,v,zeta)

 !Convert to spectral space (u -> ds and v -> gs):
call spctop3d(u,ds,0,nz)
call spctop3d(v,gs,0,nz)

 !Obtain average over theta:
usbar=zero
vsbar=zero
do iz=0,nz
  usbar=usbar+weight(iz)*ds(:,:,iz)
  vsbar=vsbar+weight(iz)*gs(:,:,iz)
enddo

 !Obtain variance over theta:
usvar=zero
vsvar=zero
do iz=0,nz
  usvar=usvar+weight(iz)*(ds(:,:,iz)-usbar)**2
  vsvar=vsvar+weight(iz)*(gs(:,:,iz)-vsbar)**2
enddo

 !Combine each into a single array for ease of processing:
usbar=usbar**2+vsbar**2
usvar=usvar+vsvar

 !Sum variance and square of the mean over wavenumber shells:
sa=zero
sb=zero

 !x and y-independent mode:
k=kmag(1,1)
sa(k)=sa(k)+f14*usvar(1,1)
sb(k)=sb(k)+f14*usbar(1,1)

 !y-independent mode:
do kx=2,ng
  k=kmag(kx,1)
  sa(k)=sa(k)+f12*usvar(kx,1)
  sb(k)=sb(k)+f12*usbar(kx,1)
enddo

 !x-independent mode:
do ky=2,ng
  k=kmag(1,ky)
  sa(k)=sa(k)+f12*usvar(1,ky)
  sb(k)=sb(k)+f12*usbar(1,ky)
enddo

 !All other modes:
do ky=2,ng
  do kx=2,ng
    k=kmag(kx,ky)
    sa(k)=sa(k)+usvar(kx,ky)
    sb(k)=sb(k)+usbar(kx,ky)
  enddo
enddo

 !Average successive values of k to remove noise:
do k=1,kmaxred
  sv(k)=f12*(sa(k-1)+sa(k))
enddo
sa(1:kmaxred)=sv(1:kmaxred)

do k=1,kmaxred
  sv(k)=f12*(sb(k-1)+sb(k))
enddo
sb(1:kmaxred)=sv(1:kmaxred)

 !Take ratio and square root to define variability:
sv(1:kmaxred)=sqrt(sa(1:kmaxred)/sb(1:kmaxred))

 !Write data:
open(50,file='var'//pind//'.asc',status='replace')
do k=1,kmaxred
  write(50,'(1x,f9.3,1x,f15.9)') dble(k)-f12,sv(k)
enddo

write(*,*)
write(*,*) ' The variability vs k is available in var'//pind//'.asc'

 !End main program
end program variability
!=======================================================================
