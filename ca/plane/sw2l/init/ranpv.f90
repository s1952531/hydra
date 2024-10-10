program ranpv
!-----------------------------------------------------------------
!    Generates a random phased baroclinic (mode 2) PV distribution
!    with spectrum S(k) = c * k^{2p-3} * exp[-(p-1)*(k/k_0)^2],
!    for integer p > 1.  Here, k_0 = <k^2 S>/<S> is the spectral
!    centroid.
!-----------------------------------------------------------------

use spectral

implicit double precision(a-h,o-z)
double precision:: qs(ng,ng),qa(ng,ng)

 !Initialise inversion constants and arrays:
call init_spectral

nw=ng/2

write(*,*) ' The baroclinic (mode 2) PV has the spectrum'
write(*,*) '  S(k) = c k^{2p-3} * exp[-(p-1)*(k/k_0)^2].'
write(*,*)
write(*,*) ' Enter p > 1 & k_0:'
read(*,*) pow,ak0

write(*,*) ' Enter the maximum baroclinic PV magnitude / f:'
read(*,*) qeddy
qeddy=cof*qeddy

write(*,*) ' Enter an integer seed for the random # generator:'
read(*,*) ngen
do i=1,ngen
  uni=rand(0)
enddo

! Generate spectral amplitude sqrt{k^{2(p-2)} * exp[-(p-1)*(k/k_0)^2]} :
efac=one/ak0**2
p1=pow-one
p2=pow-two
do ky=1,nw+1
  do kx=1,nw+1
    s=efac*(rk(kx)**2+rk(ky)**2)
    qs(kx,ky)=sqrt(efac*s**p2*exp(-p1*s))
  enddo
enddo

! Apply to generate full spectrum:
do ky=2,nw
  kyc=ng+2-ky
  do kx=2,nw
    kxc=ng+2-kx
    phix=twopi*rand(0)-pi
    phiy=twopi*rand(0)-pi
    cx=cos(phix)
    sx=sin(phix)
    cy=cos(phiy)
    sy=sin(phiy)
    amp=qs(kx,ky)
    qs(kx ,ky )=amp*cx*cy
    qs(kxc,ky )=amp*sx*cy
    qs(kx, kyc)=amp*cx*sy
    qs(kxc,kyc)=amp*sx*sy
  enddo
enddo

ky=1
do kx=2,nw
  kxc=ng+2-kx
  phix=twopi*rand(0)-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cx
  qs(kxc,ky )=amp*sx
enddo

kx=1
do ky=2,nw
  kyc=ng+2-ky
  phiy=twopi*rand(0)-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cy
  qs(kx, kyc)=amp*sy
enddo

ky=nw+1
do kx=2,nw
  kxc=ng+2-kx
  phix=twopi*rand(0)-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cx
  qs(kxc,ky )=amp*sx
enddo

kx=nw+1
do ky=2,nw
  kyc=ng+2-ky
  phiy=twopi*rand(0)-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cy
  qs(kx, kyc)=amp*sy
enddo

qs(1,1)=zero
qs(nw+1,nw+1)=zero

! Transform to physical space:
call spctop(ng,ng,qs,qa,xfactors,yfactors,xtrig,ytrig)

! Work out max/min values and scale:
qamin=qa(1,1)
qamax=qamin
do ix=1,ng
  do iy=1,ng
    qamin=min(qamin,qa(iy,ix))
    qamax=max(qamax,qa(iy,ix))
  enddo
enddo

fmult=qeddy/max(abs(qamax),abs(qamin))
do ix=1,ng
  do iy=1,ng
    qa(iy,ix)=fmult*qa(iy,ix)
  enddo
enddo

! Write PV anomaly in each layer:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
fac=-vm12/(vm11*vm22-vm12*vm21)
write(11,rec=1) zero,fac*qa
fac=-fac*vm11/vm12
write(11,rec=2) zero,fac*qa
close(11)

!-------------------------------------------------------------
! Write zero fields of divergence and acceleration divergence:
qa=zero
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qa
write(11,rec=2) zero,qa
close(11)

open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qa
write(11,rec=2) zero,qa
close(11)


end program
