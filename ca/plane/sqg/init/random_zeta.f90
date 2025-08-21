program random_zeta
!-----------------------------------------------------------------
!    Generates a random phased vorticity distribution with an
!    enstrophy spectrum ~ k^{2p-1} * exp[-(p-1)*(k/k_0)^2], p > 1.
!    The buoyancy field is obtained from this vorticity field
!    using SQG theory.  
!-----------------------------------------------------------------

use spectral

implicit none

double precision:: zs(nx,ny),zz(ny,nx)
double precision:: bs(nx,ny),bb(ny,nx)
double precision:: uni,pow,ak0,alpha,rom,ak0sqi,aksq,s,p1
double precision:: phix,phiy,cx,cy,sx,sy,amp,fmult
double precision:: bmin,bmax,zmin,zmax

integer, dimension(:), allocatable :: seed
integer:: kx,ky,kxc,kyc,k,i,ngen

! Initialise inversion constants and arrays:
call init_spectral

! Set exponent p in power spectrum:
pow=three

write(*,*) ' We assume initial enstrophy spectrum of the form'
write(*,*)
write(*,*) '   Z(k) = c*k^{2p-1}*exp[-(p-1)*(k/k_0)^2]'
write(*,*)
write(*,*) ' We take p = 3.  Enter k_0:'
read(*,*) ak0

write(*,*) ' Enter the the nominal Rossby number alpha (e.g. 0.125):'
read(*,*) alpha

! Redefine vorop to relate buoyancy to vorticity:
do ky=1,ny
  do kx=1,nx
    vorop(kx,ky)=-alpha*green(kx,ky)*(rkx(kx)**2+rky(ky)**2)
  enddo
enddo
 !Avoid division by zero below (average fields are zero):
vorop(1,1)=one

write(*,*) ' Enter the maximum |zeta|/f:'
read(*,*) rom

write(*,*) ' Enter an integer seed for the random number generator:'
read(*,*) ngen
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=ngen
do i=1,ngen
  call random_seed(put=seed)
enddo

! Generate enstrophy spectrum / k (actually, its square root):
ak0sqi=one/ak0**2
p1=pow-one
do ky=1,nwy+1
  do kx=1,nwx+1
    aksq=rkx(kx)**2+rky(ky)**2
    s=ak0sqi*aksq
    zs(kx,ky)=sqrt(s**p1*exp(-p1*s))
  enddo
enddo

! Apply to generate full spectrum:
do ky=2,nwy
  kyc=ny+2-ky
  do kx=2,nwx
    kxc=nx+2-kx
    call random_number(uni)
    phix=twopi*uni-pi
    call random_number(uni)
    phiy=twopi*uni-pi
    cx=cos(phix)
    sx=sin(phix)
    cy=cos(phiy)
    sy=sin(phiy)
    amp=zs(kx,ky)
    zs(kx ,ky )=amp*cx*cy
    zs(kxc,ky )=amp*sx*cy
    zs(kx, kyc)=amp*cx*sy
    zs(kxc,kyc)=amp*sx*sy
  enddo
enddo

ky=1
do kx=2,nwx
  kxc=nx+2-kx
  call random_number(uni)
  phix=twopi*uni-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=zs(kx,ky)
  zs(kx ,ky )=amp*cx
  zs(kxc,ky )=amp*sx
enddo

kx=1
do ky=2,nwy
  kyc=ny+2-ky
  call random_number(uni)
  phiy=twopi*uni-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=zs(kx,ky)
  zs(kx ,ky )=amp*cy
  zs(kx, kyc)=amp*sy
enddo

ky=nwy+1
do kx=2,nwx
  kxc=nx+2-kx
  call random_number(uni)
  phix=twopi*uni-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=zs(kx,ky)
  zs(kx ,ky )=amp*cx
  zs(kxc,ky )=amp*sx
enddo

kx=nwx+1
do ky=2,nwy
  kyc=ny+2-ky
  call random_number(uni)
  phiy=twopi*uni-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=zs(kx,ky)
  zs(kx ,ky )=amp*cy
  zs(kx, kyc)=amp*sy
enddo

zs(1,1)=zero
zs(nwx+1,nwy+1)=zero

! Obtain unscaled buoyancy field:
bs=zs/vorop

! Transform to physical space:
call spctop(nx,ny,zs,zz,xfactors,yfactors,xtrig,ytrig)
call spctop(nx,ny,bs,bb,xfactors,yfactors,xtrig,ytrig)

! Work out max/min values of vorticity:
zmin=minval(zz)
zmax=maxval(zz)

! Renormalise:
fmult=rom/max(abs(zmax),abs(zmin))
zz=fmult*zz
bb=fmult*bb
zmin=zmin*fmult
zmax=zmax*fmult

 !Save average buoyancy for plotting purposes:
open(44,file='average_qq.asc',status='replace')
write(44,*) zero
close(44)

! Write data:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,bb
close(11)

open(11,file='zz_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,zz
close(11)

bmin=minval(bb)
bmax=maxval(bb)

write(*,*)
write(*,'(a,f12.7,a,f11.7)') ' min b/(fN) = ',bmin,'  &  max b/(fN) = ',bmax
write(*,'(a,f12.7,a,f11.7)') ' min zeta/f = ',zmin,'  &  max zeta/f = ',zmax

end program random_zeta
