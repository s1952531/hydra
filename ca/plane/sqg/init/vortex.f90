program vortex
!-----------------------------------------------------------------
!    Superposes a random phased vorticity distribution with a 
!    spectrum Z(k) = c k^{2p-1} * exp[-(p-1)*(k/k_0)^2], p > 1,
!    and a circular vortex.
!-----------------------------------------------------------------

use spectral

implicit none
double precision:: za(ny,nx),zz(ny,nx),bb(ny,nx)
double precision:: zs(nx,ny),bs(nx,ny)
double precision:: zrat,pow,x0,y0,xfac,yfac,x,y,ssq,zavg
double precision:: uni,ak0,zpert,ak0sqi,aksq,s,p1,rom,alpha
double precision:: cx,cy,sx,sy,phix,phiy,amp
double precision:: bmin,bmax,zmin,zmax,fmult

integer, dimension(:), allocatable :: seed
integer:: i,ix,iy,k,ngen,kx,ky,kxc,kyc

! Initialise inversion constants and arrays:
call init_spectral

write(*,*)
write(*,*) ' We take zeta/f = A*(1 - s)^p + R, where s = (x/x_0)^2+(y/y_0)^2'
write(*,*) ' for s < 1, and R(x,y) a random function.'
write(*,*) ' *** Use p = 0 to instead take zeta/f = A*e^{-s} + R.'
write(*,*)
write(*,*) ' Enter the maximum zeta/f:'
read(*,*) rom
write(*,*) ' max(R)/A, p, x_0 and y_0:'
read(*,*) zrat,pow,x0,y0

xfac=one/x0
yfac=one/y0

if (pow > zero) then
   do ix=1,nx
      x=xfac*(xmin+glx*dble(ix-1))
      do iy=1,ny
         y=yfac*(ymin+gly*dble(iy-1))
         ssq=x**2+y**2
         if (ssq < one) then
            zz(iy,ix)=(one-ssq)**pow
         else
            zz(iy,ix)=zero
         endif
      enddo
   enddo
else
   do ix=1,nx
      x=xfac*(xmin+glx*dble(ix-1))
      do iy=1,ny
         y=yfac*(ymin+gly*dble(iy-1))
         zz(iy,ix)=exp(-x**2-y**2)
      enddo
   enddo
endif

 !Calculate and remove average zz; also rescale so that max(zz) = 1:
zavg=sum(zz)/dble(nx*ny)
fmult=one/(one-zavg)
zz=fmult*(zz-zavg)

write(*,*) ' We add noise, R(x,y), with a spectrum of the form'
write(*,*)
write(*,*) '   Z(k) = c*k^{2p-1}*exp[-(p-1)*(k/k_0)^2]'
write(*,*)
write(*,*) ' We take p = 3.  Enter k_0:'
read(*,*) ak0

! Set exponent p in power spectrum:
pow=three

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

! Transform to physical space as za:
call spctop(nx,ny,zs,za,xfactors,yfactors,xtrig,ytrig)
! Get max value and scale to 1:
fmult=one/maxval(za)
za=fmult*za

! Now zz (the vortex) and the perturbation have max values of 1.
! Renormalise:
zz=zz+zrat*za
fmult=rom/maxval(zz)
zz=fmult*zz

! Obtain spectral buoyancy field:
bb=zz
call ptospc(nx,ny,bb,bs,xfactors,yfactors,xtrig,ytrig)
bs=bs/vorop
bs(1,1)=zero
call spctop(nx,ny,bs,bb,xfactors,yfactors,xtrig,ytrig)

! Get all min & max values:
zmin=minval(zz)
zmax=maxval(zz)
bmin=minval(bb)
bmax=maxval(bb)

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

write(*,*)
write(*,'(a,f12.7,a,f11.7)') ' min b/(fN) = ',bmin,'  &  max b/(fN) = ',bmax
write(*,'(a,f12.7,a,f11.7)') ' min zeta/f = ',zmin,'  &  max zeta/f = ',zmax

end program vortex
