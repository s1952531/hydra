program vortex
!-----------------------------------------------------------------
!    Superposes a random phased buoyancy distribution with a 
!    spectrum Q(k) = c k^{2p-1} * exp[-(p-1)*(k/k_0)^2], p > 1,
!    and a circular vortex with zero mean buoyancy.
!-----------------------------------------------------------------

use spectral

implicit none
double precision:: qs(nx,ny),qa(ny,nx),qq(ny,nx)
double precision:: bmax,pow,x0,y0,xfac,yfac,x,y,ssq,qavg
double precision:: uni,ak0,qpert,ak0sqi,aksq,s,p1
double precision:: cx,cy,sx,sy,phix,phiy,amp
double precision:: fmult,qmin,qmax

integer, dimension(:), allocatable :: seed
integer:: i,ix,iy,k,ngen,kx,ky,kxc,kyc

! Initialise inversion constants and arrays:
call init_spectral

write(*,*)
write(*,*) ' We take b_0/N = b_max*(1 - s)^p where s = (x/x_0)^2 + (y/y_0)^2'
write(*,*) ' for s < 1, and b_0/N = 0 otherwise. *** Use p = 0 to instead'
write(*,*) ' take b_0/N = b_max*e^{-s}.'
write(*,*)
write(*,*) ' Enter b_max, p, x_0 and y_0:'
read(*,*) bmax,pow,x0,y0

xfac=one/x0
yfac=one/y0

if (pow > zero) then
   do ix=1,nx
      x=xfac*(xmin+glx*dble(ix-1))
      do iy=1,ny
         y=yfac*(ymin+gly*dble(iy-1))
         ssq=x**2+y**2
         if (ssq < one) then
            qq(iy,ix)=(one-ssq)**pow
         else
            qq(iy,ix)=zero
         endif
      enddo
   enddo
else
   do ix=1,nx
      x=xfac*(xmin+glx*dble(ix-1))
      do iy=1,ny
         y=yfac*(ymin+gly*dble(iy-1))
         qq(iy,ix)=exp(-x**2-y**2)
      enddo
   enddo
endif

qq=bmax*qq

 !Calculate and remove average qq:
qavg=sum(qq)/dble(nx*ny)

 !Save average for plotting purposes:
open(44,file='average_qq.asc',status='replace')
write(44,'(1x,f15.11)') qavg
close(44)

qq=qq-qavg

! Set exponent p in power spectrum:
pow=three

write(*,*) ' We add noise with a spectrum of the form'
write(*,*)
write(*,*) '   Q(k) = c*k^{2p-1}*exp[-(p-1)*(k/k_0)^2]'
write(*,*)
write(*,*) ' We take p = 3.  Enter k_0:'
read(*,*) ak0

write(*,*) ' Enter the maximum noise amplitude (in terms of b_0/N):'
read(*,*) qpert

write(*,*) ' Enter an integer seed for the random number generator:'
read(*,*) ngen
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=ngen
do i=1,ngen
  call random_seed(put=seed)
enddo

! Generate potential enstrophy spectrum / k (actually, its square root):
ak0sqi=one/ak0**2
p1=pow-one
do ky=1,nwy+1
  do kx=1,nwx+1
    aksq=rkx(kx)**2+rky(ky)**2
    s=ak0sqi*aksq
    qs(kx,ky)=sqrt(s**p1*exp(-p1*s))
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
    amp=qs(kx,ky)
    qs(kx ,ky )=amp*cx*cy
    qs(kxc,ky )=amp*sx*cy
    qs(kx, kyc)=amp*cx*sy
    qs(kxc,kyc)=amp*sx*sy
  enddo
enddo

ky=1
do kx=2,nwx
  kxc=nx+2-kx
  call random_number(uni)
  phix=twopi*uni-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cx
  qs(kxc,ky )=amp*sx
enddo

kx=1
do ky=2,nwy
  kyc=ny+2-ky
  call random_number(uni)
  phiy=twopi*uni-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cy
  qs(kx, kyc)=amp*sy
enddo

ky=nwy+1
do kx=2,nwx
  kxc=nx+2-kx
  call random_number(uni)
  phix=twopi*uni-pi
  cx=cos(phix)
  sx=sin(phix)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cx
  qs(kxc,ky )=amp*sx
enddo

kx=nwx+1
do ky=2,nwy
  kyc=ny+2-ky
  call random_number(uni)
  phiy=twopi*uni-pi
  cy=cos(phiy)
  sy=sin(phiy)
  amp=qs(kx,ky)
  qs(kx ,ky )=amp*cy
  qs(kx, kyc)=amp*sy
enddo

qs(1,1)=zero
qs(nwx+1,nwy+1)=zero

! Transform to physical space:
call spctop(nx,ny,qs,qa,xfactors,yfactors,xtrig,ytrig)

! Scale:
fmult=qpert/maxval(abs(qa))
qa=fmult*qa

! Combine with vortex:
qq=qq+qa

! Get new max/min values:
qmin=minval(qq)
qmax=maxval(qq)

! Write data:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

write(*,*)
write(*,'(a,f12.7,a,f11.7)') ' min b_0/N = ',qmin,'  &  max b_0/N = ',qmax

end program vortex
