program bci2l

use constants

! Computes the growth rates of all unstable modes in a two-layer
! vertically-sheared flow.

! WARNING, nz = 2 is assumed!

implicit none

double precision:: hhat(nz),lambda(nz),kdsq(nz)
double precision:: vl2m(nz,nz),vm2l(nz,nz)
double precision:: uuha(nz),dqdy(nz)
double precision:: amat11,amat21,fac
double precision:: scx,scy,rkx,rky,rksq
double precision:: w1,w2,b11,b12,b21,b22
double precision:: sig,sigmax,rkxmax,rkymax
integer:: m,iz,kx,ky

!-----------------------------------------------------------------
 !Read vertical structure file:
open(60,file='vertical.asc',status='old')
do iz=1,nz
   read(60,*) hhat(iz),kdsq(iz)
enddo
close(60)
 !hhat = mean layer depth / total mean depth (sum(hhat) = 1).
 !kdsq = f^2/(b'*H) where f = Coriolis frequency, b' = buoyancy
 !       difference between layer iz and iz+1, H = total mean depth.
 !Note: kdsq(nz) is unused

 !Ensure hhat sums to 1:
fac=one/sum(hhat)
hhat=fac*hhat

 !Define A_{11} and A_{21} following vertical.f90:
amat11= kdsq(1)/hhat(1)
amat21=-kdsq(1)/hhat(2)

!-----------------------------------------------------------------
 !Read vertical eigenvalues and eigenmodes:
open(60,file='modes.asc',status='old')
do m=1,nz
   read(60,*) lambda(m)
enddo
 !Note, lambda(m) corresponds to lambda_m in the notes.
do m=1,nz
   do iz=1,nz
      read(60,*) vl2m(iz,m),vm2l(iz,m)
   enddo
enddo
close(60)
 !vl2m converts layer quantities to  mode quantities
 !vm2l converts  mode quantities to layer quantities

!--------------------------------------------------------------------
 !Read layer-mean zonal flow:
open(60,file='meanflow.asc',status='old')
do iz=1,nz
   read(60,*) uuha(iz)
enddo
 !WARNING: we assume below only uuha(1) is non-zero!

 !Define mean PV gradients in each layer:
dqdy(1)=beta+amat11*uuha(1)
dqdy(2)=beta+amat21*uuha(1)

!--------------------------------------------------------------------
 !Determine the maximum growth rate over all possible horizontal
 !wavenumbers:
sigmax=zero
rkxmax=zero
rkymax=zero
scx=twopi/ellx
scy=   pi/elly
do ky=0,ny
   rky=scy*dble(ky)
   do kx=1,nx/2
      rkx=scx*dble(kx)
      rksq=rkx**2+rky**2
      w1=one/rksq
      w2=one/(rksq+lambda(2))
      b11=uuha(1)-dqdy(1)*(w1*vm2l(1,1)*vl2m(1,1)+w2*vm2l(1,2)*vl2m(1,2))
      b22=uuha(2)-dqdy(2)*(w1*vm2l(2,1)*vl2m(2,1)+w2*vm2l(2,2)*vl2m(2,2))
      b12=       -dqdy(1)*(w1*vm2l(1,1)*vl2m(2,1)+w2*vm2l(1,2)*vl2m(2,2))
      b21=       -dqdy(2)*(w1*vm2l(2,1)*vl2m(1,1)+w2*vm2l(2,2)*vl2m(1,2))
      sig=b11*b22-b12*b21-f14*(b11+b22)**2
      if (sig > zero) then
         sig=rkx*sqrt(sig)
         if (sig > sigmax) then
            rkxmax=rkx
            rkymax=rky
            sigmax=sig
         endif
      endif
   enddo
enddo

if (sigmax == zero) then
   write(*,*) ' The flow is stable.'
else
   write(*,*) ' The maximum growth rate is ',sigmax
   write(*,*) ' and this occurs when k_x = ',rkxmax,' and k_y = ',rkymax
endif

end program bci2l
