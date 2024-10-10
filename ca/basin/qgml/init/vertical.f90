program vertical

! This routine computes the vertical eigenvalues and eigenvectors
! of the linear system arising in the layer formulation.
! *** Uses lapack to solve the eigensystem.

use parameters

implicit none

integer,parameter:: nzm1=nz-1

 !Declarations
double precision:: hhat(nz),kdsq(nz),bmat(nz),binv(nz)
double precision:: amat(nz,nz),ainv(nz,nz),eval(nz),work(4*nz)
double precision:: fac

integer:: ipiv(nz)
integer:: iz,info,j

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
 !Note: kdsq(nz) is unused; zero can be written for this entry.

 !Ensure hhat sums to 1:
fac=1.d0/sum(hhat)
hhat=fac*hhat

!-----------------------------------------------------------------
 !Set up (tri-diagonal) matrix A:
amat=0.d0
do iz=1,nzm1
   amat(iz,iz+1)=-kdsq(iz)/hhat(iz)
   amat(iz+1,iz)=-kdsq(iz)/hhat(iz+1)
enddo
amat(1,1)=-amat(1,2)
do iz=2,nzm1
   amat(iz,iz)=-amat(iz,iz+1)-amat(iz,iz-1)
enddo
amat(nz,nz)=-amat(nz,nzm1)

 !Symmetrise:
bmat=sqrt(hhat(1)/hhat)
binv=1.d0/bmat
amat(1,2)=binv(1)*amat(1,2)*bmat(2)
amat(nz,nzm1)=binv(nz)*amat(nz,nzm1)*bmat(nzm1)
do iz=2,nzm1
   do j=-1,1
      amat(iz,iz+j)=binv(iz)*amat(iz,iz+j)*bmat(iz+j)
   enddo
enddo

!-----------------------------------------------------------------
 !Solve symmetric eigensystem with sorted eigenvalues (info = 0):
info=0
call dsyev('V','U',nz,amat,nz,eval,work,4*nz,info)

 !Barotropic mode is always present; ensure eigenvalue is zero:
eval(1)=0.d0

 !Normalise so that sum of barotropic vector is 1:
fac=0.d0
do iz=1,nz
   fac=fac+binv(iz)*amat(iz,1)
enddo
fac=1.d0/fac
do j=1,nz
   do iz=1,nz
      amat(iz,j)=fac*binv(iz)*amat(iz,j)
   enddo
enddo

 !Find inverse matrix "ainv":
ainv=amat
call dgetrf(nz,nz,ainv,nz,ipiv,info)
call dgetri(nz,ainv,nz,ipiv,work,4*nz,info)

!-----------------------------------------------------------------
 !Write data:
open(60,file='modes.asc',status='replace')
 !Eigenvalues:
do j=1,nz
   write(60,'(1x,f17.11)') eval(j)
enddo
 !Matrices for converting mode quantities to layer quantities & vice versa:
do j=1,nz
   do iz=1,nz
      write(60,'(2(1x,f17.11))') amat(iz,j),ainv(j,iz)
   enddo
enddo
 !Note transpose is necessary in ainv.
close(60)
 !These are vm2l & vl2m in spectral.f90

end program vertical
