program tetra

! Sets up any number of tetrahedral configurations (4 vortices each)
! having zero total angular momentum (only on a spherical surface).

! Written 26 July 2025 by D G Dritschel @ Wildwood Crest

use constants 

 !Declarations:
implicit none

double precision:: qc(ng,nt),xg(ng,nt),yg(ng,nt),zg(ng,nt),d(ng,nt)
double precision:: clat(ng),slat(ng)
double precision:: clon(nt),slon(nt)
double precision:: kap(4),x(4),y(4),z(4),adjust(4),scalefac(4)
double precision:: a,uni,rlat,rlon,r,phi
double precision:: dImag,dIx,dIy,dIz,dsum,qavg
integer:: ntetra,i,j,k,m

!--------------------------------------------------------------------------
write(*,*) ' We set up vortex tetrahedra each having zero angular momentum'
write(*,*) ' on a sphere. Each tetrahedron consists of vortices with'
write(*,*) ' strengths kappa_j, for j = 1, 2, 3 and 4.'
write(*,*)
write(*,*) ' Enter kappa_2/kappa_1, kappa_3/kappa_1 & kappa_4/kappa_1:'
read(*,*) kap(2),kap(3),kap(4)

write(*,*)
write(*,*) ' Enter the number of tetrahedra, N:'
read(*,*) ntetra
write(*,*) ' We take kappa_1 = 4/sqrt{N} below.'
kap(1)=four/sqrt(dble(ntetra))
kap(2)=kap(1)*kap(2)
kap(3)=kap(1)*kap(3)
kap(4)=kap(1)*kap(4)
adjust=kap/sum(kap**2)

write(*,*)
write(*,*) ' The vorticity field of each vortex has the form'
write(*,*) '     omega(x,y,z) = a*kappa_j*exp(-a*s), where'
write(*,*) ' s = 2*(1-d)/(1+d) with d = x*x_j+y*y_j+z*z_j.'
write(*,*) ' (Both (x,y,z) and (x_j,y_j,z_j) are unit vectors.)'
write(*,*) ' Enter the inverse-square length a:'
read(*,*) a

!--------------------------------------------------------------------------
 !Initialize random # generator:
do i=1,iseed
   uni=rand(0)
enddo

 !Define x,y,z at all grid points:
do j=1,ng
   rlat=(dble(j)-f12)*dl-hpi
   clat(j)=cos(rlat)
   slat(j)=sin(rlat)
enddo

do i=1,nt
   rlon=dble(i-1)*dl-pi
   clon(i)=cos(rlon)
   slon(i)=sin(rlon)
enddo

do i=1,nt
   do j=1,ng
      xg(j,i)=clat(j)*clon(i)
      yg(j,i)=clat(j)*slon(i)
      zg(j,i)=slat(j)
   enddo
enddo

!--------------------------------------------------------------------------
 !Loop over tetrahedra and accumulate vorticity field:
qc=zero
do m=1,ntetra
   !Randomly pick vortex z and longitude values:
   do k=1,4
      z(k)=two*rand(0)-one
      r=sqrt(one-z(k)**2)
      phi=pi*(two*rand(0)-one)
      x(k)=r*cos(phi)
      y(k)=r*sin(phi)
   enddo

   !Adjust to have zero angular momentum:
   dImag=1.0
   do while (dImag > small)
      dIx=sum(kap*x)
      dIy=sum(kap*y)
      dIz=sum(kap*z)
      x=x-adjust*dIx
      y=y-adjust*dIy
      z=z-adjust*dIz
      scalefac=one/sqrt(x**2+y**2+z**2)
      x=scalefac*x
      y=scalefac*y
      z=scalefac*z
      dImag=sqrt(dIx**2+dIy**2+dIz**2)
   enddo

   do k=1,4
      d=x(k)*xg+y(k)*yg+z(k)*zg
      qc=qc+a*kap(k)*exp(-2.0*a*(one-d)/(one+d))
   enddo
enddo

 !Compute and remove average:
dsum=(f1112*(clat(1)+clat(ng))+sum(clat(2:ngm1)))*dble(nt)
qavg=zero
do i=1,nt
  qavg=qavg+f1112*(clat(1)*qc(1,i)+clat(ng)*qc(ng,i))
  do j=2,ngm1
    qavg=qavg+clat(j)*qc(j,i)
  enddo
enddo
qavg=qavg/dsum
qc=qc-qavg

 !Write initial absolute vorticity field:
open(20,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qc
close(20)

write(*,*) ' PV jump based on 80 contour levels:'
write(*,*) (maxval(qc)-minval(qc))/80.d0

end program tetra
