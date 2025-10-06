program n5

! Sets up 5 vortices on a sphere, using data in 5vortex_ICs.

! Written 1 August 2025 by D G Dritschel @ Wildwood Crest

use constants 

 !Declarations:
implicit none

double precision:: qc(ng,nt),xg(ng,nt),yg(ng,nt),zg(ng,nt),d(ng,nt)
double precision:: clat(ng),slat(ng)
double precision:: clon(nt),slon(nt)
double precision:: kap(5),x(5),y(5),z(5)
double precision:: a,rlat,rlon,dsum,qavg
integer:: i,j,k

!--------------------------------------------------------------------------
open(20,file='5vortex_ICs',status='old')
do k=1,5
   read(20,*) x(k),y(k),z(k),kap(k)
enddo
close(20)

write(*,*)
write(*,*) ' The vorticity field of each vortex has the form'
write(*,*) '     omega(x,y,z) = a*kappa_j*exp(-a*s), where'
write(*,*) ' s = 2*(1-d)/(1+d) with d = x*x_j+y*y_j+z*z_j.'
write(*,*) ' (Both (x,y,z) and (x_j,y_j,z_j) are unit vectors.)'
write(*,*) ' Enter the inverse-square length a:'
read(*,*) a

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
 !Loop over vortices and accumulate vorticity field:
qc=zero
do k=1,5
  d=x(k)*xg+y(k)*yg+z(k)*zg
  qc=qc+a*kap(k)*exp(-2.0*a*(one-d)/(one+d))
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

end program n5
