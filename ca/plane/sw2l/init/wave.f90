program wavy
! Initialises a wavy height field and a state of rest for testing
! energy conservation.

use constants

implicit double precision(a-h,o-z)

double precision:: h1(ng,ng),q1(ng,ng),d1(ng,ng),g1(ng,ng)
double precision:: h2(ng,ng),q2(ng,ng),d2(ng,ng),g2(ng,ng)

write(*,*) ' We take h_tilde_1 = A_1*cos(k_1*x)*cos(l_1*y).'
write(*,*) ' Enter A_1, k_1 & l_1:'
read(*,*) amp1,kx1,ky1
dkx1=dble(kx1)
dky1=dble(ky1)

write(*,*) ' We take h_tilde_2 = A_2*cos(k_2*x)*cos(l_2*y).'
write(*,*) ' Enter A_2, k_2 & l_2:'
read(*,*) amp2,kx2,ky2
dkx2=dble(kx2)
dky2=dble(ky2)

do ix=1,ng
  x=gl*dble(ix-1)-pi
  do iy=1,ng
    y=gl*dble(iy-1)-pi
    h1(iy,ix)=amp1*cos(dkx1*x)*cos(dky1*y)
    h2(iy,ix)=amp2*cos(dkx2*x)*cos(dky2*y)
  enddo
enddo
fac1=csq1*(dkx1**2+dky1**2)
fac2=csq2*(dkx2**2+dky2**2)
g1=fac1*h1+alpha*fac2*h2
g2=fac1*h1+fac2*h2
q1=cof/(one+h1)
q1avg=dsumi*sum(q1)
q1=q1-q1avg
q2=cof/(one+h2)
q2avg=dsumi*sum(q2)
q2=q2-q2avg
d1=zero
d2=zero

open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,q1
write(11,rec=2) zero,q2
close(11)

open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,d1
write(11,rec=2) zero,d2
close(11)

open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,g1
write(11,rec=2) zero,g2
close(11)

end program wavy
