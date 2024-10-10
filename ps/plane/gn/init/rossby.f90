program rossby
! Initialises a Rossby wave

use constants

implicit double precision(a-h,o-z)

double precision:: qq(ng,ng)

write(*,*) ' The PV is given by q = beta[y+a*sin(n*x+m*y)]'
write(*,*) ' Enter a:'
read(*,*) a
write(*,*) ' Enter n & m (integer):'
read(*,*) n,m
dn=dble(n)
dm=dble(m)

do ix=1,ng
  x=gl*dble(ix-1)-pi
  do iy=1,ng
    y=gl*dble(iy-1)-pi
    qq(iy,ix)=beta*(y+a*sin(dn*x+dm*y))
  enddo
enddo

open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program
