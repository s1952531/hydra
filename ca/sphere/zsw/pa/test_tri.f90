program testtri

implicit none

integer, parameter:: n=40, n1=n-1, n2=n-2
double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0
double precision, parameter:: f12=1.d0/2.d0, f13=1.d0/3.d0
double precision:: r(n1),a0(n1),ap(n1),etd(n1),htd(n1)
double precision:: x(0:n),x3(0:n),y(0:n),ye(0:n),dydxe(0:n)
double precision:: dx(n),ybar(n),rhs(n1)
double precision:: alp(0:n),bet(0:n1),dydx(0:n)
double precision:: eps,x0,fni,fac
integer:: j,m

!-----------------------------------------------------------------
write(*,*) ' We choose the discrete x values according to'
write(*,*) '      x(j) = j/n + eps*sin(j*m*pi/n).'
write(*,*)
write(*,*) ' Enter eps and m:'
read(*,*) eps,m

fni=1.d0/dble(n)
fac=pi*dble(m)
do j=0,n
  x0=fni*dble(j)
  x(j)=x0+eps*sin(x0*fac)
enddo

x3=x**3                 !x^3
ye=x**2*(5.d0-2.d0*x3)  !y
dydxe=10.d0*x*(1.d0-x3) !dy/dx
y=f13*x3*(5.d0-x3)      !holds int_0^x{y} temporarily

!-----------------------------------------------------------------
do j=1,n
  dx(j)=x(j)-x(j-1)
  ybar(j)=(y(j)-y(j-1))/dx(j) !Average y over (x_{j-1},x_j)
enddo

do j=1,n1
  r(j)=dx(j)/dx(j+1)
  a0(j)=2.d0*(1.d0+r(j))
enddo

do j=2,n1
  ap(j-1)=r(j)
enddo

htd(1)=1.d0/a0(1)
etd(1)=-ap(1)*htd(1)

do j=2,n2
  htd(j)=1.d0/(a0(j)+etd(j-1))
  etd(j)=-ap(j)*htd(j)
enddo
htd(n1)=1.d0/(a0(n1)+etd(n2))
!-----------------------------------------------------------------

do j=1,n1
  rhs(j)=6.d0*(ybar(j+1)-ybar(j))
enddo

alp(1)=rhs(1)*htd(1)

do j=2,n1
  alp(j)=(rhs(j)-alp(j-1))*htd(j)
enddo

do j=n2,1,-1
  alp(j)=etd(j)*alp(j+1)+alp(j)
enddo

!-----------------------------------------------------------------
alp(0)=0.d0
alp(n)=0.d0

do j=0,n2
  bet(j)=f12*(r(j+1)*alp(j+1)-alp(j))
enddo
bet(n1)=-f12*alp(n1)

dydx(0)=0.d0
do j=1,n1
  dydx(j)=alp(j)/dx(j+1)
enddo
dydx(n)=0.d0

do j=0,n1
  y(j)=ybar(j+1)-f12*alp(j)-f13*bet(j)
enddo
y(n)=y(n1)+alp(n1)+bet(n1)

!-----------------------------------------------------------------
write(*,*)

open(77,file='dydx.asc',status='replace')
do j=0,n
  write(77,*) x(j),dydx(j),dydxe(j)
enddo
close(77)
write(*,*) ' R.m.s. error in dy/dx = ',sqrt(fni*sum((dydx-dydxe)**2))

open(77,file='y.asc',status='replace')
do j=0,n
  write(77,*) x(j),y(j),ye(j)
enddo
close(77)
write(*,*) ' R.m.s. error in   y   = ',sqrt(fni*sum((y-ye)**2))

end program testtri
