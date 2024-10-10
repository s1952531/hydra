!#####################################################################
!                  Point Mass dynamics on a Sphere
!#####################################################################

!      Code adapted from pvs.F by David Dritschel @ Wildwood Crest
!              ***Version 1.0 completed 30 July 2013***

!      This code solves:

!      du_j/dt = A_j*x_j + sum_{i .ne. j} s_i(x_i - x_j)/(1 - x_j*x_i)
!      dx_j/dt = u_j

!      where x_i (a 3D vector of unit length) is the position of the
!      ith mass in Cartesian coordinates on the sphere, s_i is its
!      "strength", here mass/4*pi, x_j*x_i means a dot product, and
!      A_j = S - s_j - |u_j|^2, with S = sum_i s_i is the total mass
!      divided by 4*pi.

!      The strengths (s_i) of the masses are read in initially from
!      the file "strengths.asc", the point mass positions (x_i) from
!      the file "coordinates.asc", and their velocities (u_i) from
!      the file "speeds.asc". 

!      Output is appended to these files.  They can be displayed 
!      using the python scripy animate.py.

!      These files can be generated, for example, using diy
!      (compiled in the makefile).

!      The full code consists of the following modules:

!      pms.f90        : This source - main program evolution loop;
!      parameters.f90 : User defined parameters for a simulation;
!      constants.f90  : Fixed constants used throughout;
!      variables.f90  : Global quantities that may change in time.

!      The makefile can be used to compile these modules.
!--------------------------------------------------------------------------
program pms

 ! Import parameters, constants and common arrays:
use constants
use variables

implicit none

 !Local variables:
double precision:: dtsave,tmax
integer:: nsteps,isave

!----------------------------------------------------------
 ! Define fixed arrays and constants and read initial data:
call initialise

 ! Used for saving data quasi-regularly:
dtsave=zero
nsteps=0
isave=0

tmax=t+tsim
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 ! Start the time loop:
do while (t .lt. tmax)

   ! Evolve points from t to t + dt (where dt is adapted):
  call evolve

  nsteps=nsteps+1
  dtsave=dtsave+dt
  if (dtsave .ge. tsave) then
     ! Save the point positions and speeds:
    call dump(nsteps)
    dtsave=dtsave-tsave
    nsteps=0
    isave=1
  endif

enddo

if (isave .eq. 0) call dump(nsteps)

! End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! Close files:
close(11)
close(22)
close(44)

!===============================================================

 ! Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none
double precision:: eke,ape,ene,fac
integer:: i,k,iniframe

!-----------------------------------------------------------------
 ! Read mass strengths (masses / 4*pi):
open(20,file='strengths.asc',status='old')
do i=1,n
  read(20,*) s(i)
enddo
close(20)
 ! Define total strength:
stot=zero
do i=1,n
  stot=stot+s(i)
enddo

 ! Used in invert:
do i=1,n
  unit(i)=one
enddo

 ! Read mass positions & speeds:
open(11,file='coordinates.asc',status='old')
open(22,file='speeds.asc',status='old')

if (iniframe .eq. 0) then
   ! This is a new run, starting from t = 0:
  read(11,*) t
  read(22,*) t
  do i=1,n
    read(11,*) x(i),y(i),z(i)
    read(22,*) u(i),v(i),w(i)
  enddo

   ! Open file to output energy:
  open(44,file='energy.asc',status='replace')

   ! Write initial data:
  call dump(0)

else
   ! This is a continuation of an old run:
  open(44,file='energy.asc',status='old')
  do loop=1,iniframe+1
    read(11,*) t
    read(22,*) t
    read(44,*) t,eke,ape,ene
    do i=1,n
      read(11,*) x(i),y(i),z(i)
      read(22,*) u(i),v(i),w(i)
    enddo
  enddo
endif

 ! Ensure points lie on the sphere:
do i=1,n
  fac=one/sqrt(x(i)**2+y(i)**2+z(i)**2)
  x(i)=fac*x(i)
  y(i)=fac*y(i)
  z(i)=fac*z(i)
enddo

 ! Ensure velocity is tangent to the sphere:
do i=1,n
  fac=x(i)*u(i)+y(i)*v(i)+z(i)*w(i)
  u(i)=u(i)-fac*x(i)
  v(i)=v(i)-fac*y(i)
  w(i)=w(i)-fac*z(i)
enddo

return
end subroutine

!=======================================================================

subroutine evolve
! Integrates the equations from time t to time t + dt using a 4th-order
! Runge-Kutta integration scheme.

implicit none
double precision:: xi(n),yi(n),zi(n),ui(n),vi(n),wi(n)
double precision:: xf(n),yf(n),zf(n),uf(n),vf(n),wf(n)
double precision:: xx,yy,zz,uu,vv,ww,fac
integer:: i

 !-------------------------------------
 ! RK4 predictor step to time t + dt/2:
call invert(0)
 ! invert returns the acceleration (a,b,c) at each point.

do i=1,n
  xi(i)=x(i)
  yi(i)=y(i)
  zi(i)=z(i)
  xx=xi(i)+dt2*u(i)
  yy=yi(i)+dt2*v(i)
  zz=zi(i)+dt2*w(i)
  fac=one/sqrt(xx**2+yy**2+zz**2)
  x(i)=fac*xx
  y(i)=fac*yy
  z(i)=fac*zz
  xf(i)=xi(i)+dt6*u(i)
  yf(i)=yi(i)+dt6*v(i)
  zf(i)=zi(i)+dt6*w(i)

  ui(i)=u(i)
  vi(i)=v(i)
  wi(i)=w(i)
  uu=ui(i)+dt2*a(i)
  vv=vi(i)+dt2*b(i)
  ww=wi(i)+dt2*c(i)
  fac=uu*x(i)+vv*y(i)+ww*z(i)
  u(i)=uu-fac*x(i)
  v(i)=vv-fac*y(i)
  w(i)=ww-fac*z(i)
  uf(i)=ui(i)+dt6*a(i)
  vf(i)=vi(i)+dt6*b(i)
  wf(i)=wi(i)+dt6*c(i)
enddo
 ! The use of fac above ensures all masses remain on the unit sphere
 ! and all velocities are tangent to the sphere.

 !-------------------------------------
 ! RK4 corrector step at time t + dt/2:
t=t+dt2
call invert(1)

do i=1,n
  xx=xi(i)+dt2*u(i)
  yy=yi(i)+dt2*v(i)
  zz=zi(i)+dt2*w(i)
  fac=one/sqrt(xx**2+yy**2+zz**2)
  x(i)=fac*xx
  y(i)=fac*yy
  z(i)=fac*zz
  xf(i)=xf(i)+dt3*u(i)
  yf(i)=yf(i)+dt3*v(i)
  zf(i)=zf(i)+dt3*w(i)

  uu=ui(i)+dt2*a(i)
  vv=vi(i)+dt2*b(i)
  ww=wi(i)+dt2*c(i)
  fac=uu*x(i)+vv*y(i)+ww*z(i)
  u(i)=uu-fac*x(i)
  v(i)=vv-fac*y(i)
  w(i)=ww-fac*z(i)
  uf(i)=uf(i)+dt3*a(i)
  vf(i)=vf(i)+dt3*b(i)
  wf(i)=wf(i)+dt3*c(i)
enddo

 !-----------------------------------
 ! RK4 predictor step at time t + dt:
call invert(1)

do i=1,n
  xx=xi(i)+dt*u(i)
  yy=yi(i)+dt*v(i)
  zz=zi(i)+dt*w(i)
  fac=one/sqrt(xx**2+yy**2+zz**2)
  x(i)=fac*xx
  y(i)=fac*yy
  z(i)=fac*zz
  xf(i)=xf(i)+dt3*u(i)
  yf(i)=yf(i)+dt3*v(i)
  zf(i)=zf(i)+dt3*w(i)

  uu=ui(i)+dt*a(i)
  vv=vi(i)+dt*b(i)
  ww=wi(i)+dt*c(i)
  fac=uu*x(i)+vv*y(i)+ww*z(i)
  u(i)=uu-fac*x(i)
  v(i)=vv-fac*y(i)
  w(i)=ww-fac*z(i)
  uf(i)=uf(i)+dt3*a(i)
  vf(i)=vf(i)+dt3*b(i)
  wf(i)=wf(i)+dt3*c(i)
enddo

 !-----------------------------------
 ! RK4 corrector step at time t + dt:
t=t+dt2
call invert(1)

do i=1,n
  xx=xf(i)+dt6*u(i)
  yy=yf(i)+dt6*v(i)
  zz=zf(i)+dt6*w(i)
  fac=one/sqrt(xx**2+yy**2+zz**2)
  x(i)=fac*xx
  y(i)=fac*yy
  z(i)=fac*zz

  uu=uf(i)+dt6*a(i)
  vv=vf(i)+dt6*b(i)
  ww=wf(i)+dt6*c(i)
  fac=uu*x(i)+vv*y(i)+ww*z(i)
  u(i)=uu-fac*x(i)
  v(i)=vv-fac*y(i)
  w(i)=ww-fac*z(i)
enddo

return 
end subroutine

!=======================================================================

subroutine invert(lev)
! Finds the acceleration (a,b,c) at each point.  If lev = 0, 
! this routine also computes a new time step.

implicit none
double precision,parameter:: epsa=0.002d0, epsv=0.0005d0
 ! epsa, epsv: for adapting the time step

double precision:: fac,velm,accm,dtacc,trem
integer:: lev,i,j

 !----------------------------------------------
do j=1,n
  fac=stot-s(j)-u(j)**2-v(j)**2-w(j)**2
  a(j)=fac*x(j)
  b(j)=fac*y(j)
  c(j)=fac*z(j)
enddo

do i=1,n
  unit(i-1)=one
  unit(i)=big
  do j=1,n
    fac=s(i)/(unit(j)-x(i)*x(j)-y(i)*y(j)-z(i)*z(j))
    a(j)=a(j)+(x(i)-x(j))*fac
    b(j)=b(j)+(y(i)-y(j))*fac
    c(j)=c(j)+(z(i)-z(j))*fac
  enddo
enddo
unit(n)=one

if (lev .eq. 0) then
   ! Adapt the time step to ensure sqrt{|a|_max}*dt < epsa
   ! and |u|_max*dt < epsv, where epsa & epsv are small 
   ! dimensionless numbers in the parameter statement above.

  velm=zero
  accm=zero
  do j=1,n
    velm=max(velm,u(j)**2+v(j)**2+w(j)**2)
    accm=max(accm,a(j)**2+b(j)**2+c(j)**2)
  enddo
  velm=sqrt(velm)
  accm=sqrt(accm)

   ! time step to maintain accuracy:
  dtacc=min(epsa/sqrt(accm),epsv/velm)

   ! See how close we are to the end of the simulation:
  trem=tmax-t

   !Choose the smallest time step of the above:
  dt=min(dtacc,trem)

   ! Define time interval constants:
  dt2=dt/two
  dt3=dt/three
  dt6=dt/six

endif

return 
end subroutine

!=======================================================================

subroutine dump(nsteps)
! Writes data to coordinates.asc and speeds.asc. Computes the kinetic, 
! potential and total energy and writes these to energy.asc:

implicit none
double precision:: eke,ape,ene
integer:: nsteps,i,j

loop=loop+1

 !-----------------------------------------------------------
 ! Append positions & speeds to coordinates.asc & speeds.asc:
write(11,'(f12.5)') t
write(22,'(f12.5)') t
do i=1,n
  write(11,'(3(1x,f15.12))') x(i),y(i),z(i)
  write(22,'(3(1x,f15.12))') u(i),v(i),w(i)
enddo

 !--------------------------------------------
 ! Compute the energy and write to energy.asc:

 ! Kinetic part
eke=zero
do j=1,n
  eke=eke+s(j)*(u(j)**2+v(j)**2+w(j)**2)
enddo
eke=half*eke

ape=zero
do j=2,n
  do i=1,j-1
    ape=ape+s(i)*s(j)*log(one-x(i)*x(j)-y(i)*y(j)-z(i)*z(j))
  enddo
enddo

ene=eke+ape

 ! Write energy to energy.asc:
write(44,'(f12.5,3(1x,f15.9))') t,eke,ape,ene

write(*,'(a,f12.5,a,f15.9,a,i9)') ' t = ',t,' E = ',ene, &
                                   &  ' # time steps = ',nsteps

return 
end subroutine

 !End main program
end program
!=======================================================================
