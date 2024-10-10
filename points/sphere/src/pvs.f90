!#####################################################################
!                 Point Vortex dynamics on a Sphere
!#####################################################################

!        Code adapted from pvs.F by David Dritschel @ St Andrews
!              ***Version 1.0 completed 29 October 2012***

!        This code solves:

!            dx_j/dt = sum_{i .ne. j} s_i(x_j X x_i)/(1 - x_j*x_i)

!        where x_i (a 3D vector of unit length) is the position of 
!        the ith vortex in Cartesian coordinates on the sphere, s_i
!        is its "strength", here circulation/4*pi, X means a cross
!        product, and * means a dot product.

!        The strengths (s_i) of the vortices are read in initially
!        from the file "strengths.dat", the point vortex positions
!        (x_i) from the file "points.dat", and the initial energy
!        from the file "energy.dat".  The latter two are appended
!        to over the course of the simulation. 

!        These files can be generated, for example, using trihex
!        (compiled in the makefile).

!        The full code consists of the following modules:

!        pvs.f90        : This source - main program evolution loop;
!        parameters.f90 : User defined parameters for a simulation;
!        constants.f90  : Fixed constants used throughout;
!        variables.f90  : Global quantities that may change in time.

!        The makefile can be used to compile these modules.
!--------------------------------------------------------------------------
program pvs

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
     ! Save the point vortices to points.dat:
    call dump(nsteps)
    dtsave=dtsave-tsave
    nsteps=0
    isave=1
  endif

enddo

if (isave .eq. 0) call dump(nsteps)

! End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

!===============================================================

 ! Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none
double precision:: ene
integer:: i,k

!-----------------------------------------------------------------
 ! Read vortex strengths:
open(20,file='strengths.dat',status='old')
do i=1,n
  read(20,*) s(i)
enddo
close(20)

 ! Read vortex positions:
open(11,file='points.dat',status='old')
do k=0,iniframe
  read(11,*) t
  do i=1,n
    read(11,*) x(i),y(i),z(i)
  enddo
enddo

 ! Open file containing initial energy:
open(44,file='energy.dat',status='old')
rewind 44
do k=0,iniframe
  read(44,*) t,ene
enddo

 ! Used in invert:
do i=1,n
  unit(i)=one
enddo

return
end subroutine

!=======================================================================

subroutine evolve
! Integrates the contours from time t to time t + dt using a 4th-order
! 4th-order Runge-Kutta integration scheme.

implicit none
double precision:: xi(n),yi(n),zi(n)
double precision:: xf(n),yf(n),zf(n)
double precision:: xx,yy,zz,fac
integer:: i

 !-------------------------------------
 ! RK4 predictor step to time t + dt/2:
call invert(0)
 ! invert returns the velocity at each point.

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
enddo
 ! The use of fac above ensures all vortices remain on the unit sphere.

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
enddo

return 
end subroutine

!=======================================================================

subroutine invert(lev)
! Finds the velocity at each point.  If lev = 0, this routine
! also computes a new time step.

implicit none
double precision,parameter:: epsa=0.02d0,epsv=0.005d0
 ! epsa, epsv: for adapting the time step

double precision:: fac,velom,angvm,dtacc,trem
double precision:: dx,dy,dz,du,dv,dw
integer:: lev,i,j

 !----------------------------------------------
do j=1,n
  u(j)=zero
  v(j)=zero
  w(j)=zero
enddo

do i=1,n
  unit(i-1)=one
  unit(i)=big
  do j=1,n
    fac=s(i)/(unit(j)-x(i)*x(j)-y(i)*y(j)-z(i)*z(j))
    u(j)=u(j)+(y(i)*z(j)-y(j)*z(i))*fac
    v(j)=v(j)+(z(i)*x(j)-z(j)*x(i))*fac
    w(j)=w(j)+(x(i)*y(j)-x(j)*y(i))*fac
  enddo
enddo
unit(n)=one

if (lev .eq. 0) then
   ! Adapt the time step to ensure max angular rotation of a
   ! pair of vortices multiplied by the time step dt < epsa,
   ! and the maximum velocity (on a sphere of radius 1) multiplied
   ! by the time step dt < epsv, where epsa & epsv are small 
   ! dimensionless numbers in the parameter statement above.

  angvm=zero
  do i=1,n
    unit(i-1)=one
    unit(i)=big
    do j=1,n
      dx=x(i)-x(j)
      dy=y(i)-y(j)
      dz=z(i)-z(j)
      du=u(i)-u(j)
      dv=v(i)-v(j)
      dw=w(i)-w(j)
      angvm=max(angvm,((dy*dw-dz*dv)**2+(dz*du-dx*dw)**2+(dx*dv-dy*du)**2)/ &
                   &  (unit(j)-x(i)*x(j)-y(i)*y(j)-z(i)*z(j))**2)
    enddo
  enddo
  unit(n)=one
  angvm=half*sqrt(angvm)
   ! angvm: the maximum angular rotation of any pair of vortices.

   ! Also check the max velocity (= angular velocity for a sphere
   ! of unit radius):
  velom=zero
  do i=1,n
    velom=max(velom,u(i)**2+v(i)**2+w(i)**2)
  enddo
  velom=sqrt(velom)

   ! time step to maintain accuracy:
  dtacc=min(epsa/angvm,epsv/velom)

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
! Finds the velocity at each point.  If lev = 0, this routine
! also computes a new time step.

implicit none
double precision:: ene
integer:: nsteps,i,j

 !---------------------------------------------
 ! Append point vortex positions to points.dat:
write(11,'(f12.5)') t
do i=1,n
  write(11,'(3(1x,f12.9))') x(i),y(i),z(i)
enddo

 !----------------------------------------------
 ! Compute the energy and write to standard out:
ene=zero
do i=1,n
  unit(i-1)=one
  unit(i)=two
  do j=1,n
    ene=ene-s(i)*s(j)*log(unit(j)-x(i)*x(j)-y(i)*y(j)-z(i)*z(j))
  enddo
enddo
unit(n)=one

 ! Also write energy to energy.dat:
write(44,'(1x,f12.5,1x,1p,e14.7)') t,ene

write(*,'(a,f12.5,a,f17.8,a,i9)') ' t = ',t,'   Energy = ',ene, &
                                   &  '   # time steps = ',nsteps

return 
end subroutine


 !End main program
end program
!=======================================================================
