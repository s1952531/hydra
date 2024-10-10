!#####################################################################
!           Point Mass dynamics on the two-dimensional plane
!#####################################################################

!          Code adapted from hydra/masses/sphere/src/pms.F by 
!                   David Dritschel @ Wildwood Crest
!              ***Version 1.0 completed 4 August 2013***

!      This code solves:

!      du_j/dt = 2*sum_{i .ne. j} s_i(x_i - x_j)/|x_i - x_j|^2
!      dx_j/dt = u_j

!      where x_i is the position of the ith mass, and s_i is its
!      "strength", here mass/4*pi.

!      The strengths (s_i) of the masses are read in initially
!      from the file "strengths.dat", the point mass positions
!      (x_i) from the file "points.dat", and their velocities (u_i)
!      from "speeds.dat". 
!      Output is appended to energy.dat and written unformatted 
!      and direct access to points.r4 & speeds.r4.

!      These files can be generated, for example, using ring
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
     ! Save the point positions and speeds to points.r4 & speeds.r4:
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
double precision:: eke,ape,ene,fac
real tr4,xr4(n),yr4(n)
integer:: i,k,iniframe

!-----------------------------------------------------------------
 ! Read mass strengths (masses / 4*pi):
open(20,file='strengths.dat',status='old')
do i=1,n
  read(20,*) s(i)
enddo
close(20)

 ! Used in invert:
do i=1,n
  unit(i)=zero
enddo

 ! Read mass positions & speeds:
if (iniframe .eq. 0) then
   ! This is a new run, starting from t = 0:
  open(11,file='points.dat',status='old')
  read(11,*) t
  do i=1,n
    read(11,*) x(i),y(i)
  enddo
  close(11)

  open(22,file='speeds.dat',status='old')
  read(22,*) t
  do i=1,n
    read(22,*) u(i),v(i)
  enddo
  close(22)

   ! Open output files:
  open(11,file='points.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
  open(22,file='speeds.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
  open(44,file='energy.dat',status='unknown')

   ! Write initial data:
  loop=0
  call dump(0)

else
   ! This is a continuation of an old run:
  loop=iniframe+1
   !Open input/output files:
  open(11,file='points.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
  read(11,rec=loop) tr4,xr4,yr4
  t=tr4
  do i=1,n
    x(i)=xr4(i)
    y(i)=yr4(i)
  enddo

  open(22,file='speeds.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
  read(22,rec=loop) tr4,xr4,yr4
  t=tr4
  do i=1,n
    u(i)=xr4(i)
    v(i)=yr4(i)
  enddo

  open(44,file='energy.dat',status='old')
  rewind 44
  do k=0,iniframe
    read(44,*) t,eke,ape,ene
  enddo
endif

return
end subroutine

!=======================================================================

subroutine evolve
! Integrates the equations from time t to time t + dt using a 4th-order
! Runge-Kutta integration scheme.

implicit none
double precision:: xi(n),yi(n),ui(n),vi(n)
double precision:: xf(n),yf(n),uf(n),vf(n)
double precision:: xx,yy,uu,vv,fac
integer:: i

 !-------------------------------------
 ! RK4 predictor step to time t + dt/2:
call invert(0)
 ! invert returns the acceleration (a,b) at each point.

do i=1,n
  xi(i)=x(i)
  yi(i)=y(i)
  x(i)=xi(i)+dt2*u(i)
  y(i)=yi(i)+dt2*v(i)
  xf(i)=xi(i)+dt6*u(i)
  yf(i)=yi(i)+dt6*v(i)

  ui(i)=u(i)
  vi(i)=v(i)
  u(i)=ui(i)+dt2*a(i)
  v(i)=vi(i)+dt2*b(i)
  uf(i)=ui(i)+dt6*a(i)
  vf(i)=vi(i)+dt6*b(i)
enddo

 !-------------------------------------
 ! RK4 corrector step at time t + dt/2:
t=t+dt2
call invert(1)

do i=1,n
  x(i)=xi(i)+dt2*u(i)
  y(i)=yi(i)+dt2*v(i)
  xf(i)=xf(i)+dt3*u(i)
  yf(i)=yf(i)+dt3*v(i)

  u(i)=ui(i)+dt2*a(i)
  v(i)=vi(i)+dt2*b(i)
  uf(i)=uf(i)+dt3*a(i)
  vf(i)=vf(i)+dt3*b(i)
enddo

 !-----------------------------------
 ! RK4 predictor step at time t + dt:
call invert(1)

do i=1,n
  x(i)=xi(i)+dt*u(i)
  y(i)=yi(i)+dt*v(i)
  xf(i)=xf(i)+dt3*u(i)
  yf(i)=yf(i)+dt3*v(i)

  u(i)=ui(i)+dt*a(i)
  v(i)=vi(i)+dt*b(i)
  uf(i)=uf(i)+dt3*a(i)
  vf(i)=vf(i)+dt3*b(i)
enddo

 !-----------------------------------
 ! RK4 corrector step at time t + dt:
t=t+dt2
call invert(1)

do i=1,n
  x(i)=xf(i)+dt6*u(i)
  y(i)=yf(i)+dt6*v(i)

  u(i)=uf(i)+dt6*a(i)
  v(i)=vf(i)+dt6*b(i)
enddo

return 
end subroutine

!=======================================================================

subroutine invert(lev)
! Finds the acceleration (a,b) at each point.  If lev = 0, 
! this routine also computes a new time step.

implicit none
double precision,parameter:: epsa=0.002d0, epsv=0.0005d0
 ! epsa, epsv: for adapting the time step

double precision:: fac,velm,accm,dtacc,trem
integer:: lev,i,j

 !----------------------------------------------
do j=1,n
  a(j)=zero
  b(j)=zero
enddo

do i=1,n
  unit(i-1)=zero
  unit(i)=big
  do j=1,n
    fac=s(i)/((x(i)-x(j))**2+(y(i)-y(j))**2+unit(j))
    a(j)=a(j)+(x(i)-x(j))*fac
    b(j)=b(j)+(y(i)-y(j))*fac
  enddo
enddo
unit(n)=zero

do j=1,n
  a(j)=two*a(j)
  b(j)=two*b(j)
enddo

if (lev .eq. 0) then
   ! Adapt the time step to ensure sqrt{|a|_max}*dt < epsa
   ! and |u|_max*dt < epsv, where epsa & epsv are small 
   ! dimensionless numbers in the parameter statement above.

  velm=zero
  accm=zero
  do j=1,n
    velm=max(velm,u(j)**2+v(j)**2)
    accm=max(accm,a(j)**2+b(j)**2)
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
! Writes data to points.r4 & speeds.r4. computes the kinetic, potential
! and total energy and writes these to energy.dat:

implicit none
double precision:: eke,ape,ene
real tr4,xr4(n),yr4(n)
integer:: nsteps,i,j

loop=loop+1
tr4=real(t)

 !---------------------------------------------
 ! Append point mass positions to points.dat:
do i=1,n
  xr4(i)=real(x(i))
  yr4(i)=real(y(i))
enddo
write(11,rec=loop) tr4,xr4,yr4

 !---------------------------------------------
 ! Append point mass velocities to speeds.dat:
do i=1,n
  xr4(i)=real(u(i))
  yr4(i)=real(v(i))
enddo
write(22,rec=loop) tr4,xr4,yr4

 !----------------------------------------------
 ! Compute the energy and write to standard out:

 ! Kinetic part
eke=zero
do j=1,n
  eke=eke+s(j)*(u(j)**2+v(j)**2)
enddo
eke=half*eke

ape=zero
do j=2,n
  do i=1,j-1
    ape=ape+s(i)*s(j)*log((x(i)-x(j))**2+(y(i)-y(j))**2)
  enddo
enddo

ene=eke+ape

 ! Also write energy to energy.dat:
write(44,'(f12.5,3(1x,f15.9))') t,eke,ape,ene

write(*,'(a,f12.5,a,f15.9,a,i9)') ' t = ',t,' E = ',ene, &
                                   &  ' # time steps = ',nsteps

return 
end subroutine

 !End main program
end program
!=======================================================================
