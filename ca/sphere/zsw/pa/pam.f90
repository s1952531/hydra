!#########################################################################
!          The Zonally-Averaged Spherical Single-Layer Shallow-Water
!                        Point Advection Method (PAM)
!#########################################################################

!       Code developed in October 2020 by D G Dritschel @ St Andrews.

!       This code simulates the unforced Zonally-Averaged Shallow-Water
!       Equations (SWE) in variables (q,h,v), where q is the potential
!       vorticity, h is the dimensionless height (scaled on the constant
!       mean height), and v is the meridional velocity.

!       The code conserves PV by preserving q on each moving latitude
!       phi_j, and conserves mass m_j in each phi interval.  The latter
!       determines the mean height h in each phi interval at any time.

!       The system of equations is

!       dphi_j/dt = v(phi_j,t)     --- PV conservation
!       hbar_j = m_j/dz_j          --- mass conservation (dz = d(sin(phi)))
!       dv_j/dt = -c^2*(dh/dphi)_j - u_j*(f_j+u_j*tan(phi_j))
!       u_j=(1/r_j)*int_{-pi/2}^{phi_j} r*(f - q*h)*dphi

!       where hbar is the mean height in latitude interval j, c^2 = g*H
!       is the squared short-scale gravity wave speed, f = 2*Omega*z
!       is the Coriolis frequency, z = sin(phi), r = cos(phi), q is the
!       PV, u is the zonal velocity component and v is the meridional
!       velocity component.

!       Note that this code exactly conserves angular momentum by
!       taking the PV to be piecewise constant in each phi interval,
!       q = q_j, then using this to determine the absolute zonal
!       velocity from
!                u_{aj} = u_j + omega*r_j = U_j/r_j
!       on each moving latitude phi_j(t).  Here U_j is the angular
!       momentum, determined from -sum_{i=1}^{j} q_j*m_j (constant).
!       This allows one to find the acceleration from
!       dv_j/dt = [(omega*r_j)^2 - (U_j/r_j)^2]*z_j/r_j - c^2*(dh/dphi)_j.

!       The derivative (dh/dphi)_j at each latitude phi = phi_j is found
!       by a novel quadratic spline interpolation which ensures continuity
!       of h and dh/dphi at each phi_j, and satisfies int{h*dz} = m_j
!       over each phi interval.

!       The full algorithm consists of the following modules:

!       pam.f90       : This source - main time evolution module;
!       parameters.f90: User defined parameters for a simulation;
!       constants.f90 : Fixed constants used throughout the other modules;
!----------------------------------------------------------------------------
program pam

use constants

implicit none

 !Latitudes plus cos & sin latitudes (r & z):
double precision:: phi(0:n),r(0:n),z(0:n)

 !Mass and mean PV times mass in each interval (phi_{i-1},phi_i) & -sum:
double precision:: mass(n),qm(n),uu(n1)

 !Difference in phi, r and z in each interval:
double precision:: dphi(n),dr(n),dz(n)

 !Height h at each latitude and dh/dphi:
double precision:: h(0:n),dhdp(0:n)

 !Mean height over each latitude interval:
double precision:: hbar(n)

 !Zonal and meridional velocity components:
double precision:: u(0:n),v(0:n)

 !Initial latitudes used for computing Lagrangian displacement:
double precision:: phie(0:n)

double precision:: eini,t

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

 !Evolve phi and v until end:
call evolve

 !Close files:
call finalise(0)

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none

! Local variables:
double precision:: tmp
integer:: j

!----------------------------------------------------------------------
 !Read in initial latitudes phi, mass, PV times mass & v:
open(11,file='init.asc',status='old')
phi(0)=-hpi !This value at the South Pole is not read in below
v(0)=zero   !v = 0 at both poles
do j=1,n
  read(11,*) phi(j),mass(j),qm(j),v(j)
enddo
close(11)
 !Note: mass and qm never change over the course of the simulation

 !Initialise time:
t=zero

 !Ensure that sum(qm) = zero:
tmp=sum(qm)/dble(n)
qm=qm-tmp
 !This ensures the global average vorticity is zero.

 !Accumulate partial sums of -qm for use in computing u:
uu(1)=-qm(1)
do j=2,n1
  uu(j)=uu(j-1)-qm(j)
enddo
 !Note uu(n) = zero; this is not needed.

 !Set end point values of u to zero (where they remain):
u(0)=zero
u(n)=zero

!--------------------------------------
 !Open all plain text diagnostic files:
open(15,file='ecomp.asc',status='replace')
open(16,file='cfl.asc',status='replace')
open(17,file='froude.asc',status='replace')
open(18,file='hminmax.asc',status='replace')
open(19,file='dvrms.asc',status='replace')

 !Open files for coarse grid saves:
open(21,file='phi.r8',form='unformatted', access='direct', &
                  status='replace',recl=nbytes)
open(22,file='h.r8',form='unformatted', access='direct', &
                  status='replace',recl=nbytes)
open(23,file='u.r8',form='unformatted', access='direct', &
                  status='replace',recl=nbytes)
open(24,file='v.r8',form='unformatted', access='direct', &
                  status='replace',recl=nbytes)

return
end subroutine initialise

!=======================================================================

subroutine evolve
 !Evolves fields and points
  
implicit none

integer:: step,loop

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do loop=1,nloops
   !Save the fields every nsteps time steps:
  call savedata(loop)

  do step=1,nsteps
     !Evolve flow from time t to t + dt:
    call advance
  enddo

enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Save final data:
call savedata(nloops+1)

return
end subroutine evolve

!=======================================================================

subroutine advance

! Advances phi & v from time t to t+dt.

! Uses the 4th-order Runge-Kutta method.

implicit none

 !Local variables:
double precision:: phi0(1:n1),v0(1:n1)
double precision:: phif(1:n1),vf(1:n1)
double precision:: dvdt(1:n1)

!------------------------------------------------------------------
 !Calculate dv/dt:
call accel(dvdt)

 !Compute and save various diagnostics each time step:
call diagnose

phi0=phi(1:n1)
phi(1:n1)=phi0+dt2*v(1:n1)
phif=phi0+dt6*v(1:n1)

v0=v(1:n1)
v(1:n1)=v0+dt2*dvdt
vf=v0+dt6*dvdt

 !RK4 corrector step at time t0 + dt/2:
t=t+dt2

call accel(dvdt)

phi(1:n1)=phi0+dt2*v(1:n1)
phif=phif+dt3*v(1:n1)

v(1:n1)=v0+dt2*dvdt
vf=vf+dt3*dvdt

 !RK4 predictor step at time t0 + dt:
call accel(dvdt)

phi(1:n1)=phi0+dt*v(1:n1)
phif=phif+dt3*v(1:n1)

v(1:n1)=v0+dt*dvdt
vf=vf+dt3*dvdt

 !RK4 corrector step at time t0 + dt:
t=t+dt2
call accel(dvdt)

phi(1:n1)=phif+dt6*v(1:n1)

v(1:n1)=vf+dt6*dvdt

return
end subroutine advance

!=======================================================================

subroutine accel(dvdt)

! Computes dv/dt at current time

implicit none

 !Passed variable:
double precision:: dvdt(n1)

!---------------------------------------------------------------
! Find dhdp = dh/dphi (among other things):
call convert

! Compute dv/dt:
dvdt=((omega*r(1:n1))**2-(uu/r(1:n1))**2)*z(1:n1)/r(1:n1)-csq*dhdp(1:n1)

return
end subroutine accel

!=======================================================================

subroutine convert

! Finds h and dh/dphi from the mass in each phi interval

implicit none

! Local variables:
double precision:: rho(n1),am(n1),a0(n1),ap(n1),etd(n1),htd(n1)
double precision:: dphic(n),dzc(n)
double precision:: vtd(n),wtd(n),rhs(n1)
double precision:: alp(0:n),bet(0:n1)
integer:: j

!-----------------------------------------------------------------
! Define cosine and sine of latitude (r & z):
r=cos(phi)
z=sin(phi)

! Compute differences in phi, r and z:
do j=1,n
  dphi(j)=phi(j)-phi(j-1)
  dr(j)=r(j)-r(j-1)
  dz(j)=z(j)-z(j-1)
enddo
dphic=one/dphi
dzc=one/dz

! Average value of h over interval (phi_{j-1},phi_j):
hbar=mass*dzc
! Integral of xi*dz where xi = (phi-phi_{j-1})/dphi:
vtd=(z(1:n)+dr*dphic)/dz
! Integral of xi^2*dz:
wtd=(z(1:n)+two*dphic*(r(1:n)-dz*dphic))/dz

do j=1,n1
  rho(j)=dphi(j)/dphi(j+1)
  a0(j)=two*vtd(j+1)+(one-wtd(j))*rho(j)-wtd(j+1)
enddo

do j=2,n1
  am(j)=two*(one-vtd(j))+wtd(j)-one
  ap(j-1)=wtd(j)*rho(j)
enddo

htd(1)=1.d0/a0(1)
etd(1)=-ap(1)*htd(1)

do j=2,n2
  htd(j)=1.d0/(a0(j)+am(j)*etd(j-1))
  etd(j)=-ap(j)*htd(j)
enddo
htd(n1)=1.d0/(a0(n1)+am(n1)*etd(n2))

do j=1,n1
  rhs(j)=two*(hbar(j+1)-hbar(j))
enddo

alp(1)=rhs(1)*htd(1)

do j=2,n1
  alp(j)=(rhs(j)-am(j)*alp(j-1))*htd(j)
enddo

do j=n2,1,-1
  alp(j)=etd(j)*alp(j+1)+alp(j)
enddo

!-----------------------------------------------------------------
alp(0)=0.d0
alp(n)=0.d0

do j=0,n2
  bet(j)=f12*(rho(j+1)*alp(j+1)-alp(j))
enddo
bet(n1)=-f12*alp(n1)

dhdp(0)=0.d0
do j=1,n1
  dhdp(j)=alp(j)/dphi(j+1)
enddo
dhdp(n)=0.d0

do j=0,n1
  h(j)=hbar(j+1)-vtd(j+1)*alp(j)-wtd(j+1)*bet(j)
enddo
h(n)=h(n1)+alp(n1)+bet(n1)

return
end subroutine convert

!=======================================================================

subroutine diagnose

! Diagnoses energy components and other quantities every time step.

implicit none

 !Local variables:
double precision:: wk(0:n)
double precision:: ekin,epot,etot
double precision:: cfl,romax,frmax,hmax,hmin
double precision:: prms,vrms

!-----------------------------------------------------------------
! Find hbar and h:
call convert

! Compute zonal velocity, u:
u(1:n1)=uu/r(1:n1)-omega*r(1:n1)

!-----------------------------------------------------------------
 !Compute energy components:
wk=u**2+v**2
ekin=f14*sum(mass*(wk(0:n1)+wk(1:n)))
epot=f12*csq*sum(mass*(hbar-one))
etot=ekin+epot

 !Write energies to ecomp.asc:
write(15,'(f12.5,3(1x,f16.14))') t,ekin,epot,etot
write(*,'(a,f12.5,a,f14.12)') ' t = ',t,'   E_tot = ',etot

 !Maximum Froude number:
wk=abs(v)/(cgw*h)
frmax=maxval(wk)
 !CFL based on gravity wave speed * time step / phi interval:
cfl=dt*cgw*maxval(sqrt(hbar)/dphi)
 !Maximum and minimum hbar:
hmax=maxval(hbar)
hmin=minval(hbar)

 !Rms displacement and meridional velocity:
if (t .gt. zero) then
  prms=sqrt(sum((phi(1:n1)-phie(1:n1))**2)/dble(n1))
else
   !Save initial particle locations to measure displacement for t > 0:
  phie=phi
  prms=zero
endif
vrms=sqrt(sum(v(1:n1)**2)/dble(n1))

 !Record diagnostics:
write(16,'(1x,f12.5,1x,f12.9)') t,cfl
write(17,'(1x,f12.5,1x,f12.9)') t,frmax
write(18,'(1x,f12.5,2(1x,f12.9))') t,hmin,hmax
write(19,'(1x,f12.5,2(1x,f12.9))') t,prms,vrms

if (t .eq. zero) then
   !Store initial value of energy:
  eini=etot
else
  if (abs(etot-eini) .gt. 0.05d0*eini) then
    write(*,*) &
         ' Energy error exceeds 5% of the initial energy! *** Stopping ***'
    call finalise(1)
    stop
  endif
endif

return
end subroutine diagnose

!=======================================================================

subroutine savedata(igrids)

! Saves mean phi, h, u & v at the desired save time

implicit none

 !Passed variable:
integer:: igrids

!-----------------------------------------------------------------
! Find hbar and h:
call convert

! Compute zonal velocity for output only:
r=cos(phi)
u(1:n1)=uu/r(1:n1)-omega*r(1:n1)

! Write data:
write(21,rec=igrids) t,f12*(phi(0:n1)+phi(1:n))
write(22,rec=igrids) t,hbar-one
write(23,rec=igrids) t,f12*(u(0:n1)+u(1:n))
write(24,rec=igrids) t,f12*(v(0:n1)+v(1:n))

return
end subroutine savedata

!=======================================================================

subroutine finalise(iopt)

implicit none

 !Passed parameter:
integer:: iopt  

if (iopt .eq. 0) write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(15)
close(16)
close(17)
close(18)
close(19)
close(21)
close(22)
close(23)
close(24)

return
end subroutine finalise

 !End main program
end program pam
!=======================================================================
