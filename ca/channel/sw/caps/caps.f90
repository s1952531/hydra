!#########################################################################
!               The Single-Layer Singly-Periodic Channel 
!       Shallow-Water Combined Lagrangian Advection Method (CLAM)
!#########################################################################

!       Code developed in July 2023 by D G Dritschel @ St Andrews.

!       This code simulates the unforced rotating Shallow-Water Equations
!       in variables (q,v',b), where q is the potential vorticity, v' is
!       the imbalanced y velocity component and b = Dv/Dt is the y
!       acceleration component.  Additionally, we must evolve the 
!       boundary horizontal velocities u^{+/-} for consistency.

!       Note, q is arbitrary except assumed positive. On the other hand,
!       but v & b can have either sign but must vanish at each boundary.

!       Below, qq is PV, (uu,vv) is the velocity field, bs is b or Dv/Dt,
!       hh is the dimensionless height anomaly, and zz = qq*(1+hh) - cof
!       is the relative vorticity, where cof is the (constant) Coriolis
!       frequency.  The boundary velocities u^{+} and u^{-} are uup & uum
!       at y = y_max and y_min, respectively.

!       Contour advection and generation are done internally now.
!       For details of the method, see Dritschel and Fontane,
!       J. Comput. Phys. 229, pp. 5408--5417 (2010).

!       The full algorithm consists of the following modules:

!       caps.f90      : This source - main program loop, repeats
!                       successive calls to evolve fields and recontour;
!       parameters.f90: User defined parameters for a simulation;
!       constants.f90 : Fixed constants used throughout the other modules;
!       common.f90    : Common data preserved throughout simulation 
!                       (through the recontouring--evolution cycle);
!       spectral.f90  : Fourier transform common storage and routines;
!       contours.f90  : Contour advection common storage and routines;
!       congen.f90    : Source code for contour-to-grid conversion;
!       generic.f90   : Service routines (averaging, norms, etc);
!       evolution.f90 : Main time evolution module - advects gridded 
!                       fields using a PS method along with contours.
!----------------------------------------------------------------------------
program caps

 !Import all modules:
use common

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Obtain new PV contours:
  call recont

   !Advect PV and other fields until next recontouring or end:
  call evolve

enddo

 !End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

call finalise

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit none

!----------------------------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral

 !Initialise constants and arrays for contour advection:
call init_contours

!------------------------------------------------------------------------
 !Read in y velocity and y acceleration components:
open(11,file='vv_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,vv
close(11)

open(11,file='bb_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,bs
close(11)

 !Ensure that vv = 0 on the boundaries:
if (maxval(abs(vv(0,:))) > zero .or. maxval(abs(vv(ny,:))) > zero) then
  write(*,*) ' ***The y velocity is not zero on the boundaries! Stopping!'
  stop
endif

 !Ensure that bs = 0 on the boundaries:
if (maxval(abs(bs(0,:))) > zero .or. maxval(abs(bs(ny,:))) > zero) then
  write(*,*) ' ***The y acceleration is not zero on the boundaries! Stopping!'
  stop
endif

!Convert bs to semi-spectral space:
call forfft(nyp1,nx,bs,xtrig,xfactors)
 !Apply de-aliasing filter:
call dealias(bs,1)

!----------------------------------------------------------------------
 !Read in gridded PV as qq and convert to semi-spectral space as qs:
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(11,rec=1) t,qq
close(11)

 !Copy qq to qs:
qs=qq
 !Convert qs to semi-spectral space:
call forfft(nyp1,nx,qs,xtrig,xfactors)
 !Apply de-aliasing filter:
call dealias(qs,1)

!----------------------------------------------------------------------
 !Read in zonal velocity at each boundary:
open(11,file='uum_init.r8',form='unformatted', &
      access='direct',status='old',recl=nxbytes)
read(11,rec=1) uum
close(11)

open(11,file='uup_init.r8',form='unformatted', &
      access='direct',status='old',recl=nxbytes)
read(11,rec=1) uup
close(11)

 !Calculate mean values of boundary zonal velocity (time invariant):
uumbar=sum(uum)*dnxi
uupbar=sum(uup)*dnxi
 !dnxi = 1/nx here.

 !Constant used to correct mean PV:
qoff=cof-(uupbar-uumbar)/elly
 !Note: the mean relative vorticity is -(uupbar-uumbar)/elly.

 !Convert to semi-spectral space:
call forfft(1,nx,uum,xtrig,xfactors)
call forfft(1,nx,uup,xtrig,xfactors)

 !Apply de-aliasing filter:
uum=filt*uum
uup=filt*uup

!----------------------------------------------------------------------
 !Obtain initial dimensionless height anomaly (hh), x velocity (uu),
 !relative vorticity (zz) and adjusted PV (qq) consistent with zero
 !domain-averaged zz (note: need to start with hh = 0 or a guess):
hh=zero
call invert(qq,vv,bs,uum,uup,qoff,hh,uu,zz,qt,vb)
 !Note: qq, hh, uu, zz, qt & vb are in physical space after this call.
 !      vv is already available in physical space.

!----------------------------------------------------------------------
 !Set qr = qq to build initial PV contours in recontour below:
qr=qq

 !Determine PV contour interval:
dq=(maxval(qq)-minval(qq))/dble(ncontq)

!Initially there are no contours (they are built from the gridded PV):
nq=0
nptq=0
t=zero

!----------------------------------------------------------------------
 !Open all plain text diagnostic files:
open(14,file='contours/complexity.asc',status='replace')
open(15,file='evolution/ecomp.asc',status='replace')
open(16,file='evolution/ro-fr-hm.asc',status='replace')

 !Open files for coarse grid saves:
open(32,file='evolution/qq.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(33,file='evolution/bb.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(34,file='evolution/hh.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(35,file='evolution/zz.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(36,file='evolution/uu.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)
open(37,file='evolution/vv.r8',form='unformatted',access='direct', &
                             status='replace',recl=nbytes)

 !Open files for contour writes:
open(80,file='contours/qqsynopsis.asc',status='replace')
open(83,file='contours/qqresi.r8',form='unformatted',access='direct', &
                                status='replace',recl=nbytes)

 !Define number of time steps between grid and contour saves:
ngsave=nint(tgsave/dt)
ncsave=nint(tcsave/dt)
 !*** WARNING: tgsave and tcsave should be an integer multiple of dt

return
end subroutine initialise

!=======================================================================

subroutine evolve

use evolution

implicit none

 !Advect PV and other fields until next recontouring or end:
write(*,*) 'Evolving contours and fields ...'
call advect

return 
end subroutine evolve

!=======================================================================

subroutine recont

use congen

implicit none

 !Obtain new PV contours:
write(*,*) 'Recontouring PV ...'
call recontour(qr)
write(*,'(a,i8,a,i9)') '   nq = ',nq,'   nptq = ',nptq

return 
end subroutine recont

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(14)
close(15)
close(16)

close(32)
close(33)
close(34)
close(35)
close(36)
close(37)

close(80)
close(83)

return
end subroutine finalise

 !End main program
end program caps
!=======================================================================
