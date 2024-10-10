!#########################################################################
!        The Combined Lagrangian Advection Method for a two-layer
!          quasi-geostrophic flow in periodic channel geometry
!#########################################################################

!        Code written by David Dritschel & Stuart King @ St Andrews
!              ***Version 1.0 completed 22 October 2013***
!          ***Updated for general stratification in late 2014***
!                ***Version 2.0 completed 4 March 2015***

!        Principally adapted from the stratified codes in 
!        hydra/ca/strat/sper; this code uses a semi-Lagrangian 
!        method for the gridded fields and a semi-spectral/4th-order 
!        difference method for PV inversion.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!        This code solves
!            Dq_1/Dt = therm_1*d'_1 - r*zeta_1  *** lower layer
!            Dq_2/Dt = F + therm_2*(d'2 - d'1)  *** upper layer
!            q_1 = beta*y + Lap(psi_1) - d_1/h_1 - f_0*H_b/H_1
!            q_2 = beta*y + Lap(psi_2) - (d_2-d_1)/h_2
!        in the periodic channel xmin < x < xmax ; ymin < y < ymax
!        (with free slip boundary conditions in y, periodic in x)
!        where q_j is the QG PV in layer j, psi_j is the streamfunction
!        in layer j, d_1 is the displacement (times f_0/(H_1+H_2)) of
!        the interface between the lower layer (1) and the upper layer (2),
!        d_2 is the displacement (times f_0/(H_1+H_2)) of the upper free
!        surface (zero when there is a barotropic mode), d'_j = d_j - deq_j
!        where deq_j is the equilibrium displacement (times f_0/(H_1+H_2)),
!        f_0 is the constant part of the Coriolis frequency, beta is the
!        planetary vorticity gradient, H_j is the depth of layer j, 
!        h_j = H_j/(H_1+H_2), H_b is the height of the bottom topography,
!        zeta_1 = Lap(psi_1) is the vorticity in layer 1, F is wind stress
!        forcing, therm_j is the thermal damping rate bringing interface j
!        back to a given equilibrium form deq_j, and r is the Ekman damping
!        rate applied only to the lower layer.

!        Note, d_1 & d_2 are related to psi_1 & psi_2 through
!               d_1=h_1*h_2*kd_bar^2*(psi_1-alpha*psi_2)  & 
!               d_2=h_1*h_2*(1-alpha)*kd_bar^2*psi_2
!        where alpha = rho_2/rho_1 is the density ratio and kd_bar 
!        is the prescribed mean baroclinic Rossby deformation wavenumber.
!        See constants.f90 for the construction of the vertical modes
!        from alpha, h1, h2 & kd_bar.

!        At the y boundaries, the mean x velocity takes prescribed 
!        values (see parameters.f90 and spectral.f90), except for the
!        special case of a purely barotropic mode (alpha = 1).  In this
!        case, the mean zonal velocity of this mode is set to zero.
!        The boundary velocities are found by specifying the PV gradient
!        in the lower layer, epsilon*beta.  The dimensionless parameter
!        epsilon determines the mean flow in each layer under the 
!        additional constraint that the global momentum is zero.  See
!        the script "flow-setup" for details.  Note: an additional mean
!        x velocity may be included to study e.g. the influence of 
!        bottom topography.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!        The full algorithm consists of the following modules:

!        casl.f90      : This source - main program loop, repeats successive 
!                        calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation 
!                        (through recontouring--evolution cycle);
!        evolution.f90 : Main time evolution module;
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        congen.f90    : Source code for contour-to-grid conversion;
!        generic.f90   : Generic service routines for CASL.
!----------------------------------------------------------------------------
program casl

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

   !Advect PV until next recontouring or end:
   call evolve

enddo
 !End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Close all files and end: 
call finalise

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, opens and reads
! input files, then opens output files.

implicit none

double precision:: qqmin(nz),qqmax(nz)
double precision:: td,wfac,tbase
integer:: ix,iy,iz,iread,nreads,k,m

!-----------------------------------------------------------------
 !Allocate memory for non-conservative terms if present:
if (damping) then
   allocate(sq(0:ny,0:nxm1,nz),sqpre(0:ny,0:nxm1,nz))
   allocate(wforce(0:ny,0:nxm1))
endif

 !Allocate memory for a variable upper surface displacement:
if (.not. barot) & 
     allocate(dd2eq(0:ny,0:nxm1),dd2(0:ny,0:nxm1),dd2pre(0:ny,0:nxm1))

!--------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral

 !Initialise constants and arrays for contour advection:
call init_contours

!-----------------------------------------------------------------
 !Offsets used to determine trajectory location in the grid during
 !semi-Lagrangian advection:
do ix=0,nxm1
   xig(ix)=dble(ix)
enddo

do iy=0,ny
   yig(iy)=dble(iy)
enddo

 !Read in equilibrium middle interface displacement (if present):
if (heating) then
   open(12,file='disp1eq.r8',form='unformatted', &
         access='direct',status='old',recl=2*nbytes)
   read(12,rec=1) td,dd1eq
   close(12)

   if (.not. barot) then
      open(12,file='disp2eq.r8',form='unformatted', &
            access='direct',status='old',recl=2*nbytes)
      read(12,rec=1) td,dd2eq
      close(12)
   endif
endif

 !Read in scaled bottom topography f_0*H_b/H_1 (if present):
if (topogr) then
   open(12,file='topo.r8',form='unformatted', &
         access='direct',status='old',recl=2*nbytes)
   read(12,rec=1) td,fhb
   close(12)
else
   fhb=zero
endif

if (wind) then
   !Prescribe wind stress forcing fwind*sin(2*pi*(y-y_min)/L_y)
   !(can generalise to 2D functions):
  wfac=twopi/dble(ny)
  do ix=0,nxm1
    do iy=0,ny
      wforce(iy,ix)=fwind*sin(wfac*dble(iy))
    enddo
  enddo
endif

!-----------------------------------------------------------------
 !Decide what kind of simulation this is:
if (loopinit .eq. 0) then
   !This is a new simulation starting from data in qq_init.r8

   !Start from t = 0:
   tinit=zero

   !Read in initial PV (into qs) and compute its domain average:
   open(11,file='qq_init.r8',form='unformatted', &
         access='direct',status='old',recl=2*nbytes)
   do iz=1,nz
      read(11,rec=iz) t,qs(:,:,iz)
      call average(qs(0,0,iz),qavg(iz))
   enddo
   close(11)

   !Compute contour interval for PV in each layer:
   call contint(qs,ncontq,qjump,qqmin,qqmax)
   !Write information to log file:
   do iz=1,nz
      write(*,'(a,i1,a,1x,f13.5,1x,f11.5,1x,f9.5)') ' layer',iz, &
              ' q_min, q_max, qjump = ',qqmin(iz),qqmax(iz),qjump(iz)
   enddo

   !Copy qs into qd for recontouring:
   qd=qs

   !Initially there are no contours; set relevant counters to zero:
   nq=0
   nptq=0
   do iz=1,nz
      jl2q(iz)=0
   enddo

else
   !In this case, we are restarting from existing data.  First read in
   !the contours and residual PV (see common.f90):
   call readcont(loopinit,1)
   !This also defines the entire PV field qs for restart and defines
   !the initial time, t.  It sets the value of tinit = t for use below.

endif

!--------------------------------------------------------------
 !Open files (depending on type of simulation):
if ((loopinit .eq. 0) .or. replace) then
   !This is either a new simulation starting from t = 0 (loopinit = 0)
   !or one which starts from a specified time in a simulation already
   !performed, and we only want to use the data at this time for the
   !initial conditions (no previous data is kept from the old simulation):

   if (replace) then
      !Obtain the last time in ene.asc BEFORE tinit found from contour files:
      open(15,file='evolution/ene.asc',status='old')
      td=zero
      do
         iread=0
         read(15,*,iostat=iread) td
         if (iread .ne. 0) exit 
         if (td .ge. tinit) exit
      enddo

      backspace(15)
      backspace(15)
      read(15,*) td
      write(*,*) ' Read ene.asc up to time t = ',td
      close(15)
      !Initialise counter and offset for writing direct access grid files:
      igrec=0
   else
      td=zero
      igrec=0
   endif

   !Set the correct initial time based on last entry in ene.asc:
   tinit=td
   t=tinit

   !Open all plain text diagnostic files (ending in ".asc"):
   open(12,file='evolution/monitor.asc',status='replace')
   open(13,file='evolution/norms.asc',status='replace')
   open(14,file='contours/complexity.asc',status='replace')
   open(15,file='evolution/ene.asc',status='replace')
   !Open files for coarse grid saves (unformatted files ending in ".r4"):
   open(31,file='evolution/qq1.r4',form='unformatted', &
         access='direct',status='replace',recl=nbytes)
   open(32,file='evolution/qq2.r4',form='unformatted', &
         access='direct',status='replace',recl=nbytes)
   open(33,file='evolution/dd1.r4',form='unformatted', &
         access='direct',status='replace',recl=nbytes)
   if (.not. barot) open(34,file='evolution/dd2.r4',form='unformatted', &
         access='direct',status='replace',recl=nbytes)
   !Open file for 1d PV spectra (mode 1 and mode 2; ascii file):
   open(51,file='spectra/qspec.asc',status='replace')
   !Open file for contour writes (synopsis file):
   open(80,file='contours/qqsynopsis.asc',status='replace')

   !Initialise counter for writing direct access contour files:
   icrec=1

else
   !Here, we continue the simulation from time t = tinit, reading 
   !from the data in the cont subdirectory.  All previous data is
   !kept; we append more data on to existing files after making 
   !sure that they are synchronised.  Note: the base time for the
   !simulation is read from the first line of ene.asc.
   open(12,file='evolution/monitor.asc',status='old')
   td=zero
   do
      iread=0
      read(12,*,iostat=iread) td
      if (iread .ne. 0) exit 
      if (td .ge. tinit) exit
   enddo
   backspace(12)
   if (iread .eq. 0) then
      write(*,*) ' Read monitor.asc up to time before t = ',td
   else
      write(*,*) ' Read monitor.asc up to time t = ',td
   endif

   open(13,file='evolution/norms.asc',status='old')
   td=zero
   do
      iread=0
      read(13,*,iostat=iread) td
      if (iread .ne. 0) exit 
      if (td .ge. tinit) exit
   enddo
   backspace(13)
   if (iread .eq. 0) then
      write(*,*) ' Read norms.asc up to time before t = ',td
   else
      write(*,*) ' Read norms.asc up to time t = ',td
   endif

   open(14,file='contours/complexity.asc',status='old')
   td=zero
   do
      iread=0
      read(14,*,iostat=iread) td
      if (iread .ne. 0) exit 
      if (td .ge. tinit) exit
   enddo
   backspace(14)
   if (iread .eq. 0) then
      write(*,*) ' Read complexity.asc up to time before t = ',td
   else
      write(*,*) ' Read complexity.asc up to time t = ',td
   endif

   open(15,file='evolution/ene.asc',status='old')
   !Obtain initial "base" time for simulation:
   read(15,*) tbase
   backspace(15)

   td=zero
   nreads=0
   do
      iread=0
      read(15,*,iostat=iread) td
      if (iread .ne. 0) exit 
      if (td .ge. tinit) exit
      nreads=nreads+1
   enddo
   backspace(15)
   if (iread .eq. 0) then
      write(*,*) ' Read ene.asc up to time before t = ',td
   else
      write(*,*) ' Read ene.asc up to time t = ',td
   endif

   !Open files for coarse grid saves:
   open(31,file='evolution/qq1.r4',form='unformatted', &
         access='direct',status='old',recl=nbytes)
   open(32,file='evolution/qq2.r4',form='unformatted', &
         access='direct',status='old',recl=nbytes)
   open(33,file='evolution/dd1.r4',form='unformatted', &
         access='direct',status='old',recl=nbytes)
   if (.not. barot) open(34,file='evolution/dd2.r4',form='unformatted', &
          access='direct',status='old',recl=nbytes)

   !Open file for 1d PV spectra (mode 1 & mode 2):
   open(51,file='spectra/qspec.asc',status='old')
   do m=1,nreads
      read(51,*) td
      do k=1,kmaxred
         read(51,*)
      enddo
   enddo
   write(*,*) ' Read qspec.asc up to t = ',td

   !Open file for contour writes:
   open(80,file='contours/qqsynopsis.asc',status='old')
   do m=1,loopinit
      read(80,*)
   enddo

   !Initialise counters for writing direct access files:
   igrec=nreads
   icrec=loopinit+1

   write(*,*) ' igrec = ',igrec
   write(*,*) ' icrec = ',icrec

   !Set initial time to the base time in ene.asc:
   tinit=tbase
   write(*,*) ' tinit = ',tinit

endif

return
end subroutine initialise

!=======================================================================

subroutine evolve

use evolution

implicit none

!Advect PV until next recontouring or end:
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
call recontour(qd)
write(*,'(a,i8,a,i9,a,f13.5)') ' nq = ',nq,'  nptq = ',nptq,'  t = ',t

return 
end subroutine recont

!=======================================================================

subroutine finalise

implicit none

write(*,*) ' Code completed normally'

 !Close output files (opened in subroutine initialise):
close(12) 
close(13)
close(14)
close(15)
close(31)
close(32)
close(33)
if (.not. barot) close(34)
close(51)
close(80)

return
end subroutine finalise

 !End main program
end program casl
!=======================================================================
