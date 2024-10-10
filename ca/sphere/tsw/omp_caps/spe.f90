!###########################################################################
!  The Spherical Shallow-Water Combined Lagrangian Advection Method (CLAM)
!               D.G. Dritschel, 3-22 June 2010, St Andrews
!###########################################################################

!---------------------------------------------------------------------------
!    Major revision on 3 February 2013 to use height and divergence as
!    prognostic variables, and to apply bi-harmonic divergence damping.  
!    Simplified bi-harmonic damping is also applied to residual PV.
!    De-aliasing is used on all nonlinear products.

!    f90 revision 14 February 2013 by Stuart King @ St Andrews
!    file extensions modified 24 June 2013 - all files now unformatted 
!    with direct access 
!      .r4  --> single precision file
!      .r8  --> double precision file
!      .asc --> human-readable formatted ascii file

!---------------------------------------------------------------------------

!       Time integration scheme: Iterative implicit trapezoidal rule, i.e.
!                                x^{n+1}-x^{n} = 0.5*dt*(u^{n+1}+u^{n})
!                                where u = dx/dt.  For height and divergence,
!                                a semi-implicit scheme is used as well,
!                                using the above rule with the equations
!                                delta_t + Lh = S_delta; h_t + delta = S_h
!                                where L = c^2*Lap - f^2 and S_delta & S_h
!                                are the remaining tendencies.

!       Spatial discretization:  Nodes are redistributed at a given time
!                                interval according to a weighted sum of
!                                nearby curvature values.  A semi-spectral
!                                method is used for the PV fields and for
!                                height and divergence.  Fourth-order 
!                                compact differencing is used in latitude.

!       Input data files:
!       =================
!        pv_init.r8      *total* initial gridded PV 
!         depdiv.r8      initial depth (or height) and divergence fields.
!        (hequil.r8      thermal equilibrium height field)                     
!        (topogr.r8      gridded Laplacian of topography h_b/H)
!                        

!       Output data files:
!       =================
!       Every approximate tsave units of time:
!        qq.r4           PV
!        zz.r4           relative vorticity
!        hh.r4           dimensionless height anomaly
!        uu.r4           zonal velocity
!        vv.r4           meridional velocity
!        dd.r4           divergence
!        ene-ang.asc     Kinetic, potential & total energy & angular momentum

!     Every tsim units of time (every "period"):
!      cont/synopsis.asc  number of contours, nodes, etc 
!       "   indexnnn      contour counters, starting points, etc at period "nnn" 
!       "   nodesnnn      contour nodes at period "nnn"
!      resi/pnnn          residual PV field at period "nnn"
!      grid/qq.r4         PV 
!        "  zz.r4         relative vorticity
!        "  hh.r4         dimensionless height anomaly
!        "  uu.r4         zonal velocity
!        "  vv.r4         meridional velocity
!        "  dd.r4         divergence

!       Every contour regularisation (surgery):
!        complexity.asc  number of PV contours, nodes and time
 
!       Every time step:
!        cfl.asc         time, cfl, enstrophy, max(zeta)/f_pole, min(h),
!                        max(h), max(Fr) and "twist" parameter

!==========================================================================


!     The full algorithm consists of the following modules:
!        spe.f90       : This source - main program loop, repeats successive 
!                        calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation 
!                        (through recontouring--evolution cycle);
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        congen.f90    : Source code for contour-to-grid conversion;
!        evolution.f90 : Main time evolution module - advects gridded fields 
!                        using PS method along with contours.
!----------------------------------------------------------------------------
program spe

use common

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .lt. tfin)

   !Obtain new PV contours:
  call recont

   !Advect fields until next recontouring or end:
  call evolve

enddo

!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

call finalise

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit double precision(a-h,o-z)
implicit integer(i-n)

!--------------------------------------------------------------------
call init_spectral
call init_contours

!-----------------------------------------------------------------
if (topogr) then
   !If there is a topography (h_b), read the Lap(h_b)/H, where H
   !is the mean fluid depth, and multipy it by c^2:
  open(20,file='topogr.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(20,rec=1) dum,ttdd
  close(20)

  do i=1,nt
    do j=1,ng
      ttdd(j,i)=csq*dum
    enddo
  enddo
endif

!---------------------------------------------------------------------
if (stoch) then
   !A random vorticity field with enstrophy esr*dt and
   !spectrum proportional to k^5*exp(-2k^2/k0^2) is added
   !at the end of each time step.

   !Initialize random # generator on first call:
  do i=1,iseed
    uni=rand(0)
  enddo

   !Generate squared wavenumber arrays used in the forcing:
  rksri=one/dble(ksr)
  do m=1,nt
    do j=1,ng
      plon(j,m)=(rksri*wave(m)*clati(j))**2
      glon(j,m)=exp(-plon(j,m))
    enddo
  enddo
  do k=1,nt
    plat(k)=(rksri*wave(k))**2
    glat(k)=exp(-plat(k))
  enddo

endif

!---------------------------------------------------------------------
if (thermal) then
   !Read the thermal equilibrium height anomaly in hequil.dat
  open(20,file='hequil.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(20,rec=1) dum,hhe
  close(20)
endif

!--------------------------------------------------------------------
 !Read initial gridded PV into qd for initial contouring:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(20,rec=1) t,qd
close(20)

do i=1,nt
  do j=1,ng
    qs(j,i)=qd(j,i)
  enddo
enddo

!--------------------------------------------------------
 !See if PV ramping will be used for initialisation:
ramping=((t .eq. zero) .and. (tramp .gt. zero))

if (ramping) then
   !The PV anomaly is ramped up from 0 to a fixed form 
   !over a time period of tramp; thereafter, the simulation
   !restarts with freely-evolving PV (t is reset to 0)
  do i=1,nt
    do j=1,ng
      dd(j,i)=zero
      hh(j,i)=zero
    enddo
  enddo
   !Do this to avoid writing fields during ramping:
  tsavori=tsave
  tmaxori=tmax

  tdsave=two*tramp
  tmax=tramp

  iref=0
  nper=nint(tramp/tsim)+1
else

   !Read initial height field:
  open(20,file='hh_init.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(20,rec=1) dum,hh
  close(20)

   !Read initial divergence field:
  open(20,file='dd_init.r8',form='unformatted', &
      & access='direct',status='old',recl=2*nbytes)
  read(20,rec=1) dum,dd
  close(20)

   !De-aliase fields and convert fields to semi-spectral space:
  call dealiase(hh)
  call dealiase(dd)

   !Save a copy of hh & dd in physical space:
  do m=1,nt
    do j=1,ng
      hhp(j,m)=hh(j,m)
      ddp(j,m)=dd(j,m)
    enddo
  enddo
  call revfft(ng,nt,hhp,trig,factors) 
  call revfft(ng,nt,ddp,trig,factors) 

  iref=-1
  tdsave=tsave
  nper=nperiod
endif

!--------------------------------------------------------------------
 !Open various diagnostic files:

 !Energy (kinetic, potential, total) and angular momentum:
open(17,file='ene-ang.asc',status='unknown')

 !CFL parameter & other diagnostics:
open(18,file='cfl.asc',status='unknown')

 !Contour complexity (n, npt & t):
open(19,file='complexity.asc',status='unknown')

 !h, u, v, zeta, PV & delta every tsim time units:
open(31,file='grid/hh.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(32,file='grid/uu.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(33,file='grid/vv.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(34,file='grid/zz.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(35,file='grid/qq.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(36,file='grid/dd.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

 !h, u, v, zeta, PV & delta approximately every tdsave time units:
open(41,file='hh.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(42,file='uu.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(43,file='vv.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(44,file='zz.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(45,file='qq.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(46,file='dd.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

 !Longitudinal spectra for delta & residual PV every tdsave time units: 
open(50,file='long-spec.asc',status='unknown')

 !Contour summary & residual PV files for fine-grid PV reconstruction:
open(80,file='cont/synopsis.asc',status='unknown')
open(83,file='cont/qd.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)

 !Set dump counter for writing to correct record:
idump=0

return
end subroutine

!=======================================================================

subroutine evolve

use evolution

implicit none

!Advect fields until next recontouring or end:
write(*,*) 'Evolving contours and fields ...'
call advect

return 
end subroutine

!=======================================================================

subroutine recont

use congen

implicit none

!Obtain new PV contours:
write(*,*) 'Recontouring PV ...'
call recontour
write(*,'(a,i9,a,i10,a,f12.5)') '   n = ',n,'   npt = ',npt,'   t = ',t
write(19,'(1x,f12.5,1x,i9,1x,i10)') t,n,npt

return 
end subroutine

!=======================================================================

subroutine finalise

implicit none

write(*,*) 'clam completed normally' 

 !Close output files (opened in subroutine initialise):
close(17)
close(18)
close(19)
close(31)
close(32)
close(33)
close(34)
close(35)
close(36)
close(41)
close(42)
close(43)
close(44)
close(45)
close(46)
close(50)
close(80)
close(83)

return 
end subroutine

 !End main program
end program
!=======================================================================


