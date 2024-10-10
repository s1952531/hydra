!############################################################################
!  The Spherical Thermal Shallow-Water Combined Lagrangian Advection Method
!              D.G. Dritschel, 6-8 August 2019, New York
!############################################################################

!----------------------------------------------------------------------------
!    Extension of the spherical SW model in ~/hydra/ca/sphere/sw/caps
!    allowing a more realistic thermal structure - using Salby's model.
!----------------------------------------------------------------------------

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

!       Prognostic Variables:
!       =====================
!        qq   potential vorticity (PV) = (zeta+f)/(1+eta): uses CLAM approach
!        ee   eta in definition of PV; related to height H by H/H_char = h =
!             (1+ee)^kappa where kappa=2/7 and H_char is a characteristic depth
!        tt   log dimensionless potential temperature = log(theta/theta_char)
!             where theta_char is a characteristic potential temperature
!             --- only if thermal forcing is present
!        dd   divergence of the horizontal velocity field

!       Auxiliary Variables:
!       ====================
!        zz   relative vorticity = zeta = curl of the velocity field in
!             the radial (upward) direction
!        uu   zonal velocity component
!        vv   meridional velocity component

!       Input data files:
!       =================
!        qq_init.r8    initial potential vorticity (PV)
!        hh_init.r8    initial thickness h = H/H_char
!       (tt_init.r8    initial log dimensionless potential temperature, tau)
!        dd_init.r8    initial divergence, delta
!        (hequil.r8    thermal equilibrium dimensionless height anomaly field)
!        (topogr.r8    dimensionless bottom topography h_b = H_b/H_char)

!       Output data files:
!       =================
!       Every approximate tsave units of time:
!        qq.r4           PV
!        zz.r4           relative vorticity, zeta
!        hh.r4           thickness, h
!       (tt.r4           log dimensionless potential temperature, tau)
!        uu.r4           zonal velocity, u
!        vv.r4           meridional velocity, v
!        dd.r4           divergence, delta
!        ene-ang.asc     Kinetic, potential & total energy & angular momentum

!     Every tsim units of time (every "period"):
!      cont/synopsis.asc  number of contours, nodes, etc 
!       "   indexnnn      contour counters, etc at period "nnn" 
!       "   nodesnnn      contour nodes at period "nnn"
!      resi/pnnn          residual PV field at period "nnn"
!      grid/qq.r4         PV 
!        "  zz.r4         zeta
!        "  hh.r4         h
!       ("  tt.r4         tau)
!        "  uu.r4         u
!        "  vv.r4         v
!        "  dd.r4         delta

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

implicit none

double precision:: dum,uni,rksri
integer:: i,j,m,ierr

!---------------------------------------------------------------------
 !Initialise spectral and contours modules:
call init_spectral
call init_contours

!---------------------------------------------------------------------
 !If there is topography (H_b), read h_b=H_b/H_char where H_char is a
 !characteristic fluid depth, and compute c^2*Lap(h_b) where c^2=g*H_char:
open(20,file='topogr.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes,iostat=ierr)
topogr=(ierr .eq. 0)
if (topogr) then
  read(20,rec=1) dum,hhb
  close(20)

   !Calculate Laplacian of c^2*h_b:
  topo=csq*hhb
  call laplace(topo)
endif

!---------------------------------------------------------------------
 !Initialise stochastic forcing if used:
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
  plat=(rksri*wave)**2
  glat=exp(-plat)
endif

!--------------------------------------------------------------------
 !Read initial gridded PV into qd for initial contouring:
open(20,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(20,rec=1) t,qd
close(20)

 !Copy into qs for initialisation in module evolution (see init):
qs=qd

 !Read initial thickness h (hhp) and define eta = h^{1/kappa} - 1 (eep):
open(20,file='hh_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(20,rec=1) dum,hhp
close(20)
eep=hhp**kappai-one
 !Make a spectral copy of eta for time stepping (ee):
ee=eep
call forfft(ng,nt,ee,trig,factors) 

 !Read initial log (dimensionless) potential temperature field tau (ttp)
 !only if thermal forcing is present:
if (thermal) then
  open(20,file='tt_init.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(20,rec=1) dum,ttp
  close(20)
   !Make a spectral copy of tau for time stepping (tt):
  tt=ttp
  call forfft(ng,nt,tt,trig,factors) 
endif

 !Read initial divergence field, delta (ddp):
open(20,file='dd_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(20,rec=1) dum,ddp
close(20)
 !Make a spectral copy of delta for time stepping (dd):
dd=ddp
call forfft(ng,nt,dd,trig,factors) 

 !De-aliase fields and convert fields to semi-spectral space:
call dealiase(ee)
if (thermal) call dealiase(tt)
call dealiase(dd)

 !Save copies of eta, h and possibly tau in physical space:
eep=ee
call revfft(ng,nt,eep,trig,factors)
hhp=(one+eep)**kappa
if (thermal) then
  ttp=tt
  call revfft(ng,nt,ttp,trig,factors)
endif

iref=-1

!---------------------------------------------------------------------
if (thermal) then
   !Read the thermal equilibrium height anomaly (hhe) in hequil.dat
  open(20,file='hequil.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(20,rec=1) dum,hhe
  close(20)
endif

!--------------------------------------------------------------------
 !Open various diagnostic files:

 !Energy (kinetic, potential, total) and angular momentum:
open(17,file='ene-ang.asc',status='unknown')

 !CFL parameter & other diagnostics:
open(18,file='cfl.asc',status='unknown')

 !Contour complexity (n, npt & t):
open(19,file='complexity.asc',status='unknown')

 !h, possibly tau, u, v, zeta, PV & delta every tsim time units:
open(30,file='grid/hh.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
if (thermal) then
  open(31,file='grid/tt.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
endif
open(32,file='grid/uu.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(33,file='grid/vv.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(34,file='grid/zz.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(35,file='grid/qq.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(36,file='grid/dd.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)

 !h, possibly tau, u, v, zeta, PV & delta approximately every tsave time units:
open(40,file='hh.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
if (thermal) then
  open(41,file='tt.r4',form='unformatted', &
       access='direct',status='replace',recl=nbytes)
endif
open(42,file='uu.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(43,file='vv.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(44,file='zz.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(45,file='qq.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)
open(46,file='dd.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)

 !Longitudinal spectra for h, delta, PV & residual PV every tsave time units: 
open(50,file='long-spec.asc',status='unknown')

 !Contour summary & residual PV files for fine-grid PV reconstruction:
open(80,file='cont/synopsis.asc',status='unknown')
open(83,file='cont/qd.r4',form='unformatted', &
      access='direct',status='replace',recl=nbytes)

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

write(*,*)
write(*,*) '   spe completed normally' 

 !Close output files (opened in subroutine initialise):
close(17)
close(18)
close(19)
close(30)
if (thermal) close(31)
close(32)
close(33)
close(34)
close(35)
close(36)
close(40)
if (thermal) close(41)
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
