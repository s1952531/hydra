!#####################################################################
!            The Combined Lagrangian Advection Method for
!        2D non-rotating Boussinesq flow in channel geometry
!        after a conformal transformation of the original domain. 
!#####################################################################

!     Code originally developed by Stuart King & David Dritschel @ 
!     St Andrews on 7 September 2012.

!     Code extended by David & Heidi Dritschel @ St Andrews to allow
!     for any domain geometry on 18 January 2019.

!     This code solves 
!                        Dzeta/Dt = db/dx = S_zz  
!                           Db/Dt = 0
!                          div(u) = 0
!     in a general domain having free slip boundary conditions.

!     We use contour advection for buoyancy, whereas for vorticity we
!     use both contour advection and the pseudo-spectral method as
!     described in Dritschel & Fontane, J. Comput. Phys. 229, 5408-5417
!     (2010).

!     For time integration, we use an implicit (iterated) trapezoidal 
!     time-stepping method.

!     The source S_zz is obtained directly from the buoyancy contours,
!     not the gridded field, essentially by integrating over the delta
!     function source and interpolating to nearby grid points - see
!     getzzsrc in the contours.f90 module.

!     Incompressibility is handled via the introduction of a stream-
!     function such that:
!             Lap(psi) = zeta
!     Determining psi (hence u,v) at each time step is termed the 
!     'inversion' problem.  Here this is done using a spectral method
!     together with the exact solutions of Lap(psi) = 0 to impose
!     the boundary conditions.

!     All of this is modified by the conformal transformation; see
!     uses of confac, confaci, yox and yoy to see the changes.

!     The full algorithm consists of the following modules:
!        casl.f90      : This source - main program loop, repeats successive 
!                        calls to evolve fields and recontour;
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        variables.f90 : Global quantities that may change in time;
!        common.f90    : Common data preserved throughout simulation 
!                        (through recontouring--evolution cycle);
!                      : and general service routines;
!        spectral.f90  : Fourier transform common storage and routines;
!        contours.f90  : Contour advection common storage and routines;
!        congen.f90    : Source code for contour-to-grid conversion;
!        evolution.f90 : Main time evolution module - advects gridded fields 
!                        using SL method along with contours.
!----------------------------------------------------------------------------
program caps

use common

implicit none

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tfin)

   !Obtain new buoyancy & vorticity contours:
  call recont

   !Advect buoyancy & vorticity until next recontouring or end:
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

!-----------------------------------------------------------------
 !Read in initial vorticity (zs) and buoyancy (bb):
open(11,file='zz_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) dum,zs
close(11)
open(12,file='bb_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(12,rec=1) t,bb
close(12)

 !Compute contour interval for buoyancy:
bbmax=zero
bbmin=zero
do ix=0,nx
  do iy=0,ny
    bbmax=max(bbmax,bb(iy,ix))
    bbmin=min(bbmin,bb(iy,ix))
  enddo
enddo
bjump=(bbmax-bbmin)/dble(ncontb)

 !Define a maximum time step based on the buoyancy range and grid scale:
dtmax=f12*sqrt(min(glx,gly)/(bbmax-bbmin))

 !Write information to log file:
write(*,*)
write(*,'(a,3(1x,f9.5))') ' b_min, b_max, bjump = ',bbmin,bbmax,bjump

 !Copy zs into zd for recontouring:
do ix=0,nx
  do iy=0,ny
    zd(iy,ix)=zs(iy,ix)
  enddo
enddo

 !Compute contour interval for vorticity:
call l1norm(zs,zzl1)
call l2norm(zs,zzl2)

if (zzl1 .gt. zero) then 
  zjump=(zzl2/zzl1)/dble(ncontz)
else
  zjump=zero
endif

!------------------------------------------------------------
 !Initially there are no contours:
nb=0
nptb=0
nz=0
nptz=0

 !Initialise time step so that subroutine adapt chooses a suitable one:
dt=zero

 !Set final time for simulation end:
itime=int((t+small)/tgsave)
tgrid=tgsave*dble(itime)
tfin=tgrid+tsim

!--------------------------------------------------
 !Call initialisation routines from modules:

 !Initialise inversion constants and arrays:
call init_spectral(bbmax-bbmin)
 !Initialise constants and arrays for contour advection:
call init_contours

 !Compute domain area (domarea):
qavg=zero
do ix=1,nxm1
  qavg=qavg+confac(0,ix)+confac(ny,ix)
enddo
do iy=1,nym1
  qavg=qavg+confac(iy,0)+confac(iy,nx)
enddo
qavg=f12*qavg+f14*(confac(0,0)+confac(ny,0)+confac(0,nx)+confac(ny,nx))

do ix=1,nxm1
  do iy=1,nym1
    qavg=qavg+confac(iy,ix)
  enddo
enddo

domarea=garea*qavg
 !Note: garea is the grid box area, glx*gly, of the conformal domain.

 !Useful factor for computing domain averages of quantities:
domsumi=one/qavg

 !Compute reference potential energy (for simplicity, initial PE):
eperef=zero
call potential(bb,epe)
eperef=epe

!--------------------------------------
 !Open all diagnostic plain text  files:     
open(12,file='monitor.asc',status='unknown')
open(13,file='norms.asc',status='unknown')
open(14,file='complexity.asc',status='unknown')
open(15,file='ene.asc',status='unknown')
 !Open files for vorticity and buoyancy coarse
 !grid direct writes:
open(31,file='zz.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(32,file='bb.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
 !Open files for contour writes:
open(80,file='cont/zzsynopsis.asc',status='unknown')
open(83,file='cont/zzresi.r4',form='unformatted', &
    & access='direct',status='replace',recl=nbytes)
open(90,file='cont/bbsynopsis.asc',status='unknown')

 !Initialise counter for writing direct files to the correct counter:
igrids=0

return
end subroutine

!=======================================================================

subroutine evolve

use evolution

implicit none

!Advect buoyancy & vorticity until next recontouring or end:
write(*,*) 'Evolving contours and fields ...'
call advect

return 
end subroutine

!=======================================================================

subroutine recont

use congen

implicit none

!Obtain new buoyancy contours:
write(*,*) 'Recontouring buoyancy ...'
call recontour(bb,xb,yb,bjump,bavg,nextb,indb,npb,i1b,i2b,nb,nptb,0)
write(*,'(a,i8,a,i9,a,f9.5)') '   nb = ',nb,'   nptb = ',nptb,'   db = ',bjump

!Obtain new vorticity contours:
write(*,*) 'Recontouring vorticity ...'
call recontour(zd,xz,yz,zjump,zavg,nextz,indz,npz,i1z,i2z,nz,nptz,1)
write(*,'(a,i8,a,i9,a,f9.5)') '   nz = ',nz,'   nptz = ',nptz,'   dz = ',zjump

return 
end subroutine

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
close(80)
close(83)
close(90)

return
end subroutine

 !End main program
end program
!=======================================================================

