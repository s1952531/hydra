!###########################################################################
!   The Spherical Barotropic Combined Lagrangian Advection Method (CLAM)
!  Adapted from the shallow-water version on 30 May 2014 by D.G. Dritschel
!###########################################################################

!       Uses a 4th-order Runge-Kutta scheme with an adapted time step
!       together with 4th-order compact differencing is used in latitude.

!       Input data file:
!       =================
!        qq_init.r8         initial gridded PV (absolute vorticity)

!       Output data files:
!       =================
!       Every approximate tsave units of time:
!        zz.r4              relative vorticity

!       Every tsim units of time (every "period"):
!        cont/synopsis.asc  number of contours, nodes, etc 
!         "   indexnnn      contour counters, etc at period "nnn" 
!         "   nodesnnn      contour nodes at period "nnn"
!        resi/pnnn          residual PV field at period "nnn"
!        grid/zz.r4         relative vorticity

!       Every contour regularisation (surgery):
!        complexity.asc  number of PV contours, nodes and time

!       Every time step:
!        monitor.asc        time, energy, angular momentum, enstrophy, etc.

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

real:: zzr4(ng,nt),tr4
character(len=3):: pind

!--------------------------------------------------------------------
 !Initialise OpenMP:
call omp_set_num_threads(nthreads)
call omp_set_dynamic(.false.)

!--------------------------------------------------------------------
 !Initialise modules:
call init_spectral
call init_contours

!--------------------------------------------------------------------
 !Read latest relative vorticity in grid subdirectory:
open(35,file='grid/zz.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

loop=0
do  
  loop=loop+1
  iread=0
  read(35,rec=loop,iostat=iread) tr4,zzr4
  if (iread .ne. 0) exit 
  write(*,'(a,f9.2)') ' Read field data at t = ',tr4
enddo
loop=loop-2
write(*,*) ' *** loop = ',loop

 !Copy field into double precision array:
do i=1,nt
  do j=1,ng
    qs(j,i)=dble(zzr4(j,i))
  enddo
enddo

 !FFT qs in longitude (where it remains as a semi-spectral array):
call forfft(ng,nt,qs,trig,factors) 

 !Set dump counter for writing to correct record in main job directory:
idump=loop*nint(tsim/tsave)

 !Read contours at this time:
open(80,file='cont/synopsis.asc',status='old')
do k=0,loop
  read(80,*) t,n,npt
enddo

pind='000'
write(pind(1:3),'(i3.3)') loop

open(81,file='cont/index'//pind,form='unformatted',status='old')
read(81) np(1:n),i1(1:n),ind(1:n)
close(81)

open(82,file='cont/nodes'//pind,form='unformatted',status='old')
read(82) x(1:npt),y(1:npt),z(1:npt)
close(82)

 !Generate next() array for use in congen module:
do i=1,npt
  next(i)=i+1
enddo
do j=1,n
  next(i1(j)+np(j)-1)=i1(j)
enddo

write(*,*) ' Read PV contours at t = ',t
write(*,*) ' n = ',n,' npt = ',npt

 !Read residual PV:
open(83,file='cont/qd.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)
read(83,rec=loop+1) tr4,zzr4
write(*,*) ' Read residual PV at t = ',tr4

 !Copy residual PV into double precision array:
do i=1,nt
  do j=1,ng
    qd(j,i)=dble(zzr4(j,i))
  enddo
enddo

 !Set data save flag:
iref=-1

!--------------------------------------------------------------------
 !Open various diagnostic files:

 !Energy & other diagnostics:
open(18,file='monitor.asc',status='old')
td=zero
do while (td .lt. t+1.d-5)
  read(18,*) td
enddo
backspace(18)
write(*,*) ' Read monitor.asc up to time before t = ',td

 !Contour complexity (n, npt & t):
open(19,file='complexity.asc',status='old')
td=zero
do while (td .lt. t+1.d-5)
  read(19,*) td
enddo
backspace(19)
write(*,*) ' Read complexity.asc up to time before t = ',td

 !relative vorticity approximately every tsave time units:
open(45,file='zz.r4',form='unformatted', &
    & access='direct',status='old',recl=nbytes)

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

write(*,*) 'spe completed normally' 

 !Close output files (opened in subroutine initialise):
close(18)
close(19)
close(35)
close(45)
close(80)
close(83)

return 
end subroutine

 !End main program
end program
!=======================================================================
