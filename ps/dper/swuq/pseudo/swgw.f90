! Pseudo-spectral solver for the shallow water problem with constant PV.

!==========================================================================

!     The full algorithm consists of the following modules:
!        swgw.f90      : Main code for time-stepping using a PS scheme with
!                        trapezoidal semi-implicit time-stepping.
!        parameters.f90: User defined parameters for a simulation;
!        constants.f90 : Fixed constants used throughout the other modules;
!        common.f90    : Common data preserved throughout simulation 
!        spectral.f90  : Fourier transform common storage and routines;
!----------------------------------------------------------------------------
program swgw

use common
use spectral

implicit double precision(a-h,o-z)
implicit integer(i-n)

!---------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call initialise

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .lt. tsim)
  do while (t .lt. tsave)

     !Compute max height and other diagnostics:
    hhmin=one
    hhmax=one
    frmax=zero
    ddrms=zero
    do ix=1,nx
      do iy=1,ny
        usq=uu(iy,ix)**2+vv(iy,ix)**2
        hht=one+hhp(iy,ix)
        hhmin=min(hhmin,hht)
        hhmax=max(hhmax,hht)
        frmax=max(frmax,usq/hht)
        ddrms=ddrms+ddp(iy,ix)**2
      enddo
    enddo
    ddrms=sqrt(ddrms*dsumi)
    frmax=sqrt(frmax)/cgw
     !Record diagnostics to norms.dat
    write(18,'(f12.5,1x,4(1x,f9.6))') t,hhmin,hhmax,frmax,ddrms

     !Advect fields:
    call timestep

  enddo

   !Write data 
  call writedata
  tsave=tsave+tdsave
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Close output files (opened in subroutine initialise):
close(17)
close(18)
close(41)
close(42)
close(80)

write(*,*) 'code completed normally' 

!===============================================================

 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine initialise

! Routine initialises fixed constants and arrays, and reads in
! input files, opens output files ready for writing to. 

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: wka(ny,nx),wkb(ny,nx)
!--------------------------------------------------------------------
call init_spectral

!-----------------------------------------------------------------
 !If topography exists construct the source it induces in the
 !divergence equation:
if (topo) then
  open(14,file='topo.dat',status='old')
  do ix=1,nx
    do iy=1,ny
      read(14,*) hhb(iy,ix)
    enddo 
  enddo
  close(14)
endif

!--------------------------------------------------------------------
 !Read initial gridded relative height 
 !(possibly removing topography; see above):
open(12,file='hh_init.dat',status='old')
hbar=zero
do ix=1,nx
  do iy=1,ny
    read(12,*) hhp(iy,ix)
    wka(iy,ix)=hhp(iy,ix)
    hbar=hbar+hhp(iy,ix)
  enddo
enddo
close(12)
hbar=hbar*dsumi

 !Read initial divergence field:
open(13,file='dd_init.dat',status='old')
dbar=zero
do ix=1,nx
  do iy=1,ny
    read(13,*) ddp(iy,ix)
    wkb(iy,ix)=ddp(iy,ix)
    dbar=dbar+ddp(iy,ix)
  enddo
enddo
close(13)
dbar=dbar*dsumi

 !Remove global means of hh & dd:
do ix=1,nx
  do iy=1,ny
    hhp(iy,ix)=hhp(iy,ix)-hbar
    ddp(iy,ix)=ddp(iy,ix)-dbar
  enddo
enddo
!--------------------------------------------------------

 !Define spectral fields for time-stepping:
call ptospc(nx,ny,wka,hh,xfactors,yfactors,xtrig,ytrig)
call ptospc(nx,ny,wkb,dd,xfactors,yfactors,xtrig,ytrig)

 !De-aliase fields:
do ky=1,ny
  do kx=1,nx
    hh(kx,ky)=hh(kx,ky)*dafilt(kx,ky)
    dd(kx,ky)=dd(kx,ky)*dafilt(kx,ky)
  enddo
enddo

!--------------------------------------------------------------------
 !Open various diagnostic files:

 !Energy (kinetic, potential, total):
open(17,file='energy.dat',status='unknown')

 !Various diagnostics:
open(18,file='norms.dat',status='unknown')
open(80,file='spectra.dat',status='unknown')

 !h & delta every tsave time units:
open(41,file='hh.dat',status='unknown')
open(42,file='dd.dat',status='unknown')

!--------------------------------------------------------------------

 !Write initial conditions for later post-processing:
call invert(hh,dd,uu,vv)
call writedata

tsave=tdsave

return
end subroutine

!=======================================================================

subroutine timestep
!   Integrates the equations of motion from time t to time t + dt.
!   *** Uses implicit trapezoidal integration ***

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Define local parameters and arrays:
integer,parameter:: niter=2
!   niter+1: number of interations of the implicit scheme

double precision:: hho(nx,ny),ddo(nx,ny)
double precision:: rhh(nx,ny),rdd(nx,ny),rhho(nx,ny),rddo(nx,ny)

!------------------------------------------------------------------------
 !Find gridded velocity:
call invert(hh,dd,uu,vv)

 !Find tendencies:
call tends

!--------------------------------------------------------------------
 !Store quantities at the beginning of a time step needed in the
 !iteration below and take first forward time step:
do ky=1,ny
  do kx=1,nx
    hho(kx,ky)=hh(kx,ky)
    hshh=f12*shh(kx,ky)
    rhho(kx,ky)=hshh+alp*hh(kx,ky)
    rhh(kx,ky)=hshh+rhho(kx,ky)
     !rhh = R_h in sistep below

    ddo(kx,ky)=dd(kx,ky)
    hsdd=f12*sdd(kx,ky)
    rddo(kx,ky)=hsdd+alp*dd(kx,ky)
    rdd(kx,ky)=hsdd+rddo(kx,ky)
     !rdd = R_delta in sistep below
  enddo
enddo

 !Get approximate height & divergence at t + dt:
call sistep(hho,ddo,rhh,rdd)

!--------------------------------------------------------
 !Iterate now to solve the implicit problem in time:
do iter=1,niter

   !Update velocity at this time and interpolate it to the contours:
  call invert(hh,dd,uu,vv)

   !Update tendencies at this time:
  call tends

   !Define R_h & R_d (use rhh & rdd):
  do ky=1,ny
    do kx=1,nx
      rhh(kx,ky)=f12*shh(kx,ky)+rhho(kx,ky)
      rdd(kx,ky)=f12*sdd(kx,ky)+rddo(kx,ky)
    enddo
  enddo

   !Update height & divergence at t + dt:
  call sistep(hho,ddo,rhh,rdd)
enddo

 !Update height & divergence in physical space:
do ky=1,ny
  do kx=1,nx
    rhh(kx,ky)=hh(kx,ky)
    rdd(kx,ky)=dd(kx,ky)
  enddo
enddo
call spctop(nx,ny,rhh,hhp,xfactors,yfactors,xtrig,ytrig)
call spctop(nx,ny,rdd,ddp,xfactors,yfactors,xtrig,ytrig)

 !Increment time:
t=t+dt

return
end subroutine

!==========================================================================

subroutine tends
! Computes the tendencies for height & divergence (shh & sdd, 
! omitting the linear part used in semi-implicit time stepping).

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: wkap(ny,nx),wkbp(ny,nx),wkcp(ny,nx),wkdp(ny,nx)
double precision:: wkas(nx,ny),wkbs(nx,ny),wkcs(nx,ny),wkds(nx,ny)
!--------------------------------------------------------------------------

 !Compute nonlinear height tendency:
if (topo) then
  do ix=1,nx
    do iy=1,ny
      hdif=hhp(iy,ix)-hhb(iy,ix)
      wkap(iy,ix)=hdif*uu(iy,ix)
      wkbp(iy,ix)=hdif*vv(iy,ix)
    enddo
  enddo
else
  do ix=1,nx
    do iy=1,ny
      wkap(iy,ix)=hhp(iy,ix)*uu(iy,ix)
      wkbp(iy,ix)=hhp(iy,ix)*vv(iy,ix)
    enddo
  enddo
endif

call ptospc(nx,ny,wkap,wkas,xfactors,yfactors,xtrig,ytrig)
call ptospc(nx,ny,wkbp,wkbs,xfactors,yfactors,xtrig,ytrig)

call xderiv(nx,ny,rkx,wkas,wkcs)
call yderiv(nx,ny,rky,wkbs,wkds)

do ky=1,ny
  do kx=1,nx
    shh(kx,ky)=-dafilt(kx,ky)*(wkcs(kx,ky)+wkds(kx,ky))
  enddo
enddo

 !Compute nonlinear divergence tendency:
do ix=1,nx
  do iy=1,ny
    wkap(iy,ix)=ddp(iy,ix)*uu(iy,ix)
    wkbp(iy,ix)=ddp(iy,ix)*vv(iy,ix)
  enddo
enddo

call ptospc(nx,ny,wkap,wkas,xfactors,yfactors,xtrig,ytrig)
call ptospc(nx,ny,wkbp,wkbs,xfactors,yfactors,xtrig,ytrig)

call xderiv(nx,ny,rkx,wkas,wkcs)
call yderiv(nx,ny,rky,wkbs,wkds)

do ky=1,ny
  do kx=1,nx
    sdd(kx,ky)=wkcs(kx,ky)+wkds(kx,ky)
  enddo
enddo

do ix=1,nx
  do iy=1,ny
    wkap(iy,ix)=uu(iy,ix)
    wkbp(iy,ix)=vv(iy,ix)
  enddo
enddo

 !Obtain u_x & u_y: 
call ptospc(nx,ny,wkap,wkas,xfactors,yfactors,xtrig,ytrig)
call xderiv(nx,ny,rkx,wkas,wkcs)
call yderiv(nx,ny,rky,wkas,wkds)
call spctop(nx,ny,wkcs,wkcp,xfactors,yfactors,xtrig,ytrig)
call spctop(nx,ny,wkds,wkdp,xfactors,yfactors,xtrig,ytrig)

 !Obtain v_x & v_y:
call ptospc(nx,ny,wkbp,wkbs,xfactors,yfactors,xtrig,ytrig)
call xderiv(nx,ny,rkx,wkbs,wkcs)
call yderiv(nx,ny,rky,wkbs,wkds)
call spctop(nx,ny,wkcs,wkap,xfactors,yfactors,xtrig,ytrig)
call spctop(nx,ny,wkds,wkbp,xfactors,yfactors,xtrig,ytrig)

 !Form 2*J(u,v):
do ix=1,nx
  do iy=1,ny
    wkap(iy,ix)=two*(wkcp(iy,ix)*wkbp(iy,ix)-wkdp(iy,ix)*wkap(iy,ix))
  enddo
enddo

call ptospc(nx,ny,wkap,wkas,xfactors,yfactors,xtrig,ytrig)

do ky=1,ny
  do kx=1,nx
    sdd(kx,ky)=dafilt(kx,ky)*(wkas(kx,ky)-sdd(kx,ky))
  enddo
enddo

return
end subroutine

!========================================================================

subroutine sistep(hho,ddo,rhh,rdd)
! Updates height & divergence given R_h & R_delta (in the arrays rhh & rdd) 
! as well as the original values of height & divergence (in the arrays hho
! & ddo).  Here, R_h = (shh^n+ssh^{n+1})/2+(2/dt)*hh^n where n
! refers to the current time level and n+1 refers to the next one.

implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Arguments declarations:
double precision:: hho(nx,ny),ddo(nx,ny),rhh(nx,ny),rdd(nx,ny)

!------------------------------------------------------------------------

do ky=1,ny
  do kx=1,nx
    dda=helmi(kx,ky)*(omsq(kx,ky)*rhh(kx,ky)+alp*rdd(kx,ky))
    hha=dt2*(rhh(kx,ky)-dda)
    hh(kx,ky)=two*hha-hho(kx,ky)
    dd(kx,ky)=two*dda-ddo(kx,ky)
  enddo
enddo

 !Update height & divergence in physical space:
do ky=1,ny
  do kx=1,nx
    rhh(kx,ky)=hh(kx,ky)
    rdd(kx,ky)=dd(kx,ky)
  enddo
enddo
call spctop(nx,ny,rhh,hhp,xfactors,yfactors,xtrig,ytrig)
call spctop(nx,ny,rdd,ddp,xfactors,yfactors,xtrig,ytrig)

return
end subroutine

!=======================================================================

subroutine writedata
! Writes height anomaly and divergence to hh.dat and dd.dat.  
! Also computes and writes the total energy to energy.dat.

implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: hspec(0:max(nx,ny)),dspec(0:max(nx,ny)),uspec(0:max(nx,ny))
!------------------------------------------------------------
 !Write fields for later imaging or diagnostics:
write(41,'(f12.5)') t
write(42,'(f12.5)') t
do ix=1,nx
  do iy=1,ny
    write(41,'(1x,e14.7)')  hhp(iy,ix)
    write(42,'(1x,e14.7)')  ddp(iy,ix)
  enddo
enddo

!------------------------------------------------------------------
 !Compute and write energy per unit area:
eke=zero
epe=zero
do ix=1,nx
  do iy=1,ny
    eke=eke+(one+hhp(iy,ix))*(uu(iy,ix)**2+vv(iy,ix)**2)
    epe=epe+csq*hhp(iy,ix)**2
  enddo
enddo
eke=f12*eke*dsumi
epe=f12*epe*dsumi

write(17,'(f12.5,3(1x,f15.11))') t,eke,epe,eke+epe

 !Compute 1d vorticity spectrum:
call spec1d(hh,hspec)
call spec1d(dd,dspec)

 !Normalise to take into account uneven sampling of wavenumbers 
 !in each shell [k-1/2,k+1/2]:
do k=1,kmaxred
  hspec(k)=spmf(k)*hspec(k)+1.d-32
  dspec(k)=spmf(k)*dspec(k)+1.d-32
  uspec(k)=ksqi(k)*(fsq*hspec(k)+dspec(k))
enddo

 !Write out spectra to file:
write(80,'(1x,f12.5,1x,i6)') t,kmaxred
 !kmaxred = kmax/sqrt(2) to avoid shells in the upper corner of the
 !          kx,ky plane which are not fully populated
do k=1,kmaxred
  write(80,'(4(1x,f12.8))') alk(k),log10(hspec(k)),log10(dspec(k)),log10(uspec(k))
enddo
 !Note: alk(k) = log_10(k)

return
end subroutine

!========================================================================

 !End main program
end program
!=======================================================================
