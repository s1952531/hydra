!#########################################################################
!  Finds balanced fields given only the PV anomaly in qq_init.r8 or
!  in q1.r4 & q2.r4, using Optimal PV balance with fixed PV contours
!  (see Viudez & Dritschel, J. Fluid Mech. 521, 343-352, 2004).

!  The shape of the PV distribution is held fixed but its amplitude is
!  ramped up according to q_j(x,y,t)=qb_j(x,y)*R(t) where qb_j(x,y) is
!  the "base" configuration (read in from either qq_init.r8 or q1.r4
!  q2.r4), while R(t) = 0.5*(1 - cos(pi*t/T)) is a ramp function going
!  from 0 at t = 0 to 1 at t = T, the specified ramping period.

!  This can be used both for initialisation, and for post-processing.

!       Written 3 December 2020 by D G Dritschel @ St Andrews
!#########################################################################

program fopvbal

 !Import constants and parameters:
use constants
 !Import spectral module:
use spectral

implicit none
double precision:: h1(ng,ng),d1(ng,ng),g1(ng,ng),q1(ng,ng)
double precision:: u1(ng,ng),v1(ng,ng),z1(ng,ng)
double precision:: h2(ng,ng),d2(ng,ng),g2(ng,ng),q2(ng,ng)
double precision:: u2(ng,ng),v2(ng,ng),z2(ng,ng)
double precision:: qb1(ng,ng),qs1(ng,ng),ds1(ng,ng),gs1(ng,ng)
double precision:: qb2(ng,ng),qs2(ng,ng),ds2(ng,ng),gs2(ng,ng)
double precision:: ekih,ekiv,epot,etot,t,tramp
real:: qr4(ng,ng),tr4
integer:: iopt,loop,iread,nstep

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

write(*,*) ' Choose one of the following options:'
write(*,*)
write(*,*) ' (0) initialise from data in qq_init.r8'
write(*,*) ' (1) post-process from data in q1.r4 & q2.r4'
write(*,*)
write(*,*) ' Option?'
read(*,*) iopt
write(*,*)
write(*,*) ' Ramping period, T?'
read(*,*) tramp
 !Ensure tramp is a multiple of dt:
nstep=nint(tramp/dt)
tramp=dt*dble(nstep)

if (iopt .eq. 0) then
   !Only find balance at t = 0 to initialise:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(11,rec=1) t,qb1
  read(11,rec=2) t,qb2
  close(11)

   !Find balanced fields:
  call balance

   !Save data:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,q1
  write(11,rec=2) zero,q2
  close(11)

  open(11,file='dd_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,d1
  write(11,rec=2) zero,d2
  close(11)

  open(11,file='hh_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,h1
  write(11,rec=2) zero,h2
  close(11)

  open(11,file='gg_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,g1
  write(11,rec=2) zero,g2
  close(11)

  open(11,file='zz_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,z1
  write(11,rec=2) zero,z2
  close(11)

  write(*,*)
  write(*,*) ' The balanced fields are available in xx_init.r8 where'
  write(*,*) ' xx = dd, hh, gg, zz & qq (adjusted only by a constant).'
  write(*,*)
  
else
   !Find balance at all times in q1.r4 & q2.r4:
  open(30,file='evolution/q1.r4',form='unformatted', &
        access='direct',status='old',recl=nbytes)
  open(40,file='evolution/q2.r4',form='unformatted', &
        access='direct',status='old',recl=nbytes)

  open(31,file='evolution/bq1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(41,file='evolution/bq2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(32,file='evolution/bd1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(42,file='evolution/bd2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(33,file='evolution/bh1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(43,file='evolution/bh2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(34,file='evolution/bg1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(44,file='evolution/bg2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(35,file='evolution/bz1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(45,file='evolution/bz2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(51,file='evolution/becomp.asc',status='replace')

   !Read data at all times and balance:
  loop=0
  do
    loop=loop+1
    iread=0
    read(30,rec=loop,iostat=iread) tr4,qr4
    if (iread .ne. 0) exit
    qb1=dble(qr4)
    read(40,rec=loop,iostat=iread) tr4,qr4
    qb2=dble(qr4)

     !Balance using only the PV anomaly at time t = tr4:
    write(*,'(a,f9.2)') ' *** Processing t = ',tr4
    call balance

     !Write data:
    write(31,rec=loop) tr4,real(q1)
    write(41,rec=loop) tr4,real(q2)

    write(32,rec=loop) tr4,real(d1)
    write(42,rec=loop) tr4,real(d2)

    write(33,rec=loop) tr4,real(h1)
    write(43,rec=loop) tr4,real(h2)

    write(34,rec=loop) tr4,real(g1)
    write(44,rec=loop) tr4,real(g2)

    write(35,rec=loop) tr4,real(z1)
    write(45,rec=loop) tr4,real(z2)

    write(51,'(f9.2,5(1x,f16.9))') tr4,ekiv,ekih,ekih+ekiv,epot,etot
  enddo

  close(30)
  close(40)
  close(31)
  close(41)
  close(32)
  close(42)
  close(33)
  close(43)
  close(34)
  close(44)
  close(35)
  close(45)
  close(51)

  write(*,*)
  write(*,*) ' The balanced fields are available in bx1.r4 & bx2.r4 where'
  write(*,*) ' x = d, h, g, z & q.  Also, becomp.asc contains the energies.'
  write(*,*)

endif


 !Internal subroutine definitions (inherit global variables):

contains 

!=======================================================================

subroutine balance

implicit none

 !Local variables:
double precision:: wka(ng,ng),wkb(ng,ng)
integer:: istep

!-----------------------------------------------------------------
! Initialise:
ds1=zero
ds2=zero
gs1=zero
gs2=zero
qs1=zero
qs2=zero

!-----------------------------------------------------------------
! Perform the entire time integration:
do istep=1,nstep
  t=dt*dble(istep-1)
  write(*,'(a,f8.6)') ' t/T = ',t/tramp
   !Evolve the flow from time t to t + dt:
  call advance
enddo
write(*,'(a,f8.6)') ' t/T = ',one

!-----------------------------------------------------------------
! Invert PV at final time:
call ptospc(ng,ng,qb1,qs1,xfactors,yfactors,xtrig,ytrig)
call ptospc(ng,ng,qb2,qs2,xfactors,yfactors,xtrig,ytrig)
call main_invert(qs1,qs2,ds1,ds2,gs1,gs2,h1,h2,u1,u2,v1,v2,q1,q2,z1,z2)
call spctop(ng,ng,ds1,d1,xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,ds2,d2,xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,gs1,g1,xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,gs2,g2,xfactors,yfactors,xtrig,ytrig)

! Compute energy components, divided by rho_1*H where H = H_1+H_2:
wka=hbar1*(one+h1)
wkb=hbar2*(one+h2)

! Horizontal and vertical parts of the kinetic energy:
ekih=ekmf*sum(wka*(u1**2+v1**2)+alpha*wkb*(u2**2+v2**2))
! ekmf = 0.5*glx*gly/(hbar_1+hbar_2)
ekiv=zero

! Potential energy:
epot=epmf*sum(wka**2+alpha*wkb*(two*wka+wkb)-dape)
! epmf = g*ekmf (see constants.f90)

! Total energy:
etot=ekih+epot+ekiv

return
end subroutine balance

!=======================================================================

subroutine advance

! Evolves delta and gamma by an iterative implicit method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}-F^n)/2] + N[(F^{n+1}-F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms, and N refers to the nonlinear source
! terms.  We start with a guess for F^{n+1} in N and iterate 
! niter times (see parameter statement below).

implicit none

! Local variables:
integer,parameter:: niter=2

! Spectral fields needed in time stepping:
double precision:: ds1i(ng,ng),sds1(ng,ng),nds1(ng,ng)
double precision:: ds2i(ng,ng),sds2(ng,ng),nds2(ng,ng)
double precision:: gs1i(ng,ng),sgs1(ng,ng),ngs1(ng,ng)
double precision:: gs2i(ng,ng),sgs2(ng,ng),ngs2(ng,ng)
double precision:: ttt1(ng,ng),ttt2(ng,ng)
! Other variables:
double precision:: ramp
integer:: iter

!-----------------------------------------------------------------------
! Invert PV, delta and gamma to obtain the velocity field:
call main_invert(qs1,qs2,ds1,ds2,gs1,gs2,h1,h2,u1,u2,v1,v2,q1,q2,z1,z2)
! Note: qsj, dsj & gsj are in spectral space while 
!       hj, uj, vj, qj and zj are in physical space (j = 1 & 2).

!------------------------------------------------------------------
! Start with a guess for F^{n+1} for all fields:

! Calculate the source terms (sds,sgs) for delta and h (ds,gs):
call source(sds1,sds2,sgs1,sgs2)

! Update divergence and acceleration divergence:
ds1i=ds1
gs1i=gs1
nds1=sds1+dt4i*ds1i
ngs1=sgs1+dt4i*gs1i
sds1=nds1+sds1             !2*N_tilde_delta_1
sgs1=ngs1+sgs1             !2*N_tilde_gamma_1

ds2i=ds2
gs2i=gs2
nds2=sds2+dt4i*ds2i
ngs2=sgs2+dt4i*gs2i
sds2=nds2+sds2             !2*N_tilde_delta_2
sgs2=ngs2+sgs2             !2*N_tilde_gamma_2

ttt1=rdis*sds1+sgs1
ttt2=rdis*sds2+sgs2
sds1=bb22*ttt1-bb12*ttt2   !2*bar_delta_1
sds2=bb11*ttt2-bb21*ttt1   !2*bar_delta_2

ds1=sds1-ds1i
ds2=sds2-ds2i
gs1=rdisi*(sgs1-gw11*sds1-gw12*sds2)-gs1i
gs2=rdisi*(sgs2-gw21*sds1-gw22*sds2)-gs2i

! Update current PV at t+dt:
ramp=ramp_fun(t+dt)
q1=ramp*qb1
call ptospc(ng,ng,q1,qs1,xfactors,yfactors,xtrig,ytrig)
q2=ramp*qb2
call ptospc(ng,ng,q2,qs2,xfactors,yfactors,xtrig,ytrig)

!------------------------------------------------------------------
! Iterate to improve estimates of F^{n+1}:
do iter=1,niter
  ! Perform inversion at t^{n+1} from estimated quantities:
  call main_invert(qs1,qs2,ds1,ds2,gs1,gs2,h1,h2,u1,u2,v1,v2,q1,q2,z1,z2)

  ! Calculate the source terms:
  call source(sds1,sds2,sgs1,sgs2)

  ! Update divergence and acceleration divergence:
  sds1=nds1+sds1
  sgs1=ngs1+sgs1

  sds2=nds2+sds2
  sgs2=ngs2+sgs2

  ttt1=rdis*sds1+sgs1
  ttt2=rdis*sds2+sgs2
  sds1=bb22*ttt1-bb12*ttt2
  sds2=bb11*ttt2-bb21*ttt1

  ds1=sds1-ds1i
  ds2=sds2-ds2i
  gs1=rdisi*(sgs1-gw11*sds1-gw12*sds2)-gs1i
  gs2=rdisi*(sgs2-gw21*sds1-gw22*sds2)-gs2i
enddo

return
end subroutine advance

!=======================================================================

subroutine source(sds1,sds2,sgs1,sgs2)

! Gets the source terms (sds,sgs) for delta and gamma (ds,gs) ---
! all in spectral space and in both layers (1 & 2).

! Note that (sds,sgs) only include the nonlinear terms for a 
! semi-implicit treatment, closely analogous to that described in 
! the appendix of Dritschel & Jalali (JFM, 2020).

! The spectral fields ds, gs, qd and qs are all spectrally truncated.
! Note, h, q, u & v obtained by main_invert before calling this routine
! are all spectrally truncated.

! The sources are *not* spectrally truncated as that occurs when they
! are used in subroutine advance.
  
implicit none

! Passed variables:
double precision:: sds1(ng,ng),sgs1(ng,ng)
double precision:: sds2(ng,ng),sgs2(ng,ng)

! Local variables:
double precision:: wka(ng,ng),wkb(ng,ng),wkp(ng,ng),wkq(ng,ng)
double precision:: wkf(ng,ng),wkg(ng,ng),wks(ng,ng)

!------------------------------------------------------------------
! ds source:
wka=ds1
call spctop(ng,ng,wka,wkf,xfactors,yfactors,xtrig,ytrig)
! wkf contains delta_1 in physical space
wkp=wkf*u1
wkq=wkf*v1
call divs(wkp,wkq,sds1)
! sds1 contains div(delta_1*u_1,delta_1*v_1) in spectral space
call jacob(u1,v1,wkp)
call ptospc(ng,ng,wkp,wkq,xfactors,yfactors,xtrig,ytrig)
! wkq contains J(u_1,v_1) in spectral space
sds1=two*wkq-sds1
! sds1 is now the full nonlinear ds source in layer 1
sds1(1,1)=zero ! ensures zero domain average

wka=ds2
call spctop(ng,ng,wka,wkg,xfactors,yfactors,xtrig,ytrig)
! wkg contains delta_2 in physical space
wkp=wkg*u2
wkq=wkg*v2
call divs(wkp,wkq,sds2)
! sds2 contains div(delta_2*u_2,delta_2*v_2) in spectral space
call jacob(u2,v2,wkp)
call ptospc(ng,ng,wkp,wkq,xfactors,yfactors,xtrig,ytrig)
! wkq contains J(u_2,v_2) in spectral space
sds2=two*wkq-sds2
! sds2 is now the full nonlinear ds source in layer 2
sds2(1,1)=zero ! ensures zero domain average

!-----------------------------------------------------------------
! gs source:
wkf=csq1*h1 ! Hydrostatic pressure in layer 1
wkg=csq2*h2 ! Hydrostatic pressure in layer 2

wka=wkf*u1  ! c_1^2*h_1*u_1
wkb=wkg*u2  ! c_2^2*h_2*u_2
wkf=wkf*v1  ! c_1^2*h_1*v_1
wkg=wkg*v2  ! c_2^2*h_2*v_2

! Compute div(zeta_1*u_1,zeta_1*v_1) and store in sgs1 (spectral):
wkp=z1*u1
wkq=z1*v1
call divs(wkp,wkq,sgs1)

! Form terms involving the layer depths:
wkp=wka+alpha*wkb
wkq=wkf+alpha*wkg
call divs(wkp,wkq,wks)
! wks = div(c_1^2*h_1*(u_1,v_1)+alpha*c_2^2*h_2*(u_2,v_2)) in spectral space

! Finalise gs source in layer 1:
sgs1=-cof*sgs1-rksq*wks
! Here -rksq = -k^2 = filtered Laplacian in spectral space
sgs1(1,1)=zero ! ensures zero domain average

! Compute div(zeta_2*u_2,zeta_2*v_2) and store in sgs2 (spectral):
wkp=z2*u2
wkq=z2*v2
call divs(wkp,wkq,sgs2)

! Form terms involving the layer depths:
wkp=wka+wkb
wkq=wkf+wkg
call divs(wkp,wkq,wks)
! wks = div(c_1^2*h_1*(u_1,v_1)+c_2^2*h_2*(u_2,v_2)) in spectral space

! Finalise gs source in layer 2:
sgs2=-cof*sgs2-rksq*wks
! Here -rksq = -k^2 = filtered Laplacian in spectral space
sgs2(1,1)=zero ! ensures zero domain average

return
end subroutine source

!=======================================================================

double precision function ramp_fun(time)

implicit none

double precision:: time

ramp_fun=f12-f12*cos(pi*time/tramp)

return
end function ramp_fun


 !End main program
end program fopvbal
!=======================================================================
