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
double precision:: u1(ng,ng),v1(ng,ng),z1(ng,ng),pn1(ng,ng)
double precision:: h2(ng,ng),d2(ng,ng),g2(ng,ng),q2(ng,ng)
double precision:: u2(ng,ng),v2(ng,ng),z2(ng,ng),pn2(ng,ng)
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
  open(30,file='q1.r4',form='unformatted', &
        access='direct',status='old',recl=nbytes)
  open(40,file='q2.r4',form='unformatted', &
        access='direct',status='old',recl=nbytes)

  open(31,file='bq1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(41,file='bq2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(32,file='bd1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(42,file='bd2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(33,file='bh1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(43,file='bh2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(34,file='bg1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(44,file='bg2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(35,file='bz1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(45,file='bz2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(36,file='bp1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(46,file='bp2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(51,file='becomp.asc',status='replace')

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

    write(36,rec=loop) tr4,real(pn1)
    write(46,rec=loop) tr4,real(pn2)

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
  close(36)
  close(46)
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
double precision:: htot1(ng,ng),htot2(ng,ng)
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

pn1=zero
pn2=zero
d1=zero
d2=zero

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
htot1=hbar1*(one+h1)
wka=d1*htot1

htot2=hbar2*(one+h2)
wkb=d2*htot2

! Horizontal and vertical parts of the kinetic energy:
ekih=ekmf**sum(htot1*(u1**2+v1**2)+alpha*htot2*(u2**2+v2**2))
ekiv=ekmf**sum(htot1*f13*wka**2+alpha*htot2*(wka**2+wka*wkb+f13*wkb**2))
! ekmf = 0.5*glx*gly/(hbar_1+hbar_2)

! Potential energy:
epot=epmf*sum(htot1**2+alpha*htot2*(two*htot1+htot2)-dape)
! epmf = g*ekmf (see constants.f90)

! Total energy:
etot=ekih+epot+ekiv

return
end subroutine balance

!=======================================================================

subroutine advance

! Evolves delta & gamma_l by an iterative implicit method of the form
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
if (t .gt. zero) then
  call source(sds1,sds2,sgs1,sgs2)
else
  sds1=zero
  sds2=zero
  sgs1=zero
  sgs2=zero
endif

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
t=t+dt
ramp=ramp_fun(t)
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
double precision:: htot1(ng,ng),htot2(ng,ng)
double precision:: hinv1(ng,ng),hinv2(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wkf(ng,ng),wkg(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng)

!----------------------------------------------------------------
! Find the scaled vertically-averaged non-hydrostatic pressure
! (pn1,pn2) in each layer (note: pnj = p_tilde_j/H_j^2):
call nhpsolve(wkf,wkg,wka,wkb)
! On return, d1,d2 = divergence in layers 1,2 in physical space,
! wkf,wkg = -xi_1,-xi_2 where xi_j = D(delta_j)/Dt-delta_j^2,
! while wka,wkb = J_1,J_2 (non-hydrostatic part of zeta tendency).

!---------------------------------------------------------------
! Compute next the nonlinear part of delta source in each layer.

! ==>  Layer 1  <==
! Compute div(delta_1*u_1,delta_1*v_1) and store in sds1 (spectral):
wkp=d1*u1
wkq=d1*v1
call divs(wkp,wkq,sds1)

! Compute 2*delta_1^2 + xi_1 and store in spectral space as wkq:
wkp=two*d1**2-wkf
call ptospc(ng,ng,wkp,wkq,xfactors,yfactors,xtrig,ytrig)
! *** wkf is now free to re-use

! Finalise delta_1 source, N_delta_1 (de-aliased):
sds1=wkq-sds1-pm11*gs1-pm12*gs2
! Note, the linear terms are subtracted here but they largely cancel
! similar terms in wkq-sds1 on rhs.
sds1(1,1)=zero ! ensures zero domain average

! ==>  Layer 2  <==
! Compute div(delta_2*u_2,delta_2*v_2) and store in sds2 (spectral):
wkp=d2*u2
wkq=d2*v2
call divs(wkp,wkq,sds2)

! Compute 2*delta_2^2 + xi_2 and store in spectral space as wkq:
wkp=two*d2**2-wkg
call ptospc(ng,ng,wkp,wkq,xfactors,yfactors,xtrig,ytrig)
! *** wkg is now free to re-use

! Finalise delta_2 source, N_delta_2 (de-aliased):
sds2=wkq-sds2-pm21*gs1-pm22*gs2
! Note, the linear terms are subtracted here but they largely cancel
! similar terms in wkq-sds2 on rhs.
sds2(1,1)=zero ! ensures zero domain average

!-----------------------------------------------------------------
! Compute next the nonlinear part of gamma_l source in each layer.

wkf=csq1*h1 ! Hydrostatic pressure in layer 1
wkg=csq2*h2 ! Hydrostatic pressure in layer 2

! Re-use the arrays htotj & hinvj for the fluxes needed below:
htot1=wkf*u1
htot2=wkg*u2
hinv1=wkf*v1
hinv2=wkg*v2

! ==>  Layer 1  <==
! Compute div(zeta_1*u_1,zeta_1*v_1) and store in sgs1 (spectral):
wkp=z1*u1
wkq=z1*v1
call divs(wkp,wkq,sgs1)

! Combine with J_1 (in wka) from nhpsolve above:
call ptospc(ng,ng,wka,wkq,xfactors,yfactors,xtrig,ytrig)
sgs1=cof*(wkq-sgs1)
! *** wka is now free to re-use

! Form terms involving the layer depths:
wkp=htot1+alpha*htot2
wkq=hinv1+alpha*hinv2
call divs(wkp,wkq,wka)
! wka = div(c_1^2*h_1*(u_1,v_1)+alpha*c_2^2*h_2*(u_2,v_2)) in spectral space

! Finalise gamma_l source in layer 1:
sgs1=sgs1-rksq*wka
! Here -rksq = -k^2 = filtered Laplacian in spectral space
sgs1(1,1)=zero ! ensures zero domain average

! ==>  Layer 2  <==
! Compute div(zeta_2*u_2,zeta_2*v_2) and store in sgs2 (spectral):
wkp=z2*u2
wkq=z2*v2
call divs(wkp,wkq,sgs2)

! Combine with J_2 (in wkb) from nhpsolve above:
call ptospc(ng,ng,wkb,wkq,xfactors,yfactors,xtrig,ytrig)
sgs2=cof*(wkq-sgs2)
! *** wkb is now free to re-use

! Form terms involving the layer depths:
wkp=htot1+htot2
wkq=hinv1+hinv2
call divs(wkp,wkq,wkb)
! wkb = div(c_1^2*h_1*(u_1,v_1)+c_2^2*h_2*(u_2,v_2)) in spectral space

! Finalise gamma_l source in layer 2:
sgs2=sgs2-rksq*wkb
! Here -rksq = -k^2 = filtered Laplacian in spectral space
sgs2(1,1)=zero ! ensures zero domain average

return
end subroutine source

!=======================================================================
subroutine nhpsolve(s1,s2,wka,wkb)

! Finds the scaled vertically-averaged non-hydrostatic pressure
! pnj = H_j^{-2}*bar{P}_{nj}/h_j, in layers j = 1 & 2.
! On return, dj = divergence in layer j, sj = -xi_j =
! -[D(delta_j)/Dt - delta_j^2], while (wka,wkb) = (J_1,J_2) =
! non-hydrostatic part of zeta_j tendency.

! *** All passed variables are in physical space ***
  
implicit none

! Passed variables:
double precision:: s1(ng,ng),s2(ng,ng),wka(ng,ng),wkb(ng,ng)

! Local variables:
double precision,parameter:: ptol=1.d-7
! ptol: maximum relative rms NH pressure error

double precision:: htot1(ng,ng),htot2(ng,ng),hinv1(ng,ng),hinv2(ng,ng)
double precision:: h1x(ng,ng),h1y(ng,ng),tt(ng,ng)
double precision:: r1x(ng,ng),r1y(ng,ng),t1x(ng,ng),t1y(ng,ng)
double precision:: r2x(ng,ng),r2y(ng,ng),sgx(ng,ng),sgy(ng,ng)
double precision:: ht1x(ng,ng),ht1y(ng,ng)
double precision:: ch1(ng,ng),ch2(ng,ng),ch3(ng,ng)
double precision:: wkc(ng,ng),wkd(ng,ng),wkp(ng,ng)
double precision:: perr

!---------------------------------------------------------------
! Get total dimensionless layer thicknesses (1 + h_j):
htot1=one+h1
htot2=one+h2

! Get their inverses:
hinv1=one/htot1
call dealias(hinv1)
hinv2=one/htot2
call dealias(hinv2)

! Define T = 6/(4 + 3*mu) where mu = alpha*(H_2/H_1)*(1+h_2)/(1+h_1):
wkp=mubar*htot2*hinv1
call dealias(wkp)
tt=six/(four+three*wkp)
call dealias(tt)

! Define c_hat_1,2,3:
wkp=three*tt*hinv1
call dealias(wkp)
ch1=wkp*hinv1-cona1
call dealias(ch1)
ch2=hrati*wkp*hinv2-cona2
call dealias(ch2)
wkd=tt*hinv2          !wkd = T/(1+h_2): *** preserve to define sgx & sgy
call dealias(wkd)
ch3=two*wkd*(hinv2+mubar3*hinv1)-cona3   !mubar3 = 3*mubar
call dealias(ch3)

! Calculate (h1x,h1y) = grad{h_1} and (r1x,r1y) = H_1^2*grad{h_1}/(1+h_1):
wkp=h1
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,h1x,h1y)
r1x=hbsq1*hinv1*h1x
call dealias(r1x)
r1y=hbsq1*hinv1*h1y
call dealias(r1y)

! Calculate (t1x,t1y) = tau = T*(r1x,r1y) = H_1^2*T*grad{h_1}/(1+h_1):
t1x=tt*r1x
call dealias(t1x)
t1y=tt*r1y
call dealias(t1y)

! Define tau/2 for efficiency in pressure iteration below:
ht1x=f12*t1x
ht1y=f12*t1y

! Calculate (r2x,r2y) = H_2^2*grad{h_2}/(1+h_2):
wkp=hbsq2*h2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,r2x,r2y)
r2x=hinv2*r2x
call dealias(r2x)
r2y=hinv2*r2y
call dealias(r2y)

! Calculate (sgx,sgy) = (r2x,r2y) + H_1*H_2*T*grad{h_1}/(1+h_2):
wka=wkd*h1x
call dealias(wka)
sgx=r2x+hb1hb2*wka
wka=wkd*h1y
call dealias(wka)
sgy=r2y+hb1hb2*wka
! Note: here wkd contains T/(1+h_2) from above; wkd can now be re-used

!--------------------------------------------------------------------
! Form the fixed parts (s1,s2) of the sources needed to calculate the
! non-hydrostatic pressure:

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute gamma_tilde = gamma_l + 2*(J(u,v) - delta^2) in layer 1:
wka=ds1
call spctop(ng,ng,wka,d1,xfactors,yfactors,xtrig,ytrig)
! d1 contains delta_1 in physical space
call jacob(u1,v1,wkp)
! wkp contains J(u_1,v_1) in physical space
wkp=wkp-d1**2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
wka=gs1+two*filt*wka
call spctop(ng,ng,wka,s1,xfactors,yfactors,xtrig,ytrig)
! s1 now contains gamma_tilde_1 in physical space (de-aliased)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute gamma_tilde = gamma_l + 2*(J(u,v) - delta^2) in layer 2:
wka=ds2
call spctop(ng,ng,wka,d2,xfactors,yfactors,xtrig,ytrig)
! d2 contains delta_2 in physical space
call jacob(u2,v2,wkp)
! wkp contains J(u_2,v_2) in physical space
wkp=wkp-d2**2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
wka=gs2+two*filt*wka
call spctop(ng,ng,wka,s2,xfactors,yfactors,xtrig,ytrig)
! s2 now contains gamma_tilde_2 in physical space (de-aliased)

!======================================================================
! Iterate to find approximate solution:

perr=one
do while (perr .gt. ptol)

  ! Find full S_1 (rhs of 1st pressure equation) in spectral space (wka):
  wkp=f23*pn1-ahrsq*pn2     !P = (2/3)*pn1 - alpha*(H_2/H_1)^2*pn2
  wkc=wkp*t1x
  wkd=wkp*t1y
  call divs(wkc,wkd,wka)    !wka = div{P*tau} (spectral)
  wkp=s1+ch1*wkp            !wkp = gamma_tilde_1 + c_hat_1*P
  call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=wkc-wka               !wka = gamma_tilde_1 + c_hat_1*P - div{P*tau}

  ! Find full S_2 (rhs of 2nd pressure equation) in spectral space (wkb):
  wkc=pn1*ht1x+pn2*sgx
  wkd=pn1*ht1y+pn2*sgy
  call divs(wkc,wkd,wkb)    !wkb = div{(1/2)*pn1*tau+pn2*sigma} (spectral)
  wkp=s2-ch2*pn1+ch3*pn2    !wkp = gamma_tilde_2 - c_hat_2*pn1 + c_hat_3*pn2
  call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
  wkb=wkc-wkb               !wkb = gamma_tilde_2 - c_hat_2*pn1 + c_hat_3*pn2
                            !     -div{(1/2)*pn1*tau+pn2*sigma}

  ! Obtain next approximation for pn1 & pn2 (hold in wkc & wkd temporarily):
  wkp=pi11*wka+pi12*wkb
  call spctop(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)

  wkp=pi21*wka+pi22*wkb
  call spctop(ng,ng,wkp,wkd,xfactors,yfactors,xtrig,ytrig)

  ! Compute relative rms difference error:
  perr=sqrt(sum((wkc-pn1)**2+(wkd-pn2)**2)/sum(pn1**2+pn2**2))

! Copy wkc & wkd into pn1 & pn2 with relaxation (essential!):
  pn1=f12*(pn1+wkc)
  pn2=f12*(pn2+wkd)
  ! From tests, this 50/50 blend of old and new solutions leads to
  ! exponential convergence (50% error reduction each iteration).
enddo

!----------------------------------------------------------
! Compute -xi_1 = -[D(delta_1)/Dt - delta_1^2] -> s1
!     and -xi_2 = -[D(delta_2)/Dt - delta_2^2] -> s2:
wkp=f23*pn1-ahrsq*pn2
s1=(ch1+cona1)*wkp
s2=(ch3+cona3)*pn2-(ch2+cona2)*pn1

! Compute non-hydrostatic parts of zeta_1 tendency (-> wka):
wkp=tt*wkp
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call gradient(wkb,wkc,wkd)
wka=r1x*wkd-r1y*wkc                     ! wka = J_1

! Compute non-hydrostatic parts of zeta_2 tendency (-> wkb):
wkp=pn2
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
call gradient(wkb,sgx,sgy)

wkp=pn2*hinv2+hhrati*pn1*hinv1          ! pn2/(1+h_2)+(H_1/2*H_2)*pn1/(1+h_1)
call dealias(wkp)
wkp=tt*wkp                              ! wkp = kappa in the notes
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call gradient(wkb,wkc,wkd)              ! (wkc,wkd) = grad{kappa}
wkb=r2x*sgy-r2y*sgx+hb1hb2*(h1x*wkd-h1y*wkc)   ! wkb = J_2

!------------------------------------------------------------------
wkp=pn1*(one+h1)

return
end subroutine

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
