module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

 !Various PV fields:
double precision:: qt1(numLatGridPts,numLongGridPts),qc1(numLatGridPts,numLongGridPts),qd1(numLatGridPts,numLongGridPts)
double precision:: qt2(numLatGridPts,numLongGridPts),qc2(numLatGridPts,numLongGridPts),qd2(numLatGridPts,numLongGridPts)

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================
subroutine advect

! Main subroutine for advecting fields and contours

implicit none

 !Local variables:
double precision,parameter:: twistmax=2.5d0
!      twistmax: the maximum value of the time integral of |zeta|_max
!                between regularisations of the contours.
integer,parameter:: nregmax=20
!      Every nregmax contour regularisations, the code stops
!      to rebuild the contours in a separate memory space.
integer:: ireg,itime,jtime
logical:: ggen

!-----------------------------------------------------------------------
 !Define fixed arrays and constants and read initial data:
call init

 !Used for regularising contours:
twist=zero

 !Counter used for counting number of contour regularisations done:
ireg=0

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Start the time loop:
do while (t .le. tsim)

   !Save data periodically:
  itime=nint(t/dt)
  jtime=itime/gridSaveInterval
  if (gridSaveInterval*jtime .eq. itime) then
    call inversion
    call savegrid(jtime+1)
    ggen=.false.
  else
    ggen=.true.
  endif
   !ggen is used to indicate if calling inversion is needed in advance below
  jtime=itime/contourSaveInterval
  if (contourSaveInterval*jtime .eq. itime) call savecont(jtime+1)

   !Perform contour surgery or recontouring when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1

     !Don't continue if maximum number of regularisations reached:
    if (ireg .eq. nregmax) then
       !Prepare PV residual qResid for recontouring (and preserve qs):
      call prepare
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the PV contours (surgery + node redistribution):
    call surgery
     !Record active contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i9,1x,i10)') t,numContoursTot,numNodesTot

    twist=twist-twistmax
  endif

   !Advect flow from time t to t + dt:
  call advance(ggen)
  
enddo
!End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Possibly save final data:
itime=nint(t/dt)
jtime=itime/gridSaveInterval
if (gridSaveInterval*jtime .eq. itime) then
  call inversion
  call savegrid(jtime+1)
endif
jtime=itime/contourSaveInterval
if (contourSaveInterval*jtime .eq. itime) call savecont(jtime+1)

return
end subroutine advect

!=======================================================================

subroutine init

! Initialises residual PV for normal time integration following 
! contour regeneration. Note qSpec1 & qSpec2 contain the semi-spectral
! PV anomaly in each layer.

implicit none

! Local variable:
double precision:: qa(numLatGridPts,numLongGridPts,numLayers)

!------------------------------------------------------------------
 !Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,numContoursTot,numNodesTot

!------------------------------------------------------------
 !Convert PV contours (x,y,z) to gridded PV anomaly as qa:
call con2grid(qa)

 !Convert each component to semi-spectral space as qc1 & qc2:
qc1=qa(:,:,1)
call forfft(numLatGridPts,numLongGridPts,qc1,trig,factors) 
qc2=qa(:,:,2)
call forfft(numLatGridPts,numLongGridPts,qc2,trig,factors) 

 !Define (semi-spectral) residual PV qd = (1-F)*(qs-qc) in each layer:
qd1=qSpec1-qc1
qd2=qSpec2-qc2

call hipass(qd1)
call hipass(qd2)

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current PV anomaly field is stored in qs (semi-spectral), and
! the residual PV field needed in congen.f90 will be stored in qResid.

implicit none

! Local variables:
double precision:: zpf(numLatGridPts,numLongGridPts),avgzpf
integer:: i,iz

!------------------------------------------------------------
 !Re-initialise qs & qd before recontouring:
 !          Reset  qs = F*(qs-qc) + qc + qd
 !            and  qd = (1-F)*(qs-qc)
 !Here F & (1-F) are low & high pass filter functions (see spectral.f90)
call reset
 !Note, qc is obtained in reset by a call to con2grid.

 !Obtain full PV anomaly on the grid as qResid:
qResid(:,:,1)=qSpec1
call revfft(numLatGridPts,numLongGridPts,qResid(:,:,1),trig,factors)
qResid(:,:,2)=qSpec2
call revfft(numLatGridPts,numLongGridPts,qResid(:,:,2),trig,factors)

 !Add f to define the total PV, including the planetary part:
do iz=1,numLayers
  do i=1,numLongGridPts
    qResid(:,i,iz)=qResid(:,i,iz)+cof
  enddo
enddo

 !Correct average PV in each layer by enforcing zero average vorticity:
zpf=qResid(:,:,1)*(one+h1)
avgzpf=average(zpf)
qResid(:,:,1)=qResid(:,:,1)-avgzpf

zpf=qResid(:,:,2)*(one+h2)
avgzpf=average(zpf)
qResid(:,:,2)=qResid(:,:,2)-avgzpf

 !Obtain unique polar values of qResid for use in recontouring:
do iz=1,numLayers
  qSPole(iz) = sum(qResid( 1,:,iz))/dble(numLongGridPts)
  qNPole(iz) = sum(qResid(numLatGridPts,:,iz))/dble(numLongGridPts)
enddo
 !These are averages of points half a latitude spacing from the poles

 !Obtain PV anomaly due to contours (qc) on the grid:
call revfft(numLatGridPts,numLongGridPts,qc1,trig,factors)
call revfft(numLatGridPts,numLongGridPts,qc2,trig,factors)
 !Note: qc1 & qc2 were obtained in reset and are in this module's commons.

 !Define residual PV (qResid) needed for recontouring:
do i=1,numLongGridPts
  qResid(:,i,1)=qResid(:,i,1)-qc1(:,i)-cof
  qResid(:,i,2)=qResid(:,i,2)-qc2(:,i)-cof
enddo

 !Note: We are leaving this module after this; qd will be redefined
 !      upon re-entry in subroutine init.

return
end subroutine prepare

!=======================================================================

subroutine advance(ggen)

! Advances PV from time t to t+dt by a combination of contour 
! advection (for PV contours) and the pseudo-spectral method 
! (for all semi-spectral fields, namely qs, qd, ds & gs).

! Uses an iterative implicit method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}-F^n)/2] + N[(F^{n+1}-F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms, and N refers to the nonlinear source
! terms.  We start with a guess for F^{n+1} in N and iterate 
! niter times (see parameter statement below).

implicit none

 !Passed variable:
logical:: ggen

 !Local variables:
integer,parameter:: niter=2

 !Semi-spectral fields needed in time stepping (layers 1 and 2):
double precision:: qsi1(numLatGridPts,numLongGridPts),qsi2(numLatGridPts,numLongGridPts)
double precision:: sqs1(numLatGridPts,numLongGridPts),sqs2(numLatGridPts,numLongGridPts)

double precision:: qdi1(numLatGridPts,numLongGridPts),qdi2(numLatGridPts,numLongGridPts)
double precision:: qdm1(numLatGridPts,numLongGridPts),qdm2(numLatGridPts,numLongGridPts)
double precision:: sqd1(numLatGridPts,numLongGridPts),sqd2(numLatGridPts,numLongGridPts)

double precision:: dsi1(numLatGridPts,numLongGridPts),dsi2(numLatGridPts,numLongGridPts)
double precision:: nds1(numLatGridPts,numLongGridPts),nds2(numLatGridPts,numLongGridPts)
double precision:: sds1(numLatGridPts,numLongGridPts),sds2(numLatGridPts,numLongGridPts)

double precision:: gsi1(numLatGridPts,numLongGridPts),gsi2(numLatGridPts,numLongGridPts)
double precision:: ngs1(numLatGridPts,numLongGridPts),ngs2(numLatGridPts,numLongGridPts)
double precision:: sgs1(numLatGridPts,numLongGridPts),sgs2(numLatGridPts,numLongGridPts)

 !Initial contour nodes and velocities needed in time stepping:
double precision:: xi(numNodesTot),yi(numNodesTot),zi(numNodesTot),u(numNodesTot),v(numNodesTot),w(numNodesTot)

 !Other local quantities:
double precision:: xtmp,ytmp,ztmp,fac
integer:: i,iter

!------------------------------------------------------------------
 !Invert PV and compute velocity at current time level, say t=t^n:
if (ggen) call inversion
 !If ggen is false, inversion was called previously at this time level.
 !A zero argument above means qc does not need to be regenerated.

 !Compute twist parameter and save various diagnostics each time step:
call diagnose

 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset  qs = F*(qs-qc) + qc + qd
 !            and  qd = (1-F)*(qs-qc)
 !Here F & (1-F) are low & high pass filter functions (see spectral.f90)
call reset

!--------------------------------------------------------------------
 !Start by obtaining a guess for F^{n+1} for all contours and fields:

 !Contours:
call velint(u1,u2,v1,v2,u,v,w)
 !Here (u,v,w) stands for the Cartesian velocity on the contours
 !at time level n, e.g. u(x^n,t^n).

 !Loop over nodes i:
do i=1,numNodesTot
   !Store x^n+0.5*dt*u(x^n,t^n) for efficiency in the iteration below:
  xtmp=x(i)+dt2*u(i)
  ytmp=y(i)+dt2*v(i)
  ztmp=z(i)+dt2*w(i)
   !Ensure point remains on the unit sphere:
  fac=one/sqrt(xtmp**2+ytmp**2+ztmp**2)
  xi(i)=fac*xtmp
  yi(i)=fac*ytmp
  zi(i)=fac*ztmp
   !Initial guess for x^{n+1}:
  xtmp=x(i)+dt*u(i)
  ytmp=y(i)+dt*v(i)
  ztmp=z(i)+dt*w(i)
   !Ensure point remains on the unit sphere:
  fac=one/sqrt(xtmp**2+ytmp**2+ztmp**2)
  x(i)=fac*xtmp
  y(i)=fac*ytmp
  z(i)=fac*ztmp
enddo

 !Store delta and gamma at the current time level "n" after projection
 !onto vertical modes in dsi & gsi:
dsi1=deltaSpec1*vm11+deltaSpec2*vm12
dsi2=deltaSpec1*vm21+deltaSpec2*vm22

gsi1=gammaSpec1*vm11+gammaSpec2*vm12
gsi2=gammaSpec1*vm21+gammaSpec2*vm22
 !dsi & gsi remain fixed throughout the iteration below (do not re-use).

 !Calculate the source terms for PV (sqs,sqd) as well as those for
 !divergence and acceleration divergence (sds,sgs) for each *layer*:
call source(sqs1,sqs2,sqd1,sqd2,sds1,sds2,sgs1,sgs2)
 !It is now safe to re-use ds & gs below. At the end of each iteration
 !step, ds & gs contain improved estimates for delta^{n+1} & gamma^{n+1}
 !for each *layer*.

 !Store PV field quantities that remain fixed (do not re-use):
qdi1=qd1
qdi2=qd2
qdm1=qd1+dt4*sqd1
qdm2=qd2+dt4*sqd2

qsi1=qSpec1+dt2*sqs1
qsi2=qSpec2+dt2*sqs2 
 !Note: the PV evolution does not couple the layers.

 !Initial guess for qd & qs at time level "n+1":
qd1=diss*(qdm1+dt4*sqd1)-qdi1
qd2=diss*(qdm2+dt4*sqd2)-qdi2

qSpec1=qsi1+dt2*sqs1
qSpec2=qsi2+dt2*sqs2

 !Store arrays needed for updating delta and gamma (in mode space):
nds1=sds1*vm11+sds2*vm12+dt4i*dsi1 !N_delta^n + (4/dt)*delta^n (mode 1)
nds2=sds1*vm21+sds2*vm22+dt4i*dsi2 !N_delta^n + (4/dt)*delta^n (mode 2)

ngs1=sgs1*vm11+sgs2*vm12+dt4i*gsi1 !N_gamma^n + (4/dt)*gamma^n (mode 1)
ngs2=sgs1*vm21+sgs2*vm22+dt4i*gsi2 !N_gamma^n + (4/dt)*gamma^n (mode 2)
 !These remain fixed throughout the iteration below (do not re-use).

 !Obtain estimates for ds (delta) and gs (gamma) at time level "n+1":
call dgstep(dsi1,dsi2,gsi1,gsi2,nds1,nds2,ngs1,ngs2,sds1,sds2,sgs1,sgs2)
 !The first 8 arrays are *fixed*, the last 4 are modified in dgstep.
 !Here, on entry, sds = N_delta^n & sgs = N_gamma^n for each *layer*.
 !On return, ds & gs (in commons) contain first estimates for
 !delta^{n+1} & gamma^{n+1} for each *layer*.

 !If there is baroclinic PV forcing (qForcing), update its form:
if (forcing) call pvforce
 !Note: the forcing stays constant in time within a time step.

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call inversion

   !Calculate the source terms for PV (sqs,sqd) as well as those for
   !divergence and acceleration divergence (sds,sgs) for each *layer*:
  call source(sqs1,sqs2,sqd1,sqd2,sds1,sds2,sgs1,sgs2)
   !It is now safe to re-use ds & gs below. At the end of each iteration
   !step, ds & gs contain improved estimates for delta^{n+1} & gamma^{n+1}
   !for each *layer*.

   !Update contours:
  call velint(u1,u2,v1,v2,u,v,w)
  do i=1,numNodesTot
    xtmp=xi(i)+dt2*u(i)
    ytmp=yi(i)+dt2*v(i)
    ztmp=zi(i)+dt2*w(i)
    fac=one/sqrt(xtmp**2+ytmp**2+ztmp**2)
    x(i)=fac*xtmp
    y(i)=fac*ytmp
    z(i)=fac*ztmp
  enddo
   !Now (x,y,z) contain a new guess for x^{n+1}, etc.

   !Update PV fields:
  qd1=diss*(qdm1+dt4*sqd1)-qdi1
  qd2=diss*(qdm2+dt4*sqd2)-qdi2

  qSpec1=qsi1+dt*sqs1
  qSpec2=qsi2+dt*sqs2

   !Update delta and gamma:
  call dgstep(dsi1,dsi2,gsi1,gsi2,nds1,nds2,ngs1,ngs2,sds1,sds2,sgs1,sgs2)
   !The first 8 arrays are *fixed*, the last 4 are modified in dgstep.
   !On entry, sds = N_delta^{n+1} & sgs = N_gamma^{n+1} for each *layer*.
   !On return, ds & gs (in commons) contain improved estimates for
   !delta^{n+1} & gamma^{n+1} for each *layer*.
enddo

 !Apply latitudinal hyperviscous damping to qd, ds & gs:
call latdamp(qd1)
call latdamp(qd2)

call latdamp(deltaSpec1)
call latdamp(deltaSpec2)

call latdamp(gammaSpec1)
call latdamp(gammaSpec2)
 !The longitudinal part is incorporated above in diss, rdis and simp
 !for ds & gs (see dgstep and subroutine init_spectral in spectral.f90).

 !Advance time:
t=t+dt

return
end subroutine advance

!=======================================================================

subroutine dgstep(dsi1,dsi2,gsi1,gsi2,nds1,nds2,ngs1,ngs2,sds1,sds2,sgs1,sgs2)
 !Updates the divergence (delta) and the acceleration divergence (gamma)
 !at the next time level (n+1) given their values dsi & gsi at the
 !current time level (n), the fixed arrays nds & ngs containing
 !N_delta^n + (4/dt)*delta^n & N_gamma^n + (4/dt)*gamma^n, and the
 !source terms sds = N_delta^{n+1} & sgs = N_gamma^{n+1}.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !***Note: On entry, sds & sgs contain the sources for each *layer* and
 !   are converted to modes here. The arrays dsi & gsi contain delta^n
 !   & gamma^n for each *mode*, while nds & ngs are the fixed parts of
 !   the source terms for delta and gamma at time level n, again for
 !   each *mode*.  ===> These 8 arrays should *not* be modified. <===
 !   On return, ds & gs (in commons) contain delta^{n+1) & gamma^{n+1}
 !   for each *layer*.
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 !Uses a semi-implicit treatment of the equations, which in mode space
 !are delta_t = gamma + N_delta and gamma_t = G(delta) + N_gamma, where
 !G = c^2*Lap-f^2 is the gravity-wave operator (Lap = the Laplacian).
 !For any variable F, we replace F_t by (F^{n+1} - F^n)/dt, and the
 !terms on the right-hand side by their averages, (F^{n+1} + F^n)/2.

implicit none

 !Passed variables:
double precision:: dsi1(numLatGridPts,numLongGridPts),dsi2(numLatGridPts,numLongGridPts)
double precision:: nds1(numLatGridPts,numLongGridPts),nds2(numLatGridPts,numLongGridPts)
double precision:: sds1(numLatGridPts,numLongGridPts),sds2(numLatGridPts,numLongGridPts)

double precision:: gsi1(numLatGridPts,numLongGridPts),gsi2(numLatGridPts,numLongGridPts)
double precision:: ngs1(numLatGridPts,numLongGridPts),ngs2(numLatGridPts,numLongGridPts)
double precision:: sgs1(numLatGridPts,numLongGridPts),sgs2(numLatGridPts,numLongGridPts)

!---------------------------------------------------------------
 !First convert layer sources in (sds,sgs) to modes as (ds,gs)
 !in commons; it is safe to do this as ds & gs are available for
 !use, and they are overwritten here on return with the solution
 !for delta^{n+1} and gamma^{n+1} for each *layer*:
deltaSpec1=sds1*vm11+sds2*vm12
deltaSpec2=sds1*vm21+sds2*vm22

gammaSpec1=sgs1*vm11+sgs2*vm12
gammaSpec2=sgs1*vm21+sgs2*vm22

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !Update delta^{n+1} & gamma^{n+1} for mode 1:
sds1=nds1+deltaSpec1            !N_delta^{n+1} + N_delta^n + (4/dt)*delta^n
sgs1=ngs1+gammaSpec1            !N_gamma^{n+1} + N_gamma^n + (4/dt)*gamma^n
gammaSpec1=rdis*sds1+sgs1       !Source for delta^n + delta^{n+1}
call simp1(gammaSpec1,sgs1,deltaSpec1) !sgs = delta^{n} + delta^{n+1}; ds is not used
deltaSpec1=sgs1-dsi1            !delta^{n+1}
gammaSpec1=rdis*sgs1-sds1-gsi1  !gamma^{n+1}

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !Update delta^{n+1} & gamma^{n+1} for mode 2:
sds2=nds2+deltaSpec2            !N_delta^{n+1} + N_delta^n + (4/dt)*delta^n
sgs2=ngs2+gammaSpec2            !N_gamma^{n+1} + N_gamma^n + (4/dt)*gamma^n
gammaSpec2=rdis*sds2+sgs2       !Source for delta^n + delta^{n+1}
call simp2(gammaSpec2,sgs2,deltaSpec2) !sgs = delta^{n} + delta^{n+1}; ds is not used
deltaSpec2=sgs2-dsi2            !delta^{n+1}
gammaSpec2=rdis*sgs2-sds2-gsi2  !gamma^{n+1}

!-------------------------------------------------------------
 !Finally convert mode quantities (deltaSpec1,deltaSpec2,gammaSpec1,gammaSpec2) to layers:
sds1=deltaSpec1*vt11+deltaSpec2*vt12
deltaSpec2 =deltaSpec1*vt21+deltaSpec2*vt22
deltaSpec1=sds1

sgs1=gammaSpec1*vt11+gammaSpec2*vt12
gammaSpec2 =gammaSpec1*vt21+gammaSpec2*vt22
gammaSpec1=sgs1
 !Now gs & gs contain delta^{n+1} & gamma^{n+1} for each *layer*

 !Above, subroutines simpj (for each mode j) invert R^2-G to obtain
 !delta^{n+1}+delta^n, where R = 2/dt + D (= rdis in the code; here
 !D is the longitudinal part of the hyperdiffusion) and
 !G = c^2*Lap - f^2.

return
end subroutine dgstep

!=======================================================================

subroutine source(sqs1,sqs2,sqd1,sqd2,sds1,sds2,sgs1,sgs2)

! Gets the source terms (sqs,sqd) for the PV (qs,qd) as well as those
! (sds,sgs) for divergence and acceleration divergence (ds,gs) ---
! --- all in semi-spectral space and for each *layer*.

! Note that (sds,sgs) only include the nonlinear terms needed for a
! semi-implicit treatment, analogous to that described in the appendix
! of Mohebalhojeh & Dritschel (2004). See also the comments in
! subroutine dgstep.

! The semi-spectral fields ds, gs, qd and qs are all spectrally truncated.
! Note, h, u, v & z obtained by main_invert before calling this routine
! are also all spectrally truncated.

implicit none

 !Passed variables:
double precision:: sqs1(numLatGridPts,numLongGridPts),sqs2(numLatGridPts,numLongGridPts),sqd1(numLatGridPts,numLongGridPts),sqd2(numLatGridPts,numLongGridPts)
double precision:: sds1(numLatGridPts,numLongGridPts),sds2(numLatGridPts,numLongGridPts),sgs1(numLatGridPts,numLongGridPts),sgs2(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: wka(numLatGridPts,numLongGridPts),wkb(numLatGridPts,numLongGridPts),wkc(numLatGridPts,numLongGridPts),wkd(numLatGridPts,numLongGridPts)
double precision:: wke(numLatGridPts,numLongGridPts),wkf(numLatGridPts,numLongGridPts),wkg(numLatGridPts,numLongGridPts),wkh(numLatGridPts,numLongGridPts)
double precision:: om1(numLatGridPts,numLongGridPts),om2(numLatGridPts,numLongGridPts)
double precision:: usq1(numLatGridPts,numLongGridPts),usq2(numLatGridPts,numLongGridPts)
double precision:: d1(numLatGridPts,numLongGridPts),d2(numLatGridPts,numLongGridPts),rhs(numLatGridPts),avgval
integer:: i

!--------------------------------------------------------------
 !First define Omega = u/r (where r = cos(lat)) for each layer:
do i=1,numLongGridPts
  om1(:,i)=secantLat*u1(:,i)
  om2(:,i)=secantLat*u2(:,i)
enddo
 !*** Do not re-use om1 & om2 below ***

!---------------------------------------------------------------
 !qd source (can be computed in each layer independently):

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain qd source for layer 1 (sqd1):
call deriv(numLatGridPts,numLongGridPts,rk,qd1,wkg)
call revfft(numLatGridPts,numLongGridPts,wkg,trig,factors) ! wkg = d(qd1)/d(lon)
wka=qd1
call revfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
call latder(wka,wkh)                ! wkh = d(qd1)/d(lat)
sqd1=-om1*wkg-v1*wkh

 !Add baroclinic PV forcing (qForcing) if present:
if (forcing) sqd1=sqd1-mrat2*qForcing
 !mrat2 = m_2/(m_1+m_2) where m_j is the mass of layer j (see constants.f90)

 !Add non-conservative term associated with thermal damping if present:
if (thermal) then
  do i=1,numLongGridPts
    sqd1(:,i)=sqd1(:,i)+thermDampRate*(qAnom1(:,i)+cof)*(one-one/(one+h1(:,i)))
  enddo
endif

 !Convert to semi-spectral space and de-alias:
call dealias(sqd1)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain qd source for layer 2 (sqd2):
call deriv(numLatGridPts,numLongGridPts,rk,qd2,wkg)
call revfft(numLatGridPts,numLongGridPts,wkg,trig,factors) ! wkg = d(qd2)/d(lon)
wka=qd2
call revfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
call latder(wka,wkh)                ! wkh = d(qd2)/d(lat)
sqd2=-om2*wkg-v2*wkh

 !Add baroclinic PV forcing (qForcing) if present:
if (forcing) sqd2=sqd2+mrat1*qForcing
 !mrat1 = m_1/(m_2+m_2) where m_j is the mass of layer j (see constants.f90)

 !Add non-conservative term associated with thermal damping if present:
if (thermal) then
  do i=1,numLongGridPts
    sqd2(:,i)=sqd2(:,i)+thermDampRate*(qAnom2(:,i)+cof)*(one-one/(one+h2(:,i)))
  enddo
endif

 !Convert to semi-spectral space and de-alias:
call dealias(sqd2)

!---------------------------------------------------------------
 !qs source (can be computed in each layer independently):

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain qs source for layer 1 (sqs1):
call deriv(numLatGridPts,numLongGridPts,rk,qSpec1,wkg)
call revfft(numLatGridPts,numLongGridPts,wkg,trig,factors) ! wkg = d(qSpec1)/d(lon)
wka=qSpec1
call revfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
call latder(wka,wkh)                ! wkh = d(qSpec1)/d(lat)
do i=1,numLongGridPts
  sqs1(:,i)=-om1(:,i)*wkg(:,i)-v1(:,i)*(wkh(:,i)+beta) !include beta*v1
enddo
 !Convert to semi-spectral space and de-alias:
call dealias(sqs1)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain qs source for layer 2 (sqs2):
call deriv(numLatGridPts,numLongGridPts,rk,qSpec2,wkg)
call revfft(numLatGridPts,numLongGridPts,wkg,trig,factors) ! wkg = d(qSpec2)/d(lon)
wka=qSpec2
call revfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
call latder(wka,wkh)                ! wkh = d(qSpec2)/d(lat)
do i=1,numLongGridPts
  sqs2(:,i)=-om2(:,i)*wkg(:,i)-v2(:,i)*(wkh(:,i)+beta) !include beta*v2
enddo
 !Convert to semi-spectral space and de-alias:
call dealias(sqs2)

!---------------------------------------------------------------
 !delta source (can be computed in each layer independently):

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain delta source for layer 1 (sds1):

 !Compute A = z*u/r + (1/r)*dv/d(lon) -> wka
 !    and B = z*v/r - (1/r)*du/d(lon) -> wkb:
 !where z = sin(lat) and r = cos(lat):

 !Convert v1 & om1 (= u_1/r) to semi-spectral space as wka & wkb:
wka=v1
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
wkb=om1
call forfft(numLatGridPts,numLongGridPts,wkb,trig,factors) 

 !Compute longitudinal derivatives of wka & wkb:
call deriv(numLatGridPts,numLongGridPts,rk,wka,wkc)
call deriv(numLatGridPts,numLongGridPts,rk,wkb,wkd)

 !Recover physical fields (inverse FFT):
call revfft(numLatGridPts,numLongGridPts,wkc,trig,factors) 
call revfft(numLatGridPts,numLongGridPts,wkd,trig,factors) 

 !Complete definition of A & B:
do i=1,numLongGridPts
  wka(:,i)=sinLat*om1(:,i)+secantLat*wkc(:,i)
  wkb(:,i)=tanLat*v1(:,i)-wkd(:,i)
enddo
 !tlat = tan(lat) above; *** Do not re-use wka or wkb ***

 !Get physical space velocity divergence -> d1:
d1=deltaSpec1
call revfft(numLatGridPts,numLongGridPts,d1,trig,factors) 
 !*** Do not re-use d1 ***

 !Get physical space derivatives of divergence:
call deriv(numLatGridPts,numLongGridPts,rk,deltaSpec1,wkg)
call revfft(numLatGridPts,numLongGridPts,wkg,trig,factors) ! wkg = d(delta)/d(lon)
call latder(d1,wkc)                 ! wkc = d(delta)/d(lat)

 !Compute the nonlinear part of the velocity divergence tendency:
usq1=u1**2+v1**2 ! *** Do not re-use usq1 ***
sds1=two*(wka*(z1-wka)-wkb*(wkb+d1))-usq1-d1**2-om1*wkg-v1*wkc
 !om1 = u_1/r on rhs above; *** wka & wkb now safe to re-use below ***

 !Add effect of fixed topography (if present) --- layer 1 only:
if (topoFlag) sds1=sds1-specDeltaTopoSrc

 !Convert to semi-spectral space and de-alias:
call dealias(sds1)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain delta source for layer 2 (sds2):

 !Compute A = z*u/r + (1/r)*dv/d(lon) -> wka
 !    and B = z*v/r - (1/r)*du/d(lon) -> wkb:
 !where z = sin(lat) and r = cos(lat):

 !Convert v2 & om2 (= u_2/r) to semi-spectral space as wka & wkb:
wka=v2
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
wkb=om2
call forfft(numLatGridPts,numLongGridPts,wkb,trig,factors) 

 !Compute longitudinal derivatives of wka & wkb:
call deriv(numLatGridPts,numLongGridPts,rk,wka,wkc)
call deriv(numLatGridPts,numLongGridPts,rk,wkb,wkd)

 !Recover physical fields (inverse FFT):
call revfft(numLatGridPts,numLongGridPts,wkc,trig,factors) 
call revfft(numLatGridPts,numLongGridPts,wkd,trig,factors) 

 !Complete definition of A & B:
do i=1,numLongGridPts
  wka(:,i)=sinLat*om2(:,i)+secantLat*wkc(:,i)
  wkb(:,i)=tanLat*v2(:,i)-wkd(:,i)
enddo
 !tlat = tan(lat) above; *** Do not re-use wka or wkb ***

 !Get physical space velocity divergence -> d2:
d2=deltaSpec2
call revfft(numLatGridPts,numLongGridPts,d2,trig,factors) 
 !*** Do not re-use d2 ***

 !Get physical space derivatives of divergence:
call deriv(numLatGridPts,numLongGridPts,rk,deltaSpec2,wkg)
call revfft(numLatGridPts,numLongGridPts,wkg,trig,factors) ! wkg = d(delta)/d(lon)
call latder(d2,wkc)                 ! wkc = d(delta)/d(lat)

 !Compute the nonlinear part of the velocity divergence tendency:
usq2=u2**2+v2**2 ! *** Do not re-use usq2 ***
sds2=two*(wka*(z2-wka)-wkb*(wkb+d2))-usq2-d2**2-om2*wkg-v2*wkc
 !om2 = u_2/r on rhs above; *** wka & wkb now safe to re-use below ***

 !Convert to semi-spectral space and de-alias:
call dealias(sds2)

!---------------------------------------------------------------
 !gamma source (these are coupled through pressure effects)

 !First obtain c^2*Lap(div(h*(u,v))) in each layer as wkg & wkh;
 !both of these terms appear in each layer's gamma source term.

!::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Compute c_1^2*div((u_1,v_1)*h_1) -> wke:
wka=h1
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
call deriv(numLatGridPts,numLongGridPts,rk,wka,wkd)
call revfft(numLatGridPts,numLongGridPts,wkd,trig,factors) ! wkd = dh_1/d(lon)
call latder(h1,wkc)                 ! wkc = dh_1/d(lat)
wke=d1*h1+om1*wkd+v1*wkc
 !Note om1 = u_1/r on rhs above

 !If thermal damping is present, add its contribution:
if (thermal) wke=wke+thermDampRate*h1

 !Compute Laplacian of wke = div(h_1*(u_1,v_1)) after de-aliasing:
call dealias(wke)
call laplace(wke,wkg)
wkg=csq1*wkg
 !wkg = c_1^2*Lap(div(h_1*(u_1,v_1))) in semi-spectral space

!::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Compute c_2^2*div((u_2,v_2)*h_2) -> wke:
wka=h2
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
call deriv(numLatGridPts,numLongGridPts,rk,wka,wkd)
call revfft(numLatGridPts,numLongGridPts,wkd,trig,factors) ! wkd = dh_2/d(lon)
call latder(h2,wkc)                 ! wkc = dh_2/d(lat)
wke=d2*h2+om2*wkd+v2*wkc
 !Note om2 = u_2/r on rhs above

 !If thermal damping is present, add its contribution:
if (thermal) wke=wke+thermDampRate*h2

 !Compute Laplacian of wke after de-aliasing:
call dealias(wke)
call laplace(wke,wkh)
wkh=csq2*wkh
 !wkh = c_2^2*Lap(div(h_2*(u_2,v_2))) in semi-spectral space

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain gamma source for layer 1 (sgs1):

 !Define wkf = f*zeta_1 & wkd = 2*f*beta*v_1:
do i=1,numLongGridPts
  wkf(:,i)=cof*z1(:,i)            !cof = f
  wkd(:,i)=dfb*v1(:,i)            !dfb = 2*f*beta
enddo

 !Compute div((u_1,v_1)*Z_1) where Z_1 = f*zeta_1:
wka=wkf
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors)
call deriv(numLatGridPts,numLongGridPts,rk,wka,wke)
call revfft(numLatGridPts,numLongGridPts,wke,trig,factors) ! wke = dZ/d(lon)
call latder(wkf,wkc)                ! wkc = dZ/d(lat)

 !Store div((u_1,v_1)*Z_1) + 2*f*beta*v_1 in sgs1 temporarily:
sgs1=d1*wkf+om1*wke+v1*wkc+wkd
 !De-alias and convert to semi-spectral space:
call dealias(sgs1)

 !Define B = c_1^2*h_1 + (rho_2/rho_1)*c_2^2*h_2 + (u_1^2 + v_1^2)/2 -> wkb:
wkb=csq1*h1+densityRatio*csq2*h2+f12*usq1
 !De-alias and convert to semi-spectral space:
call dealias(wkb)
 !Add effect of fixed topography (if present) --- layer 1 only:
if (topoFlag) wkb=wkb+specGammaTopoSrc
 !Calculate dB/d(lon) -> wkf (keep this semi-spectral):
call deriv(numLatGridPts,numLongGridPts,rk,wkb,wkf)

 !Complete calculation of (de-aliased) nonlinear part of gamma tendency:
sgs1=wkg+densityRatio*wkh+fpole*wkf-sgs1
 !wkg & wkh are defined above as the Laplacian of c^2*div(h*(u,v))
 !in layers 1 and 2, respectively.

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain gamma source for layer 2 (sgs2):

 !Define wkf = f*zeta_2 & wkd = 2*f*beta*v_2:
do i=1,numLongGridPts
  wkf(:,i)=cof*z2(:,i)            !cof = f
  wkd(:,i)=dfb*v2(:,i)            !dfb = 2*f*beta
enddo

 !Compute div((u_2,v_2)*Z_2) where Z_2 = f*zeta_2:
wka=wkf
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors)
call deriv(numLatGridPts,numLongGridPts,rk,wka,wke)
call revfft(numLatGridPts,numLongGridPts,wke,trig,factors) ! wke = dZ/d(lon)
call latder(wkf,wkc)                ! wkc = dZ/d(lat)

 !Store div((u_2,v_2)*Z_2) + 2*f*beta*v_2 in sgs2 temporarily:
sgs2=d2*wkf+om2*wke+v2*wkc+wkd
 !De-alias and convert to semi-spectral space:
call dealias(sgs2)

 !Define B = c_1^2*h_1 + c_2^2*h_2 + (u_2^2 + v_2^2)/2 -> wkb:
wkb=csq1*h1+csq2*h2+f12*usq2
 !De-alias and convert to semi-spectral space:
call dealias(wkb)
 !Calculate dB/d(lon) -> wkf (keep this semi-spectral):
call deriv(numLatGridPts,numLongGridPts,rk,wkb,wkf)

 !Complete calculation of (de-aliased) nonlinear part of gamma tendency:
sgs2=wkg+wkh+fpole*wkf-sgs2
 !wkg & wkh are defined above as the Laplacian of c^2*div(h*(u,v))
 !in layers 1 and 2, respectively.

!-----------------------------------------------------
 !Remove global mean values of sds & sgs:
rhs=sds1(:,1)*cosLat
avgval=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sds1(:,1)=sds1(:,1)-avgval

rhs=sds2(:,1)*cosLat
avgval=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sds2(:,1)=sds2(:,1)-avgval

rhs=sgs1(:,1)*cosLat
avgval=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sgs1(:,1)=sgs1(:,1)-avgval

rhs=sgs2(:,1)*cosLat
avgval=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sgs2(:,1)=sgs2(:,1)-avgval

return
end subroutine source

!=======================================================================

subroutine inversion

! Finds the gridded dimensionless height (thickness) anomaly (h) and
! velocity field (u,v) from the PV contours and the (semi-spectral)
! divergence (ds) and acceleration divergence (gs).
  
implicit none

! Local variable:
double precision:: qa(numLatGridPts,numLongGridPts,numLayers)

!------------------------------------------------------------
 !Convert PV contours (x,y,z) to gridded PV anomaly as qa:
call con2grid(qa)

 !Convert each component to semi-spectral space as qc1 & qc2:
qc1=qa(:,:,1)
call forfft(numLatGridPts,numLongGridPts,qc1,trig,factors) 
qc2=qa(:,:,2)
call forfft(numLatGridPts,numLongGridPts,qc2,trig,factors) 

 !Combine fields to update qTot with full (semi-spectral) field,
 !qTot = F*(qs-qc) + qc + qd, where F is a low pass filter:
qt1=qSpec1-qc1
qt2=qSpec2-qc2

call lopass(qt1)
call lopass(qt2)

qt1=qt1+qc1+qd1
qt2=qt2+qc2+qd2

 !Invert PV, divergence and acceleration divergence to obtain the
 !dimensionless height anomaly and velocity field, as well as the
 !gridded PV anomaly and relative vorticity (see spectral.f90):
call main_invert(qt1,qt2,deltaSpec1,deltaSpec2,gammaSpec1,gammaSpec2,h1,h2,u1,u2,v1,v2,qAnom1,qAnom2,z1,z2)
 !Note: qTot, ds & gs (for each layer) are in semi-spectral space while 
 !      h, u, v, q and z are in physical space.

return
end subroutine inversion

!=======================================================================

subroutine reset
! Resets qs & qd at the beginning of each time step, or
! just before recontouring at the end of a run.

implicit none

 !Local variable:
double precision:: qa(numLatGridPts,numLongGridPts,numLayers)

!-------------------------------------------------------------------
 !Reset  qs = F*(qs-qc) + qc + qd  and  qd = (1-F)*(qs-qc)
 !where F & (1-F) are low & high pass filter functions
 !(see spectral.f90)

 !Convert PV contours (x,y,z) to gridded PV anomaly as qa:
call con2grid(qa)

 !Convert each component to semi-spectral space as qc1 & qc2:
qc1=qa(:,:,1)
call forfft(numLatGridPts,numLongGridPts,qc1,trig,factors) 
qc2=qa(:,:,2)
call forfft(numLatGridPts,numLongGridPts,qc2,trig,factors) 

 !Reset qs & qd in each layer:
qSpec1=qSpec1-qc1
qSpec2=qSpec2-qc2

call lopass(qSpec1)
call lopass(qSpec2)

qSpec1=qSpec1+qc1+qd1
qSpec2=qSpec2+qc2+qd2

qd1=qSpec1-qc1
qd2=qSpec2-qc2

call hipass(qd1)
call hipass(qd2)

return
end subroutine reset

!=======================================================================

subroutine pvforce
! Updates the PV forcing by blending a new random field into the
! existing one in qForcing.

implicit none

! Local variables:
double precision:: wkb(numLatGridPts,numLongGridPts)
double precision:: fac

!-------------------------------------------------------------------
 !Generate forcing at t+dt:
call generate_forcing(wkb,rmsPvForcing)

 !Blend with existing forcing (Markovian process):
qForcing=wold*qForcing+wnew*wkb
 !See constants.f90 for wold and wnew

 !Restore rms value of blended field:
call getrms(qForcing,fac)
fac=rmsPvForcing/fac
qForcing=fac*qForcing

return
end subroutine pvforce

!=======================================================================

subroutine diagnose

! Computes the twist parameter, the time integral of |zeta|_max, and
! various quantities at every time step to monitor the flow evolution.

implicit none

 !Local variables:
double precision:: wka(numLatGridPts,numLongGridPts),wkb(numLatGridPts,numLongGridPts)
double precision:: zmax,frmax,zrms

!----------------------------------------------------------------------
! Compute diagnostics:
zmax=max(maxval(abs(z1)),maxval(abs(z2)))

wka=(u1**2+v1**2)/(one+h1)
wkb=(u2**2+v2**2)/(one+h2)
frmax=sqrt(max(csq1i*maxval(wka),csq2i*maxval(wkb)))

wka=mrat1*(z1**2+mubar*z2**2)
zrms=sqrt(average(wka))
! above, mubar = (rho_2*H_2)/(rho_1*H_1) is the layer mass ratio
! and mrat1 = 1/(1+mubar)

!Increment the integral of |zeta|_max:
twist=twist+dt*zmax

 !Record |zeta|_max/f_pole, Fr_max and zeta_rms/fpole:
write(16,'(1x,f12.5,3(1x,f14.8))') t,zmax/fpole,frmax,zrms/fpole

return
end subroutine diagnose

!=======================================================================

subroutine savegrid(igrids)

! Saves fields, energy components, angular momenta and spectra
! of various fields at the desired save time

implicit none

 !Passed variable:
integer:: igrids

 !Local variables:
double precision:: wkp(numLatGridPts,numLongGridPts),wka(numLatGridPts,numLongGridPts),wkb(numLatGridPts,numLongGridPts)
double precision:: spec1(numLatGridPts),spec2(numLatGridPts)
double precision:: ekin,epot,etot,angm
integer:: i,m

!-----------------------------------------------------------------
! Compute energy components, divided by rho_1*H where H = H_1+H_2:
wka=hAvgLayer1*(one+h1) !Physical thickness of layer 1
wkb=hAvgLayer2*(one+h2) !Physical thickness of layer 2
! *** Do not re-use wka & wkb ***

! Horizontal and vertical parts of the kinetic energy:
wkp=wka*(u1**2+v1**2)+densityRatio*wkb*(u2**2+v2**2) !kinetic energy density
ekin=twoPi*average(wkp)
! This is the area integral of the kinetic energy density
! densityRatio = rho_2/rho_1 (see constants.f90)

! Potential energy:
if (topoFlag) then
   wkp=hAvgLayer1*botTopo+wka !Physical height of the top of layer 1
   wkp=wka*(two*wkp-wka)+densityRatio*wkb*(two*wkp+wkb)-dape
else
   wkp=wka**2+densityRatio*wkb*(two*wka+wkb)-dape
endif
epot=twoPi*gravity*average(wkp)
! wkp: potential energy density
! dape = hbar_1**2+(rho_2/rho_1)*hbar_2*(2*hbar_1+hbar_2): ensures that
!        a state of rest has zero potential energy (see constants.f90)
! *** wka & wkb are now free to re-use ***

 !Total energy (sum of the above):
etot=ekin+epot

 !Compute angular momentum:
do i=1,numLongGridPts
  wka(:,i)=cosLat*((one+h1(:,i))*u1(:,i)+omega*cosLat*h1(:,i))
  wkb(:,i)=cosLat*((one+h2(:,i))*u2(:,i)+omega*cosLat*h2(:,i))
enddo
 !wka & wkb: angular momentum densities in layers 1 and 2
 !Omega = planetary rotation rate (= 2*pi, see constants.f90)
wkp=wka+densityRatio*wkb
angm=fourpi*average(wkp)

 !Write energies & angular momentum to ecomp.asc:
write(15,'(f13.6,4(1x,f16.9))') t,ekin,epot,etot,angm
write(*,'(a,f13.6,a,f14.7,a,f14.7)') ' t = ',t,'     E = ',etot, &
                                     '     A = ',angm

!-----------------------------------------------------------------
 !Compute 1d longitudinal power spectra for various fields and
 !for each vertical mode. Layer quantities F_j are projected
 !onto the barotropic mode using F_1*mrat1 + F_2*mrat2 and
 !onto the baroclinic mode using F_2 - F_1, where
 !mrat1 = m_1/(m_1+m_2) and mrat2 = m_2/(m_1+m_2),
 !with m_j equal to the mean mass of layer j.

! Relative vorticity:
wka=z1
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors)
wkb=z2
call forfft(numLatGridPts,numLongGridPts,wkb,trig,factors)
wkp=wka*mrat1+wkb*mrat2    !Barotropic
call longspec(wkp,spec1)
spec1=log10(spec1+1.d-32)
wkp=wkb-wka                !Baroclinic
call longspec(wkp,spec2)
spec2=log10(spec2+1.d-32)
write(51,'(f13.6,1x,i5)') t,numLatGridPts
do m=1,numLatGridPts
  write(51,'(3(1x,f12.8))') alm(m),spec1(m),spec2(m)
enddo

! Divergence:
wkp=deltaSpec1*mrat1+deltaSpec2*mrat2    !Barotropic
call longspec(wkp,spec1)
spec1=log10(spec1+1.d-32)
wkp=deltaSpec2-deltaSpec1                !Baroclinic
call longspec(wkp,spec2)
spec2=log10(spec2+1.d-32)
write(52,'(f13.6,1x,i5)') t,numLatGridPts
do m=1,numLatGridPts
  write(52,'(2(1x,f12.8))') alm(m),spec1(m),spec2(m)
enddo

! Acceleration divergence:
wkp=gammaSpec1*mrat1+gammaSpec2*mrat2    !Barotropic
call longspec(wkp,spec1)
spec1=log10(spec1+1.d-32)
wkp=gammaSpec2-gammaSpec1                !Baroclinic
call longspec(wkp,spec2)
spec2=log10(spec2+1.d-32)
write(53,'(f13.6,1x,i5)') t,numLatGridPts
do m=1,numLatGridPts
  write(53,'(2(1x,f12.8))') alm(m),spec1(m),spec2(m)
enddo

! Dimensionless height (thickness) anomaly:
wka=h1
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors)
wkb=h2
call forfft(numLatGridPts,numLongGridPts,wkb,trig,factors)
wkp=wka*mrat1+wkb*mrat2    !Barotropic
call longspec(wkp,spec1)
spec1=log10(spec1+1.d-32)
wkp=wkb-wka                !Baroclinic
call longspec(wkp,spec2)
spec2=log10(spec2+1.d-32)
write(54,'(f13.6,1x,i5)') t,numLatGridPts
do m=1,numLatGridPts
  write(54,'(3(1x,f12.8))') alm(m),spec1(m),spec2(m)
enddo

!-----------------------------------------------------------------
! Write various gridded fields to direct access files:

! PV anomaly:
write(31,rec=igrids) real(t),real(qAnom1)
write(41,rec=igrids) real(t),real(qAnom2)

! Velocity divergence:
wkp=deltaSpec1
call revfft(numLatGridPts,numLongGridPts,wkp,trig,factors)
write(32,rec=igrids) real(t),real(wkp)
wkp=deltaSpec2
call revfft(numLatGridPts,numLongGridPts,wkp,trig,factors)
write(42,rec=igrids) real(t),real(wkp)

! Acceleration divergence:
wkp=gammaSpec1
call revfft(numLatGridPts,numLongGridPts,wkp,trig,factors)
write(33,rec=igrids) real(t),real(wkp)
wkp=gammaSpec2
call revfft(numLatGridPts,numLongGridPts,wkp,trig,factors)
write(43,rec=igrids) real(t),real(wkp)

! Dimensionless height (thickness) anomaly:
write(34,rec=igrids) real(t),real(h1)
write(44,rec=igrids) real(t),real(h2)

! Relative vorticity:
write(35,rec=igrids) real(t),real(z1)
write(45,rec=igrids) real(t),real(z2)

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont(irec)

! Saves PV contours for post-processing and imaging

implicit none

 !Passed variable:
integer:: irec

 !Local variables:
double precision:: qa(numLatGridPts,numLongGridPts)
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') irec

 !Write contours to the contours subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') numContoursTot,numNodesTot,t,dq

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
qa=qt1-qc1
call revfft(numLatGridPts,numLongGridPts,qa,trig,factors)
write(83,rec=2*irec-1) real(t),real(qa)
qa=qt2-qc2
call revfft(numLatGridPts,numLongGridPts,qa,trig,factors)
write(83,rec=2*irec) real(t),real(qa)
! Note: both layers are saved above

 !Save PV contours:
open(81,file='contours/qqindex'//pind,form='unformatted',status='replace')
write(81) numNodesOnContourList(1:numContoursTot),firstNodeIdxList(1:numContoursTot),contsQIdx(1:numContoursTot),contourIdxToLayerIdx(1:numContoursTot)
close(81)

open(82,file='contours/qqnodes'//pind,form='unformatted',status='replace')
write(82) x(1:numNodesTot),y(1:numNodesTot),z(1:numNodesTot)
close(82)

return
end subroutine savecont

!=======================================================================

 !Main end module
end module evolution
