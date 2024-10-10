module evolution

! Module containing subroutines to evolve PV contours and all fields.

use common

implicit none

 !PV field after conversion from contours and residual as in CLAM:
double precision:: qc(0:ny,0:nxm1),qd(0:ny,0:nxm1)
 !Gradient of the geopotential anomaly:
double precision:: phix(0:ny,0:nxm1),phiy(0:ny,0:nxm1)

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
integer:: ireg
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
  call savedata(ggen)
   
   !Perform contour surgery or recontouring when twist is large enough:
  if (twist .gt. twistmax) then
    ireg=ireg+1

     !Don't continue if maximum number of regularisations reached:
    if (ireg .eq. nregmax) then
       !Prepare PV residual qr for recontouring (and preserve qs):
      call prepare
       !Exit module and go to recontouring:
      return
    endif

     !Regularise the PV contours (surgery + node redistribution):
    call surgery
     !Record active contour complexity to complexity.asc:
    write(14,'(1x,f12.5,1x,i9,1x,i10)') t,nq,nptq

    twist=twist-twistmax
  endif

   !Advect flow from time t to t + dt:
  call advance(ggen)

enddo
 !End of time loop
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

 !Possibly save final data:
call savedata(ggen)

return
end subroutine advect

!=======================================================================

subroutine init

! Initialises residual PV for normal time integration following 
! contour regeneration

implicit none

!------------------------------------------------------------------
 !Record active contour complexity to complexity.asc:
write(14,'(1x,f12.5,1x,i8,1x,i9)') t,nq,nptq

!-------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)
 !Convert qc to semi-spectral space:
call forfft(nyp1,nx,qc,xtrig,xfactors)

 !Define (semi-spectral) residual PV qd = (1-F)[qs-qc]:
qd=qs-qc
qc=qd !Re-use qc as a work array:
call lowpass(qc)
 !qc now contains F[qs-qc] in semi-spectral space.
qd=qd-qc

return
end subroutine init

!=======================================================================

subroutine prepare

! This routine is called just before exiting to contour regeneration.
! The current (semi-spectral) PV field is stored in qs, and the
! residual PV needed in congen.f90 is stored in qr.

implicit none

!-----------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)

 !Convert qc to semi-spectral space:
call forfft(nyp1,nx,qc,xtrig,xfactors)

 !Define spectral PV (qs) and PV residual (qr):
qr=qs-qc
call lowpass(qr)
qr=qr+qd           !qr = F(qs-qc) + qd
qs=qr+qc           !qs = F(qs-qc) + qd + qc

 !Convert the PV residual to physical space for contouring purposes:
call revfft(nyp1,nx,qr,xtrig,xfactors)
 !Note: we are leaving this module next and qd will be redefined
 !upon re-entry in subroutine init.

return
end subroutine prepare

!=======================================================================

subroutine savedata(ggen)

! Saves data periodically

implicit none

 !Passed variable:
logical:: ggen

 !Local variables:
integer:: itime,jtime

!-----------------------------------------------------
itime=nint(t/dt)
jtime=itime/ngsave

if (ngsave*jtime .eq. itime) then
  call inversion
  call savegrid(jtime+1)
  ggen=.false.
else
  ggen=.true.
endif
 !ggen is used to work out if calling inversion is needed in advance

jtime=itime/ncsave
if (ncsave*jtime .eq. itime) call savecont(jtime+1)

return
end subroutine savedata

!=======================================================================

subroutine advance(ggen)

! Advances PV from time t to t+dt by a combination of contour
! advection (for PV contours) and the pseudo-semi-spectral method 
! (for all semi-spectral fields, namely qs, qd, aa, bb, uum & uup).

! Uses an iterative implicit (trapezoidal) method of the form
!
! (F^{n+1}-F^n)/dt = L[(F^{n+1}+F^n)/2] + N[(F^{n+1}+F^n)/2]
!
! for a field F, where n refers to the time level, L refers to
! the linear source terms (if present), and N refers to the
! nonlinear source terms.  We start with a guess for F^{n+1}
! in N and iterate niter times (see parameter statement below).

implicit none

 !Passed variable:
logical:: ggen

 !Local variables:
integer,parameter:: niter=2

 !Semi-spectral fields needed in PV time stepping:
double precision:: qsi(0:ny,0:nxm1),sqs(0:ny,0:nxm1)
double precision:: qdi(0:ny,0:nxm1),qdm(0:ny,0:nxm1),sqd(0:ny,0:nxm1)

 !Semi-spectral fields needed in acceleration time stepping:
double precision:: aai(0:ny,0:nxm1),saa(0:ny,0:nxm1),faa(0:ny,0:nxm1)
double precision:: bbi(0:ny,0:nxm1),sbb(0:ny,0:nxm1),fbb(0:ny,0:nxm1)
double precision:: taa(0:ny,0:nxm1),gaa(0:ny,0:nxm1)
double precision:: tbb(0:ny,0:nxm1),gbb(0:ny,0:nxm1)

 !Spectral fields needed in boundary zonal velocity time stepping:
double precision:: uumi(0:nxm1),suum(0:nxm1),fuum(0:nxm1),duum(0:nxm1)
double precision:: uupi(0:nxm1),suup(0:nxm1),fuup(0:nxm1),duup(0:nxm1)
double precision:: aami(0:nxm1),aam(0:nxm1)
double precision:: aapi(0:nxm1),aap(0:nxm1)

 !Contour positions and velocities needed in time stepping:
double precision:: xqi(nptq),yqi(nptq),uq(nptq),vq(nptq)

 !Other local quantities:
double precision:: xx,yy
integer:: i,j,iter

!-------------------------------------------------------------------
 !Compute height and velocity at current time level, say t=t^n:
if (ggen) call inversion
 !If ggen is false, inversion was called previously at this time level.
 !Note, inversion returns qc in semi-spectral space, needed just below.

 !Re-initialise qs & qd at the beginning of the time step:
 !          Reset qs = F*(qs-qc) + qc + qd
 !            and qd = (1-F)*(qs-qc)
 !Here F is a low pass filter (see spectral.f90)
qsi=qs-qc
call lowpass(qsi)
qdi=qsi+qd          !qdi = F*(qs-qc) + qd
qs=qdi+qc           !qs now reset; qdi above is now the new qs-qc
qd=qdi
call lowpass(qdi)   !This is F(qs-qc)
qd=qd-qdi           !qd now reset (high pass filter of qs-qc)

 !Compute twist parameter and save various diagnostics each time step:
call diagnose

!------------------------------------------------------------------
 !Start with a guess for F^{n+1} for all contours and fields:

 !Contours:
call velint(uu,vv,uq,vq)
 !Here (uq,vq) is the contour velocity at time level n, i.e. u(x^n,t^n)
do i=1,nptq
   !Store x^n+0.5*dt*u(x^n,t^n) for efficiency in the iteration below:
  xx=xq(i)+dt2*uq(i)
  yy=yq(i)+dt2*vq(i)
  xqi(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  yqi(i)=min(ymax,max(ymin,yy))
   !Preliminary guess for x^{n+1}:
  xx=xq(i)+dt*uq(i)
  yy=yq(i)+dt*vq(i)
  xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
  yq(i)=min(ymax,max(ymin,yy))
enddo

 !Calculate the source terms for the PV fields, the acceleration and
 !and the boundary zonal velocity:
call source(sqs,sqd,saa,sbb,taa,tbb,suum,suup)

 !Update PV fields:
qsi=qs+dt2*sqs
qs=qs+dt*sqs
qdi=qd
qdm=qd+dt4*sqd
do j=0,ny
  qd(j,:)=diss*(qdm(j,:)+dt4*sqd(j,:))-qdi(j,:)
enddo
 !diss is a hyperdiffusion operator (see spectral.f90)

 !Update zonal velocity at y boundaries (below dt4i = 4/dt):
uumi=uum
uupi=uup
fuum=suum+dt4i*uumi
fuup=suup+dt4i*uupi
suum=fuum+suum           !2*S^{-}_check
call xderiv1d(suum,duum) !x derivative of the above
suup=fuup+suup           !2*S^{+}_check
call xderiv1d(suup,duup) !x derivative of the above
uum=divu*(rdis*suum+betu*duup-alpu*duum)-uumi
uup=divu*(rdis*suup+alpu*duup-betu*duum)-uupi

 !Compute the acceleration a = -phi_x at the boundaries (spectral space):
aami=-phix(0,:)
call forfft(1,nx,aami,xtrig,xfactors)
aapi=-phix(ny,:)
call forfft(1,nx,aapi,xtrig,xfactors)

 !Update interior acceleration field (aa,bb):
aai=aa
bbi=bb
faa=saa+dt4i*aai
fbb=sbb+dt4i*bbi
saa=faa+saa          !2*\check{N}_a
sbb=fbb+sbb          !2*\check{N}_b
gaa=taa+dt4i*aai
gbb=tbb+dt4i*bbi
taa=gaa+taa          !2*\check{T}_a
tbb=gbb+tbb          !2*\check{T}_b
 !Put (twice) boundary acceleration in edge values of sbb array:
sbb(0 ,:)=two*aami
sbb(ny,:)=two*aapi
 !Solve semi-implicit system:
call sistep(saa,sbb,taa,tbb)
 !On return, saa = aa^n + aa^{n+1} & sbb = bb^n + bb^{n+1}
aa=saa-aai
bb=sbb-bbi

!------------------------------------------------------------------
 !Iterate to improve estimates of F^{n+1}:
do iter=1,niter
   !Perform inversion at t^{n+1} from estimated quantities:
  call inversion

   !Calculate the source terms for the PV fields, the acceleration and
   !and the boundary zonal velocity:
  call source(sqs,sqd,saa,sbb,taa,tbb,suum,suup)

   !Interpolate gridded velocity (uu,vv) at contour nodes as (uq,vq):
  call velint(uu,vv,uq,vq)
  do i=1,nptq
    xx=xqi(i)+dt2*uq(i)
    yy=yqi(i)+dt2*vq(i)
    xq(i)=oms*(xx-ellx*dble(int(xx*hlxi)))
    yq(i)=min(ymax,max(ymin,yy))
  enddo
   !Now (xq,yq) contain a new guess for x^{n+1}.

   !Update PV fields:
  qs=qsi+dt2*sqs
  do j=0,ny
    qd(j,:)=diss*(qdm(j,:)+dt4*sqd(j,:))-qdi(j,:)
  enddo

   !Update zonal velocity at y boundaries:
  suum=fuum+suum           !2*S^{-}_check
  call xderiv1d(suum,duum) !x derivative of the above
  suup=fuup+suup           !2*S^{+}_check
  call xderiv1d(suup,duup) !x derivative of the above
  uum=divu*(rdis*suum+betu*duup-alpu*duum)-uumi
  uup=divu*(rdis*suup+alpu*duup-betu*duum)-uupi

   !Compute the acceleration a = -phi_x at the boundaries (spectral space):
  aam=-phix(0,:)
  call forfft(1,nx,aam,xtrig,xfactors)
  aap=-phix(ny,:)
  call forfft(1,nx,aap,xtrig,xfactors)

   !Update interior acceleration field (aa,bb):
  saa=faa+saa          !2*\check{N}_a
  sbb=fbb+sbb          !2*\check{N}_b
  taa=gaa+taa          !2*\check{T}_a
  tbb=gbb+tbb          !2*\check{T}_b
   !Put (summed) boundary acceleration in edge values of sbb array:
  sbb(0 ,:)=aami+aam
  sbb(ny,:)=aapi+aap
   !Solve semi-implicit system:
  call sistep(saa,sbb,taa,tbb)
   !On return, saa = aa^n + aa^{n+1} & sbb = bb^n + bb^{n+1}
  aa=saa-aai
  bb=sbb-bbi
enddo

 !Advance time:
t=t+dt

return
end subroutine advance

!====================================================================
subroutine sistep(saa,sbb,taa,tbb)

! Solves the semi-implicit system for the acceleration (a,b). Returns
! the sums of these quantities at time levels n and n+1, from which
! the values at level n+1 are easily computed in the calling routine.

! Here saa = 2*\check{N}_a, sbb = 2*\check{N}_b, taa = 2*\check{T}_a
! and  tbb = 2*\check{T}_b, (twice) the source terms needed to find
! 2*\bar{a} = a^n + a^{n+1} and 2*\bar{b} = b^n + b^{n+1}.

! Important: the edge values sbb, which are strictly zero, here must
! instead contain the sums of the edge values of a (= -phi_x) at time
! levels n and n+1. The edge values of saa are not used.

! All quantities are in semi-spectral space.

 !Declarations:
implicit none

 !Passed variables:
double precision:: saa(0:ny,0:nxm1),sbb(0:ny,0:nxm1)
double precision:: taa(0:ny,0:nxm1),tbb(0:ny,0:nxm1)

 !Local variables:
double precision:: wka(0:ny,0:nxm1),wkb(0:ny,0:nxm1)
integer:: j,m

!---------------------------------------------------------------------
 !Form r.h.s. of the equation for a^n + a^{n+1}:
call xderiv(taa,wka)
call yderiv(wka,wkb)
do j=1,nym1
  wkb(j,:)=-cofi*(ksq*tbb(j,:)+wkb(j,:)+kdsq*sbb(j,:))-csqi*rdis*saa(j,:)
enddo
 !wkb = 2*\check{S}_a in the notes.

 !Solve the tri-diagonal system for a^n + a^{n+1}:
 !(Below ap = 1/dy^2)
do m=0,nxm1
  aa(0,m)=sbb(0,m)   !sbb(0,m)  = a^n + a^{n+1} at the lower boundary (known)
  do j=1,nym1
    aa(j,m)=(wkb(j,m)-ap*aa(j-1,m))*htdt(j,m)
  enddo
  aa(ny,m)=sbb(ny,m) !sbb(ny,m) = a^n + a^{n+1} at the upper boundary (known)
  do j=nym1,1,-1
    aa(j,m)=etdt(j,m)*aa(j+1,m)+aa(j,m)
  enddo
enddo

 !Obtain b^n + b^{n+1}:
call xderiv(aa,wka)
call yderiv(wka,bb) !bb = zero on boundaries after this call
do j=1,nym1
  bb(j,:)=kksqi*(geoi*(rdis*aa(j,:)-saa(j,:))-bb(j,:))
enddo
 !Above, geoi = f/c^2 and kksqi = 1/K^2.

return
end subroutine sistep

!=======================================================================

subroutine inversion

! Finds the gridded dimensionless height anomaly (hh), velocity
! field (uu,vv) and relative vorticity (zz) from the PV contours (qc)
! and auxillary PV fields (qs,qd), the acceleration (aa,bb) and the
! zonal velocities (uum & uup) at y = y_min & y_max.

! hh, uu, vv and zz are all returned in physical space along with
! the gradient of phi = c^2*h_tilde in (phix,phiy).

implicit none

 !Local work array:
double precision:: wk(0:ny,0:nxm1)

!---------------------------------------------------------------
 !Convert PV contours (xq,yq) to gridded values as qc:
call con2grid(qc)

 !Convert qc to semi-spectral space:
call forfft(nyp1,nx,qc,xtrig,xfactors)

 !Combine fields to update qq with full (semi-spectral) field,
 !qq = F[qs-qc]+qc+qd, where F is a low pass filter:
qq=qs-qc
call lowpass(qq)
qq=qq+qc+qd

 !Convert qq to physical space:
call revfft(nyp1,nx,qq,xtrig,xfactors)

!---------------------------------------------------------------
 !Obtain the dimensionless height anomaly (hh) and velocity
 !field (uu,vv) as well as the relative vorticity (zz) by
 !inversion (see spectral.f90 for full details):
call invert(qq,aa,bb,uum,uup,qoff,hh,uu,vv,zz)

 !Calculate the gradient of phi = c^2*h_tilde -> (phix,phiy):
wk=csq*hh
call forfft(nyp1,nx,wk,xtrig,xfactors)
call gradient(wk,phix,phiy)
 !Add boundary values of phiy (set to zero by gradient):
phiy(0 ,:)=-cof*uu(0 ,:)
phiy(ny,:)=-cof*uu(ny,:)

return
end subroutine inversion

!=======================================================================

subroutine source(sqs,sqd,saa,sbb,taa,tbb,suum,suup)

! Gets the source terms (sqs,sqd) for the PV fields (qs,qd), as well
! as (saa,sbb,taa,tbb) for the acceleration divergence (aa,bb) and
! finally (suum,suup) for the boundary zonal velocities (uum,uup) ---
! all in semi-spectral space.  Note that saa ... suup only include
! the nonlinear terms for a semi-implicit treatment, closely analogous
! to that described in the appendix of Mohebalhojeh & Dritschel (2004).

! In the channel SW notes (chsw.pdf), saa & sbb correspond to N_a & N_b
! while taa & tbb correspond to T_a & T_b. Note, N_a = T_a - P_x and
! N_b = T_b - P_y (in semi-spectral space), where P = phi*(b_x-a_y)/f.
! We need N_a, N_b, T_a & T_b and P separately for the semi-implicit
! time stepping used in subroutine sistep.
  
! The spectral fields qs, qd, aa, bb, uum & uup are all de-aliased.
! Furthermore, hh, phix, phiy, uu, vv & zz, obtained by calling
! subroutine inversion before calling this routine, are all de-aliased
! but are in physical space.
  
implicit none

 !Passed variables:
double precision:: sqs(0:ny,0:nxm1),sqd(0:ny,0:nxm1)
double precision:: saa(0:ny,0:nxm1),sbb(0:ny,0:nxm1)
double precision:: taa(0:ny,0:nxm1),tbb(0:ny,0:nxm1)
double precision:: suum(0:nxm1),suup(0:nxm1)

 !x derivative of velocity field (u_x,v_x):
double precision:: ux(0:ny,0:nxm1),vx(0:ny,0:nxm1)
 !y derivative of velocity field (u_y,v_y):
double precision:: uy(0:ny,0:nxm1),vy(0:ny,0:nxm1)

 !x derivative of acceleration field (a_x,b_x):
double precision:: ax(0:ny,0:nxm1),bx(0:ny,0:nxm1)
 !y derivative of acceleration field (a_y,a_y):
double precision:: ay(0:ny,0:nxm1),by(0:ny,0:nxm1)

 !Geopotential anomaly:
double precision:: phi(0:ny,0:nxm1)

 !x derivative of phi_l on the boundaries:
double precision:: plmx(0:nxm1),plpx(0:nxm1)

 !Index needed for loops over y:
integer:: j

!---------------------------------------------------------------
 !qs source --- only nonlinear advection term is needed:
call gradient(qs,taa,tbb)
sqs=-uu*taa-vv*tbb
 !Convert to spectral space:
call forfft(nyp1,nx,sqs,xtrig,xfactors)
 !Apply de-aliasing filter:
call dealias(sqs,1)

!---------------------------------------------------------------
 !qd source --- only nonlinear advection term is needed:
call gradient(qd,taa,tbb)
sqd=-uu*taa-vv*tbb
 !Convert to spectral space:
call forfft(nyp1,nx,sqd,xtrig,xfactors)
 !Apply de-aliasing filter:
call dealias(sqd,1)

!---------------------------------------------------------------
 !Nonlinear part of aa & bb sources (N_a & N_b in the notes):

 !Compute geopotential anomaly:
phi=csq*hh

 !Compute x derivative of velocity field -> (ux,vx):
taa=uu
call forfft(nyp1,nx,taa,xtrig,xfactors)
call xderiv(taa,ux)
call revfft(nyp1,nx,ux,xtrig,xfactors)
taa=vv
call forfft(nyp1,nx,taa,xtrig,xfactors)
call xderiv(taa,vx)
call revfft(nyp1,nx,vx,xtrig,xfactors)

 !Compute y derivative of velocity field -> (uy,vy):
uy=vx-zz
call yderiv(vv,vy)

 !Compute x derivative of acceleration field -> (ax,bx):
call xderiv(aa,ax)
call revfft(nyp1,nx,ax,xtrig,xfactors)
call xderiv(bb,bx)
call revfft(nyp1,nx,bx,xtrig,xfactors)

 !Compute y derivative of acceleration field -> (ay,by):
call yderiv(aa,ay)
call yderiv(bb,by)
 !Here, we use central differencing; boundary values not needed.
 !Also, aa(0,:) and aa(ny,:) must contain the correct boundary
 !acceleration, namely -phi_x, while bb(0,:) = bb(ny,:) = 0.
ay(0,:)=zero
by(0,:)=zero
ay(ny,:)=zero
by(ny,:)=zero
call revfft(nyp1,nx,ay,xtrig,xfactors)
call revfft(nyp1,nx,by,xtrig,xfactors)

 !Compute P = phi*(b_x - a_y)/f (valid only for interior grid points):
taa=cofi*phi*(bx-ay)
 !Transform to semi-spectral space and de-alias:
call forfft(nyp1,nx,taa,xtrig,xfactors)
call dealias(taa,1)

 !Calculate P_x and store in saa temporarily:
call xderiv(taa,saa)

 !Calculate P_y explicitly using only known boundary values:
do j=1,nym1
  taa(j,:)=ap*(aa(j+1,:)-two*aa(j,:)+aa(j-1,:))
enddo
 !Only do this for interior grid points (note: ap = 1/dy^2)
taa(0,:)=zero
taa(ny,:)=zero
call revfft(nyp1,nx,taa,xtrig,xfactors)
 !taa = a_yy in physical space
call yderiv(bx,tbb)
 !tbb = b_xy in physical space
 !Now compute P_y = (phi_y*(b_x-a_y) + phi*(b_xy-ayy))/f:
sbb=cofi*(phiy*(bx-ay)+phi*(tbb-taa))
 !Transform to semi-spectral space and de-alias:
call forfft(nyp1,nx,sbb,xtrig,xfactors)
call dealias(sbb,1)

 !Compute remaining nonlinear terms, T_a & T_b (-> taa & tbb):
taa=ux*phix+vx*phiy-uu*ax-vv*ay
call forfft(nyp1,nx,taa,xtrig,xfactors)
call dealias(taa,1)
tbb=uy*phix+vy*phiy-uu*bx-vv*by
call forfft(nyp1,nx,tbb,xtrig,xfactors)
call dealias(tbb,1)

 !Complete aa & bb source term calculation:
saa=taa-saa
sbb=tbb-sbb
 !Note: these are only defined at the interior grid points.

!---------------------------------------------------------------
 !Source terms for the boundary zonal velocities:
suum=alpu*uum-betu*uup   !phi^{-}_l in semi-spectral space
call xderiv1d(suum,plmx) !x derivative of the above
suum=phix(0,:)+uu(0,:)*ux(0,:)
call forfft(1,nx,suum,xtrig,xfactors)
suum=plmx-filt*suum

suup=betu*uum-alpu*uup   !phi^{+}_l in semi-spectral space
call xderiv1d(suup,plpx) !x derivative of the above
suup=phix(ny,:)+uu(ny,:)*ux(ny,:)
call forfft(1,nx,suup,xtrig,xfactors)
suup=plpx-filt*suup

return
end subroutine source

!=======================================================================

subroutine diagnose

! Computes the twist parameter, the time integral of |zeta|_max, and
! various quantities at every time step to monitor the flow evolution.

implicit none

 !Local variables:
double precision:: wk(0:ny,0:nxm1)
double precision:: ro,fr,hmin,hmax

!----------------------------------------------------------------------
 !Compute Rossby number:
ro=maxval(abs(zz))/cof

 !Compute Froude number:
wk=(uu**2+vv**2)/(one+hh)
fr=sqrt(maxval(wk))/cgw

 !Compute h_min and h_max:
hmin=minval(hh)
hmax=maxval(hh)

 !Write data to ro-fr-hm.asc in the evolution subdirectory:
write(16,'(1x,f12.5,4(1x,f12.8))') t,ro,fr,hmin,hmax

 !Increment the integral of |zeta|_max:
twist=twist+dt*maxval(abs(zz))

return
end subroutine diagnose

!=======================================================================

subroutine savegrid(igrids)

! Saves various fields and energy components at the desired save time

implicit none

 !Passed variable:
integer:: igrids

 !Local variables:
double precision:: wk(0:ny,0:nxm1)
double precision:: ekin,epot,etot

!---------------------------------------------------------------
 !Compute kinetic energy:
wk=(one+hh)*(uu**2+vv**2)
call average(wk,ekin)
ekin=f12*domarea*ekin

 !Compute potential energy:
wk=hh**2
call average(wk,epot)
epot=f12*domarea*csq*epot

 !Compute total energy:
etot=ekin+epot

 !Write data to ecomp.asc in evolution subdirectory:
write(15,'(f13.6,3(1x,f16.9))') t,ekin,epot,etot
write(*,'(a,f13.6,a,f13.6)') ' t = ',t,'  E_tot = ',etot

!---------------------------------------------------------------
 !Write various gridded fields to direct access files:
write(31,rec=igrids) t,qq
wk=aa
call revfft(nyp1,nx,wk,xtrig,xfactors)
write(32,rec=igrids) t,wk
wk=bb
call revfft(nyp1,nx,wk,xtrig,xfactors)
write(33,rec=igrids) t,wk
write(34,rec=igrids) t,hh
write(35,rec=igrids) t,zz
write(36,rec=igrids) t,uu
write(37,rec=igrids) t,vv

return
end subroutine savegrid

!=======================================================================
      
subroutine savecont(irec)

! Saves PV contours for post-processing and imaging

implicit none

 !Passed variable:
integer:: irec

 !Local variables:
double precision:: wk(0:ny,0:nxm1)
character(len=3):: pind

!---------------------------------------------------------------
write(*,'(a,f12.5)') ' Saving contours at t = ',t
write(pind(1:3),'(i3.3)') irec

 !Write contours to the contours subdirectory:
write(80,'(i8,1x,i9,1x,f12.5,1x,f16.12)') nq,nptq,t,dq

 !Save residual needed to build ultra-fine-grid PV for plotting purposes:
call con2grid(qc)
wk=qq-qc
write(83,rec=irec) t,wk

 !Save PV contours:
open(81,file='contours/qqindex'//pind,form='unformatted',status='replace')
write(81) npq(1:nq),i1q(1:nq),indq(1:nq)
close(81)

open(82,file='contours/qqnodes'//pind,form='unformatted',status='replace')
write(82) xq(1:nptq),yq(1:nptq)
close(82)

return
end subroutine savecont

!=======================================================================

 !Main end module
end module evolution
