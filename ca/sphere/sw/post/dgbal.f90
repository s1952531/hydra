!######################################################################
!  Finds balanced fields given only the PV anomaly in qq_init.r8 or
!  in qq.r4, using first-order delta-gamma balance in which delta_t
!  and gamma_t are set to zero to diagnose all fields.  Here gamma
!  is the linearised acceleration divergence.

!  This can be used both for initialisation and for post-processing.

! Adapted from planar code on 11 Feb 2020 by D G Dritschel @ St Andrews
!######################################################################

program dgbal

 ! Import constants and parameters:
use constants
 ! Import spectral module:
use spectral

implicit none
double precision:: hh(nLatGridPts,nLongGridPts),dd(nLatGridPts,nLongGridPts),gg(nLatGridPts,nLongGridPts),uu(nLatGridPts,nLongGridPts),vv(nLatGridPts,nLongGridPts)
double precision:: qb(nLatGridPts,nLongGridPts),qq(nLatGridPts,nLongGridPts),zz(nLatGridPts,nLongGridPts),bb(nLatGridPts,nLongGridPts)
double precision:: qs(nLatGridPts,nLongGridPts),ds(nLatGridPts,nLongGridPts),gs(nLatGridPts,nLongGridPts)
double precision:: bs(nLatGridPts,nLongGridPts)
double precision:: ekin,epot,etot,angm,t,berro
real:: qqr4(nLatGridPts,nLongGridPts),tr4
integer:: iopt,loop,iread,nstep,i

!---------------------------------------------------------
 ! Initialise inversion constants and arrays:
call init_spectral

write(*,*) ' Choose one of the following options:'
write(*,*)
write(*,*) ' (0) initialise from data in qq_init.r8'
write(*,*) ' (1) post-process from data in qq.r4'
write(*,*)
write(*,*) ' Option?'
read(*,*) iopt
write(*,*)

if (iopt .eq. 0) then
   ! Only find balance at t = 0 to initialise:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(11,rec=1) t,qb
  close(11)

   ! If topographic forcing is present, read in initial state:
  if (isTopoForcing) then
    open(11,file='bb_init.r8',form='unformatted', &
          access='direct',status='old',recl=2*nbytes)
    read(11,rec=1) t,bb
    close(11)
  endif
  
   ! Define PV anomaly:
  do i=1,nLongGridPts
    qb(:,i)=qb(:,i)-corFreq
  enddo

   ! Find balanced fields:
  call balance(berro)

   ! Save data:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,qq
  close(11)

  open(11,file='dd_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,dd
  close(11)

  open(11,file='hh_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,hh
  close(11)

  open(11,file='gg_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,gg
  close(11)

  open(11,file='zz_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,zz
  close(11)

  write(*,*)
  write(*,*) ' The balanced fields are available in xx_init.r8 where'
  write(*,*) ' xx = dd, hh, gg, zz & qq (adjusted only by a constant).'
  write(*,*)
  
else
   ! Find balance at all times in qq.r4:
  open(30,file='evolution/qq.r4',form='unformatted', &
        access='direct',status='old',recl=nbytes)

  if (isTopoForcing) then
     !Also read in the current toporaphic forcing:
    open(40,file='evolution/bb.r4',form='unformatted', &
          access='direct',status='old',recl=nbytes)
  endif

  open(31,file='evolution/bqq.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(32,file='evolution/bdd.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(33,file='evolution/bhh.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(34,file='evolution/bgg.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(35,file='evolution/bzz.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(51,file='evolution/becomp.asc',status='replace')

  open(61,file='evolution/benorm.asc',status='replace')

   ! Read data at all times and balance:
  loop=0
  do
    loop=loop+1
    iread=0
    read(30,rec=loop,iostat=iread) tr4,qqr4
    if (iread .ne. 0) exit
    qb=dble(qqr4)
    t=dble(tr4)

    if (isTopoForcing) then
       ! Read in current topographic forcing:
      read(40,rec=loop) tr4,qqr4
      bb=dble(qqr4)
    endif
    
     ! Define PV anomaly:
    do i=1,nLongGridPts
      qb(:,i)=qb(:,i)-corFreq
    enddo

     ! Balance using only the PV anomaly at time t = t:
    write(*,'(a,f9.2)') ' *** Processing t = ',t
    call balance(berro)

     ! Write data:
    write(31,rec=loop) tr4,real(qq)
    write(32,rec=loop) tr4,real(dd)
    write(33,rec=loop) tr4,real(hh)
    write(34,rec=loop) tr4,real(gg)
    write(35,rec=loop) tr4,real(zz)

    write(51,'(f13.6,4(1x,f16.9))') t,ekin,epot,etot,angm

    write(61,'(f9.2,1x,f11.7)') t,log10(berro)
  enddo

  close(30)
  if (isTopoForcing) close(40)
  close(31)
  close(32)
  close(33)
  close(34)
  close(35)
  close(51)
  close(61)

  write(*,*)
  write(*,*) ' The balanced fields are available in bxx.r4 where xx = qq,'
  write(*,*) ' dd, hh, gg & zz.'
  write(*,*)
  write(*,*) ' Also, becomp.asc contains the energy components, while'
  write(*,*) ' benorm.asc contains log10(relative energy norm error).'
  write(*,*)

endif


 ! Internal subroutine definitions (inherit global variables):

contains

!=======================================================================

subroutine balance(berro)

implicit none

! Passed variable (returned):
double precision:: berro

! Local variables:
double precision:: wkp(nLatGridPts,nLongGridPts)
integer:: i

!-----------------------------------------------------------------
! Iterate to find balanced fields:
call iterate(berro)

!-----------------------------------------------------------------
! Invert PV again after convergence:

! Convert qb to semi-spectral space as qs:
qs=qb
call forfft(nLatGridPts,nLongGridPts,qs,trig,factors) 
! Invert:
call main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
 !Note: qs, ds & gs are in semi-spectral space while 
 !      hh, uu, vv, qq and zz are in physical space.
dd=ds
call revfft(nLatGridPts,nLongGridPts,dd,trig,factors)
gg=gs
call revfft(nLatGridPts,nLongGridPts,gg,trig,factors)

! Compute energy components, divided by mean depth H:
wkp=(one+hh)*(uu**2+vv**2)
ekin=twopi*average(wkp)
wkp=hh**2
epot=twopi*csq*average(wkp)
etot=ekin+epot
do i=1,nLongGridPts
  wkp(:,i)=clat*((one+hh(:,i))*uu(:,i)+omega*clat*hh(:,i))
enddo
angm=fourpi*average(wkp)

return
end subroutine balance

!=======================================================================

subroutine iterate(berro)

! Solves the equations resulting from setting delta_t = gamma_t = 0.
!
! Exploits the semi-implicit formulation, here modified for an
! infinite time step and no dissipation.

implicit none

! Number of steps to reach full PV anomaly field:
integer,parameter:: nsteps=10
double precision,parameter:: dstep=one/dble(nsteps)

! Error in energy norm below which solution is said to be converged:
double precision,parameter:: tole=1.d-10

! Passed variable (returned):
double precision:: berro

! Local variables:
double precision:: hho(nLatGridPts,nLongGridPts),uuo(nLatGridPts,nLongGridPts),vvo(nLatGridPts,nLongGridPts)
double precision:: dso(nLatGridPts,nLongGridPts),gso(nLatGridPts,nLongGridPts)
double precision:: dsp(nLatGridPts,nLongGridPts),gsp(nLatGridPts,nLongGridPts),hhp(nLatGridPts,nLongGridPts)
double precision:: sds(nLatGridPts,nLongGridPts),sgs(nLatGridPts,nLongGridPts)
double precision:: berr,qmult
integer:: step

!-----------------------------------------------------------------------
! Initialise:
ds=zero
gs=zero

dsp=zero
gsp=zero
hhp=zero

dso=zero
gso=zero

hho=zero
uuo=zero
vvo=zero

do step=1,nsteps
  qmult=dstep*dble(step)
  qs=qmult*qb
  call forfft(nLatGridPts,nLongGridPts,qs,trig,factors) 

   !Extrapolate previous solutions:
  ds=two*dso-dsp
  gs=two*gso-gsp
  hh=two*hho-hhp
  dsp=dso
  gsp=gso
  hhp=hho
  
  ! Iterate to find balance (minimise error between two iterates):
  berro=one
  do
     ! Invert PV, delta and gamma to obtain the velocity field:
    call main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
     ! Note: qs, ds & gs are in spectral space while 
     !       hh, uu, vv, qq and zz are in physical space.

     ! Compute energy norm error:
    berr=average((uu-uuo)**2+(vv-vvo)**2+csq*(hh-hho)**2)

    if (berro .lt. tole .and. berr .lt. tole) then
       ! Error below minimum specified:
      berro=berr
       ! Exit the loop
      exit
    endif

     ! Copy new fields into old ones for next iteration:
    uuo=uu
    vvo=vv
    hho=hh
    dso=ds
    gso=gs
    berro=berr

     ! Compute NL source terms for delta and gamma:
    call source(sds,sgs)

     ! Update divergence and acceleration divergence:
    call helminv(sgs,ds,gs)
    ds=-csqi*ds
    gs=-sds
  enddo
enddo

write(*,*) ' Converged!  Relative energy error = ',berro

return
end subroutine iterate

!=======================================================================

subroutine source(sds,sgs)

! Gets the source terms (sds,sgs) for delta and gamma (ds,gs) ---
! all in semi-spectral space.

! Note that (sds,sgs) only include the nonlinear terms for a 
! semi-implicit treatment, closely analogous to that described in 
! the appendix of Dritschel & Jalali (JFM, 2020).
  
implicit none

 ! Passed variables:
double precision:: sds(nLatGridPts,nLongGridPts),sgs(nLatGridPts,nLongGridPts)

 ! Local variables:
double precision:: wka(nLatGridPts,nLongGridPts),wkb(nLatGridPts,nLongGridPts),wkc(nLatGridPts,nLongGridPts),wkd(nLatGridPts,nLongGridPts)
double precision:: wke(nLatGridPts,nLongGridPts),wkf(nLatGridPts,nLongGridPts),wkp(nLatGridPts,nLongGridPts),wkq(nLatGridPts,nLongGridPts)
double precision:: wku(nLatGridPts,nLongGridPts),wkv(nLatGridPts,nLongGridPts)
double precision:: dd(nLatGridPts,nLongGridPts),rhs(nLatGridPts),avgval
integer:: i

!---------------------------------------------------------------
 !Sources for velocity & acceleration divergence:

 !Compute A = z*U + dV/d(lon) -> wka & B = z*V - dU/d(lon) -> wkb,
 !where U = u/r & V = v/r, while z = sin(lat) and r = cos(lat):
do i=1,nLongGridPts
  wku(:,i)=clati*uu(:,i)
  wkv(:,i)=clati*vv(:,i)
  wka(:,i)=wkv(:,i)
  wkb(:,i)=wku(:,i)
enddo
 !*** Do not re-use wku ***

 !Convert wka & wkb to semi-spectral space:
call forfft(nLatGridPts,nLongGridPts,wka,trig,factors) 
call forfft(nLatGridPts,nLongGridPts,wkb,trig,factors) 

 !Compute longitudinal derivatives:
call deriv(nLatGridPts,nLongGridPts,rk,wka,wkc)
call deriv(nLatGridPts,nLongGridPts,rk,wkb,wkd)

 !Recover physical fields:
call revfft(nLatGridPts,nLongGridPts,wkc,trig,factors) 
call revfft(nLatGridPts,nLongGridPts,wkd,trig,factors) 

 !Complete definition of A & B and define wkf = f*zeta & wkb = 2*f*beta*v:
do i=1,nLongGridPts
  wka(:,i)=slat*wku(:,i)+wkc(:,i) !slat = z = sin(lat)
  wkb(:,i)=slat*wkv(:,i)-wkd(:,i)
  wkf(:,i)=corFreq*zz(:,i)            !cof = f
  wkd(:,i)=dfb*vv(:,i)            !dfb = 2*f*beta
enddo
 !*** Do not re-use wka, wkb, wkf or wkd ***

 !Get physical space velocity divergence -> dd:
dd=ds
call revfft(nLatGridPts,nLongGridPts,dd,trig,factors) 
 !*** Do not re-use dd ***

 !Get physical space derivatives of divergence:
call deriv(nLatGridPts,nLongGridPts,rk,ds,wkv)
call revfft(nLatGridPts,nLongGridPts,wkv,trig,factors) ! wkv = d(delta)/d(lon)
call latder(dd,wkc)                 ! wkc = d(delta)/d(lat)

 !Compute the nonlinear part of the velocity divergence tendency (sds):
wkp=uu**2+vv**2 ! *** Do not re-use wkp ***
sds=two*(wka*(zz-wka)-wkb*(wkb+dd))-wkp-dd**2-wku*wkv-vv*wkc
 !wku = u/r on rhs above; *** wka & wkb now safe to re-use below ***

 !If Equivalent Barotropic, deal with nonlinear height term:
if (eqbarot) then
  bs=Rocpi*(one+hh)**Rocp-hh
   !Add effect of topographic forcing if present:
  if (isTopoForcing) bs=bs+bb
  call forfft(nLatGridPts,nLongGridPts,bs,trig,factors) ! *** Do not re-use bs ***
  call laplace(bs,wka)
  call revfft(nLatGridPts,nLongGridPts,wka,trig,factors)
  sds=sds-csq*wka
else
   !Standard SW case with R/c_p = 1; just add forcing if present:
  if (isTopoForcing) then
    bs=bb
    call forfft(nLatGridPts,nLongGridPts,bs,trig,factors)
    call laplace(bs,wka)
    call revfft(nLatGridPts,nLongGridPts,wka,trig,factors)
    sds=sds-csq*wka
  endif
endif

 !Convert to semi-spectral space and de-alias:
call dealias(sds)

 !Compute the nonlinear part of the acceleration divergence tendency (sgg):

 !First compute div((u,v)*Z) (immediately below wkf = Z = f*zeta):
wka=wkf
call forfft(nLatGridPts,nLongGridPts,wka,trig,factors)
call deriv(nLatGridPts,nLongGridPts,rk,wka,wkv)
call revfft(nLatGridPts,nLongGridPts,wkv,trig,factors) ! wkv = dZ/d(lon)
call latder(wkf,wkc)                ! wkc = dZ/d(lat)

 !Store div((u,v)*Z) + 2*f*beta*v in sgs temporarily:
sgs=dd*wkf+wku*wkv+vv*wkc+wkd
 !De-alias and convert to semi-spectral space:
call dealias(sgs)
 !*** wkf can now be re-used ***

 !Define B = c^2*h + (u^2 + v^2)/2 -> wkb:
wkb=csq*hh+f12*wkp
 !De-alias and convert to semi-spectral space:
call dealias(wkb)
 !Add effect of topographic forcing and/or eq. bt. nonlinearity if present:
if (isTopoForcing .or. eqbarot) wkb=wkb+csq*bs
 !Calculate dB/d(lon) -> wkf (keep this semi-spectral):
call deriv(nLatGridPts,nLongGridPts,rk,wkb,wkf)

 !Compute div((u,v)*h) -> wke:
wka=hh
call forfft(nLatGridPts,nLongGridPts,wka,trig,factors) 
call deriv(nLatGridPts,nLongGridPts,rk,wka,wkv)
call revfft(nLatGridPts,nLongGridPts,wkv,trig,factors) ! wkv = dh/d(lon)
call latder(hh,wkc)                 ! wkc = dh/d(lat)
wke=dd*hh+wku*wkv+vv*wkc
 !Note wku = u/r on rhs above

 !If thermal damping is present, add contribution:
if (thermal) wke=wke+rth*hh

 !Compute Laplacian of div((u,v)*h) in wke after de-aliasing:
call dealias(wke)     ! Now wke is in semi-spectral space
call laplace(wke,wka) ! wka is returned in semi-spectral space

 !Complete calculation of nonlinear part of gamma tendency:
sgs=csq*wka+fpole*wkf-sgs
 !Note: each part has been separately de-aliased.

!-----------------------------------------------------
 !Remove global mean values of sds & sgs:
rhs=sds(:,1)*clat
avgval=(f1112*(rhs(1)+rhs(nLatGridPts))+sum(rhs(2:nLatGridPtsMin1)))*rsumi
sds(:,1)=sds(:,1)-avgval

rhs=sgs(:,1)*clat
avgval=(f1112*(rhs(1)+rhs(nLatGridPts))+sum(rhs(2:nLatGridPtsMin1)))*rsumi
sgs(:,1)=sgs(:,1)-avgval

end subroutine source

 ! End main program
end program dgbal
!=======================================================================
