!#########################################################################
!  Re-initialises a flow with balanced fields obtained from the conditions 
!  delta_t=gamma_t=0 using data previously set up with a data generation
!  routine.  Assumes the previous data has delta = gamma = 0.

!  Starts with SW balance to get a smooth starting approximation.

!           Written 6/4/2018 by D G Dritschel @ St Andrews
!#########################################################################

program dgbalini

 !Import contants, parameters and common arrays:
use constants
use spectral

implicit none

 !Physical fields:
double precision:: qq(ng,ng),hh(ng,ng),dd(ng,ng),gg(ng,ng)
double precision:: uu(ng,ng),vv(ng,ng),zz(ng,ng)
double precision:: aa(ng,ng),bb(ng,ng),cc(ng,ng)
double precision:: htot(ng,ng),hx(ng,ng),hy(ng,ng),ht(ng,ng)
double precision:: ut(ng,ng),vt(ng,ng),ztn(ng,ng)
double precision:: wkg(ng,ng),wkh(ng,ng),wkm(ng,ng)
double precision:: wkp(ng,ng),wkq(ng,ng),wkz(ng,ng)
double precision:: badd(ng,ng)
double precision:: hhpre(ng,ng),ddpre(ng,ng),ggpre(ng,ng),zzpre(ng,ng)

 !Spectral fields:
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng),uds(ng,ng),vds(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng),wke(ng,ng)
double precision:: sds(ng,ng),sgs(ng,ng)

 !Other constants:
double precision,parameter:: tole=1.d-10
 !tole: relative energy norm error in successive iterates when finding
 !      hh, uu & vv from qq, dd & gg.  The energy error is computed from 
 !      <(u-u0)^2+(v-v0)^2+c^2*(h-h0)^2>/<u^2+v^2+c^2*h^2>
 !      where <:> means a domain average and (u0,v0,h0) is the previous
 !      guess for (u,v,h).
double precision:: qadd,qbar,fqbar,t,uio,vio
double precision:: ddrmserr,ggrmserr,toterr,toterrpre

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!----------------------------------------------------------------------
 !Read in gridded PV anomaly and convert to spectral space as qs:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)
 !Note: zz typically has zero domain average, whereas the actual
 !      PV anomaly may not since this is determined by the 
 !      requirement that the mean relative vorticity is zero;
 !      qr is corrected upon calling main_invert in spectral.f90

 !Convert to spectral space (zz is overwritten; the PV is recovered below):
call ptospc(ng,ng,zz,qs,xfactors,yfactors,xtrig,ytrig)

 !Ensure domain average qs is zero (this does not matter):
qs(1,1)=zero
 !Spectrally-truncate for use in de-aliasing:
qs=filt*qs

!----------------------------------------------------------------------
 !Read in gridded divergence and convert to spectral space as ds:
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,dd
close(11)

zz=dd
call ptospc(ng,ng,zz,ds,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
ds(1,1)=zero

!----------------------------------------------------------------------
 !Read in gridded acceleration divergence and convert to spectral space
 !as gs:
open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,gg
close(11)

zz=gg
call ptospc(ng,ng,zz,gs,xfactors,yfactors,xtrig,ytrig)
 !Domain average must be zero:
gs(1,1)=zero

 !Obtain initial dimensionless height anomaly (hh), velocity (uu,vv),
 !relative vorticity (zz) and adjusted PV anomaly consistent with 
 !zero domain averaged zz:
call main_invert(qs,ds,gs,hh,uu,vv,qq,zz)
 !Note: qs, ds & gs are in spectral space while 
 !      hh, uu, vv, qq and zz are in physical space.

!-----------------------------------------------------------------
!Iterate to find the SW balanced fields (corresponding to H = 0):

!Energy norm error (must be > tole to start):
toterrpre=one/small**2
toterr=f12
hhpre=hh
ddpre=dd
ggpre=gg
zzpre=zz
do while (toterr .gt. tole)
  !Obtain balanced estimate for gamma (gg):
  call jacob(uu,vv,wkp)
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  wkp=dd*uu
  wkq=dd*vv
  call divs(wkp,wkq,wka)
  gs=filt*(wka-two*wkb)
  wka=gs
  call spctop(ng,ng,wka,gg,xfactors,yfactors,xtrig,ytrig)
  ggrmserr=sum((gg-ggpre)**2)/(sum(ggpre**2)+small)

  !Obtain balanced estimate for delta (dd):
  wkp=hh*uu
  wkq=hh*vv
  call divs(wkp,wkq,wka)
  wkp=zz*uu
  wkq=zz*vv
  call divs(wkp,wkq,wkb)
  ds=helm*(cof*wkb-c2g2*wka)
  wka=ds
  call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)
  ddrmserr=sum((dd-ddpre)**2)/(sum(ddpre**2)+small)

  !Find height anomaly field (hh):
  htot=one+hh
  qadd=-dsumi*sum(qq*htot)
  qq=qq+qadd
  qbar=dsumi*sum(qq)
  wkp=cof*(qq+hh*(qq-qbar))-gg
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  fqbar=cof*qbar
  wka=filt*wkb/(opak-fqbar)
  call spctop(ng,ng,wka,hh,xfactors,yfactors,xtrig,ytrig)
   !wkp: corrected de-aliased height field (to be hh below)
  htot=one+hh

  !Obtain relative vorticity field (zz):
  wkp=qq*htot
  qadd=-dsumi*sum(wkp)
  qq=qq+qadd
  zz=htot*(cof+qq)-cof

  !Obtain velocity field (uu,vv):
  call ptospc(ng,ng,zz,wkb,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkb
  wkb=filt*wkb
  call spctop(ng,ng,wkb,zz,xfactors,yfactors,xtrig,ytrig)
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)
  wke=rlap*ds
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkc)
  wkb=wka-wkb
  wkd=wkc+wkd
  call spctop(ng,ng,wkb,uu,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,vv,xfactors,yfactors,xtrig,ytrig)

  !Add mean flow (uio,vio):
  uio=-sum(hh*uu)*dsumi
  vio=-sum(hh*vv)*dsumi
  uu=uu+uio
  vv=vv+vio

  !Compute overall error:
  toterr=f12*(ddrmserr+ggrmserr)

  write(*,*) ' relative delta error = ',ddrmserr
  write(*,*) ' relative gamma error = ',ggrmserr

  !If error is going up again, stop and save fields:
  if (toterrpre .lt. toterr) exit

  !Otherwise continue with another iteration:
  hhpre=hh
  ddpre=dd
  ggpre=gg
  zzpre=zz
  toterrpre=toterr
enddo

write(*,*) ' Minimum error = ',toterrpre

!-----------------------------------------------------------------
!Iterate to find the GN balanced fields from SW balanced fields:

!Energy norm error (must be > tole to start):
toterrpre=one/small**2
toterr=f12
hhpre=hh
ddpre=dd
ggpre=gg
zzpre=zz
do while (toterr .gt. tole)
  !Obtain balanced estimate for gamma (gg):
  call jacob(uu,vv,wkp)
  call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
  hx=dd*uu
  hy=dd*vv
  call divs(hx,hy,wka)
  gs=filt*(wka-two*wkc)
  wka=gs
  call spctop(ng,ng,wka,gg,xfactors,yfactors,xtrig,ytrig)
  ggrmserr=sum((gg-ggpre)**2)/(sum(ggpre**2)+small)

  !Obtain balanced estimate for delta (dd):
  wka=gs+two*wkc
  call spctop(ng,ng,wka,aa,xfactors,yfactors,xtrig,ytrig)
  aa=aa-two*dd**2
  call dealias(aa)
  wkp=hh
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  call xderiv(ng,ng,hrkx,wka,wkd)
  call spctop(ng,ng,wkd,hx,xfactors,yfactors,xtrig,ytrig)
  call yderiv(ng,ng,hrky,wka,wkd)
  call spctop(ng,ng,wkd,hy,xfactors,yfactors,xtrig,ytrig)
  htot=one+hh
  wkp=hh*uu
  wkq=hh*vv
  call divs(wkp,wkq,wkd)
  wkd=filt*wkd
  wke=c2g2*wkd
  call spctop(ng,ng,wkd,ht,xfactors,yfactors,xtrig,ytrig)
  ht=-ht-dd
  wkz=zz+cof
  wkp=zz*uu
  wkq=zz*vv
  call divs(wkp,wkq,wkd)
  call spctop(ng,ng,wkd,ztn,xfactors,yfactors,xtrig,ytrig)
  call jacob(aa,hh,wkp)
  call dealias(wkp)
  ztn=hbsq3*wkp*htot-ztn
  wkp=f12*(uu**2+vv**2)
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  wka=filt*wka
  call xderiv(ng,ng,hrkx,wka,wkd)
  call spctop(ng,ng,wkd,wkp,xfactors,yfactors,xtrig,ytrig)
  call yderiv(ng,ng,hrky,wka,wkd)
  call spctop(ng,ng,wkd,wkq,xfactors,yfactors,xtrig,ytrig)
  wkh=htot**2
  call dealias(wkh)
  wkm=htot*aa
  call dealias(wkm)
  wkh=wkh*wkm
  call ptospc(ng,ng,wkh,wka,xfactors,yfactors,xtrig,ytrig)
  wka=filt*wka
  call xderiv(ng,ng,hrkx,wka,wkd)
  call spctop(ng,ng,wkd,wkg,xfactors,yfactors,xtrig,ytrig)
  call yderiv(ng,ng,hrky,wka,wkd)
  call spctop(ng,ng,wkd,wkh,xfactors,yfactors,xtrig,ytrig)
  wkm=hbsq3/htot
  call dealias(wkm)
  ut=wkm*wkg-wkp+wkz*vv-csq*hx
  call dealias(ut)
  vt=wkm*wkh-wkq-wkz*uu-csq*hy
  call dealias(vt)
  call jacob(ut,vv,wkp)
  call jacob(uu,vt,wkq)
  wkp=wkp+wkq
  call dealias(wkp)
  bb=two*htot*wkp
  call dealias(bb)
  call ptospc(ng,ng,ztn,wkb,xfactors,yfactors,xtrig,ytrig)
  wkp=ht*aa
  call dealias(wkp)
  wkp=htot*(two*wkp+bb)
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  wka=cof*wkb+wke-hbsq3*rksq*wka
  wkp=htot*ht
  call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
  wkc=filt*wkc
  call xderiv(ng,ng,hrkx,wkc,wkd)
  call spctop(ng,ng,wkd,wkg,xfactors,yfactors,xtrig,ytrig)
  call yderiv(ng,ng,hrky,wkc,wkd)
  call spctop(ng,ng,wkd,wkh,xfactors,yfactors,xtrig,ytrig)
  wkp=aa*wkg+bb*hx
  wkq=aa*wkh+bb*hy
  call divs(wkp,wkq,wkb)
  ds=-helm*(wka+hbsq3*wkb)
  wka=ds
  call spctop(ng,ng,wka,dd,xfactors,yfactors,xtrig,ytrig)
  ddrmserr=sum((dd-ddpre)**2)/(sum(ddpre**2)+small)

  !Find height anomaly field (hh):
  wkp=dd**2
  call dealias(wkp)
  badd=gg-two*wkp
  call jacob(hh,dd,cc)
  cc=hbsq3*cc
  call dealias(cc)
  qadd=-dsumi*sum(qq*htot)
  qq=qq+qadd
  qbar=dsumi*sum(qq)
  call jacob(uu,vv,wkp)
  call dealias(wkp)
  bb=hbsq3*htot*(badd+two*wkp)
  call dealias(bb)
  hx=bb*hx
  hy=bb*hy
  call divs(hx,hy,wka)
  wkp=htot*bb
  call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
  wkc=wka-rksq*wkc
  wkp=cof*(qq+hh*(qq-qbar)-htot*cc)-gg
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  fqbar=cof*qbar
  wka=filt*(wkb+wkc)/(opak-fqbar)
  call spctop(ng,ng,wka,hh,xfactors,yfactors,xtrig,ytrig)
  htot=one+hh

  !Obtain relative vorticity field (zz):
  call jacob(hh,dd,cc)
  cc=hbsq3*cc
  call dealias(cc)
  wkp=qq*htot
  qadd=-dsumi*sum(wkp)
  qq=qq+qadd
  zz=htot*(cof+qq-cc)-cof

  !Obtain velocity field (uu,vv):
  call ptospc(ng,ng,zz,wkb,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkb
  wkb=filt*wkb
  call spctop(ng,ng,wkb,zz,xfactors,yfactors,xtrig,ytrig)
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)
  wke=rlap*ds
  call xderiv(ng,ng,hrkx,wke,wka)
  call yderiv(ng,ng,hrky,wke,wkc)
  wkb=wka-wkb
  wkd=wkc+wkd
  call spctop(ng,ng,wkb,uu,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,vv,xfactors,yfactors,xtrig,ytrig)

  !Add mean flow (uio,vio):
  uio=-sum(hh*uu)*dsumi
  vio=-sum(hh*vv)*dsumi
  uu=uu+uio
  vv=vv+vio

  !Compute overall error:
  toterr=f12*(ddrmserr+ggrmserr)

  write(*,*) ' relative delta error = ',ddrmserr
  write(*,*) ' relative gamma error = ',ggrmserr

  !If error is going up again, stop and save fields:
  if (toterrpre .lt. toterr) exit

  !Otherwise continue with another iteration:
  hhpre=hh
  ddpre=dd
  ggpre=gg
  zzpre=zz
  toterrpre=toterr
enddo

write(*,*) ' Minimum error = ',toterrpre

!-----------------------------------------------------------------
!Write data:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,dd
close(11)

open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,gg
close(11)

write(*,*)
write(*,*) ' Initial fields balanced and re-written.'

 !End main program
end program
!=======================================================================
