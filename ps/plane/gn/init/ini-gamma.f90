!#########################################################################
!  Generates gamma (acceleration divergence) from the fields of 
!  q_l = zeta - f*h_tilde, delta and gamma_l = f*zeta - c^2*Lap(h_tilde).  
!  Re-writes gg_init.r8, replacing gamma_l by gamma.

!           Written 18/6/2019 by D G Dritschel @ St Andrews
!#########################################################################

program gamma0

 !Import spectral module:
use spectral

implicit none
 !Tolerance used for solving for gamma_tilde:
double precision,parameter:: toler=1.d-11

 !Physical arrays:
double precision:: uu(ng,ng),vv(ng,ng),hh(ng,ng),htot(ng,ng)
double precision:: dd(ng,ng),qq(ng,ng),gg(ng,ng),aa(ng,ng)
double precision:: bb(ng,ng),hx(ng,ng),hy(ng,ng),wkp(ng,ng)

 !Spectral arrays:
double precision:: qs(ng,ng),ds(ng,ng),gs(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: wkd(ng,ng),wke(ng,ng),wkf(ng,ng),wkg(ng,ng)

 !Other local variables:
double precision:: t,uio,vio,gerr

!----------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Read q_l:
open(31,file='ql_init.r8',form='unformatted',access='direct', &
                 & status='old',recl=2*nbytes)
read(31,rec=1) t,qq
close(31)

 !Read delta:
open(31,file='dd_init.r8',form='unformatted',access='direct', &
                 & status='old',recl=2*nbytes)
read(31,rec=1) t,dd
close(31)

 !Read gamma_l:
open(31,file='gg_init.r8',form='unformatted',access='direct', &
                 & status='old',recl=2*nbytes)
read(31,rec=1) t,gg
close(31)

 !Define spectral quantities:
aa=qq
call ptospc(ng,ng,aa,qs,xfactors,yfactors,xtrig,ytrig)
aa=dd
call ptospc(ng,ng,aa,ds,xfactors,yfactors,xtrig,ytrig)
aa=gg
call ptospc(ng,ng,aa,gs,xfactors,yfactors,xtrig,ytrig)

 !Invert linearised PV for the depth anomaly in spectral space:
wka=helm*(cof*qs-gs)

 !Define relative vorticity, zeta:
wkb=qs+cof*wka

 !Invert Laplace operator on zeta & delta to define velocity:
wkc=rlap*wkb
wkd=rlap*ds

 !Calculate derivatives spectrally:
call xderiv(ng,ng,hrkx,wkd,wke)
call yderiv(ng,ng,hrky,wkd,wkf)
call xderiv(ng,ng,hrkx,wkc,wkd)
call yderiv(ng,ng,hrkx,wkc,wkg)

 !Define velocity components:
wke=wke-wkg
wkf=wkf+wkd

 !Bring quantities back to physical space:
call spctop(ng,ng,wka,hh,xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,wke,uu,xfactors,yfactors,xtrig,ytrig)
call spctop(ng,ng,wkf,vv,xfactors,yfactors,xtrig,ytrig)

 !Add mean flow (uio,vio):
uio=-sum(hh*uu)*dsumi
vio=-sum(hh*vv)*dsumi
uu=uu+uio
vv=vv+vio

 !Obtain fixed part of gamma_tilde:
call jacob(uu,vv,qq)
aa=gg+two*(qq-dd**2)

 !De-alias and keep in spectral space as wkg:
call ptospc(ng,ng,aa,wkg,xfactors,yfactors,xtrig,ytrig)
wkg=filt*wkg
 !Put physical copy in initial guess for gamma_tilde:
wka=adop*wkg
call spctop(ng,ng,wka,gg,xfactors,yfactors,xtrig,ytrig)

 !Store de-aliased term h_tilde*(2+h_tilde) for use below:
wkb=hh*(two+hh)
call dealias(wkb)

 !Iterate to find gamma_tilde (in the variable gg, re-used here):
write(*,*) ' Iterating to find gamma_tilde:'
gerr=one
do while (gerr .gt. toler)
   !Calculate B = (H^2/3)*(1+h)*gamma_tilde (store in bb and de-alias):
  bb=hbsq3*(one+hh)*gg
  call dealias(bb)

   !Calculate h_x and h_y:
  wkp=hh
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
   !wka is hh in spectral space
  call xderiv(ng,ng,hrkx,wka,wkd)
  call spctop(ng,ng,wkd,hx,xfactors,yfactors,xtrig,ytrig)
   !hx is dh/dx in physical space
  call yderiv(ng,ng,hrky,wka,wkd)
  call spctop(ng,ng,wkd,hy,xfactors,yfactors,xtrig,ytrig)
   !hy is dh/dy in physical space

   !Compute div(hx*B,hy*B) (to be put into wka, in spectral space):
  hx=bb*hx
  hy=bb*hy
  call divs(hx,hy,wka)

   !Compute (H^2/3)*h_tilde*(2+h_tilde)*gamma_tilde:
  wkp=hbsq3*wkb*gg
  call ptospc(ng,ng,wkp,wke,xfactors,yfactors,xtrig,ytrig)

   !Invert P operator to find new guess for gamma_tilde:
  wkc=adop*(wkg-rksq*wke+wka)

   !Convert to physical space and calculate error:
  call spctop(ng,ng,wkc,wkp,xfactors,yfactors,xtrig,ytrig)
  gerr=sqrt(sum((wkp-gg)**2)/sum(gg**2))

  write(*,'(a,f16.12)') ' Error in gamma_tilde = ',gerr

   !Re-store gamma_tilde:
  gg=wkp
enddo
write(*,*) ' Converged!'

 !Finish definition of gamma by subtracting 2*(J(u,v)-delta^2):
gg=gg-two*(qq-dd**2)
call dealias(gg)

 !Re-write gamma:
open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,gg
close(11)

write(*,*)
write(*,*) ' Acceleration divergence re-written to gg_init.r8'

 !End main program
end program gamma0
!=======================================================================
