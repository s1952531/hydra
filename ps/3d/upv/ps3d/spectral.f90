module spectral

! Module containing subroutines for spectral operations, inversion, etc.

use constants
use sta2dfft

 !Common arrays, constants:

 !Spectral operators:
double precision:: rksq(ng,ng),hlapi(ng,ng)
double precision:: filt(ng,ng),diss(ng,ng)

 !Tridiagonal arrays for inverting the scaled QG operator with Neumann BCs:
double precision:: etdv(ng,ng,0:nz),htdv(ng,ng,0:nz)
double precision:: ap(ng,ng),apb(ng,ng)
double precision:: zlin(0:nz)

 !Tridiagonal arrays for z differentiation (Neumann BCs):
double precision:: etd1(1:nzm1),htd1(1:nzm1)

 !For 2D FFTs:
double precision:: hrkx(ng),hrky(ng),rk(ng)
double precision:: xtrig(2*ng),ytrig(2*ng)
integer:: xfactors(5),yfactors(5)

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=============================================================

subroutine init_spectral
! Initialises this module

implicit none

!Local variables:
double precision:: a0(ng,ng),a0b(ng,ng)
double precision:: rkmax,anu,rkfsqi
integer:: kx,ky,k,iz

!----------------------------------------------------------------------
 !Set up 2D FFTs:
call init2dfft(ng,ng,twopi,twopi,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define wavenumbers:
rk(1)=zero
do k=1,ng/2-1
  rk(k+1)   =hrkx(2*k)
  rk(ng+1-k)=hrkx(2*k)
enddo
rk(ng/2+1)=hrkx(ng)

!-----------------------------------------------------------------------
 !Define a variety of spectral operators:

 !Squared wavenumber array (rksq) and Gaussian (de-aliasing) filter (filt):
rkfsqi=one/(eps*dble(ng))**2
 !The filter assumes vertical to horizontal scales are in the ratio f/N
do ky=1,ng
   do kx=1,ng
      rksq(kx,ky)=rk(kx)**2+rk(ky)**2
      filt(kx,ky)=exp(-rkfsqi*rksq(kx,ky))
   enddo
enddo
 !Hyperviscosity coefficient (Dritschel, Gottwald & Oliver, JFM (2017)):
anu=cdamp*cof/rkmax**(2*nnu)
 !Hyperviscous operator (diss):
diss=anu*rksq**nnu

 !Inverse horizontal Laplacian (hlapi):
rksq(1,1)=one
hlapi=-one/rksq
rksq(1,1)=zero

!-----------------------------------------------------------------------
 !Tridiagonal coefficients depending only on kx and ky:
a0=-two*dzisq-f56*rksq
a0b=-dzisq-f13*rksq
ap=dzisq-f112*rksq
apb=dzisq-f16*rksq

 !Tridiagonal arrays for inverting the QG stretched Laplace operator
 !with Neumann boundary conditions (4th-order compact scheme):
htdv(:,:,0)=one/a0b
etdv(:,:,0)=-apb*htdv(:,:,0)
do iz=1,nzm1
   htdv(:,:,iz)=one/(a0+ap*etdv(:,:,iz-1))
   etdv(:,:,iz)=-ap*htdv(:,:,iz)
enddo
htdv(:,:,nz)=one/(a0b+apb*etdv(:,:,nzm1))

 !Horizontal mean part must be solved for separately:
htdv(1,1,:)=zero
etdv(1,1,:)=zero

 !Tridiagonal arrays for computing phi_z with given boundary values
 !(4th-order compact scheme):
htd1(1)=one/f23
etd1(1)=-f16*htd1(1)
do iz=2,nzm1
   htd1(iz)=one/(f23+f16*etd1(iz-1))
   etd1(iz)=-f16*htd1(iz)
enddo

 !Define zlin = z + D (vanishes at bottom, iz = 0):
do iz=0,nz
   zlin(iz)=dz*dble(iz)
enddo

return 
end subroutine init_spectral

!=================================================================

subroutine getpsi(b0,ss,sx,sy)
! Checks the converge of the Monge-Ampere problem for psi

implicit none

double precision,parameter:: toler=1.d-12

 !Passed arrays:
double precision:: b0(ng,ng)
double precision:: ss(ng,ng,0:nz),sx(ng,ng,0:nz),sy(ng,ng,0:nz)

 !Local variables:
double precision:: pp(ng,ng,0:nz)
double precision:: ssa(ng,ng,0:nz),sxa(ng,ng,0:nz),sya(ng,ng,0:nz)
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: sxx(ng,ng),syy(ng,ng),sxy(ng,ng)
double precision:: pz0(ng,ng)
double precision:: sserr,sxerr,syerr
integer:: iz,iter

!---------------------------------------------------------
! Convert b_0 to spectral space as wkb:
wkc=b0
call ptospc(ng,ng,wkc,wkb,xfactors,yfactors,xtrig,ytrig)

! Store (spectrally-truncated) b0:
pz0=filt*wkb
! This is the value of phi_z (spectral) at z = 0

!-----------------------------------------------------------------
! QG solution: source is zero only at the top boundary:
pp(:,:,nz)=-dzi*pz0*htdv(:,:,nz)
! This is a partial tridiagonal solve (dzi = 1/dz)
do iz=nzm1,0,-1
   pp(:,:,iz)=etdv(:,:,iz)*pp(:,:,iz+1)+pp(:,:,iz)
enddo

! Initially psi = phi:
ss=pp

! Compute grad_h{psi} -> (sx,sy):
do iz=0,nz
   call xderiv(ng,ng,hrkx,pp(:,:,iz),sx(:,:,iz))
   call yderiv(ng,ng,hrky,pp(:,:,iz),sy(:,:,iz))
enddo

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Iterate to converge on psi, psi_x and psi_y:
sserr=one
sxerr=one
syerr=one

! Save physical space copies for computing error below:
call spctop3d(ss,ssa,0,nz) 
call spctop3d(sx,sxa,0,nz) 
call spctop3d(sy,sya,0,nz) 

do while (sserr > toler .or. sxerr > toler .or. syerr > toler)
   sserr=zero
   sxerr=zero
   syerr=zero

   do iz=0,nz
      !Compute psi_xx:
      call xderiv(ng,ng,hrkx,sx(:,:,iz),wkb)
      call spctop(ng,ng,wkb,sxx,xfactors,yfactors,xtrig,ytrig)

      !Compute psi_yy:
      call yderiv(ng,ng,hrky,sy(:,:,iz),wkb)
      call spctop(ng,ng,wkb,syy,xfactors,yfactors,xtrig,ytrig)

      !Compute psi_xy:
      call yderiv(ng,ng,hrky,sx(:,:,iz),wkb)
      call spctop(ng,ng,wkb,sxy,xfactors,yfactors,xtrig,ytrig)

      !Compute Xi = 2*(psi_xx*psi_yy - psi_xy^2):
      wkc=two*(sxx*syy-sxy**2)
      call ptospc(ng,ng,wkc,wkb,xfactors,yfactors,xtrig,ytrig)

      !Finish calculation of psi = phi - nabla_h^{-2}{Xi}:
      ss(:,:,iz)=pp(:,:,iz)-hlapi*wkb

      !Compute psi_x and psi_y:
      call xderiv(ng,ng,hrkx,ss(:,:,iz),sx(:,:,iz))
      call yderiv(ng,ng,hrky,ss(:,:,iz),sy(:,:,iz))

      !Check convergence:
      wkb=ss(:,:,iz)
      call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
      sserr=max(sserr,maxval(abs(wkc-ssa(:,:,iz))))
      ssa(:,:,iz)=wkc

      wkb=sx(:,:,iz)
      call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
      sxerr=max(sxerr,maxval(abs(wkc-sxa(:,:,iz))))
      sxa(:,:,iz)=wkc

      wkb=sy(:,:,iz)
      call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
      syerr=max(syerr,maxval(abs(wkc-sya(:,:,iz))))
      sya(:,:,iz)=wkc
   enddo
   write(*,'(3(a,1p,e14.7))') ' E_psi = ',sserr, &
      '  E_psi_x = ',sxerr,'  E_psi_y = ',syerr
enddo

do iz=0,nz
   call xderiv(ng,ng,hrkx,sx(:,:,iz),wkb)
   call spctop(ng,ng,wkb,sxx,xfactors,yfactors,xtrig,ytrig) !psi_xx
   call yderiv(ng,ng,hrky,sy(:,:,iz),wkb)
   call spctop(ng,ng,wkb,syy,xfactors,yfactors,xtrig,ytrig) !psi_yy
   call yderiv(ng,ng,hrky,sx(:,:,iz),wkb)
   call spctop(ng,ng,wkb,sxy,xfactors,yfactors,xtrig,ytrig) !psi_xy
   wkc=two*(sxx*syy-sxy**2) !Xi
   call ptospc(ng,ng,wkc,wkb,xfactors,yfactors,xtrig,ytrig)
   wkb=filt*wkb
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   
   wkb=rksq*(ss(:,:,iz)-pp(:,:,iz))
   call spctop(ng,ng,wkb,wka,xfactors,yfactors,xtrig,ytrig)
   sserr=maxval(abs(wka-wkc))
   write(*,*) ' iz = ',iz,' max error in MA equation = ',sserr
enddo

! Converged:
ss=ssa
sx=sxa
sy=sya

return 
end subroutine getpsi

!=================================================================

subroutine qgbal(b0,bb,xxi,eta,zz,pp,uu,vv,gg)
! Assuming QG balance (hydrostatic and geostrophic), and given
! only the surface distribution of buoyancy (b0), this routine
! finds the scaled buoyancy anomaly b'/(f*N), and the scaled vorticity
! components xi/N, eta/N and zeta/f. Returns these quantities in bb, xxi,
! eta and zz, respectively. Also returns the perturbation pressure / f^2
! in pp, as well as u/f & v/f in uu & vv, and the static stability (Gamma)
! in gg.

! Note, b0 is the surface buoyancy scaled by f*N.
  
implicit none

 !Passed arrays:
double precision:: b0(ng,ng)
double precision:: pp(ng,ng,0:nz),zz(ng,ng,0:nz),bb(ng,ng,0:nz)
double precision:: sx(ng,ng,0:nz),sy(ng,ng,0:nz)
double precision:: uu(ng,ng,0:nz),vv(ng,ng,0:nz)
double precision:: gg(ng,ng,0:nz)

 !Local variables:
double precision:: pz(ng,ng,0:nz)
double precision:: xxi(ng,ng,0:nz),eta(ng,ng,0:nz)
double precision:: pz0(ng,ng),wkb(ng,ng),wkc(ng,ng)
double precision:: b0bar,q0,sfac
integer:: iz

!-------------------------------------------------------------------
! Horizontal mean value of b_0(x,y):
b0bar=dsumi*sum(b0)
! Initial guess for the PV anomaly, q_0:
q0=b0bar*depthi
! depthi = 1/D where D is the scaled depth (NH/f).

write(*,*)
write(*,*) ' QG solution for PV anomaly:'
write(*,'(2(a,1p,e14.7))') ' q_0 = ',q0

! Convert b_0 to spectral space as pz0:
wkc=b0
call ptospc(ng,ng,wkc,pz0,xfactors,yfactors,xtrig,ytrig)
! Factor used to convert a physical mean value to a spectral mean value:
sfac=pz0(1,1)/b0bar
! *** Do not re-use pz0 below!

!-------------------------------------------------------------------
! QG solution: source is zero only at the top boundary:
pp(:,:,nz)=-dzi*pz0*htdv(:,:,nz)
! This is a partial tridiagonal solve (dzi = 1/dz)
do iz=nzm1,0,-1
   pp(:,:,iz)=etdv(:,:,iz)*pp(:,:,iz+1)
enddo

! Compute phi_z -> pz:
call zderiv(pp,pz,pz0)
! Here, phi_z = 0 at iz = 0, and phi_z = pz0 (defined above) at iz = nz.

! Overwrite horizontal-mean part with q_0*(z+D):
do iz=0,nz
   pz(1,1,iz)=sfac*q0*zlin(iz)
enddo

! Save physical copy of pz in bb, the 3d buoyancy anomaly field:
call spctop3d(pz,bb,0,nz)

! Compute xi/N = -phi_zx, eta/N = -phi_zy and zeta/f = Lap_h(phi)
! and store in xxi, eta and zz:
do iz=0,nz
   call xderiv(ng,ng,hrkx,pz(:,:,iz),wkb)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   xxi(:,:,iz)=-wkc
   call yderiv(ng,ng,hrky,pz(:,:,iz),wkb)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   eta(:,:,iz)=-wkc
   wkb=-rksq*pp(:,:,iz)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   zz(:,:,iz)=wkc
enddo

! Compute static stability, Gamma (QG approximation of it):
gg=one-zz

! Add horizontal mean part of phi by vertically integrating \bar\phi_z:
pp(1,1,0)=zero
do iz=1,nz
   pp(1,1,iz)=pp(1,1,iz-1)+dz2*(pz(1,1,iz)+pz(1,1,iz-1))
enddo

! Compute u/f, v/f & p/f^2 and store in uu, vv & pp:
do iz=0,nz
   call yderiv(ng,ng,hrky,pp(:,:,iz),wkb)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   uu(:,:,iz)=-wkc

   call xderiv(ng,ng,hrky,pp(:,:,iz),wkb)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   vv(:,:,iz)=wkc

   wkb=pp(:,:,iz)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   pp(:,:,iz)=wkc
enddo

return 
end subroutine qgbal

!=================================================================

subroutine bcbal(b0,bb,px,py,zz,pp,rr,ss,gg)
! Assuming Charney-Bolin balance, delta = delta_t = 0 (delta = u_x + v_y)
! and given only the surface distribution of buoyancy (b0), this routine
! finds the scaled buoyancy anomaly b'/(f*N), and the scaled vorticity
! components xi/N, eta/N and zeta/f. Returns these quantities in bb, px,
! py and zz, respectively. Also returns the perturbation pressure / f^2
! in pp, as well as u/f & v/f in rr & ss, and the static stability (Gamma)
! in gg.

! Note, b0 is the surface buoyancy scaled by f*N.

implicit none

 !Maximum absolute difference in phi_z, psi or psi_z between iterates:
double precision,parameter:: toler=1.d-12

 !Passed arrays:
double precision:: b0(ng,ng)
double precision:: pp(ng,ng,0:nz),zz(ng,ng,0:nz),bb(ng,ng,0:nz)
double precision:: sx(ng,ng,0:nz),sy(ng,ng,0:nz)
double precision:: rr(ng,ng,0:nz),ss(ng,ng,0:nz)
double precision:: gg(ng,ng,0:nz)

 !Local variables:
double precision:: px(ng,ng,0:nz),py(ng,ng,0:nz),pz(ng,ng,0:nz)
double precision:: sz(ng,ng,0:nz),szx(ng,ng,0:nz),szy(ng,ng,0:nz)
double precision:: pza(ng,ng,0:nz),ssa(ng,ng,0:nz),sza(ng,ng,0:nz)
double precision:: szxa(ng,ng,0:nz),szya(ng,ng,0:nz)
double precision:: pz0(ng,ng),xxi(ng,ng)
double precision:: sxx(ng,ng),syy(ng,ng),sxy(ng,ng)
double precision:: szxx(ng,ng),szyy(ng,ng),szxy(ng,ng)
double precision:: wkb(ng,ng),wkc(ng,ng)
double precision:: zpzbar(0:nz)
double precision:: b0bar,q0,sfac
double precision:: pzerr,sserr,szerr
integer:: iz,iter

!-------------------------------------------------------------------
! Horizontal mean value of b_0(x,y):
b0bar=dsumi*sum(b0)
! Initial guess for the PV anomaly, q_0:
q0=b0bar*depthi
! depthi = 1/D where D is the scaled depth (NH/f).

write(*,*)
write(*,*) ' QG solution for PV anomaly:'
write(*,'(2(a,1p,e14.7))') ' q_0 = ',q0

! Convert b_0 to spectral space as pz0:
wkc=b0
call ptospc(ng,ng,wkc,pz0,xfactors,yfactors,xtrig,ytrig)
! Factor used to convert a physical mean value to a spectral mean value:
sfac=pz0(1,1)/b0bar

! Spectrally truncate for de-aliasing below:
pz0=filt*pz0
! *** Do not re-use pz0 below!

!-------------------------------------------------------------------
! Initial QG solution: source is zero only at the top boundary:
pp(:,:,nz)=-dzi*pz0*htdv(:,:,nz)
! This is a partial tridiagonal solve (dzi = 1/dz)
do iz=nzm1,0,-1
   pp(:,:,iz)=etdv(:,:,iz)*pp(:,:,iz+1)
enddo

! Compute grad_h{phi} -> (px,py):
do iz=0,nz
   call xderiv(ng,ng,hrkx,pp(:,:,iz),px(:,:,iz))
   call yderiv(ng,ng,hrky,pp(:,:,iz),py(:,:,iz))
enddo

! Compute phi_z -> pz:
call zderiv(pp,pz,pz0)
! Here, phi_z = 0 at iz = 0, and phi_z = pz0 (defined above) at iz = nz.

! Overwrite horizontal-mean part with q_0*(z+D):
do iz=0,nz
   pz(1,1,iz)=sfac*q0*zlin(iz)
enddo

! Save physical copy of pz in pza for testing convergence below:
call spctop3d(pz,pza,0,nz)

!-------------------------------------------------------------------
! Store psi and derivatives of psi (initially psi = phi):
ss=pp
sx=px
sy=py
sz=pz

! Save physical space copies of psi & psi_z for testing convergence:
call spctop3d(ss,ssa,0,nz) 
sza=pza

! Compute psi_zx and psi_zy:
do iz=0,nz
   call xderiv(ng,ng,hrkx,sz(:,:,iz),szx(:,:,iz))
   call yderiv(ng,ng,hrky,sz(:,:,iz),szy(:,:,iz))
enddo
! Save physical space copies in szxa & szya:
call spctop3d(szx,szxa,0,nz) 
call spctop3d(szy,szya,0,nz) 

!-------------------------------------------------------------------
! Find first correction to phi, using psi = phi:
do iz=0,nz
   !Compute Lap_h(psi) = psi_xx+psi_yy:
   wkb=-rksq*ss(:,:,iz)
   call spctop(ng,ng,wkb,szxy,xfactors,yfactors,xtrig,ytrig)

   !Compute 1/(1+Lap_h(psi)) and spectrally truncate:
   wkc=one/(one+szxy)
   call deal2d(wkc)

   !Compute psi_zx^2 + psi_zy^2 + Lap_h(psi)^2:
   wkb=szxa(:,:,iz)**2+szya(:,:,iz)**2+szxy**2
   call deal2d(wkb)

   !Add PV anomaly q0 to the above and multiply by 1/(1+Lap_h(psi)):
   wkb=(q0+wkb)*wkc

   !Convert source of Lap(phi) to spectral space as rr:
   call ptospc(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   rr(:,:,iz)=wkc

   ! Compute horizontal mean of phi_z*Lap_h(psi):
   zpzbar(iz)=dsumi*sum(pza(:,:,iz)*szxy)
   ! This is needed to update the horizontal mean value of phi_z below.
enddo

! Solve tridiagonal problem for vertical potential phi -> pp:
call qginvert(rr,pp,pz0)

! Compute grad_h{phi} -> (px,py):
do iz=0,nz
   call xderiv(ng,ng,hrkx,pp(:,:,iz),px(:,:,iz))
   call yderiv(ng,ng,hrky,pp(:,:,iz),py(:,:,iz))
enddo

! Compute phi_z -> pz:
call zderiv(pp,pz,pz0)

! Update PV anomaly:
q0=(b0bar+dsumi*sum(b0*szxy))*depthi
! Above, szxy = Lap_h(psi) at the top boundary.
  
! Overwrite horizontal-mean part:
do iz=0,nz
   pz(1,1,iz)=sfac*(q0*zlin(iz)-zpzbar(iz))
enddo

! Save physical copies of phi_z for testing convergence:
call spctop3d(pz,rr,0,nz)
pzerr=maxval(abs(rr-pza))
write(*,*) ' First nonlinear correction with psi = phi:'
write(*,'(2(a,1p,e14.7))') ' q_0 = ',q0,'  E_phi_z = ',pzerr
pza=rr

write(*,*)
write(*,*) '       q_0             E_phi_z           E_psi           E_psi_z'

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Iterate to converge on full solution:
sserr=one
szerr=one
do while (pzerr > toler .or. sserr > toler .or. szerr > toler)
   !--------------------------------------------------------------
   ! Iterate to converge on psi, psi_z and phi_z simultaneously:
   sserr=zero
   szerr=zero
   do iz=0,nz
      !Compute psi_xx:
      call xderiv(ng,ng,hrkx,sx(:,:,iz),wkb)
      call spctop(ng,ng,wkb,sxx,xfactors,yfactors,xtrig,ytrig)

      !Compute psi_yy:
      call yderiv(ng,ng,hrky,sy(:,:,iz),wkb)
      call spctop(ng,ng,wkb,syy,xfactors,yfactors,xtrig,ytrig)

      !Compute psi_xy:
      call yderiv(ng,ng,hrky,sx(:,:,iz),wkb)
      call spctop(ng,ng,wkb,sxy,xfactors,yfactors,xtrig,ytrig)

      !Compute Xi = 2*(psi_xx*psi_yy - psi_xy^2); store in xxi:
      wkc=two*(sxx*syy-sxy**2)
      xxi=wkc
      call ptospc(ng,ng,wkc,wkb,xfactors,yfactors,xtrig,ytrig)

      !Finish calculation of psi = phi - nabla_h^{-2}{Xi}:
      ss(:,:,iz)=pp(:,:,iz)-hlapi*wkb

      !Compute psi_x and psi_y:
      call xderiv(ng,ng,hrkx,ss(:,:,iz),sx(:,:,iz))
      call yderiv(ng,ng,hrky,ss(:,:,iz),sy(:,:,iz))

      !Check convergence:
      wkb=ss(:,:,iz)
      call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
      sserr=max(sserr,maxval(abs(wkc-ssa(:,:,iz))))
      ssa(:,:,iz)=wkc

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      !Compute psi_zxx:
      call xderiv(ng,ng,hrkx,szx(:,:,iz),wkb)
      call spctop(ng,ng,wkb,szxx,xfactors,yfactors,xtrig,ytrig)

      !Compute psi_zyy:
      call yderiv(ng,ng,hrky,szy(:,:,iz),wkb)
      call spctop(ng,ng,wkb,szyy,xfactors,yfactors,xtrig,ytrig)

      !Compute psi_zxy:
      call yderiv(ng,ng,hrky,szx(:,:,iz),wkb)
      call spctop(ng,ng,wkb,szxy,xfactors,yfactors,xtrig,ytrig)

      !Compute Xi_z where Xi = 2*(psi_xx*psi_yy - psi_xy^2):
      wkc=two*(szxx*syy+sxx*szyy-two*sxy*szxy)
      call ptospc(ng,ng,wkc,wkb,xfactors,yfactors,xtrig,ytrig)

      !Finish calculation of psi_z = phi_z - nabla_h^{-2}{Xi_z}:
      sz(:,:,iz)=pz(:,:,iz)-hlapi*wkb

      !Compute psi_zx and psi_zy:
      call xderiv(ng,ng,hrkx,sz(:,:,iz),szx(:,:,iz))
      call yderiv(ng,ng,hrky,sz(:,:,iz),szy(:,:,iz))

      !Check convergence:
      wkb=sz(:,:,iz)
      call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
      szerr=max(szerr,maxval(abs(wkc-sza(:,:,iz))))
      sza(:,:,iz)=wkc

      !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      !Compute phi_zx:
      call xderiv(ng,ng,hrkx,pz(:,:,iz),wkb)
      call spctop(ng,ng,wkb,sxx,xfactors,yfactors,xtrig,ytrig)

      !Compute phi_zy:
      call yderiv(ng,ng,hrky,pz(:,:,iz),wkb)
      call spctop(ng,ng,wkb,syy,xfactors,yfactors,xtrig,ytrig)
      
      !Compute Lap_h(phi) = phi_xx+phi_yy:
      wkb=-rksq*pp(:,:,iz)
      call spctop(ng,ng,wkb,sxy,xfactors,yfactors,xtrig,ytrig)

      !Compute Lap_h(psi) = psi_xx+psi_yy:
      wkb=-rksq*ss(:,:,iz)
      call spctop(ng,ng,wkb,szxy,xfactors,yfactors,xtrig,ytrig)

      !Compute 1/(1+Lap_h(psi)) and spectrally truncate:
      wkc=one/(one+szxy)
      call deal2d(wkc)

      !Compute Xi + phi_zx*psi_zx + phi_zy*psi_zy + Lap_h(phi)*Lap_h(psi):
      wkb=xxi+sxx*szxa(:,:,iz)+syy*szya(:,:,iz)+sxy*szxy
      call deal2d(wkb)

      !zz is free to store Lap_h(psi), presently in szxy:
      zz(:,:,iz)=szxy

      !Add PV anomaly q0 to the above and multiply by 1/(1+Lap_h(psi)):
      wkb=(q0+wkb)*wkc

      !Convert source of Lap(phi) to spectral space as rr:
      call ptospc(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
      rr(:,:,iz)=wkc

      ! Compute horizontal mean of phi_z*Lap_h(psi):
      zpzbar(iz)=dsumi*sum(pza(:,:,iz)*zz(:,:,iz))
      ! This is needed to update the horizontal mean value of phi_z below.
   enddo

   ! Save physical space copies of psi_zx & psi_zy in szxa & szya:
   call spctop3d(szx,szxa,0,nz)
   call spctop3d(szy,szya,0,nz)

   !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   ! Solve tridiagonal problem for vertical potential phi -> pp:
   call qginvert(rr,pp,pz0)

   ! Compute grad_h{phi} -> (px,py):
   do iz=0,nz
      call xderiv(ng,ng,hrkx,pp(:,:,iz),px(:,:,iz))
      call yderiv(ng,ng,hrky,pp(:,:,iz),py(:,:,iz))
   enddo

   ! Compute phi_z -> pz:
   call zderiv(pp,pz,pz0)

   ! Update PV anomaly:
   q0=(b0bar+dsumi*sum(b0*zz(:,:,nz)))*depthi
   ! Above, zz(:,:,nz) = Lap_h(psi) at the top boundary.
  
   ! Overwrite horizontal-mean part:
   do iz=0,nz
      pz(1,1,iz)=sfac*(q0*zlin(iz)-zpzbar(iz))
   enddo

   ! Save physical copy of pz in rr for testing convergence:
   call spctop3d(pz,rr,0,nz)
   pzerr=maxval(abs(rr-pza))
   pza=rr

   write(*,'(4(3x,1p,e14.7))') q0,pzerr,sserr,szerr
enddo
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
write(*,*)
write(*,*) ' Converged!'

! Compute static stability, Gamma:
do iz=0,nz
   !Compute phi_zx:
   call xderiv(ng,ng,hrkx,pz(:,:,iz),wkb)
   call spctop(ng,ng,wkb,sxx,xfactors,yfactors,xtrig,ytrig)

   !Compute phi_zy:
   call yderiv(ng,ng,hrky,pz(:,:,iz),wkb)
   call spctop(ng,ng,wkb,syy,xfactors,yfactors,xtrig,ytrig)

   !Compute phi_zx*psi_zx + phi_zy*psi_zy and spectrally truncate:
   wkb=sxx*szxa(:,:,iz)+syy*szya(:,:,iz)+sxy*szxy
   call deal2d(wkb)

   !Compute 1/(1+Lap_h(psi)) and spectrally truncate:
   wkc=one/(one+zz(:,:,iz))
   call deal2d(wkc)
      
   !Compute Gamma and spectrally truncate:
   wkb=(one+q0+wkb)*wkc
   call deal2d(wkb)
   gg(:,:,iz)=wkb
enddo

! pza already contains phi_z in physical space; use for b' = phi_z:
bb=pza

! Define xi/N = -psi_zx & eta/N = -psi_zy and store in px & py:
px=-szxa
py=-szya
! Note: zz already contains zeta/f.

! Add horizontal mean part of phi by vertically integrating \bar\phi_z:
pp(1,1,0)=zero
do iz=1,nz
   pp(1,1,iz)=pp(1,1,iz-1)+dz2*(pz(1,1,iz)+pz(1,1,iz-1))
enddo

! Compute p/f^2, u/f and v/f and store in (pp,rr,ss):
do iz=0,nz
   wkb=pp(:,:,iz)
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   pp(:,:,iz)=wkc

   wkb=sy(:,:,iz) !sy = dpsi/dy = -u/f
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   rr(:,:,iz)=-wkc

   wkb=sx(:,:,iz) !sx = dpsi/dx = v/f
   call spctop(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
   ss(:,:,iz)=wkc
enddo

return 
end subroutine bcbal

!=================================================================

subroutine zderiv(pp,pz,pztop)
! Computes the z derivative of a function pp given that its
! derivative is zero when iz = 0 and equal to pztop when iz = nz.
! Uses a 4th-order compact difference method.

implicit none

 !Passed arrays:
double precision:: pp(ng,ng,0:nz),pz(ng,ng,0:nz),pztop(ng,ng)

 !Local variable:
integer:: iz

!---------------------------------------------------------
pz(:,:,0)=zero
do iz=1,nzm1
   pz(:,:,iz)=hdzi*(pp(:,:,iz+1)-pp(:,:,iz-1))
enddo
pz(:,:,nz)=pztop

! Get interior derivatives by 4th-order compact difference solution:
pz(:,:,nzm1)=pz(:,:,nzm1)-f16*pztop
pz(:,:,1)=pz(:,:,1)*htd1(1)
do iz=2,nzm1
   pz(:,:,iz)=(pz(:,:,iz)-f16*pz(:,:,iz-1))*htd1(iz)
enddo
do iz=nzm2,1,-1
   pz(:,:,iz)=etd1(iz)*pz(:,:,iz+1)+pz(:,:,iz)
enddo

return
end subroutine zderiv

!=================================================================

subroutine qginvert(rr,pp,pz0)
! Inverts the QG operator on the source rr to obtain the solution
! pp in spectral space (at each height). Uses a 4th-order compact
! difference method. pz0 is psi_z (spectral) at the upper surface.

implicit none

 !Passed arrays:
double precision:: rr(ng,ng,0:nz),pp(ng,ng,0:nz),pz0(ng,ng)

 !Local variable:
integer:: iz

!---------------------------------------------------------
pp(:,:,0)=f13*rr(:,:,0)+f16*rr(:,:,1)
do iz=1,nzm1
   pp(:,:,iz)=f112*(rr(:,:,iz-1)+rr(:,:,iz+1))+f56*rr(:,:,iz)
enddo
pp(:,:,nz)=f13*rr(:,:,nz)+f16*rr(:,:,nzm1)-dzi*pz0
! dzi = 1/dz above.

pp(:,:,0)=pp(:,:,0)*htdv(:,:,0)
do iz=1,nzm1
   pp(:,:,iz)=(pp(:,:,iz)-ap*pp(:,:,iz-1))*htdv(:,:,iz)
enddo
pp(:,:,nz)=(pp(:,:,nz)-apb*pp(:,:,nzm1))*htdv(:,:,nz)

do iz=nzm1,0,-1
   pp(:,:,iz)=etdv(:,:,iz)*pp(:,:,iz+1)+pp(:,:,iz)
enddo

return
end subroutine qginvert

!=================================================================

subroutine jacob(aa,bb,cs)
! Computes the (xy) Jacobian of aa and bb and returns it in cs.
! aa and bb are in physical space while cs is in spectral space

! NOTE: aa and bb are assumed to be spectrally truncated (de-aliased).

implicit none

 !Passed arrays:
double precision:: aa(ng,ng),bb(ng,ng),cs(ng,ng)

 !Work arrays:
double precision:: ax(ng,ng),ay(ng,ng),bx(ng,ng),by(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng)

!---------------------------------------------------------
wkb=aa
call ptospc(ng,ng,wkb,wka,xfactors,yfactors,xtrig,ytrig)
 !Get derivatives of aa:
call xderiv(ng,ng,hrkx,wka,wkb)
call spctop(ng,ng,wkb,ax,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,wkb)
call spctop(ng,ng,wkb,ay,xfactors,yfactors,xtrig,ytrig)

wkb=bb
call ptospc(ng,ng,wkb,wka,xfactors,yfactors,xtrig,ytrig)
 !Get derivatives of bb:
call xderiv(ng,ng,hrkx,wka,wkb)
call spctop(ng,ng,wkb,bx,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,wkb)
call spctop(ng,ng,wkb,by,xfactors,yfactors,xtrig,ytrig)

wkb=ax*by-ay*bx
call ptospc(ng,ng,wkb,cs,xfactors,yfactors,xtrig,ytrig)
 !The output is *not* spectrally truncated!

return
end subroutine jacob

!=================================================================

subroutine divs(aa,bb,cs)
! Computes the divergence of (aa,bb) and returns it in cs.
! Both aa and bb in physical space but cs is in spectral space.

implicit none

 !Passed arrays:
double precision:: aa(ng,ng),bb(ng,ng)   !Physical
double precision:: cs(ng,ng)             !Spectral

 !Work arrays:
double precision:: wkp(ng,ng)            !Physical
double precision:: wka(ng,ng),wkb(ng,ng) !Spectral

!---------------------------------------------------------
wkp=aa
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call xderiv(ng,ng,hrkx,wka,wkb)

wkp=bb
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call yderiv(ng,ng,hrky,wka,cs)

cs=wkb+cs

return
end subroutine divs

!=================================================================

subroutine ptospc3d(fp,fs,izbeg,izend)
! Transforms a physical 3d field fp to spectral space (horizontally)
! as the array fs.

implicit none

 !Passed variables:
double precision:: fp(ng,ng,0:nz)  !Physical
double precision:: fs(ng,ng,0:nz)  !Spectral
integer:: izbeg,izend

 !Work arrays:
double precision:: wkp(ng,ng)  !Physical
double precision:: wks(ng,ng)  !Spectral
integer:: iz

!---------------------------------------------------------
do iz=izbeg,izend
   wkp=fp(:,:,iz)
   call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
   fs(:,:,iz)=wks
enddo

return
end subroutine ptospc3d

!=================================================================

subroutine spctop3d(fs,fp,izbeg,izend)
! Transforms a spectral 3d field fs to physical space (horizontally)
! as the array fp.

implicit none

 !Passed variables:
double precision:: fp(ng,ng,0:nz)  !Physical
double precision:: fs(ng,ng,0:nz)  !Spectral
integer:: izbeg,izend

 !Work arrays:
double precision:: wkp(ng,ng)  !Physical
double precision:: wks(ng,ng)  !Spectral
integer:: iz

!---------------------------------------------------------
do iz=izbeg,izend
   wks=fs(:,:,iz)
   call spctop(ng,ng,wks,wkp,xfactors,yfactors,xtrig,ytrig)
   fp(:,:,iz)=wkp
enddo

return
end subroutine spctop3d

!=================================================================

subroutine deal3d(fp)
! Filters (horizontally) a physical 3d field fp (overwrites fp).

implicit none

 !Passed variable:
double precision:: fp(ng,ng,0:nz)  !Physical

 !Local variables:
double precision:: wkp(ng,ng)  !Physical
double precision:: wks(ng,ng)  !Spectral
integer:: iz

do iz=0,nz
   wkp=fp(:,:,iz)
   call ptospc(ng,ng,wkp,wks,xfactors,yfactors,xtrig,ytrig)
   wks=filt*wks
   call spctop(ng,ng,wks,wkp,xfactors,yfactors,xtrig,ytrig)
   fp(:,:,iz)=wkp
enddo

return
end subroutine deal3d

!=================================================================

subroutine deal2d(fp)
! Filters (horizontally) a physical 2d field fp (overwrites fp).

implicit none

 !Passed variable:
double precision:: fp(ng,ng)  !Physical

 !Local variable:
double precision:: fs(ng,ng)  !Spectral

call ptospc(ng,ng,fp,fs,xfactors,yfactors,xtrig,ytrig)
fs=filt*fs
call spctop(ng,ng,fs,fp,xfactors,yfactors,xtrig,ytrig)

return
end subroutine deal2d
!===================================================================

end module spectral
