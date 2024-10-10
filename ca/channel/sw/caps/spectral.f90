module spectral

!   To initialise, use

! call init_spectral

!   To find the (gridded) dimensionless height anomaly (hh),
!   x velocity component (uu), relative vorticity (zz), linearised
!   PV anomaly zeta - f*h (qt) and balanced y velocity (vb), use

! call invert(qq,vv,bs,uum,uup,qoff,hh,uu,zz,qt,vb)

!   where qq is the PV, vv & bs are the semi-spectral fields of v &
!   Dv/Dt, while uum & uup are the zonal velocities at y = ymin & ymax.
!   Note, uum & uup are in semi-spectral space. Upon return, the PV qq
!   is adjusted by a constant offset (qoff) to be consistent with the
!   constant average value of zz (see caps.f90).

use constants
use generic
use stafft

 !Declarations:
implicit none

 !Maximum x wavenumber:
integer,parameter:: nw=nx/2,nwm1=nw-1

 !Common arrays, constants:
double precision:: rkx(0:nxm1)
double precision:: xtrig(2*nx)
integer:: xfactors(5)

 !Spectral de-aliasing & Butterworth low-pass filters (in x only):
double precision:: filt(0:nxm1),bflo(0:nxm1)

 !Tridiagonal arrays needed in inversion (w = f*u equation):
double precision:: htdu(nym1,0:nxm1),etdu(nym1,0:nxm1)

 !Tridiagonal arrays needed to recover balanced part of v:
double precision:: htdv(nym1,0:nxm1),etdv(nym1,0:nxm1)

 !Tridiagonal arrays needed in semi-implicit time stepping:
double precision:: htdt(nym1,0:nxm1),etdt(nym1,0:nxm1)

 !Hyperdiffusion arrays and squared horizontal wavenumber:
double precision:: diss(0:nxm1),rdis(0:nxm1),ksq(0:nxm1)

 !For boundary velocity evolution:
double precision:: alpu,betu,divu(0:nxm1)

contains

!====================================================================
subroutine init_spectral
 !Initialises this module.
  
 !Declarations:
implicit none

double precision:: scx,rkmax,rkf,anu,fac,a00,a0
integer:: m,l,j

!---------------------------------------------------------------------
 !Set up FFTs:
call initfft(nx,xfactors,xtrig)

 !Define x wavenumbers:
scx=two*pi/dble(ellx)
do m=0,nw
  rkx(m)=scx*dble(m)
enddo
do m=1,nwm1
  rkx(nx-m)=rkx(m)
enddo

 !Define k^2:
ksq=rkx**2

!---------------------------------------------------------------------
 !Define spectral filters (in x only):
rkmax=scx*dble(nw)
rkf=f23*rkmax
anu=cdamp*cof/rkmax**(2*nnu)
fac=9.d0/dble(nw*nw)

 !Cosine series wavenumbers:
do m=1,nw
  diss(m)=anu*ksq(m)**nnu
  if (rkx(m) .gt. rkf) then
    filt(m)=zero
    bflo(m)=zero
  else
    filt(m)=one
     !Butterworth low-pass filter:
    bflo(m)=one/(one+(fac*ksq(m))**2)
  endif
enddo
 !Do same for sine series wavenumbers:
do m=1,nwm1
  l=nx-m
  filt(l)=filt(m)
  bflo(l)=bflo(m)
  diss(l)=diss(m)
enddo
 !Constant mode (m = 0):
filt(0)=one
bflo(0)=one
diss(0)=zero

 !Used in semi-implicit time stepping of v & b = Dv/Dt:
rdis=dt2i+diss
 !Re-define diss for use in q_d, uum & uup evolution:
diss=two/(one+dt2*diss)

 !Used in semi-implicit time stepping of uum & uup:
betu= cgw/sinh(kd*elly)
alpu=betu*cosh(kd*elly)
divu=one/(rdis**2+csq*ksq)

!---------------------------------------------------------------------
 !Set up tridiagonal inversion arrays for inversion of u:
 !(Below ap = 1/dy^2 and kdsq = k_d^2 = (f/c)^2:
a00=-two*ap-kdsq
do m=0,nxm1
  htdu(1,m)=one/a00
  etdu(1,m)=-ap*htdu(1,m)
  do j=2,nym1
    htdu(j,m)=one/(a00+ap*etdu(j-1,m))
    etdu(j,m)=-ap*htdu(j,m)
  enddo
enddo

!---------------------------------------------------------------------
 !Set up tridiagonal inversion arrays to obtain balanced part of v:
do m=0,nxm1
  a0=a00-ksq(m)
  htdv(1,m)=one/a0
  etdv(1,m)=-ap*htdv(1,m)
  do j=2,nym1
    htdv(j,m)=one/(a0+ap*etdv(j-1,m))
    etdv(j,m)=-ap*htdv(j,m)
  enddo
enddo

!---------------------------------------------------------------------
 !Set up tridiagonal inversion arrays for v & Dv/Dt time stepping:
 !(Below rdis = 2/dt + diss, where diss = C*f*(k/k_max)^(2p) is the
 !hyperviscosity operator defined above in this subroutine)
do m=0,nxm1
  a0=a00-ksq(m)-csqi*rdis(m)**2
  htdt(1,m)=one/a0
  etdt(1,m)=-ap*htdt(1,m)
  do j=2,nym1
    htdt(j,m)=one/(a0+ap*etdt(j-1,m))
    etdt(j,m)=-ap*htdt(j,m)
  enddo
enddo

return
end subroutine init_spectral

!====================================================================
subroutine invert(qq,vv,bs,uum,uup,qoff,hh,uu,zz,qt,vb)
! Computes the velocity field, dimensionless height anomaly and
! other field used for computing sources and for time stepping.
  
! Input:  PV (qq), v (vv), Dv/Dt (bs), boundary zonal velocity
!         (uum & uup at y = y_min & y_max) and PV offset (qoff).
!         *** Note: bs, uum & uup are in semi-spectral space.

! Output: dimensionless height anomaly (hh), x velocity component (uu)
!         relative vorticity (zz), linearised PV anomaly zeta - f*h (qt)
!         and balanced y velocity (vb). The PV is also corrected by a
!         constant offset to be consistent with the constraint on the
!         global mean vorticity. hh, uu & zz are in physical space.

! *** Important: a guess for hh must be provided (could use hh = 0).  

 !Declarations:
implicit none

 !Passed variables:
double precision:: qq(0:ny,0:nxm1),vv(0:ny,0:nxm1),bs(0:ny,0:nxm1)
double precision:: hh(0:ny,0:nxm1),uu(0:ny,0:nxm1),zz(0:ny,0:nxm1)
double precision:: qt(0:ny,0:nxm1),vb(0:ny,0:nxm1)
double precision:: uum(0:nxm1),uup(0:nxm1),qoff

 !Local variables:
double precision,parameter:: htol=1.e-10 !error tolerance for hh
double precision:: wk(0:ny,0:nxm1),rr(0:ny,0:nxm1),pp(0:ny,0:nxm1)
double precision:: herr,avgp,avgr,pm
integer:: j,m

!---------------------------------------------------------------------
! Compute f*v_x -> wk (do not modify in iteration below):
rr=cof*vv
call forfft(nyp1,nx,rr,xtrig,xfactors)
call xderiv(rr,wk)

!---------------------------------------------------------------------
! Iterate to find zonal velocity u and dimensionless height anomaly h:
herr=one
do while (herr .gt. htol)
   !Offset PV by a constant to satisfy mean vorticity constraint:
  call adjustpv(hh,qq,qoff)

   !Define Q = f*(q-f)*(1+h) in semi-spectral space -> rr temporarily:
  rr=cof*(qq-cof)*(one+hh)
  call forfft(nyp1,nx,rr,xtrig,xfactors)
  call dealias(rr,1)

   !Define r = f*v_x - Q in semi-spectral space -> rr:
  rr=wk-rr

   !Form r.h.s. of equation for w = f*u, i.e. r_y + k_d^2*b -> uu:
  call yderiv(rr,uu)
  uu=uu+kdsq*bs
  
   !Invert Helmholtz equation for w in semi-spectral space -> uu:
   !(Below ap = 1/dy^2)
  do m=0,nxm1
    uu(0 ,:)=cof*uum    !Known boundary value at lower edge
    do j=1,nym1
      uu(j,m)=(uu(j,m)-ap*uu(j-1,m))*htdu(j,m)
    enddo
    uu(ny,:)=cof*uup    !Known boundary value at upper edge
    do j=nym1,1,-1
      uu(j,m)=etdu(j,m)*uu(j+1,m)+uu(j,m)
    enddo
  enddo

   !Obtain new estimate for geopotential, c^2*h -> pp:
  do m=0,nxm1
    pp(0,m)=zero  !To contain phi' in the notes (zero at y = y_min)
    avgp=zero     !To contain <phi'>, the y average of phi'
    avgr=zero     !To contain <r>

     !Integrate using the trapezoidal rule (uu = w = f*u and bs = b):
    do j=1,ny
      pp(j,m)=pp(j-1,m)-hgly*(uu(j,m)+uu(j-1,m)+bs(j,m)+bs(j-1,m))
      avgp=avgp+pp(j,m)+pp(j-1,m)
      avgr=avgr+rr(j,m)+rr(j-1,m)
    enddo

     !Obtain phi at y = ymin (phi^{-}):
    pm=f12*dnyi*(ldsq*avgr-avgp)-pmcc*(uup(m)-uum(m))
     !Above, dnyi = 1/ny, ldsq = c^2/f^2 & pmcc = c^2/(f*L_y)

     !Add phi^{-} to phi' to complete calculation of geopotential phi:
    pp(:,m)=pp(:,m)+pm
  enddo

   !Compute dimensionless height anomaly from h = phi/c^2 -> rr:
  rr=csqi*pp
   !FFT back to physical space:
  call revfft(nyp1,nx,rr,xtrig,xfactors)

   !Compute max abs error:
  herr=maxval(abs(rr-hh))
   !Store updated height field:
  hh=rr
enddo
 !At this point, we have converged on h and w = f*u.

 !Compute zonal velocity from u = w/f:
uu=cofi*uu
 !FFT back to physical space:
call revfft(nyp1,nx,uu,xtrig,xfactors)

 !Compute relative vorticity:
zz=(one+hh)*qq-cof
 !Potentially de-alias if needed:
!call dealias(zz,0)

 !Compute linearised PV anomaly, zeta - f*h (qt):
qt=zz-cof*hh

 !Compute the balanced y velocity (vb) by inverting G{vb} = c^2*qt_x:
wk=qt
call forfft(nyp1,nx,wk,xtrig,xfactors)
call xderiv(wk,vb)
do m=0,nxm1
  vb(0,m)=zero
  vb(1,m)=vb(1,m)*htdv(1,m)
  do j=2,nym1
    vb(j,m)=(vb(j,m)-ap*vb(j-1,m))*htdv(j,m)
  enddo
  vb(ny,m)=zero
  do j=nym2,1,-1
    vb(j,m)=etdv(j,m)*vb(j+1,m)+vb(j,m)
  enddo
enddo
call revfft(nyp1,nx,vb,xtrig,xfactors)

return
end subroutine invert

!====================================================================
subroutine adjustpv(hh,qq,qoff)

! Adjusts the PV to ensure the domain-mean vorticity has the
! correct value (see definition of qoff in caps.f90)

  !Declarations:
implicit none

 !Passed variables:
double precision:: hh(0:ny,0:nxm1),qq(0:ny,0:nxm1)
double precision:: qoff

 !Local variables:
double precision:: za(0:ny,0:nxm1)
double precision:: qb

za=(one+hh)*qq
call average(za,qb)
qb=qoff-qb
qq=qq+qb

return
end subroutine adjustpv

!====================================================================
subroutine gradient(ff,fx,fy)

! Calculates the gradient of a field f, passed as ff in semi-spectral
! space. The result, (f_x,f_y), is returned as (fx,fy) in physical space.

 !Declarations:
implicit none

 !Passed variables:
double precision:: ff(0:ny,0:nxm1),fx(0:ny,0:nxm1),fy(0:ny,0:nxm1)

!---------------------------------------------------------------------
 !Obtain f_x spectrally:
call xderiv(ff,fx)
 !Return to physical space:
call revfft(nyp1,nx,fx,xtrig,xfactors)

 !Obtain f_y by centred differences (ignoring boundaries - not needed):
call yderiv(ff,fy)
 !Return to physical space:
call revfft(nyp1,nx,fy,xtrig,xfactors)

return
end subroutine gradient

!====================================================================
subroutine xderiv1d(ff,fx)

! Calculates the x derivative of a 1D array f, passed as ff in spectral
! space. The result, f_x, is returned as fx, also in spectral space.

 !Declarations:
implicit none

 !Passed variables:
double precision:: ff(0:nxm1),fx(0:nxm1)

 !Local variables:
integer:: m,l

!---------------------------------------------------------------------
 !Obtain f_x spectrally:
fx( 0)=zero
fx(nw)=zero
do m=1,nwm1
  l=nx-m
  fx(m)=-rkx(m)*ff(l)
  fx(l)= rkx(m)*ff(m)
enddo

return
end subroutine xderiv1d

!====================================================================
subroutine xderiv(ff,fx)

! Calculates the x derivative of a field f, passed as ff in semi-spectral
! space. The result, f_x, is returned as fx, also in semi-spectral space.

 !Declarations:
implicit none

 !Passed variables:
double precision:: ff(0:ny,0:nxm1),fx(0:ny,0:nxm1)

 !Local variables:
integer:: m,l

!---------------------------------------------------------------------
 !Obtain f_x spectrally:
fx(:, 0)=zero
fx(:,nw)=zero
do m=1,nwm1
  l=nx-m
  fx(:,m)=-rkx(m)*ff(:,l)
  fx(:,l)= rkx(m)*ff(:,m)
enddo

return
end subroutine xderiv

!====================================================================
subroutine yderiv(ff,fy)

! Calculates the y derivative of a field f, passed as ff in either
! physical or semi-spectral space.  The result, f_y, is returned
! as fy in the same space.

! *** Note: the boundary values are set to zero as they are normally
!     multiplied by the meridional velocity vv, which is zero there.  

 !Declarations:
implicit none

 !Passed variables:
double precision:: ff(0:ny,0:nxm1),fy(0:ny,0:nxm1)

 !Local variable:
integer:: j

!---------------------------------------------------------------------
 !Obtain f_y by centred differences (ignoring boundaries - not needed):
fy(0,:)=zero
do j=1,ny-1
  fy(j,:)=hglyi*(ff(j+1,:)-ff(j-1,:))
enddo
fy(ny,:)=zero
 !Above, hglyi = 1/(2*dy).

return
end subroutine yderiv

!====================================================================
subroutine dealias(ff,iopt)

! Spectrally truncates a field ff (in x wavenumber space only).

! If iopt = 0, ff is passed in physical space
! If iopt = 1, ff is passed in semi-spectral space
  
 !Declarations:
implicit none

 !Passed variables:
double precision:: ff(0:ny,0:nxm1)
integer:: iopt

 !Local variable:
integer:: j

!---------------------------------------------------------------------
 !Transform in x if ff is passed in physical space: 
if (iopt .eq. 0) call forfft(nyp1,nx,ff,xtrig,xfactors)

 !Apply de-aliasing filter:
do j=0,ny
  ff(j,:)=filt*ff(j,:)
enddo

 !Transform back if ff is expected in physical space: 
if (iopt .eq. 0) call revfft(nyp1,nx,ff,xtrig,xfactors)

return
end subroutine dealias

!====================================================================
subroutine lowpass(ff)

! Low-pass filters a semi-spectral field ff (in both x wavenumber
! space and by repeated 1-2-1 averages in y, while preserving
! boundary values (except as modified by the x filtering).
  
 !Declarations:
implicit none

 !Passed variable (semi-spectral):
double precision:: ff(0:ny,0:nxm1)

 !Local variables:
double precision:: wk(0:ny,0:nxm1)
integer:: j

!---------------------------------------------------------------------
 !Apply spectral low-pass filter in x:
do j=0,ny
  ff(j,:)=bflo*ff(j,:)
enddo

 !Filter in y by two successive 1-2-1 averages preserving boundary values:
wk( 0,:)=ff( 0,:)
wk(ny,:)=ff(ny,:)
do j=1,ny-1
  wk(j,:)=f14*(ff(j-1,:)+ff(j+1,:))+f12*ff(j,:)
enddo
do j=1,ny-1
  ff(j,:)=f14*(wk(j-1,:)+wk(j+1,:))+f12*wk(j,:)
enddo

return
end subroutine lowpass

end module spectral
