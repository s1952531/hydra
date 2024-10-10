module spectral

!   To initialise, use

! call init_spectral

!   To find the (gridded) dimensionless height anomaly (hh),
!   velocity field (uu,vv) and relative vorticity (zz), use

! call invert(qq,aa,bb,uum,uup,qoff,hh,uu,vv,zz)

!   where qq is the PV, (aa,bb) is the acceleration, and uum & uup
!   are the zonal velocities at y = ymin & ymax. Note, only (aa,bb),
!   uum & uup are in semi-spectral space. Upon return the PV qq is
!   adjusted by a constant offset (qoff) to be consistent with the
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

 !Tridiagonal arrays needed in inversion (f*u equation for k = 0):
double precision:: htdu(nym1),etdu(nym1)

 !Tridiagonal arrays needed in inversion (phi equation for k > 0):
double precision:: htdp(nym1,nxm1),etdp(nym1,nxm1)

 !Tridiagonal arrays needed in semi-implicit time stepping:
double precision:: htdt(nym1,0:nxm1),etdt(nym1,0:nxm1)

 !Hyperdiffusion arrays:
double precision:: diss(0:nxm1),rdis(0:nxm1)

 !k^2, 1/k^2, k^2 + k^d^2 & 1/(k^2 + k^d^2):
double precision:: ksq(0:nxm1),ksqi(nxm1),kksq(0:nxm1),kksqi(0:nxm1)

 !For boundary velocity evolution:
double precision:: alpu(0:nxm1),betu(0:nxm1),divu(0:nxm1)

contains

!====================================================================
subroutine init_spectral
 !Initialises this module.
  
 !Declarations:
implicit none

double precision:: kk(0:nxm1),emkl(0:nxm1)
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

 !Define k^2, k^2 + k_d^2 and the latter's inverse:
ksq=rkx**2
kksq=ksq+kdsq
kksqi=one/kksq

 !Define 1/k^2 for k > 0:
ksqi=one/ksq(1:nxm1)

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

 !Used in semi-implicit time stepping of aa & bb:
rdis=dt2i+diss
 !Re-define diss for use in q_d, uum & uup evolution:
diss=two/(one+dt2*diss)

 !Used in semi-implicit time stepping of uum & uup:
kk=sqrt(kksq)
emkl=exp(-kk*elly)
fac=cof/elly
betu=two*fac*emkl/(one-emkl**2)
alpu=fac*(one+emkl**2)/(one-emkl**2)
divu=one/(rdis**2+ksq*cof**2*kksqi)

!---------------------------------------------------------------------
 !Set up tridiagonal inversion arrays for inversion of f*u when m = 0:
 !(Below ap = 1/dy^2 and kdsq = (f/c)^2)
a00=-two*ap
a0=a00-kdsq
htdu(1)=one/a0
etdu(1)=-ap*htdu(1)
do j=2,nym1
  htdu(j)=one/(a0+ap*etdu(j-1))
  etdu(j)=-ap*htdu(j)
enddo

!---------------------------------------------------------------------
 !Set up tridiagonal inversion arrays for inversion of phi when m > 0:
 !(Below kksq = (f/c)^2 + k^2)
do m=1,nxm1
  a0=a00-kksq(m)
  htdp(1,m)=one/a0
  etdp(1,m)=-ap*htdp(1,m)
  do j=2,nym1
    htdp(j,m)=one/(a0+ap*etdp(j-1,m))
    etdp(j,m)=-ap*htdp(j,m)
  enddo
enddo

!---------------------------------------------------------------------
 !Set up tridiagonal inversion arrays for acceleration time stepping:
 !(Below rdis = 2/dt + diss, where diss = C*f*(k/k_max)^(2p) is the
 !hyperviscosity operator defined above in this subroutine)
do m=0,nxm1
  a0=a00-kksq(m)-csqi*rdis(m)**2
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
subroutine invert(qq,aa,bb,uum,uup,qoff,hh,uu,vv,zz)
! Computes the velocity field and dimensionless height anomaly.
  
! Input:  PV (qq), acceleration (aa,bb), boundary zonal velocity
!         (uum at y = ymin & uup at y = ymax) and PV offset (qoff).
!         *** Note: only the input field qq is in physical space.

! Output: dimensionless height anomaly (hh), velocity field (uu,vv)
!         and relative vorticity (zz). The PV is also corrected by
!         a constant offset to be consistent with the constraint on
!         the global mean vorticity. All fields are in physical space.

! *** Important: a guess for hh must be provided (could use hh = 0).  

 !Declarations:
implicit none

 !Passed variables:
double precision:: qq(0:ny,0:nxm1),aa(0:ny,0:nxm1),bb(0:ny,0:nxm1)
double precision:: hh(0:ny,0:nxm1),uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
double precision:: zz(0:ny,0:nxm1),uum(0:nxm1),uup(0:nxm1)
double precision:: qoff

 !Local variables:
double precision,parameter:: ptol=1.e-9 !error tolerance for phi
double precision:: gg(0:ny,0:nxm1),pp(0:ny,0:nxm1)
double precision:: ppm(nxm1),ppp(nxm1)
double precision:: zuu(0:ny),zpp(0:ny),zqq(0:ny)
double precision:: perr,avgp
integer:: j,m

! Initialise by computing div(a,b) in semi-spectral space -> gg:
call xderiv(aa,pp)
call yderiv(bb,gg)
! The y derivative is set to zero at the boundaries but is not used.
gg=pp+gg
! *** Do not re-use gg in the while loop below ***

! Store boundary values of semi-spectral phi for m > 0 (i.e. k > 0):
ppm=pp(0 ,1:nxm1)*ksqi
ppp=pp(ny,1:nxm1)*ksqi
! On rhs above, pp = a_x in semi-spectral space and ksqi = 1/k^2.

!---------------------------------------------------------------------
! Iterate to find dimensionless height anomaly h and zonal velocity u:
perr=one
do while (perr .gt. ptol)
   !Offset PV by a constant to satisfy mean vorticity constraint:
  call adjustpv(hh,qq,qoff)

   !Define Q = f*(q-f)*(1+h) in semi-spectral space -> pp temporarily:
  pp=cof*(qq-cof)*(one+hh)
  call forfft(nyp1,nx,pp,xtrig,xfactors)
   !Store zonal part (for m = 0) in zqq for use below:
  zqq=pp(:,0)

   !Define r = Q - div(a,b) in semi-spectral space -> pp:
  pp=pp-gg
  
   !Solve for phi (for m > 0 only) in semi-spectral space -> pp:
   !(Below ap = 1/dy^2)
  do m=1,nxm1
    pp(0 ,m)=ppm(m)    !Known boundary value at lower edge
    do j=1,nym1
      pp(j,m)=(pp(j,m)-ap*pp(j-1,m))*htdp(j,m)
    enddo
    pp(ny,m)=ppp(m)    !Known boundary value at upper edge
    do j=nym1,1,-1
      pp(j,m)=etdp(j,m)*pp(j+1,m)+pp(j,m)
    enddo
     !Obtain f*u = -(b + phi_y) by central differences:
    do j=1,nym1
      uu(j,m)=-bb(j,m)-hglyi*(pp(j+1,m)-pp(j-1,m))
    enddo
  enddo

   !Solve for the x-independent part of f*u (i.e. m = 0) -> zuu:
   !Define source term for elliptic problem:
  do j=1,nym1
    zuu(j)=kdsq*bb(j,0)-hglyi*(zqq(j+1)-zqq(j-1))
  enddo
   !Above kdsq = (f/c)^2.

   !Solve tri-diagonal problem:
  zuu(0) =cof*uum(0)    !Known boundary value at lower edge
  do j=1,nym1
    zuu(j)=(zuu(j)-ap*zuu(j-1))*htdu(j)
  enddo
  zuu(ny)=cof*uup(0)    !Known boundary value at upper edge
  do j=nym1,1,-1
    zuu(j)=etdu(j)*zuu(j+1)+zuu(j)
  enddo

   !Solve for the x-independent part of phi -> zpp:
  zpp(0)=zero
  do j=1,ny
    zpp(j)=zpp(j-1)-hgly*(zuu(j)+zuu(j-1)+bb(j,0)+bb(j-1,0))
  enddo
   !Above, hgly = 1/(2*dy) while zqq contains Q for m = 0.
   !Compute average of zpp over y and subtract from zpp:
  avgp=dnyi*(sum(zpp(1:nym1))+f12*zpp(ny))
   !Above, dnyi = 1/ny.
  zpp=zpp-avgp

   !Compute f*u = -(b + phi_y) by central differences:
  do j=1,nym1
    zuu(j)=-bb(j,0)-hglyi*(zpp(j+1)-zpp(j-1))
  enddo

   !Insert zonal parts into pp & uu arrays:
  pp(:,0)=zpp
  uu(:,0)=zuu

   !Insert boundary values for f*u:
  uu(0 ,:)=cof*uum
  uu(ny,:)=cof*uup
  
   !FFT phi back to physical space as zz:
  zz=pp
  call revfft(nyp1,nx,zz,xtrig,xfactors)

   !Compute max abs error in phi:
  perr=maxval(abs(zz-csq*hh))

   !Store updated dimensionless height anomaly field, h:
  hh=csqi*zz
enddo
 !At this point, we have converged on h and w = f*u.

 !Compute zonal velocity from u = w/f:
uu=cofi*uu
 !FFT back to physical space:
call revfft(nyp1,nx,uu,xtrig,xfactors)

 !Compute meridional velocity, v = (a + phi_x)/f:
call xderiv(pp,vv)
 !Above, pp still contains phi in semi-spectral space.
vv=cofi*(aa+vv)
 !FFT back to physical space:
call revfft(nyp1,nx,vv,xtrig,xfactors)
 !Ensure identically zero boundary values for v:
vv(0 ,:)=zero
vv(ny,:)=zero

 !Compute relative vorticity:
zz=(one+hh)*qq-cof

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
