module spectral

!   To initialise, use

! call init_spectral

!   To find the (gridded) dimensionless height anomaly (hh), use

! call height(qq,gg,uum,uup,qoff,hh,zz)

!   where qq is the PV, gg is the acceleration divergence, while
!   uum & uup are the zonal velocities at y = ymin & ymax, all in
!   semi-spectral space.  As a bi-product, the relative vorticity zz
!   is also returned (in semi-spectral space).  Note that the PV qq
!   is adjusted by a constant offset (qoff) to be consistent with 
!   the fixed average value of the relative vorticity. 

!   To find the (gridded) velocity field (uu,vv), use

! call velocity(dd,zz,uum,uup,uu,vv)

!   where dd is the divergence in semi-spectral space, and the other
!   variables are the same as above.  On return, dd & zz are converted
!   to physical space.

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

 !Block-tridiagonal arrays needed to obtain velocity field
double precision:: a0b(0:ny,nwm1),b0b(0:ny,nwm1)
double precision:: et1(0:ny,nwm1),et2(0:ny,nwm1)

 !Tridiagonal arrays needed to obtain height field:
double precision:: hth(0:ny,0:nxm1),eth(0:ny,0:nxm1)

 !Tridiagonal arrays needed in semi-implicit time stepping:
double precision:: htd(0:ny,0:nxm1),etd(0:ny,0:nxm1)

 !Hyperdiffusion arrays:
double precision:: diss(0:nxm1),rdis(0:nxm1)

contains

!====================================================================
subroutine init_spectral

 !Declarations:
implicit none

double precision:: scx,a00,a0,b0,deti
double precision:: rkmax,rkf,anu,fac
integer:: m,l,j,nhyp

!---------------------------------------------------------------------
 !Set up FFTs:
call initfft(nx,xfactors,xtrig)

 !Define x wavenumbers:
scx=twopi/dble(ellx)
do m=0,nw
  rkx(m)=scx*dble(m)
enddo
do m=1,nwm1
  rkx(nx-m)=rkx(m)
enddo

!---------------------------------------------------------------------
 !Define spectral filters (in x only):
rkmax=scx*dble(nw)
rkf=f23*rkmax
fac=one/dble(nw/3)
nhyp=2*nnu
anu=cdamp*cof/rkmax**nhyp

 !Cosine series wavenumbers:
do m=1,nw
  diss(m)=anu*rkx(m)**nhyp
  if (rkx(m) .gt. rkf) then
    filt(m)=zero
    bflo(m)=zero
  else
    filt(m)=one
     !Butterworth low-pass filter:
    bflo(m)=one/(one+(fac*rkx(m))**4)
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

 !Used in semi-implicit time stepping of cc & gg:
rdis=dt2i+diss
 !Re-define diss for use in q_d, uum & uup evolution:
diss=two/(one+dt2*diss)

!---------------------------------------------------------------------
 !Set up block-tridiagonal inversion arrays for velocity calculation:
 !(Below ap = 1/dy^2, ape = 2/dy^2 and c0 = 2/dy)
do m=1,nwm1
  a0=-ape-rkx(m)**2
  b0=c0*rkx(m)

  deti=one/(a0**2-b0**2)
  a0b(0,m)=a0*deti
  b0b(0,m)=b0*deti
  et1(0,m)=-ape*a0b(0,m)
  et2(0,m)= ape*b0b(0,m)

  do j=1,ny-1
    a0b(j,m)=a0+et1(j-1,m)*ap
    b0b(j,m)=   et2(j-1,m)*ap
    deti=one/(a0b(j,m)**2-b0b(j,m)**2)
    a0b(j,m)=a0b(j,m)*deti
    b0b(j,m)=b0b(j,m)*deti
    et1(j,m)=-ap*a0b(j,m)
    et2(j,m)= ap*b0b(j,m)
  enddo

  a0b(ny,m)= a0+et1(ny-1,m)*ape
  b0b(ny,m)=-b0+et2(ny-1,m)*ape
  deti=one/(a0b(ny,m)**2-b0b(ny,m)**2)
  a0b(ny,m)=a0b(ny,m)*deti
  b0b(ny,m)=b0b(ny,m)*deti
enddo

!---------------------------------------------------------------------
 !Set up tridiagonal inversion arrays for c^2*height anomaly:
 !(Below ap = 1/dy^2, ape = 2/dy^2, kdsq = (f/c)^2)
a00=-ape-kdsq
do m=0,nxm1
  a0=a00-rkx(m)**2
  hth(0,m)=one/a0
  eth(0,m)=-ape*hth(0,m)
  do j=1,ny-1
    hth(j,m)=one/(a0+ap*eth(j-1,m))
    eth(j,m)=-ap*hth(j,m)
  enddo
  hth(ny,m)=one/(a0+ape*eth(ny-1,m))
enddo

!---------------------------------------------------------------------
 !Set up tridiagonal inversion arrays for divergence time stepping:
 !(Below ap = 1/dy^2, ape = 2/dy^2, kdsq = (f/c)^2, rdis = 2/dt + diss
 !and diss = C*f*(k/k_max)^(2p) is the hyperviscosity operator 
 !defined above in this subroutine)
a00=-ape-kdsq
do m=0,nxm1
  a0=a00-rkx(m)**2-csqi*rdis(m)**2
  htd(0,m)=one/a0
  etd(0,m)=-ape*htd(0,m)
  do j=1,ny-1
    htd(j,m)=one/(a0+ap*etd(j-1,m))
    etd(j,m)=-ap*htd(j,m)
  enddo
  htd(ny,m)=one/(a0+ape*etd(ny-1,m))
enddo

return
end subroutine init_spectral

!====================================================================
subroutine height(qq,gg,uum,uup,qoff,hh,zz)

! Input:  PV (qq), acceleration divergence (gg), boundary zonal
!         velocity (uum at y = ymin & uup at y = ymax) and PV
!         offset (qoff).
!         *** Note: only the input field qq is in physical space.

! Output: dimensionless height anomaly (hh) and relative vorticity
!         (zz). The PV is also corrected by a constant offset to be
!         consistent with the constraint on the global mean vorticity.
!         *** Note: only the output zz is in semi-spectral space.

! *** Important: a guess for hh must be provided (could use hh = 0).  

 !Declarations:
implicit none

 !Passed variables:
double precision:: qq(0:ny,0:nxm1),gg(0:ny,0:nxm1)
double precision:: hh(0:ny,0:nxm1),zz(0:ny,0:nxm1)
double precision:: uum(0:nxm1),uup(0:nxm1)
double precision:: qoff

 !Local variables:
double precision,parameter:: htol=1.e-10 !error tolerance for hh
double precision:: wk(0:ny,0:nxm1)
double precision:: db,qb,herr

!---------------------------------------------------------------------
! Iterate to find hh and zz (and correct qq by a constant offset):
herr=one
do while (herr .gt. htol)
   !Compute PV offset:
  wk=(one+hh)*qq
  call average(wk,qb)
  qb=qoff-qb
   !Correct PV:
  qq=qq+qb
   !Define relative vorticity (has correct domain average)
  zz=(one+hh)*qq-cof
   !Define rhs of the Helmholtz equation for hh:
  wk=cof*(zz-cof*hh)
   !Transform in x:
  call forfft(nyp1,nx,wk,xtrig,xfactors)
   !Subtract gamma (already in semi-spectral space):
  wk=wk-gg
   !Solve Helmholtz equation for hh; new estimate now in wk:
  call helmholtz(wk,uum,uup)
   !FFT dimensionless height anomaly back to physical space:
  call revfft(nyp1,nx,wk,xtrig,xfactors)
   !Compute max abs error:
  herr=maxval(abs(wk-hh))
   !Store updated height field:
  hh=wk
enddo

! Convert vorticity to semi-spectral space and de-alias:
call forfft(nyp1,nx,zz,xtrig,xfactors)
call dealias(zz,1)

return
end subroutine height

!====================================================================
subroutine velocity(dd,zz,uum,uup,uu,vv)

! Given the divergence (dd), vorticity (zz) and the boundary values
! uum & uup at y = y_min & y_max respectively in semi-spectral space,
! this routine finds the velocity field (uu,vv) on the grid.
  
! *** Important: before this routine is called, it must be ensured
! that the global mean vorticity equals the mean velocity difference
! between y = ymin and ymax, see subroutine height and use of qoff.

 !Declarations:
implicit none

 !Passed variables:
double precision:: dd(0:ny,0:nxm1),zz(0:ny,0:nxm1)
double precision:: uu(0:ny,0:nxm1),vv(0:ny,0:nxm1)
double precision:: uum(0:nxm1),uup(0:nxm1)

 !Local variables:
double precision:: psi_c(0:ny),chi_c(0:ny)
double precision:: psi_s(0:ny),chi_s(0:ny)
double precision:: zb,db
integer:: l,m,j

!---------------------------------------------------------------------
! Ensure zero global mean divergence:
db=(f12*(dd(0,0)+dd(ny,0))+sum(dd(1:nym1,0)))*dnyi
dd(:,0)=dd(:,0)-db

!---------------------------------------------------------------------
 !Solve the block tri-diagonal systems for the cosine part of psi
 !and the sine part of chi, then vice-versa:
 !(Below ap = 1/dy^2, ape = 2/dy^2 and c0 = 2/dy)
do m=1,nwm1
  l=nx-m

!---------------------------------------
!  cosine part of psi, sine part of chi:
  zb=zz(0,m)-c0*uum(m)
  db=dd(0,l) 
  psi_c(0)=zb*a0b(0,m)-db*b0b(0,m)
  chi_s(0)=db*a0b(0,m)-zb*b0b(0,m)

  do j=1,ny-1
    zb=zz(j,m)-psi_c(j-1)*ap
    db=dd(j,l)-chi_s(j-1)*ap
    psi_c(j)=zb*a0b(j,m)-db*b0b(j,m)
    chi_s(j)=db*a0b(j,m)-zb*b0b(j,m)
  enddo

  zb=zz(ny,m)+c0*uup(m)-psi_c(ny-1)*ape
  db=dd(ny,l)          -chi_s(ny-1)*ape
  psi_c(ny)=zb*a0b(ny,m)-db*b0b(ny,m)
  chi_s(ny)=db*a0b(ny,m)-zb*b0b(ny,m)

  do j=ny-1,0,-1
    psi_c(j)=et1(j,m)*psi_c(j+1)+et2(j,m)*chi_s(j+1)+psi_c(j)
    chi_s(j)=et1(j,m)*chi_s(j+1)+et2(j,m)*psi_c(j+1)+chi_s(j)
  enddo

!---------------------------------------
!  sine part of psi, cosine part of chi:
  zb=zz(0,l)-c0*uum(l)
  db=dd(0,m) 
  psi_s(0)=zb*a0b(0,m)+db*b0b(0,m)
  chi_c(0)=db*a0b(0,m)+zb*b0b(0,m)

  do j=1,ny-1
    zb=zz(j,l)-psi_s(j-1)*ap
    db=dd(j,m)-chi_c(j-1)*ap
    psi_s(j)=zb*a0b(j,m)+db*b0b(j,m)
    chi_c(j)=db*a0b(j,m)+zb*b0b(j,m)
  enddo

  zb=zz(ny,l)+c0*uup(l)-psi_s(ny-1)*ape
  db=dd(ny,m)          -chi_c(ny-1)*ape
  psi_s(ny)=zb*a0b(ny,m)+db*b0b(ny,m)
  chi_c(ny)=db*a0b(ny,m)+zb*b0b(ny,m)

  do j=ny-1,0,-1
    psi_s(j)=et1(j,m)*psi_s(j+1)-et2(j,m)*chi_c(j+1)+psi_s(j)
    chi_c(j)=et1(j,m)*chi_c(j+1)-et2(j,m)*psi_s(j+1)+chi_c(j)
  enddo

!-----------------------------------------------------------------
!  Differentiate and combine to get non-zonal velocity components:
  do j=1,ny-1
    uu(j,m)=-rkx(m)*chi_s(j)-hglyi*(psi_c(j+1)-psi_c(j-1))
    uu(j,l)= rkx(m)*chi_c(j)-hglyi*(psi_s(j+1)-psi_s(j-1))
    vv(j,m)=-rkx(m)*psi_s(j)+hglyi*(chi_c(j+1)-chi_c(j-1))
    vv(j,l)= rkx(m)*psi_c(j)+hglyi*(chi_s(j+1)-chi_s(j-1))
  enddo
enddo

! Set Nyquist wavenumber component to zero: 
uu(:,nw)=zero
vv(:,nw)=zero

! Insert known boundary values:  
uu(0,:)=uum
vv(0,:)=zero

uu(ny,:)=uup
vv(ny,:)=zero

! Add zonal flow corresponding to wavenumber m = 0:
do j=1,ny-1
  uu(j,0)=uu(j-1,0)-hgly*(zz(j,0)+zz(j-1,0))
  vv(j,0)=vv(j-1,0)+hgly*(dd(j,0)+dd(j-1,0))
enddo
! Note: the boundary values were set just above.

! FFT velocity field back to physical space:
call revfft(nyp1,nx,uu,xtrig,xfactors)
call revfft(nyp1,nx,vv,xtrig,xfactors)
! Note, uu & vv are automatically de-aliased as dd, zz, uum & uup are.

return
end subroutine velocity

!====================================================================
subroutine helmholtz(hh,uum,uup)

! Solves c^2*Lap(h) - f^2*h = S (passed as hh in semi-spectral space),
! with the boundary conditions c^2*h_y = -f*u at y = ymin and y_max.
! Note, u = uum at y = ymin and u = uup at y = ymax.  These must be
! in semi-spectral space, where the solution hh is returned.
  
 !Declarations:
implicit none

 !Passed variables:
double precision:: hh(0:ny,0:nxm1)
double precision:: uum(0:nxm1),uup(0:nxm1)

 !Local variables:
double precision:: havg
integer:: j,m

!---------------------------------------------------------------------
 !Solve the tri-diagonal system for c^2*h:
 !(Below ap = 1/dy^2, ape = 2/dy^2 and c0f = 2*f/dy)
do m=0,nxm1
  hh(0,m)=(hh(0,m)-c0f*uum(m))*hth(0,m)
  do j=1,ny-1
    hh(j,m)=(hh(j,m)-ap*hh(j-1,m))*hth(j,m)
  enddo
  hh(ny,m)=(hh(ny,m)+c0f*uup(m)-ape*hh(ny-1,m))*hth(ny,m)
  do j=ny-1,0,-1
    hh(j,m)=eth(j,m)*hh(j+1,m)+hh(j,m)
  enddo
enddo

 !Obtain h after dividing by c^2:
hh=csqi*hh

! Ensure zero global mean:
havg=(f12*(hh(0,0)+hh(ny,0))+sum(hh(1:nym1,0)))*dnyi
hh(:,0)=hh(:,0)-havg

! De-alias in x:
call dealias(hh,1)

return
end subroutine helmholtz

!====================================================================
subroutine sistep(cc,gg,suum,suup)

! Solves semi-implicit system for the sums of cc & gg over the time
! interval (t,t+dt).  Here, suum and suup are the sums of the zonal
! velocity sources on y = y_min and y_max.  On input, cc and gg
! contain 2*N_delta and 2*N_gamma, (twice) the nonlinear source
! terms in the notes (see Dritschel & Jalali, JFM 900, A33, 2020;
! equations (B12) & (B13) in particular).

 !Declarations:
implicit none

 !Passed variables:
double precision:: cc(0:ny,0:nxm1),gg(0:ny,0:nxm1)
double precision:: suum(0:nxm1),suup(0:nxm1)

 !Local variables:
double precision:: wk(0:ny,0:nxm1)
integer:: j,m

!---------------------------------------------------------------------
 !Form r.h.s. of the equation for 2*cc_bar = cc^n + cc^{n+1}:
do j=0,ny
  wk(j,:)=gg(j,:)+rdis*cc(j,:)
enddo
 !Initialise 2*gg_bar (need to add R*(2*cc_bar) after):
gg=-cc

 !Solve the tri-diagonal system for -c^2*(2*cc_bar) to be able to
 !use the tri-diagonal coefficients initialised previously:
 !(Below ap = 1/dy^2, ape = 2/dy^2 and c0f = 2*f/dy)
do m=0,nxm1
  cc(0,m)=(wk(0,m)-c0f*suum(m))*htd(0,m)
  do j=1,ny-1
    cc(j,m)=(wk(j,m)-ap*cc(j-1,m))*htd(j,m)
  enddo
  cc(ny,m)=(wk(ny,m)+c0f*suup(m)-ape*cc(ny-1,m))*htd(ny,m)
  do j=ny-1,0,-1
    cc(j,m)=etd(j,m)*cc(j+1,m)+cc(j,m)
  enddo
enddo

 !Obtain 2*cc_bar after multiplying by -1/c^2:
cc=-cc*csqi

 !Complete 2*gg_bar:
do j=0,ny
  gg(j,:)=gg(j,:)+rdis*cc(j,:)
enddo
 !On rhs, gg = -2*N_delta

return
end subroutine sistep

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
subroutine lowpass(ff,iopt)

! Low-pass filters a field ff (in both x wavenumber space and by
! repeated 1-2-1 averages in y, preserving boundary values 
! (except as modified by the x filtering).

! If iopt = 0, ff is passed in physical space
! If iopt = 1, ff is passed in semi-spectral space
  
 !Declarations:
implicit none

 !Passed variables:
double precision:: ff(0:ny,0:nxm1)
integer:: iopt

 !Local variables:
double precision:: wk(0:ny,0:nxm1)
integer:: j

!---------------------------------------------------------------------
 !Transform in x if ff is passed in physical space: 
if (iopt .eq. 0) call forfft(nyp1,nx,ff,xtrig,xfactors)

 !Apply low-pass filter in x:
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

 !Transform back if ff is expected in physical space: 
if (iopt .eq. 0) call revfft(nyp1,nx,ff,xtrig,xfactors)

return
end subroutine lowpass

end module spectral
