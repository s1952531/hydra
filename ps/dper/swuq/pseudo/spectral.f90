module spectral

use constants
use stafft
use sta2dfft

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Common arrays, constants:

 !FFT commons:
integer:: xfactors(5),yfactors(5)
double precision:: xtrig(2*nx),ytrig(2*ny)
double precision:: rkx(nx),rky(ny)

 !Green, Helmholtz inversion arrays:
double precision:: green(nx,ny),helmi(nx,ny),dafilt(nx,ny)
double precision:: omsq(nx,ny)

 !Spectrum calculation arrays:
double precision:: spmf(0:max(nx,ny)),alk(max(nx,ny)),ksqi(max(nx,ny))
integer:: kmag(0:nx-1,0:ny-1),kmax,kmaxred

contains

!===========================
subroutine init_spectral

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: akx(nx),aky(ny)
double precision:: dfilx(nx),dfily(ny)
!------------------------------------------------------------------------
 !Initialise the FFT module:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,rkx,rky)

 !Define x & y wavenumbers for use in defining Green function
 !etc below.
nwx=nx/2
akx(1)=zero
do kx=2,nx-nwx
  ktmp=rkx(2*(kx-1))
  kc=nx+2-kx
  akx(kx)=ktmp
  akx(kc)=ktmp
enddo
if (mod(nx,2) .eq. 0) akx(nwx+1)=rkx(nx)

nwy=ny/2
aky(1)=zero
do ky=2,ny-nwy
  ktmp=rky(2*(ky-1))
  kc=ny+2-ky
  aky(ky)=ktmp
  aky(kc)=ktmp
enddo
if (mod(ny,2) .eq. 0) aky(nwy+1)=rky(ny)

 !Define Green function, and Helmholtz operator:
anu=drate*omega/(rkx(nx)*rky(ny))**hyppow
do ky=1,ny
  do kx=1,nx
    rksq=akx(kx)**2+aky(ky)**2
    omsq( kx,ky)=csq*rksq+fsq
    green(kx,ky)=-one/(rksq+1.d-32)
    helmi(kx,ky)= one/(omsq(kx,ky)+alpsq+alp*anu*rksq**hyppow)
  enddo
enddo
green(1,1)=zero

 !Define de-aliasing filter:
rkyf=int(dble(ny)/three)
do ky=1,ny  
  if (aky(ky) .lt. rkyf) then
    dfily(ky)=one
  else
    dfily(ky)=zero
  endif
enddo
rkxf=int(dble(nx)/three)
do kx=1,nx  
  if (akx(kx) .lt. rkxf) then
    dfilx(kx)=one
  else
    dfilx(kx)=zero
  endif
enddo

do ky=1,ny
  do kx=1,nx
    dafilt(kx,ky)=dfilx(kx)*dfily(ky)
  enddo
enddo
 !Ensure global means are zero:
dafilt(1,1)=zero

 !Initialise arrays for computing the spectrum of any field:
scx=twopi/ellx
scy=twopi/elly
rkxmax=scx*dble(nwx)
rkymax=scy*dble(nwy)
delk=sqrt(scx*scy)
delki=one/delk
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do k=1,kmax
  spmf(k)=zero
enddo
do ky=0,ny-1
  do kx=0,nx-1
    k=nint(sqrt(akx(kx+1)**2+aky(ky+1)**2)*delki)
    kmag(kx,ky)=k
    spmf(k)=spmf(k)+one
  enddo
enddo
 !Compute spectrum multiplication factor (spmf) to account for unevenly
 !sampled shells and normalise spectra by 8/(nx*ny) so that the sum
 !of the spectrum is equal to the L2 norm of the original field:
snorm=four*pi/dble(nx*ny)
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(delk*dble(k))
  ksqi(k)=one/(delk*dble(k))**2
enddo
 !Only output shells which are fully occupied (k <= kmaxred):
kmaxred=nint(sqrt((rkxmax**2+rkymax**2)/two)*delki)
!------------------------------------------------------------------------

return 
end subroutine

!==========================================================================

subroutine invert(hh,dd,uu,vv)
! invert returns the velocity field (uu,vv) in physical space given vorticity
! and divergence (cof*hh,dd) in spectral space.

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

 !Arguments declarations:
double precision:: hh(nx,ny),dd(nx,ny)
double precision:: uu(ny,nx),vv(ny,nx)
 !Work arrays:
double precision:: wka(nx,ny),wkb(nx,ny),ax(nx,ny),ay(nx,ny),bx(nx,ny),by(nx,ny)

!------------------------------------------------------------------------

 !Obtain streanfunction and velocity potential:
do ky=1,ny
  do kx=1,nx
    wka(kx,ky)=green(kx,ky)*cof*hh(kx,ky)    
    wkb(kx,ky)=green(kx,ky)*dd(kx,ky)    
  enddo
enddo

 !Take derivatives and form uu & vv:
call xderiv(nx,ny,rkx,wka,ax)
call yderiv(nx,ny,rky,wka,ay)
call xderiv(nx,ny,rkx,wkb,bx)
call yderiv(nx,ny,rky,wkb,by)

do ky=1,ny
  do kx=1,nx
    bx(kx,ky)=bx(kx,ky)-ay(kx,ky)
    by(kx,ky)=by(kx,ky)+ax(kx,ky)    
  enddo
enddo

call spctop(nx,ny,bx,uu,xfactors,yfactors,xtrig,ytrig)
call spctop(nx,ny,by,vv,xfactors,yfactors,xtrig,ytrig)

!------------------------------------------------------------------------
return 
end subroutine

!===================================================================

subroutine spec1d(var,spec)
! Computes the 1d spectrum of a field rvar which is
! periodic in x & y, and returns the result in spec.

implicit double precision (a-h,o-z)
implicit integer(i-n)

 !Passed arrays:
double precision:: var(0:nx-1,0:ny-1),spec(0:max(nx,ny))
!---------------------------------------------------------------------

do k=0,kmax
  spec(k)=zero
enddo

 !x and y-independent mode:
k=kmag(0,0)
spec(k)=spec(k)+f14*var(0,0)**2

 !y-independent mode:
do kx=1,nx-1
  k=kmag(kx,0)
  spec(k)=spec(k)+f12*var(kx,0)**2
enddo

 !x-independent mode:
do ky=1,ny-1
  k=kmag(0,ky)
  spec(k)=spec(k)+f12*var(0,ky)**2
enddo

 !All other modes:
do ky=1,ny-1
  do kx=1,nx-1
    k=kmag(kx,ky)
    spec(k)=spec(k)+var(kx,ky)**2
  enddo
enddo

return
end subroutine

!==========================================================================

end module     
