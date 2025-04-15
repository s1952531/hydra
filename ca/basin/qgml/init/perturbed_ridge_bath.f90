program random_bath
! |-----------------------------------------------------------|
! |      This routine sets up a single SN-ridge bathymetry    |
! |      perturbed by random data with k^{-2}-spectrum.       |
! |      Optionally, one can include a linear slope in y.     |
! |-----------------------------------------------------------|

use constants
use sta2dfft

implicit none

 !Doubled domain in x & y:
integer,parameter:: nxe=2*nx,nye=2*ny
 !Maximum x & y wavenumbers in doubled domain (which is periodic):
integer,parameter:: nwx=nxe/2,nwy=nye/2

 !Common arrays, constants:
double precision:: bb(nye,nxe)
double precision:: ss(nxe,nye)
double precision:: br(0:ny,0:nx)
double precision:: hrkx(nxe),hrky(nye),rkx(nxe),rky(nye)
double precision:: xtrig(2*nxe),ytrig(2*nye)
double precision:: fac,s,rms,amp,k2,phase,slope,width,height

integer:: xfactors(5),yfactors(5)
integer:: kx,ky,k,i,ix,iy
integer, dimension(:), allocatable :: seed

!-----------------------------------------------------------------------
 !Initialise random number generator:
call random_seed(size=k)
allocate(seed(1:k))
seed(:)=iseed
do i=1,iseed
  call random_seed(put=seed)
enddo

!----------------------------------------------------------------------
 !Initialise FFTs:
call init2dfft(nxe,nye,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define x wavenumbers:
rkx(1)=zero
do kx=1,nwx-1
  rkx(kx+1)    =hrkx(2*kx)
  rkx(nxe+1-kx)=hrkx(2*kx)
enddo
rkx(nwx+1)=hrkx(nx)

 !Define y wavenumbers:
rky(1)=zero
do ky=1,nwy-1
  rky(ky+1)    =hrky(2*ky)
  rky(nye+1-ky)=hrky(2*ky)
enddo
rky(nwy+1)=hrky(nye)

!-----------------------------------------------------
 !Compute k^{-2} random field:
write(*,*) ' Generating a pertubed ridge-type topography.'
write(*,*) ' Enter the rms value of the perturbation (determines scaling factor):'
read(*,*) rms
write(*,*) ' This field is superposed on top of a linear sloping'
write(*,*) ' bathymetry in y with zero mean height and a south-north ridge. Enter the slope:'
read(*,*) slope
write(*,*) ' Enter the width of the ridge:'
read(*,*) width
write(*,*) ' Enter the height of the ridge:'
read(*,*) height

 !Generate spectrum:
do ky=1,nwy+1
  do kx=1,nwx+1
    k2=rkx(kx)**2+rky(ky)**2
    if (k2>0.0) then
       !Set amplitude proportional to k^{-2} with random phase:
      amp=rms/sqrt(k2) ! k^{-1} in amplitude -> k^{-2} in power spectrum
      call random_number(s)
      phase=twopi*s
      ss(kx,ky)=amp*cos(phase)+(0.0,1.0)*amp*sin(phase) ! Complex spectrum
    else
      ss(kx,ky)=(0.0,0.0) ! Avoid division by zero at k=0
    endif
  enddo
enddo

! Remove mean and Nyquist frequency:
ss(1,1)=zero
ss(nwx+1,nwy+1)=zero

! Transform to physical space:
call spctop(nxe,nye,ss,bb,xfactors,yfactors,xtrig,ytrig)

! Work out rms:
fac=sqrt(sum(bb**2)/dble(nxe*nye))

! Renormalise field:
fac=rms/fac
bb=fac*bb

! Add uniform sloping part in y and Gaussian ridge centered around x=0:
do ix=0,nx
  do iy=0,ny
    br(iy,ix)=bb(iy+1,ix+1)+slope*(ycen+gly*dble(iy)) &
                           +height*exp(-((glx*dble(0.5*nx-ix)-xcen)/width)**2.0)
  enddo
enddo

! Write data to a file:
open(11,file='bath.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nhbytes)
write(11,rec=1) zero,br
close(11)

end program random_bath
