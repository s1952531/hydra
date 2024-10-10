program upowerspec
!  -------------------------------------------------------------------------
!  |  Computes the 1D power spectrum for the tracer field from data        | 
!  |  in the fine grid subdirectory, previously created using genfg.f90    |
!  |                                                                       |
!  |  Output is to the formatted "spectra<nnn>", where <nnn> = 000, 001,   |
!  |  etc is the period to be analysed.                                    |
!  |                                                                       |
!  -------------------------------------------------------------------------

! Use constants and parameters from local job directory:
use constants
use sta2dfft

implicit none

 !Maximum x & y wavenumbers:
integer,parameter:: nwx=nxu/2,nwy=nyu/2

 !Common arrays, constants:
double precision:: rksq(nxu,nyu)
double precision:: rkx(nxu),hrkx(nxu)
double precision:: rky(nyu),hrky(nyu)

double precision:: xtrig(2*nxu),ytrig(2*nyu)
integer:: xfactors(5),yfactors(5)

double precision:: spmf(0:max(nxu,nyu)),alk(max(nxu,nyu))
integer:: kmag(nxu,nyu),kmax,kmaxred

integer:: period,k
real:: tr4,qqr4(nyu,nxu)
double precision:: qq(nyu,nxu),qqs(nyu,nxu),qspec(0:max(nxu,nyu))
character(len=3):: pind


!====================================================================!

call init_spectral

!---------------------------------------------------------------
 !Select time frame:
write(*,*) ' Time period (choose 0 for the initial one)?'
read(*,*) period
write(pind(1:3),'(i3.3)') period

 !Open input data file:
open(44,file='fine/qq'//pind//'.r4',form='unformatted',status='old')

 !Read data:
read(44) tr4,qqr4
close(44)
write(*,'(a,f12.5)') ' t = ',tr4

 !Convert to double precision:
qq=dble(qqr4)

call ptospc(nxu,nyu,qq,qqs,xfactors,yfactors,xtrig,ytrig)

call spec1d(qqs,qspec,0)

 !Open and write output file:
open(20,file='fine/spectra'//pind,status='unknown')
do k=1,kmaxred
  write(20,'(2(1x,f12.8))') log10(dble(k)),log10(qspec(k+1)+1.d-32)
enddo
close(20)

write(*,*)
write(*,*) ' The fine-grid vorticity spectra are in fine/spectra'//pind

contains

!===========================
subroutine init_spectral

implicit none

 !Local variables:
double precision:: rkxf(nxu),rkyf(nyu)
double precision:: scx,rkxmax,scy,rkymax
double precision:: delk,delki,fac,rksqu,snorm
integer:: kx,ky,k,iy

!----------------------------------------------------------------------
 !Set up FFTs:
call init2dfft(nxu,nyu,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)

 !Define x wavenumbers and filtered x wavenumbers:
rkx(1)=zero
do kx=1,nwx-1
  rkx(kx+1)    =hrkx(2*kx)
  rkx(nxu+1-kx)=hrkx(2*kx)
enddo
rkx(nwx+1)=hrkx(nxu)

scx=twopi/ellx
rkxmax=scx*dble(nwx)

do kx=1,nxu
  rkxf(kx)=rkx(kx)*exp(-36.d0*(rkx(kx)/rkxmax)**36.d0)
  hrkx(kx)=hrkx(kx)*exp(-36.d0*(hrkx(kx)/rkxmax)**36.d0)
enddo

 !Define y wavenumbers and filtered y wavenumbers:
rky(1)=zero
do ky=1,nwy-1
  rky(ky+1)    =hrky(2*ky)
  rky(nyu+1-ky)=hrky(2*ky)
enddo
rky(nwy+1)=hrky(nyu)

scy=twopi/elly
rkymax=scy*dble(nwy)

!-----------------------------------------------------------------------
 !Define squared total wavenumber
do ky=1,nyu
  do kx=1,nxu
    rksq(kx,ky)=rkx(kx)**2+rky(ky)**2
  enddo
enddo

!-----------------------------------------------------------------------
 !Initialise arrays for computing the spectrum of any field:
delk=sqrt((scx**2+scy**2)*f12)
delki=one/delk
kmax=nint(sqrt(rkxmax**2+rkymax**2)*delki)
do k=0,kmax
  spmf(k)=zero
enddo
do ky=1,nyu
  do kx=1,nxu
    k=nint(sqrt(rksq(kx,ky))*delki)
    kmag(kx,ky)=k
    spmf(k)=spmf(k)+one
  enddo
enddo
 !Compute spectrum multiplication factor (spmf) to account for unevenly
 !sampled shells and normalise spectra by 8/(nx*ny) so that the sum
 !of the spectrum is equal to the L2 norm of the original field:
snorm=four*pi/dble(nxu*nyu)   !!!!!!Check this with dgd....
spmf(0)=zero
do k=1,kmax
  spmf(k)=snorm*dble(k)/spmf(k)
  alk(k)=log10(delk*dble(k))
enddo
 !Only output shells which are fully occupied (k <= kmaxred):
kmaxred=nint(sqrt((rkxmax**2+rkymax**2)/two)*delki)

return 
end subroutine

!===================================================================

subroutine spec1d(ss,spec,iopt)
! Computes the 1d spectrum of a spectral field ss and returns the
! result in spec.
! If iopt = 1, the spectral field is multiplied first by K^2
!              and is returned modified!

implicit none

 !Passed variables:
double precision:: ss(nxu,nyu),spec(0:max(nxu,nyu))
integer:: iopt

 !Local variables:
integer:: kx,ky,k

!--------------------------------------------------------
if (iopt .eq. 1) then
   !Multiply spectral field by K^2 (filtered, see init_spectral):
  do ky=1,nyu
    do kx=1,nxu
      ss(kx,ky)=ss(kx,ky)*rksq(kx,ky)
    enddo
  enddo
endif

do k=0,kmax
  spec(k)=zero
enddo

 !x and y-independent mode:
k=kmag(1,1)
spec(k)=spec(k)+f14*ss(1,1)**2

 !y-independent mode:
do kx=2,nxu
  k=kmag(kx,1)
  spec(k)=spec(k)+f12*ss(kx,1)**2
enddo

 !x-independent mode:
do ky=2,nyu
  k=kmag(1,ky)
  spec(k)=spec(k)+f12*ss(1,ky)**2
enddo

 !All other modes:
do ky=2,nyu
  do kx=2,nxu
    k=kmag(kx,ky)
    spec(k)=spec(k)+ss(kx,ky)**2
  enddo
enddo

return
end subroutine

!===================================================================

end program
