program solitary
! Initialises a solitary wave (independent of y).

! Only run with zero Coriolis frequency!  This is a non-rotating solution.

! There is no need to run balinit after as this produces the initial
! fields of q_l, delta and gamma_l.

implicit double precision(a-h,o-z)

write(*,*) ' Enter the grid resolution ng:'
read(*,*) ng

call gendata(ng)

 !Internal subroutine definitions (inherit global variables):
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

contains

!=======================================================================

subroutine gendata(ng)

implicit double precision(a-h,o-z)

double precision,parameter:: pi=3.14159265358979323846264d0

double precision:: qq(ng,ng),dd(ng,ng),gg(ng,ng),hh(ng,ng)

write(*,*) ' Enter the amplitude a of the solitary wave:'
read(*,*) a
write(*,*) ' Enter the undisturbed depth h_0:'
read(*,*) h0
! Compute inverse width b:
b=0.5d0*sqrt(3.d0*a/(h0**2*(h0+a)))
write(*,'(a,f9.6)') '  This gives an inverse width b = ',b
write(*,'(a,f9.6)') '  Note: sech^2(b*pi) = ',1.d0/cosh(b*pi)**2

write(*,*) ' ng = ',ng

! Define displacement height of wave (store in hh temporarily):
havg=zero
gl=2.d0*pi/dble(ng)
do ix=1,ng
  x=gl*dble(ix-1)-pi
  eta=a/cosh(b*x)**2
  havg=havg+eta
  do iy=1,ng
    hh(iy,ix)=eta
  enddo
enddo
havg=havg/dble(ng)
hbar=h0+havg
write(*,'(a,f8.6)') '  Mean fluid depth H = ',hbar

! Add h0 to the displacement to define the total wave height:
hh=hh+h0

write(*,*) ' Enter sqrt{gH}:'
read(*,*) cgw

! Compute gravity:
g=cgw**2/hbar

! Compute the wave speed c:
c=sqrt(g*(h0+a))
write(*,'(a,f9.6)') '  Solitary wave speed c = ',c

period=2.d0*pi/c
write(*,'(a,f9.6)') '  Transit time of domain T = 2*pi/c = ',period

! Ideal time step:
dtp=1.d0/dble(ng)
! Number of periods between data saves:
nsaves=20
tsave=period/dble(nsaves)
ndt=nint(tsave/dtp)
dt=tsave/dble(ndt)
write(*,'(a,f9.6)') '  Suggested time step to use dt = ',dt

ubar=c*havg/hbar
write(*,'(a,f8.6)') '  Mean fluid speed U = ',ubar

! Define the divergence u_x and linearised acceleration divergence:
tb=2.d0*b
gfac=-g*tb*b
dfac=-tb*c*h0
ainv=1.d0/a
do ix=1,ng
  x=gl*dble(ix-1)-pi
  eta=a/cosh(b*x)**2
  thbx=tanh(b*x)
  delta=dfac*eta*thbx/(h0+eta)**2
  gammal=gfac*eta*(2.d0*thbx**2-ainv*eta)
  do iy=1,ng
    dd(iy,ix)=delta
    gg(iy,ix)=gammal
    qq(iy,ix)=0.d0
  enddo
enddo

!-------------------------------------------------------------
! Write data:
nbytes=4*(ng*ng+1)

!Linearised PV anomaly (=0):
open(11,file='ql_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) 0.d0,qq
close(11)

!Divergence:
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) 0.d0,dd
close(11)

!Acceleration divergence:
open(11,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) 0.d0,gg
close(11)

! Write a file containing parameters to insert in the parameters.f90 file:
open(20,file='params.asc',status='replace')
write(20,*) hbar,ubar,cgw
write(20,*) dt,tsave,period
close(20)

open(88,file='ecross.asc',status='replace')
do ix=1,ng
  x=gl*dble(ix-1)-pi
  eta=a/cosh(b*x)**2
  write(88,'(3(1x,f15.10))') x,h0+eta,c*(1.d0-h0/(h0+eta))
enddo

return
end subroutine gendata

 !End main program
end program solitary

!=======================================================================
