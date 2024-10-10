!#########################################################################
!  Computes an initially balanced flow (having delta = delta_t = 0) and
!  uniform PV equal to the background resting state, for a specified
!  scaled surface pressure distribution p_0(x,y) in p0_init.r8.
!  Note: p_0 = p'/f^2 where f is the Coriolis frequency.

!  The results are output to various files in the 3d subdirectory.

!  Updated 7/12/2023 by D G Dritschel @ St Andrews
!#########################################################################

program balance

 !Import spectral module:
use spectral

implicit none

 !Various quantities needed below:
double precision:: p0(ng,ng)
double precision:: bb(ng,ng,0:nz),pp(ng,ng,0:nz)
double precision:: ox(ng,ng,0:nz),oy(ng,ng,0:nz),oz(ng,ng,0:nz)
double precision:: ux(ng,ng,0:nz),uy(ng,ng,0:nz)
double precision:: gg(ng,ng,0:nz)
double precision:: sfac
double precision:: t

!---------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Read scaled surface pressure, p0 (= p'/f^2):
open(11,file='p0_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nhbytes)
read(11,rec=1) t,p0
close(11)

!---------------------------------------------------------------
 !Find balanced 3D flow:
call bcbal(p0,bb,ox,oy,oz,pp,ux,uy,gg)

 !Undo scaling of variables for direct use in PS3D or in EPIC:
sfac=cof*bvf
bb=sfac*bb
ox=bvf*ox
oy=bvf*oy
oz=cof*oz
sfac=cof**2
pp=sfac*pp
ux=cof*ux
uy=cof*uy

!------------------------------------------------------------------
 !Write data:
write(*,*)
write(*,*) ' Writing x vorticity component to 3d/ox.r4'
open(11,file='3d/ox.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(ox)
close(11)

write(*,*) ' Writing y vorticity component to 3d/oy.r4'
open(11,file='3d/oy.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(oy)
close(11)

write(*,*) ' Writing z vorticity component to 3d/oz.r4'
open(11,file='3d/oz.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(oz)
close(11)

write(*,*) ' Writing buoyancy anomaly to 3d/ba.r4'
open(11,file='3d/ba.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(bb)
close(11)

write(*,*) ' Writing perturbation pressure to 3d/pp.r4'
open(11,file='3d/pp.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(pp)
close(11)

write(*,*) ' Writing x velocity component to 3d/ux.r4'
open(11,file='3d/ux.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(ux)
close(11)

write(*,*) ' Writing y velocity component to 3d/uy.r4'
open(11,file='3d/uy.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(uy)
close(11)

write(*,*) ' Writing static stability Gamma to 3d/gg.r4'
open(11,file='3d/gg.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(gg)
close(11)

 !Compute Richardson number:
sfac=one/bvf**2
gg=sfac*(ox**2+oy**2)/gg
write(*,*) ' Writing Richardson number to 3d/ri.r4'
open(11,file='3d/ri.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(gg)
close(11)

write(*,*)
write(*,*) ' ==> See parameters.f90 for all parameters used.'
write(*,*)

 !End main program
end program balance
!=======================================================================
