!#########################################################################
!  Computes psi, psi_x and psi_y using subrouting getpsi in spectral.f90
!#########################################################################

program psi

 !Import spectral module:
use spectral

implicit none

 !Various quantities needed below:
double precision:: b0(ng,ng)
double precision:: ss(ng,ng,0:nz),sx(ng,ng,0:nz),sy(ng,ng,0:nz)
double precision:: t

!---------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Read dimensionless surface bouyancy, b0:
open(11,file='b0_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nhbytes)
read(11,rec=1) t,b0
close(11)

!---------------------------------------------------------------
 !Find psi, psi_x and psi_y:
call getpsi(b0,ss,sx,sy)

!------------------------------------------------------------------
 !Write data:
write(*,*)
write(*,*) ' Writing psi to 3d/ss.r4'
open(11,file='3d/ss.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(ss)
close(11)

write(*,*)
write(*,*) ' Writing psi_x to 3d/sx.r4'
open(11,file='3d/sx.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(sx)
close(11)

write(*,*)
write(*,*) ' Writing psi_y to 3d/sy.r4'
open(11,file='3d/sy.r4',form='unformatted',access='direct', &
     status='replace',recl=ntbytes)
write(11,rec=1) real(zero),real(sy)
close(11)

 !End main program
end program psi
!=======================================================================
