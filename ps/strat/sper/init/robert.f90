program robert

use constants

! This routine sets up initial buoyancy and vorticity fields for the
! Robert test case (see EPIC-2D paper).

implicit none

double precision,parameter:: T0=303.15d0 !Reference temperature (degrees K)
double precision,parameter:: g=9.80665 !Acceleration due to gravity (m/s^2)
double precision:: zz(0:ny,0:nxm1),bb(0:ny,0:nxm1)
double precision:: xw,yw,rw,sw,aw,bw
double precision:: xc,yc,rc,sc,ac,bc
double precision:: x,y,r
integer:: ix,iy

!---------------------------------------------------------------------
write(*,'(a,2(f6.2),a)') &
     ' We consider a domain of height and length: ',elly,ellx,'m.'
write(*,'(a,2(i6))') ' The y and x grid resolution is: ',ny,nx

write(*,*)
write(*,*) ' -----------------------------------------------------'
write(*,*) ' Enter the centre (x_w,y_w) of the warm bubble (in m):'
read(*,*) xw,yw

write(*,*) ' Enter the inner radius R_w (in m):'
read(*,*) rw

write(*,*) ' Enter the decay length sigma_w (in m):'
read(*,*) sw

write(*,*) ' Enter the temperature perturbation A (in degrees K):'
read(*,*) aw

bw=g*aw/T0

write(*,*)
write(*,*) ' -----------------------------------------------------'
write(*,*) ' Enter the centre (x_c,y_c) of the cold bubble (in m):'
read(*,*) xc,yc

write(*,*) ' Enter the inner radius R_c (in m):'
read(*,*) rc

write(*,*) ' Enter the decay length sigma_c (in m):'
read(*,*) sc

write(*,*) ' Enter the temperature perturbation A (in degrees K):'
read(*,*) ac

bc=g*ac/T0

!---------------------------------------------------------------------
 !Set up buoyancy distribution:
do ix=0,nxm1
  x=xmin+glx*dble(ix)
  do iy=0,ny
    y=ymin+gly*dble(iy)
    r=sqrt((x-xw)**2+(y-yw)**2)
    if (r <= rw) then
      bb(iy,ix)=bw
    else
      bb(iy,ix)=bw*exp(-((r-rw)/sw)**2)
    endif
    r=sqrt((x-xc)**2+(y-yc)**2)
    if (r <= rc) then
      bb(iy,ix)=bb(iy,ix)+bc
    else
      bb(iy,ix)=bb(iy,ix)+bc*exp(-((r-rc)/sc)**2)
    endif
    zz(iy,ix)=zero
  enddo
enddo

 !Write vorticity distribution to file:
open(20,file='zz_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,zz
close(20)

 !Write buoyancy distribution to file:
open(20,file='bb_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,bb
close(20)

 !Write input data to file:
open(12,file='input_for_robert',status='unknown')
write(12,'(2(1x,f5.1),a)') xw,yw,' ! centre of warm bubble (m)'
write(12,'(1x,f5.1,a)') rw,' ! inner radius of warm bubble (m)'
write(12,'(1x,f5.1,a)') sw,' ! decay length "   "     "     " '
write(12,'(1x,f5.2,a)') aw,' ! temperature anomaly (K)'
write(12,'(2(1x,f5.1),a)') xc,yc,' ! centre of cold bubble (m)'
write(12,'(1x,f5.1,a)') rc,' ! inner radius of cold bubble (m)'
write(12,'(1x,f5.1,a)') sc,' ! decay length "   "     "     " '
write(12,'(1x,f5.2,a)') ac,' ! temperature anomaly (K)'
close(12)
      
end program robert
