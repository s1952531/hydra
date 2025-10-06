program extend
!  ---------------------------------------------------
!  | Extends data in qq_init.r8 to double resolution |
!  ---------------------------------------------------

 !Import constants, parameters and arrays needed for FFTs:
use constants
use sta2dfft

implicit none
double precision:: qq(ny,nx),dqq(2*ny,2*nx)
double precision:: qs(nx,ny),dqs(2*nx,2*ny)
double precision:: xtrig(2*nx),ytrig(2*ny)
double precision:: dxtrig(4*nx),dytrig(4*ny)
double precision:: hrkx(nx),hrky(ny)
double precision:: dhrkx(nx),dhrky(ny)
double precision:: t
integer:: xfactors(5),yfactors(5)
integer:: dxfactors(5),dyfactors(5)
integer:: kx,ky

!---------------------------------------------------------
open(11,file='qq_init.r8',form='unformatted', &
      access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,qq
close(11)

 !Initialise FFTs:
 !Set up FFTs:
call init2dfft(nx,ny,ellx,elly,xfactors,yfactors,xtrig,ytrig,hrkx,hrky)
call init2dfft(2*nx,2*ny,ellx,elly,dxfactors,dyfactors,dxtrig,dytrig,dhrkx,dhrky)

call ptospc(nx,ny,qq,qs,xfactors,yfactors,xtrig,ytrig)

dqs=zero
do ky=1,ny/2+1
   do kx=1,nx/2+1
      dqs(kx,ky)=two*qs(kx,ky)
   enddo
enddo

do ky=1,ny/2+1
   do kx=nx,nx/2+2,-1
      dqs(nx+kx,ky)=two*qs(kx,ky)
   enddo
enddo

do ky=ny,ny/2+2,-1
   do kx=1,nx/2+1
      dqs(kx,ny+ky)=two*qs(kx,ky)
   enddo
enddo

do ky=ny,ny/2+2,-1
   do kx=nx,nx/2+2,-1
      dqs(nx+kx,ny+ky)=two*qs(kx,ky)
   enddo
enddo

call spctop(2*nx,2*ny,dqs,dqq,dxfactors,dyfactors,dxtrig,dytrig)

! Write data:
open(11,file='dqq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=32*ngridp+8)
write(11,rec=1) t,dqq
close(11)

 !End main program
end program extend
!=======================================================================
