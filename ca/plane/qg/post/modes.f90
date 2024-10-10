program modes
!  ---------------------------------------------------------------------------
!  |   Computes C_m + i*S_m = |q|_max^{-1}*int_{xy} q(x,y)*(x+i*y)^m dx dy   |
!  |   over the inner half of the domain, for m = m_1 to m_2 (see below),    |
!  |   from the gridded PV in qq.r4.                                         |
!  |                                                                         |
!  |   The output is t vs {C_m^2+S_m^2}^{1/2m} for m = m_1 to m_2, for all   |
!  |   times available in qq.r4.                                             |
!  |                                                                         |
!  |   Written 30 April 2024 by D G Dritschel @ St Andrews                   |
!  ---------------------------------------------------------------------------

 !Import constants and parameters:
use constants

implicit none

 !Declarations:
integer,parameter:: ix1=nx/4+1,ix2=3*nx/4
integer,parameter:: iy1=ny/4+1,iy2=3*ny/4
integer,parameter:: mmax=10

integer:: loop,iread,m,m1,m2
integer:: ix,iy

complex(kind=8):: z(iy1:iy2,ix1:ix2),zm
double precision:: xg(ix1:ix2),yg(iy1:iy2)
double precision:: cc(mmax),ss(mmax)
double precision:: qq(ny,nx),t,fac,pow
real:: qqr4(ny,nx),tr4

!---------------------------------------------------------------
 !Open file containing PV anomaly field:
open(31,file='qq.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

write(*,*) ' Modes m = m_1 to m_2 are computed.  Enter m_1 & m_2:'
read(*,*) m1,m2

do ix=ix1,ix2
   xg(ix)=xmin+glx*dble(ix-1)
enddo

do iy=iy1,iy2
   yg(iy)=ymin+gly*dble(iy-1)
enddo

do ix=ix1,ix2
   do iy=iy1,iy2
      z(iy,ix)=cmplx(xg(ix),yg(iy))
   enddo
enddo

open(22,file='modes.asc',status='replace')

!---------------------------------------------------------------
 !Read data and process:
loop=0
do
   loop=loop+1
   iread=0
   read(31,rec=loop,iostat=iread) tr4,qqr4
   if (iread .ne. 0) exit 
   write(*,'(a,f12.5)') ' Processing t = ',tr4

   !Convert PV to double precision as qq:
   qq=dble(qqr4)
   t=dble(tr4)

   fac=one/maxval(abs(qq))
   
   !Initialise
   cc(m1:m2)=zero
   ss(m1:m2)=zero

   !Compute integrals:
   do m=m1,m2
      do ix=ix1,ix2
         do iy=iy1,iy2
            zm=z(iy,ix)**m
            cc(m)=cc(m)+qq(iy,ix)*dreal(zm)
            ss(m)=ss(m)+qq(iy,ix)*dimag(zm)
         enddo
      enddo
      cc(m)=((fac*cc(m))**2+(fac*ss(m))**2)**(one/(two*dble(m)))
   enddo

   !Write data:
   write(22,'(1x,f13.5,8(1x,e14.7))') tr4,(cc(m),m=m1,m2)

enddo

 !Close output files:
close(22)
close(31)

write(*,*)
write(*,*) ' The results are ready in modes.asc'
write(*,*)

 !End main program
end program modes
!=======================================================================
