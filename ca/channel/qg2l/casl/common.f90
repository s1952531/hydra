module common

 !Import constants, parameters and common arrays needed for inversion etc:
use constants
use variables
use contours
use spectral
use generic

 !Quantities which need to be preserved between recontouring and evolution:

 !Gridded PV fields (full & residual):
double precision:: qs(0:ny,0:nxm1,nz),qd(0:ny,0:nxm1,nz)

 !Equilibrium interior interface displacement and scaled topography f_0*h_b:
double precision:: dd1eq(0:ny,0:nxm1),fhb(0:ny,0:nxm1)

 !Upper interface displacement (if no barotropic mode is present):
double precision,allocatable,dimension(:,:):: dd2eq,dd2,dd2pre

 !Non-conservative terms (if present):
double precision,allocatable,dimension(:,:,:):: sq,sqpre
double precision,allocatable,dimension(:,:):: wforce

 !For semi-Lagrangian advection: 
double precision:: xig(0:nxm1),yig(0:ny)
double precision,parameter:: xigmax=dble(nx),yigmax=dble(ny) 

 !Initial time of simulation:
double precision:: tinit

contains

!=============================================================
subroutine readcont(iind,iniqs)

! Reads time loop = iind from the contour and residual PV data
! in the cont subdirectory and, optionally, initialises qs if
! iniqs = 1 above.

! Sets the initial time, tinit, from data read from qqsynopsis.asc

implicit none

 !Passed variables:
integer:: iind,iniqs

 !Local variables:
real:: tr4,qdr4(0:ny,0:nxm1)
integer:: iop(nm)
integer:: i,j,l,ibeg,jbeg,iend,iz
character(len=3):: pind

!-----------------------------------------------------------------------
write(pind(1:3),'(i3.3)') iind

 !Open contour files:
if (replace) then
   !Here, data from an old run (in restart-*) is used to initialise a new one:
   open(40,file='contours/restart-qqsynopsis.asc',status='old')
   do i=1,iind
      read(40,*) nq,nptq,t,qjump(1),qjump(2),qavg(1),qavg(2)
   enddo
   close(40)

   open(40,file='contours/restart-qqindex'//pind,form='unformatted', &
         access='direct',status='old',recl=20*nq)
   read(40,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq),layq(1:nq),iop(1:nq)
   close(40)

   open(40,file='contours/restart-qqnodes'//pind,form='unformatted', &
         access='direct',status='old',recl=16*nptq)
   read(40,rec=1) xq(1:nptq),yq(1:nptq)
   close(40)

   open(40,file='contours/restart-qqresi'//pind,form='unformatted', &
         access='direct',status='old',recl=nbytes)
   do iz=1,nz
      read(40,rec=iz) tr4,qdr4
      qd(:,:,iz)=dble(qdr4)
   enddo
   close(40)

else
   !Here, an old run is being continued from a specific time loop, iind:
   open(40,file='contours/qqsynopsis.asc',status='old')
   do i=1,iind
      read(40,*) nq,nptq,t,qjump(1),qjump(2),qavg(1),qavg(2)
   enddo
   close(40)

   open(40,file='contours/qqindex'//pind,form='unformatted', &
         access='direct',status='old',recl=20*nq)
   read(40,rec=1) npq(1:nq),i1q(1:nq),indq(1:nq),layq(1:nq),iop(1:nq)
   close(40)

   open(40,file='contours/qqnodes'//pind,form='unformatted', &
         access='direct',status='old',recl=16*nptq)
   read(40,rec=1) xq(1:nptq),yq(1:nptq)
   close(40)

   open(40,file='contours/qqresi'//pind,form='unformatted', &
         access='direct',status='old',recl=nbytes)
   do iz=1,nz
      read(40,rec=iz) tr4,qdr4
      qd(:,:,iz)=dble(qdr4)
   enddo
   close(40)
endif

 !Time corresponding to the data read:
tinit=t

!--------------------------------------------------------------------
 !Determine beginning and ending nodes and contours in each layer and
 !construct nextq array:
do l=1,nz
   jl2q(l)=0
enddo

do j=1,nq
   ibeg=i1q(j)
   iend=ibeg+npq(j)-1
   i2q(j)=iend
   il2q(layq(j))=i2q(j)
   jl2q(layq(j))=j
   do i=ibeg,iend-1
      nextq(i)=i+1
   enddo
   if (iop(j) .eq. 0) then
      nextq(iend)=0
   else
      nextq(iend)=ibeg
   endif
enddo

ibeg=1
jbeg=1

do l=1,nz
   if (jl2q(l) .gt. 0) then
      il1q(l)=ibeg
      ibeg=il2q(l)+1
      jl1q(l)=jbeg
      jbeg=jl2q(l)+1
   endif
enddo

if (iniqs .eq. 1) then
   !Initialise entire PV field, qs, by combining contours and residual PV:
   call con2grid(qs)
   qs=qs+qd
endif
 
return
end subroutine readcont

!=======================================================================

 !Main end module
end module common
