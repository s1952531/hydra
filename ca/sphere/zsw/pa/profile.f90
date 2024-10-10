program profile
!  -------------------------------------------------------------------------
!  |   Computes the Lagrangian mean and time mean latitude displacement d, |
!  |   dimensionless height anomaly and zonal velocity u.  Also computes   |
!  |   the rms departure from the mean.                                    |
!  |                                                                       |
!  |   The quantities are output vs initial latitude in the file davg.asc  |
!  |   and uavg.asc.                                                       |
!  -------------------------------------------------------------------------

use constants

implicit none

double precision:: q(n),phie(0:n),phi(0:n),u(0:n)
double precision:: da(0:n),ua(0:n),dv(0:n),uv(0:n)
double precision:: h(n),ha(n),hv(n)
double precision:: t,t1,t2,tmax,afac
integer:: j,loop,navg,iread

!----------------------------------------------------------------------
 !Read in initial latitudes and store in phie:
open(11,file='init.asc',status='old')
phie(0)=-hpi !This value at the South Pole is not read in below
do j=1,n
  read(11,*) phie(j)
enddo
close(11)

!---------------------------------------------------------------
 !Read final time in data:
open(15,file='ecomp.asc',status='old',access='append')
backspace(15)
read(15,*) tmax
close(15)

write(*,'(a,f9.2)') ' The final time in the simulation is ',tmax
write(*,*) ' Time averaging is done over t1 <= t <= t2.  Enter t1 & t2:'
read(*,*) t1,t2

if (t2 .gt. tmax) t2=tmax

 !Initialise time averages:
da=zero
ha=zero
ua=zero
dv=zero
hv=zero
uv=zero
navg=0

 !Set edge values:
phi(0)=-hpi
phi(n)=hpi
u(0)=zero
u(n)=zero

 !Open input data files:
open(31,file='phi.r8',form='unformatted', access='direct', &
                  status='old',recl=nbytes)
open(32,file='h.r8',form='unformatted', access='direct', &
                  status='old',recl=nbytes)
open(33,file='u.r8',form='unformatted', access='direct', &
                  status='old',recl=nbytes)

 !Read data and process:
loop=0
do
  loop=loop+1
  read(31,rec=loop,iostat=iread) t,q
  if (iread .ne. 0) exit
  if (t-1.d-6 .lt. t2) then
    do j=1,n1
      phi(j)=two*q(j)-phi(j-1)
    enddo
    read(32,rec=loop) t,h
    read(33,rec=loop) t,q
    do j=1,n1
      u(j)=two*q(j)-u(j-1)
    enddo

    if (t+1.d-6 .gt. t1) then
       !Accumulate time averages:
      navg=navg+1
      da=da+phi-phie
      dv=dv+(phi-phie)**2
      ha=ha+h
      hv=hv+h**2
      ua=ua+u
      uv=uv+u**2
    endif
  endif
enddo
close(31)
close(33)

 !Finish time averaging:
afac=one/dble(navg)
da=da*afac
dv=sqrt(dv*afac)
ha=ha*afac
hv=sqrt(abs(hv*afac-ha**2))
ua=ua*afac
uv=sqrt(uv*afac)

 !Open and write output files:
open(21,file='davg.asc',status='replace')
open(22,file='havg.asc',status='replace')
open(23,file='uavg.asc',status='replace')
do j=0,n
  write(21,'(4(1x,f14.10))') phie(j),da(j),da(j)-dv(j),da(j)+dv(j)
  write(23,'(4(1x,f14.10))') phie(j),ua(j),ua(j)-uv(j),ua(j)+uv(j)
enddo
do j=1,n
  write(22,'(4(1x,f14.10))') f12*(phie(j-1)+phie(j)), &
                             ha(j),ha(j)-hv(j),ha(j)+hv(j)
enddo
close(21)
close(22)
close(23)

write(*,*)
write(*,*) ' The results are ready in davg.asc, havg.asc & uavg.asc.'
write(*,*)
write(*,*) ' Plot average displacement dphi, h and u by typing'
write(*,*)
write(*,*) ' python analyse.py'
write(*,*)

end program profile
