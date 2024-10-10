program running
!  ----------------------------------------------------------------------
!  |   Computes the running average of d_rms, v_rms and energy.         |
!  |   Output is to the files advrms.asc and aecomp.asc                 |
!  ----------------------------------------------------------------------

use constants

implicit none

double precision:: t0,t,delt,dtmid,tmax,afac
double precision:: drms,vrms,adrms,avrms
double precision:: ekin,epot,etot,aekin,aepot,aetot
integer:: loop,step,nspa,navg

!---------------------------------------------------------------
 !Read final time in data:
open(15,file='ecomp.asc',status='old',access='append')
backspace(15)
read(15,*) tmax
close(15)

write(*,'(a,f9.2)') ' The final time in the simulation is ',tmax
write(*,*) ' A running average is done over a time window delta_t.'
write(*,*) ' Enter delta_t.'
read(*,*) delt

nspa=nint(delt/dt)
dtmid=dt*dble(nspa)/two
navg=nint(tmax/dt)/nspa
afac=one/dble(nspa+1)

!Open data files and process:
open(21,file='dvrms.asc',status='old')
open(22,file='ecomp.asc',status='old')
open(31,file='advrms.asc',status='replace')
open(32,file='aecomp.asc',status='replace')

read(21,*) t0,drms,vrms
read(22,*) t0,ekin,epot,etot
do loop=1,navg
  adrms=drms
  avrms=vrms
  aekin=ekin
  aepot=epot
  aetot=etot
  do step=1,nspa
    read(21,*) t,drms,vrms
    read(22,*) t,ekin,epot,etot
    adrms=adrms+drms
    avrms=avrms+vrms
    aekin=aekin+ekin
    aepot=aepot+epot
    aetot=aetot+etot
  enddo
  adrms=adrms*afac
  avrms=avrms*afac
  aekin=aekin*afac
  aepot=aepot*afac
  aetot=aetot*afac
  write(31,'(1x,f12.5,2(1x,f12.9))') t0+dtmid,adrms,avrms
  write(32,'(f12.5,3(1x,f16.11))') t0+dtmid,aekin,aepot,aetot
  t0=t 
enddo

write(*,*)
write(*,*) ' The results are ready in advrms.asc and aecomp.asc'

end program running
