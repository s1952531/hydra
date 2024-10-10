program timestep

use constants

! This routine reads evolution/monitor.asc and computes the time-step
! estimates based on CFL, buoyancy frequency and velocity strain.
! Writes timestep.asc in the evolution subdirectory.

implicit none

 !For defining the max strain & buoyancy frequency based time step:
double precision,parameter:: alpha=0.1d0

 !For controlling numerical stability (CFL_max <= 0.8 recommended):
double precision,parameter:: cflmax=0.8d0
double precision,parameter:: cflpf=cflmax*glmin

double precision:: t,uumax,bfmax,ggmax
double precision::   dtcfl,dtbfm,dtggm

integer:: iread

!---------------------------------------------------------------------
 !Open input files:
open(31,file='evolution/zz.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)
open(32,file='evolution/bb.r4',form='unformatted',access='direct', &
                 & status='old',recl=nbytes)

open(22,file='evolution/monitor.asc',status='old')

open(55,file='evolution/timestep.asc',status='replace')

!----------------------------------------------------------------------
 !Read data and process:
do
  iread=0
  read(22,*,iostat=iread) t,uumax,bfmax,ggmax
  if (iread .ne. 0) exit

  dtcfl=cflpf/(uumax+small)
  dtbfm=alpha/(bfmax+small)
  dtggm=alpha/(ggmax+small)
  
   !Save diagnostics:
  write(55,'(1x,f13.6,3(1x,1p,e14.7))') t,dtcfl,dtbfm,dtggm
enddo

close(31)
close(32)
close(55)

write(*,*) ' t vs dt_cfl, dt_bfm and dt_ggm is available in evolution/timestep.asc'

end program timestep
