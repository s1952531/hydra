program select
! Selects one time from bb.r4 and writes new_b0_init.r8 (double precision)

use constants

implicit none

real:: qqr4(ny,nx),tr4
double precision:: qq(ny,nx)
double precision:: tp,delt,t
integer:: loop

write(*,*) ' Enter time to convert:'
read(*,*) tp

 !Find appropriate record number to read for time chosen:
delt=0.1d0*tgsave
if (tp < delt) then
   loop=1
else
   !Read ene-ens.asc to find record number:
   open(21,file='ene-ens.asc',status='old')
   t=zero
   loop=0
   do while (t+delt < tp)
      read(21,*) t
      loop=loop+1
   enddo
   close(21)
endif

 !Open file containing the scaled surface buoyancy field:
open(31,file='bb.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
read(31,rec=loop) tr4,qqr4
close(31)
write(*,'(a,f5.1)') ' Writing data at t = ',t
qq=dble(qqr4)

open(11,file='new_qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(11,rec=1) zero,qq
close(11)

end program select
