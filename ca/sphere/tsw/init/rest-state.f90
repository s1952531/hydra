program rest

! Sets up a rest state consisting of PV varying as
!                2*Omega*z

use constants 

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: slat(ng),hh(ng,nt),qq(ng,nt)

!------------------------------------------------------------
 !Define sin(latitude):
do j=1,ng
  rlat=dl*(dble(j)-f12)-hpi
  slat(j)=sin(rlat)
enddo

 !Define PV:
do i=1,nt
  do j=1,ng
    qq(j,i)=fpole*slat(j)
  enddo
enddo

 !Write initial height field (h = 1):
hh=one
open(20,file='hh_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

 !Write zero initial divergence field:
hh=zero
open(20,file='dd_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
      access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

if (thermal) then
   !Write zero thermal equilibrium height field:
  open(20,file='hequil.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(20,rec=1) zero,hh
  close(20)

   !Write zero log potential temperature field:
  open(20,file='tt.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(20,rec=1) zero,hh
  close(20)
endif

end program
