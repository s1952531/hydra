program patches

! Sets up two vortex patches in the domain

use constants 

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision:: clat(ng),slat(ng)
double precision:: clon(nt),slon(nt)
double precision:: hh(ng,nt),qq(ng,nt)
double precision:: latc1,latc2,lonc1,lonc2

!------------------------------------------------------------
 !Define latitude and longitude arrays:
do j=1,ng
  rlat=dl*(dble(j)-f12)-hpi
  clat(j)=cos(rlat)
  slat(j)=sin(rlat)
enddo

do i=1,nt
  rlon=dl*dble(i-1)-pi
  clon(i)=cos(rlon)
  slon(i)=sin(rlon)
enddo

!----------------------------------------------------------
 !Define hh (used for writing h & d) and PV:
latc1=zero
latc2=zero
lonc1=zero*cos(latc1)
lonc2=1.d0*cos(latc2)

alpha=-1.0d0
beta=0.5d0
do i=1,nt
  do j=1,ng
    rlat=dl*(dble(j)-f12)-hpi
    rlon=dl*dble(i-1)-pi
    hh(j,i)=zero
    arg1=(rlat-latc1)**2+(clat(j)*rlon-lonc1)**2
    arg2=(rlat-latc2)**2+(clat(j)*rlon-lonc2)**2

    qq(j,i)=fpole*slat(j)
    if (sqrt(arg1) .lt. beta) then 
      qq(j,i)=fpole*slat(j)+alpha
      ioff=i-ng/2
      qq(j,ioff)=fpole*slat(j)+one
    endif
  enddo
enddo

 !Write initial height field:
open(20,file='hh_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

 !Write initial divergence field:
open(20,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,hh
close(20)

 !Write initial divergence field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,qq
close(20)

end program
