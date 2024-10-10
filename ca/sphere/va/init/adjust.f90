program adjust
! Rossby adjustment test case

implicit none
integer:: ng

! Specify resolution:
write(*,*) ' Number of latitudinal divisions?'
read(*,*) ng

call gendata(ng)

contains 

!=======================================================================

subroutine gendata(ng)

 !Declarations:
implicit double precision(a-h,o-z)
implicit integer(i-n)

double precision,parameter:: pi=3.141592653589793238462643383279502884197169399375105820974944592307816d0
double precision,parameter:: hpi=pi/2.d0,twopi=2.d0*pi,pif=pi/180.d0
double precision,parameter:: f1112=11.d0/12.d0
double precision,parameter:: f=4.d0*pi

double precision:: hh(ng,2*ng),dd(ng,2*ng),qq(ng,2*ng),wka(ng,2*ng)
double precision:: clat(ng),slat(ng)
double precision:: clon(2*ng),slon(2*ng)

! Number of longitudinal divisions:
nt=2*ng

ngridp=ng*nt
nbytes=4*(ngridp+1)

!--------------------------------------------------------------
! Define perturbation used for the depth field:
write(*,*) ' We consider the dimensionless depth anomaly field'
write(*,*) ' h_tilde = A*[1-Z^2]^2 where Z = x*dot*x0, x & x0'
write(*,*) ' are vectors, and x0 is specified by (lon0,lat0).'
write(*,*) ' Enter A (suggest 0.25):'
read(*,*) apert
write(*,*) ' Enter lon0 in degrees (suggest 50):'
read(*,*) rlonp
write(*,*) ' Enter lat0 in degrees (suggest 40):'
read(*,*) rlatp

dlon=twopi/dble(nt)
! dlon: the maximum grid length (at the equator)

dlat=pi/dble(ng)
! dlat: the grid spacing in latitude

!------------------------------------------------------------
do j=1,ng
  rlat=dlat*(dble(j)-0.5d0)-hpi
  clat(j)=cos(rlat)
  slat(j)=sin(rlat)
enddo

do i=1,nt
  rlon=dlon*dble(i-1)-pi
  slon(i)=sin(rlon)
  clon(i)=cos(rlon)
enddo

r0=cos(pif*rlatp)
x0=r0*cos(pif*rlonp)
y0=r0*sin(pif*rlonp)
z0=sin(pif*rlatp)

do i=1,nt
  do j=1,ng
    x=clat(j)*clon(i)
    y=clat(j)*slon(i)
    z=slat(j)
    dot=(x*x0+y*y0+z*z0)**2
    hh(j,i)=apert*(1.d0-dot)**2
  enddo
enddo

 !Compute 4th-order average:
rsum=f1112*(clat(1)+clat(ng))+sum(clat(2:ng-1))
 !Note, rsum = 1/sin(dl/2) - sin(dl/2)/6 exactly.
rsumi=1.d0/rsum
dsumi=rsumi/dble(nt)

do i=1,nt
  wka(:,i)=clat*hh(:,i)
enddo

vsum=0.d0
do i=1,nt
  vsum=vsum+f1112*(wka(1,i)+wka(ng,i))+sum(wka(2:ng-1,i))
enddo
average=vsum*dsumi
hh=hh-average

dd=0.d0
do i=1,nt
  qq(:,i)=f*slat/(1.d0+hh(:,i))
enddo

 !Write initial height field:
open(20,file='hh_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) 0.d0,hh
close(20)

 !Write initial divergence field:
open(20,file='dd_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) 0.d0,dd
close(20)

 !Write initial PV field:
open(20,file='qq_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) 0.d0,qq
close(20)

return 
end subroutine gendata

 !End main program
end program adjust
