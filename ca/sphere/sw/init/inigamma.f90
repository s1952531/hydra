program inigamma

!=====================================================
! Initialises the acceleration divergence given the
! PV, the height and the divergence fields.  
!=====================================================

 !Import spectral module:
use spectral

 !Declarations:
implicit none

double precision:: zz(nLatGridPts,nLongGridPts),dd(nLatGridPts,nLongGridPts),hh(nLatGridPts,nLongGridPts)
double precision:: zs(nLatGridPts,nLongGridPts),ds(nLatGridPts,nLongGridPts),hs(nLatGridPts,nLongGridPts)
double precision:: ud(nLatGridPts,nLongGridPts),uu(nLatGridPts,nLongGridPts),gg(nLatGridPts,nLongGridPts)
double precision:: wka(nLatGridPts,nLongGridPts),wkb(nLatGridPts,nLongGridPts)
double precision:: t
integer:: i,m

!------------------------------------------------------------
 !Initialise spectral module:
call init_spectral

!----------------------------------------------------------------------
 !Read in the dimensionless height anomaly, hh:
open(11,file='hh_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,hh
close(11)

 !Read in the velocity divergence, dd:
open(11,file='dd_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,dd
close(11)

 !Read in gridded PV, qs = (zeta+f)/(1+h), where zeta is the relative
 !vorticity and h is the dimensionless height anomaly:
open(11,file='qq_init.r8',form='unformatted', &
    & access='direct',status='old',recl=2*nbytes)
read(11,rec=1) t,zz
close(11)

 !Define relative vorticity:
do i=1,nLongGridPts
  zz(:,i)=(one+hh(:,i))*zz(:,i)-corFreq
enddo

 !Create a spectral versions of hh, zz and dd:
zs=zz
call forfft(nLatGridPts,nLongGridPts,zs,trig,factors) 
ds=dd
call forfft(nLatGridPts,nLongGridPts,ds,trig,factors) 
hs=hh
call forfft(nLatGridPts,nLongGridPts,hs,trig,factors) 

 !Compute Laplacian of h:
call laplace(hs,gg) ! gg is returned in semi-spectral space
 !Bring gg back to physical space:
call revfft(nLatGridPts,nLongGridPts,gg,trig,factors) 

 !Invert Laplace's operator on the divergence (ds):
call laplinv(ds,wka,wkb)
 !Here the divergence potential is wka while d(wka)/dlat = wkb,
 !the divergent meridional velocity (all in semi-spectral space)

 !Compute divergent zonal velocity and store in ud:
call deriv(nLatGridPts,nLongGridPts,rk,wka,ud) 
do m=1,nLongGridPts
  ud(:,m)=clati*ud(:,m)
enddo
 !ud = (1/r)*d(wka)/dlon where r = cos(lat)

 !Invert Laplace's operator on the relative vorticity (zs):
call laplinv(zs,wka,wkb)
 !Here the streamfunction is wka while d(wka)/dlat = wkb.

 !Complete calculation of zonal velocity, uu:
uu=ud-wkb
 !Bring uu back to physical space:
call revfft(nLatGridPts,nLongGridPts,uu,trig,factors) 

 !Add everything up to define acceleration divergence, gg:
do i=1,nLongGridPts
  gg(:,i)=corFreq*zz(:,i)-bet*uu(:,i)-csq*gg(:,i)
enddo

 !Write data:
open(20,file='gg_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) 0.d0,gg
close(20)

write(*,*) ' Acceleration divergence initialised in gg_init.r8'
end program inigamma
