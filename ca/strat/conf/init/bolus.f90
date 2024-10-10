program bolus

use constants

! This routine sets up initial buoyancy and (zero) vorticity fields
! corresponding to Magda Carr's experiment involving up upper linearly
! stratified layer above a homogeneous layer, separated initially on the
! left by a homogeneous region extending through part of the depth.

! Written by D. G. Dritschel on 23 May 2019.
implicit none

 !Local vorticity and buoyancy arrays:
double precision:: zz(0:ny,0:nx),bb(0:ny,0:nx)

 !Arrays related to the conformal map:
double precision:: xo(0:ny,0:nx),yo(0:ny,0:nx)

 !Constant parameters:

 !Gravity (m/s^2)
double precision,parameter:: gravity=9.80665d0
 !Cross sectional width of tank (m)
double precision,parameter:: wtank=0.4d0 
 !Width of initial homogeneous surface layer (m)
double precision,parameter:: x1=0.6d0    
 !Depth of linearly stratified layer (m)
double precision,parameter:: ds=0.08d0   
 !Density at surface (kg/m^3)
double precision,parameter:: rho1=1025.d0
 !Density at bottom  (kg/m^3)
double precision,parameter:: rho2=1045.d0
 !Buoyancy at surface (m/s^2)
double precision,parameter:: b1=gravity*(rho2-rho1)/rho2
 !Squared buoyancy frequency in lin. strat. zone
double precision,parameter:: bfsq=b1/ds  

 !Other variables:
double precision:: xomin,xomax
double precision:: yomin,yomax
double precision:: vol,depth,y1,y2,y3
integer:: ix,iy

!-----------------------------------------------------------------------
! Read coords.r8 to get the original coordinates as a function of the
! conformal coordinates:
open(11,file='coords.r8',form='unformatted', &
     access='direct',status='old',recl=2*nbytes)
read(11,rec=1) xo
read(11,rec=2) yo
close(11)

xomin=minval(xo)
xomax=maxval(xo)
yomin=minval(yo)
yomax=maxval(yo)

!-----------------------------------------------------------------------
! Enter initialisation parameters:
write(*,'(2(a,f5.2))') ' The domain extends from x = ',xomin,' to ',xomax
write(*,'(2(a,f5.2))') '                and from y = ',yomin,' to ',yomax
write(*,*) ' with a ramped bottom in the right hand portion of the domain.'
write(*,*) 
write(*,*) ' The experimental tank is 0.4m in width.'
write(*,*) 
write(*,*) ' Initially, a homogeneous region of density 1025 kg/m^2 and of'
write(*,*) ' volume V lies above a linearly stratified zone of depth 0.08m' 
write(*,*) ' joining to a homogeneous region of density 1045 kg/m^2 below,'
write(*,*) ' for x < 0.6m from the left hand edge of the tank.'
write(*,*) 
write(*,*) ' On the right, a linearly stratified zone of depth 0.08m joins'
write(*,*) ' a homogeneous region of density 1045 kg/m^2 below.'
write(*,*) 
write(*,*) ' The density transitions across x = 0.6m over a distance of dx.'
write(*,*) 
write(*,*) ' Enter the volume V in litres:'
read(*,*) vol
depth=0.001d0*vol/(wtank*x1)
y3=yomax-ds
y2=yomax-depth
y1=y2-ds
write(*,'(a,f8.5,a)') ' This corresponds to a depth of ',depth,'m.'

!-----------------------------------------------------------------------
! Set up buoyancy profiles to the left and right of the transition:
bb=zero
zz=zero
do ix=0,nx
  do iy=0,ny
    if (xo(iy,ix) .lt. x1) then
      if (yo(iy,ix) .gt. y2) then
        bb(iy,ix)=b1
      else if (yo(iy,ix) .gt. y1) then
        bb(iy,ix)=bfsq*(yo(iy,ix)-y1)
      endif
    else
      if (yo(iy,ix) .gt. y3) then
        bb(iy,ix)=bfsq*(yo(iy,ix)-y3)
      endif
    endif
  enddo
enddo

 !Write vorticity distribution to file:
open(20,file='zz_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,zz
close(20)

 !Write buoyancy distribution to file:
open(20,file='bb_init.r8',form='unformatted', &
    & access='direct',status='replace',recl=2*nbytes)
write(20,rec=1) zero,bb
close(20)

!==========================================================================
end program bolus
