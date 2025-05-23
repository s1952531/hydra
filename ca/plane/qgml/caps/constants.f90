module constants

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! This module contains all the non-modifiable parameters as well as
! all quantities which never change throughout a simulation.
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 !Include all modifiable parameters for use below:
use parameters

 !Integer kind parameters used in contours.f90 and congen.f90:
 !  =>  dbleint represents up to 10^18
 !  =>  halfint represents up to 127
integer,parameter:: dbleint=selected_int_kind(16)
integer,parameter:: halfint=selected_int_kind(-1)

 !Grid resolution - 1:
integer,parameter:: ngm1=ng-1

!For reading & writing direct access data:
integer,parameter:: nhgridp=ng*ng,nhbytes=4*(nhgridp+1)
integer,parameter:: ngridp=nhgridp*nz,nbytes=4*(ngridp+1)

 !Fine grid used normally in contour -> grid conversion: 
integer,parameter:: mgf=4,ngf=mgf*ng
 !mgf:  fine grid/coarse grid ratio (4 is the standard value)

 !Ultra-fine grid used in contouring: 
integer,parameter:: mgu=16,ngu=mgu*ng
 !mgu:  ultra-fine grid/coarse grid ratio (16 is the standard value)

 !Maximum number of contour levels (used in surgery and congen):
integer,parameter:: nlevm=2000
 !nlevm: up to 2*nlevm contour levels are allowed (very generous)

 !Maximum number of contour nodes:
integer,parameter:: npzm=250*nhgridp,npm=npzm*nz
 !Maximum number of contours:
integer,parameter:: nm=npm/20+npm/200
 !Maximum number of nodes on any single contour:
integer,parameter:: nprm=npzm/10
 !Maximum number of nodes in any contour level:
integer,parameter:: nplm=npzm/2

 !Generic double precision numerical constants:
double precision,parameter:: zero=0.d0,one=1.d0,two=2.d0,three=3.d0
double precision,parameter:: four=4.d0,six=6.d0,f12=one/two,f13=one/three
double precision,parameter:: f23=two/three
double precision,parameter:: f14=one/four,f16=one/six,f32=three/two
double precision,parameter:: twopi=two*pi,thrpi=three*pi,pinv=one/pi
double precision,parameter:: small=1.d-12,small3=small*small*small
double precision,parameter:: oms=one-small

 !Grid constants:
double precision,parameter:: domarea=twopi*twopi
double precision,parameter:: gl=twopi/dble(ng),gli=dble(ng)/twopi
double precision,parameter:: garea=gl*gl,danorm=one/dble(ng*ng)

 !Logicals:
logical,parameter:: friction=(rekman > zero)
logical,parameter:: wind=(fwind /= zero)
logical,parameter:: beffect=(beta /= zero)

end module constants
