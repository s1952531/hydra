program euler
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Converts the output of pms.f90 to Euler angle coordinates and their 
! corresponding momenta.

!          *** Only sensible for 2 masses on a sphere ***

! Input files:
! -------------
! coordinates.asc:    positions  (x,y,z) of all masses
! speeds.asc:         velocities (u,v,w) of all masses
! strengths.asc:      "strengths" s = m/(4*pi) of all masses

! Output files:
! -------------
! ecoord.asc:         Euler angles (psi,theta,phi_1,phi_2) in radians 
! emomen.asc:         Corresponding momenta
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 ! Import parameters and constants:
use parameters
use constants

implicit none

double precision:: psi,theta,phi1,phi2
double precision:: ppsi,ptheta,pphi1,pphi2
double precision:: x1,y1,z1,x2,y2,z2
double precision:: px1,py1,pz1,px2,py2,pz2
double precision:: kap1,kap2,t
double precision:: cx,cz,ee,ff,gg,hh,aa,bb,cc
double precision:: z1sq,z2sq,rad,qp,qm
double precision:: sphi1,sphi2,cphi1,cphi2
double precision:: stheta,ctheta,spsi,cpsi
integer:: iread,nsol

!--------------------------------------------------------------------
if (n .ne. 2) then
  write(*,*) ' The number of masses must be 2 --- stopping!'
  stop
endif

! Open output files:
open(33,file='ecoord.asc',status='replace')
open(44,file='emomen.asc',status='replace')

! Open input files:
open(11,file='coordinates.asc',status='old')
open(22,file='speeds.asc',status='old')

! Read "strengths" (mass/(4*pi)):
open(20,file='strengths.asc',status='old')
read(20,*) kap1
read(20,*) kap2
close(20)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! Loop over all times in the data:
do
  iread=0
  ! Read mass positions & speeds:
  read(11,*,iostat=iread) t
  if (iread .ne. 0) exit
  read(11,*) x1,y1,z1
  read(11,*) x2,y2,z2

  read(22,*) t
  read(22,*) px1,py1,pz1
  read(22,*) px2,py2,pz2

  ! Convert speeds to momenta:
  px1=kap1*px1
  py1=kap1*py1
  pz1=kap1*pz1

  px2=kap2*px2
  py2=kap2*py2
  pz2=kap2*pz2

  ! Find the Euler angles from the Cartesian coordinates:
  cz=x1*y2-x2*y1
  cx=y1*z2-y2*z1
  ee=cx**2
  ff=ee+(x1*z2)**2+(x2*z1)**2
  z1sq=z1**2
  z2sq=z2**2
  gg=(ee+x1**2+x2**2)*z1sq
  hh=(two*x1*x2*z1)**2

  aa=ff**2-hh*z2sq
  bb=ff*gg-half*hh*(z1sq+z2sq)
  cc=gg**2-hh*z1sq

  rad=sqrt(bb**2-aa*cc)
  qp=sqrt((bb+rad)/aa)
  qm=sqrt((bb-rad)/aa)
  !qp & qm give the values of |sin(phi_1)|.

  nsol=0

  !Check each solution for consistency:
  !+++ +++ +++ +++ +++ +++ +++ +++ +++ +++
  if (qp .gt. one) then
    write(*,*) ' |sin(phi_1)| > 1 for positive root!'
  else
    !Consider only sin(theta) = z_1/sin(phi_1) >= 0:
    sphi1=sign(qp,z1)
    stheta=z1/sphi1
    sphi2=sphi1*z2/z1

    !++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++
    !Consider positive sign of cos(phi_1):
    cphi1=sqrt(one-sphi1**2)

    !+ + + + + + + + + + + + + + + + + + + 
    !Consider positive sign of cos(phi_2):
    cphi2=sqrt(one-sphi2**2)
    ctheta=z1*(x1*cphi2-x2*cphi1)/(cx*sphi1)
    if (abs(stheta) .le. one .and. abs(ctheta) .le. one) then
      spsi=(x1*cphi2-x2*cphi1)/cz
      cpsi=(y2*cphi1-y1*cphi2)/cz
      if (abs(spsi) .le. one .and. abs(cpsi) .le. one) then
        if (abs(x1*cpsi+y1*spsi-cphi1) .lt. 1.d-8 .and. &
            abs(x2*cpsi+y2*spsi-cphi2) .lt. 1.d-8 .and. &
            abs(y1*cpsi-x1*spsi-ctheta*sphi1) .lt. 1.d-8 .and. &
            abs(y2*cpsi-x2*spsi-ctheta*sphi2) .lt. 1.d-8 .and. &
            abs(stheta*sphi1-z1) .lt. 1.d-8 .and. &
            abs(stheta*sphi2-z2) .lt. 1.d-8) then
          psi=atan2(spsi,cpsi)
          theta=atan2(stheta,ctheta)
          phi1=atan2(sphi1,cphi1)
          phi2=atan2(sphi2,cphi2)
          nsol=nsol+1
          write(33,'(f12.5,4(1x,f14.11))') t,psi,theta,phi1,phi2
          write(44,'(f12.5,4(1x,f14.11))') t,cpsi,ctheta,cphi1,cphi2
        endif
      endif
    endif

    !- - - - - - - - - - - - - - - - - - -
    !Consider negative sign of cos(phi_2):
    cphi2=-sqrt(one-sphi2**2)
    ctheta=z1*(x1*cphi2-x2*cphi1)/(cx*sphi1)
    if (abs(stheta) .le. one .and. abs(ctheta) .le. one) then
      spsi=(x1*cphi2-x2*cphi1)/cz
      cpsi=(y2*cphi1-y1*cphi2)/cz
      if (abs(spsi) .le. one .and. abs(cpsi) .le. one) then
        if (abs(x1*cpsi+y1*spsi-cphi1) .lt. 1.d-8 .and. &
            abs(x2*cpsi+y2*spsi-cphi2) .lt. 1.d-8 .and. &
            abs(y1*cpsi-x1*spsi-ctheta*sphi1) .lt. 1.d-8 .and. &
            abs(y2*cpsi-x2*spsi-ctheta*sphi2) .lt. 1.d-8 .and. &
            abs(stheta*sphi1-z1) .lt. 1.d-8 .and. &
            abs(stheta*sphi2-z2) .lt. 1.d-8) then
          psi=atan2(spsi,cpsi)
          theta=atan2(stheta,ctheta)
          phi1=atan2(sphi1,cphi1)
          phi2=atan2(sphi2,cphi2)
          nsol=nsol+1
          write(33,'(f12.5,4(1x,f14.11))') t,psi,theta,phi1,phi2
          write(44,'(f12.5,4(1x,f14.11))') t,cpsi,ctheta,cphi1,cphi2
        endif
      endif
    endif

    !-- -- -- -- -- -- -- -- -- -- -- -- --
    !Consider negative sign of cos(phi_1):
    cphi1=-sqrt(one-sphi1**2)

    !+ + + + + + + + + + + + + + + + + + + 
    !Consider positive sign of cos(phi_2):
    cphi2=sqrt(one-sphi2**2)
    ctheta=z1*(x1*cphi2-x2*cphi1)/(cx*sphi1)
    if (abs(stheta) .le. one .and. abs(ctheta) .le. one) then
      spsi=(x1*cphi2-x2*cphi1)/cz
      cpsi=(y2*cphi1-y1*cphi2)/cz
      if (abs(spsi) .le. one .and. abs(cpsi) .le. one) then
        if (abs(x1*cpsi+y1*spsi-cphi1) .lt. 1.d-8 .and. &
            abs(x2*cpsi+y2*spsi-cphi2) .lt. 1.d-8 .and. &
            abs(y1*cpsi-x1*spsi-ctheta*sphi1) .lt. 1.d-8 .and. &
            abs(y2*cpsi-x2*spsi-ctheta*sphi2) .lt. 1.d-8 .and. &
            abs(stheta*sphi1-z1) .lt. 1.d-8 .and. &
            abs(stheta*sphi2-z2) .lt. 1.d-8) then
          psi=atan2(spsi,cpsi)
          theta=atan2(stheta,ctheta)
          phi1=atan2(sphi1,cphi1)
          phi2=atan2(sphi2,cphi2)
          nsol=nsol+1
          write(33,'(f12.5,4(1x,f14.11))') t,psi,theta,phi1,phi2
          write(44,'(f12.5,4(1x,f14.11))') t,cpsi,ctheta,cphi1,cphi2
        endif
      endif
    endif

    !- - - - - - - - - - - - - - - - - - -
    !Consider negative sign of cos(phi_2):
    cphi2=-sqrt(one-sphi2**2)
    ctheta=z1*(x1*cphi2-x2*cphi1)/(cx*sphi1)
    if (abs(stheta) .le. one .and. abs(ctheta) .le. one) then
      spsi=(x1*cphi2-x2*cphi1)/cz
      cpsi=(y2*cphi1-y1*cphi2)/cz
      if (abs(spsi) .le. one .and. abs(cpsi) .le. one) then
        if (abs(x1*cpsi+y1*spsi-cphi1) .lt. 1.d-8 .and. &
            abs(x2*cpsi+y2*spsi-cphi2) .lt. 1.d-8 .and. &
            abs(y1*cpsi-x1*spsi-ctheta*sphi1) .lt. 1.d-8 .and. &
            abs(y2*cpsi-x2*spsi-ctheta*sphi2) .lt. 1.d-8 .and. &
            abs(stheta*sphi1-z1) .lt. 1.d-8 .and. &
            abs(stheta*sphi2-z2) .lt. 1.d-8) then
          psi=atan2(spsi,cpsi)
          theta=atan2(stheta,ctheta)
          phi1=atan2(sphi1,cphi1)
          phi2=atan2(sphi2,cphi2)
          nsol=nsol+1
          write(33,'(f12.5,4(1x,f14.11))') t,psi,theta,phi1,phi2
          write(44,'(f12.5,4(1x,f14.11))') t,cpsi,ctheta,cphi1,cphi2
        endif
      endif
    endif
  endif

  !--- --- --- --- --- --- --- --- --- ---
  if (qm .gt. one) then
    write(*,*) ' |sin(phi_1)| > 1 for positive root!'
  else
    !Consider only sin(theta) = z_1/sin(phi_1) >= 0:
    sphi1=sign(qm,z1)
    stheta=z1/sphi1
    sphi2=sphi1*z2/z1

    !++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++ ++
    !Consider positive sign of cos(phi_1):
    cphi1=sqrt(one-sphi1**2)

    !+ + + + + + + + + + + + + + + + + + + 
    !Consider positive sign of cos(phi_2):
    cphi2=sqrt(one-sphi2**2)
    ctheta=z1*(x1*cphi2-x2*cphi1)/(cx*sphi1)
    if (abs(stheta) .le. one .and. abs(ctheta) .le. one) then
      spsi=(x1*cphi2-x2*cphi1)/cz
      cpsi=(y2*cphi1-y1*cphi2)/cz
      if (abs(spsi) .le. one .and. abs(cpsi) .le. one) then
        if (abs(x1*cpsi+y1*spsi-cphi1) .lt. 1.d-8 .and. &
            abs(x2*cpsi+y2*spsi-cphi2) .lt. 1.d-8 .and. &
            abs(y1*cpsi-x1*spsi-ctheta*sphi1) .lt. 1.d-8 .and. &
            abs(y2*cpsi-x2*spsi-ctheta*sphi2) .lt. 1.d-8 .and. &
            abs(stheta*sphi1-z1) .lt. 1.d-8 .and. &
            abs(stheta*sphi2-z2) .lt. 1.d-8) then
          psi=atan2(spsi,cpsi)
          theta=atan2(stheta,ctheta)
          phi1=atan2(sphi1,cphi1)
          phi2=atan2(sphi2,cphi2)
          nsol=nsol+1
          write(33,'(f12.5,4(1x,f14.11))') t,psi,theta,phi1,phi2
          write(44,'(f12.5,4(1x,f14.11))') t,cpsi,ctheta,cphi1,cphi2
        endif
      endif
    endif

    !- - - - - - - - - - - - - - - - - - -
    !Consider negative sign of cos(phi_2):
    cphi2=-sqrt(one-sphi2**2)
    ctheta=z1*(x1*cphi2-x2*cphi1)/(cx*sphi1)
    if (abs(stheta) .le. one .and. abs(ctheta) .le. one) then
      spsi=(x1*cphi2-x2*cphi1)/cz
      cpsi=(y2*cphi1-y1*cphi2)/cz
      if (abs(spsi) .le. one .and. abs(cpsi) .le. one) then
        if (abs(x1*cpsi+y1*spsi-cphi1) .lt. 1.d-8 .and. &
            abs(x2*cpsi+y2*spsi-cphi2) .lt. 1.d-8 .and. &
            abs(y1*cpsi-x1*spsi-ctheta*sphi1) .lt. 1.d-8 .and. &
            abs(y2*cpsi-x2*spsi-ctheta*sphi2) .lt. 1.d-8 .and. &
            abs(stheta*sphi1-z1) .lt. 1.d-8 .and. &
            abs(stheta*sphi2-z2) .lt. 1.d-8) then
          psi=atan2(spsi,cpsi)
          theta=atan2(stheta,ctheta)
          phi1=atan2(sphi1,cphi1)
          phi2=atan2(sphi2,cphi2)
          nsol=nsol+1
          write(33,'(f12.5,4(1x,f14.11))') t,psi,theta,phi1,phi2
          write(44,'(f12.5,4(1x,f14.11))') t,cpsi,ctheta,cphi1,cphi2
        endif
      endif
    endif

    !-- -- -- -- -- -- -- -- -- -- -- -- --
    !Consider negative sign of cos(phi_1):
    cphi1=-sqrt(one-sphi1**2)

    !+ + + + + + + + + + + + + + + + + + + 
    !Consider positive sign of cos(phi_2):
    cphi2=sqrt(one-sphi2**2)
    ctheta=z1*(x1*cphi2-x2*cphi1)/(cx*sphi1)
    if (abs(stheta) .le. one .and. abs(ctheta) .le. one) then
      spsi=(x1*cphi2-x2*cphi1)/cz
      cpsi=(y2*cphi1-y1*cphi2)/cz
      if (abs(spsi) .le. one .and. abs(cpsi) .le. one) then
        if (abs(x1*cpsi+y1*spsi-cphi1) .lt. 1.d-8 .and. &
            abs(x2*cpsi+y2*spsi-cphi2) .lt. 1.d-8 .and. &
            abs(y1*cpsi-x1*spsi-ctheta*sphi1) .lt. 1.d-8 .and. &
            abs(y2*cpsi-x2*spsi-ctheta*sphi2) .lt. 1.d-8 .and. &
            abs(stheta*sphi1-z1) .lt. 1.d-8 .and. &
            abs(stheta*sphi2-z2) .lt. 1.d-8) then
          psi=atan2(spsi,cpsi)
          theta=atan2(stheta,ctheta)
          phi1=atan2(sphi1,cphi1)
          phi2=atan2(sphi2,cphi2)
          nsol=nsol+1
          write(33,'(f12.5,4(1x,f14.11))') t,psi,theta,phi1,phi2
          write(44,'(f12.5,4(1x,f14.11))') t,cpsi,ctheta,cphi1,cphi2
        endif
      endif
    endif

    !- - - - - - - - - - - - - - - - - - -
    !Consider negative sign of cos(phi_2):
    cphi2=-sqrt(one-sphi2**2)
    ctheta=z1*(x1*cphi2-x2*cphi1)/(cx*sphi1)
    if (abs(stheta) .le. one .and. abs(ctheta) .le. one) then
      spsi=(x1*cphi2-x2*cphi1)/cz
      cpsi=(y2*cphi1-y1*cphi2)/cz
      if (abs(spsi) .le. one .and. abs(cpsi) .le. one) then
        if (abs(x1*cpsi+y1*spsi-cphi1) .lt. 1.d-8 .and. &
            abs(x2*cpsi+y2*spsi-cphi2) .lt. 1.d-8 .and. &
            abs(y1*cpsi-x1*spsi-ctheta*sphi1) .lt. 1.d-8 .and. &
            abs(y2*cpsi-x2*spsi-ctheta*sphi2) .lt. 1.d-8 .and. &
            abs(stheta*sphi1-z1) .lt. 1.d-8 .and. &
            abs(stheta*sphi2-z2) .lt. 1.d-8) then
          psi=atan2(spsi,cpsi)
          theta=atan2(stheta,ctheta)
          phi1=atan2(sphi1,cphi1)
          phi2=atan2(sphi2,cphi2)
          nsol=nsol+1
          write(33,'(f12.5,4(1x,f14.11))') t,psi,theta,phi1,phi2
          write(44,'(f12.5,4(1x,f14.11))') t,cpsi,ctheta,cphi1,cphi2
        endif
      endif
    endif
  endif

  if (nsol .ne. 2) write(*,'(i1,a,f12.5)') nsol,' solutions found at t = ',t

  ! Find the momenta corresponding to the Euler angles:



  ! Write out results (view using plotcol):
!  write(44,'(f12.5,4(1x,f14.9))' ) t,ppsi,ptheta,pphi1,pphi2
enddo

! Close files:
close(11)
close(22)
close(33)
close(44)

write(*,*)
write(*,*) ' The Euler angles are available in ecoord.asc'
write(*,*) ' and the corresponding momenta are available in emomen.asc'

!End main program
end program
!=======================================================================
