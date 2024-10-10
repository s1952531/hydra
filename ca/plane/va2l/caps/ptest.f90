!#########################################################################
!          Test code for the inversion of the non-hydrostatic 
!          pressure.  Compile using make ptest clean install
!#########################################################################

!  Code developed in early 2020 by D G Dritschel & M R Jalali @ St Andrews

program ptest

! Import spectral module:
use spectral

implicit none

! Declarations:
double precision,parameter:: ptol=1.d-7
! ptol: maximum relative rms NH pressure error

double precision:: epn1(ng,ng),epn2(ng,ng)
double precision:: pn1(ng,ng),pn2(ng,ng)
double precision:: h1(ng,ng),h2(ng,ng) !dimensionless depth anomalies
double precision:: htot1(ng,ng),htot2(ng,ng)
double precision:: hinv1(ng,ng),hinv2(ng,ng)
double precision:: mu(ng,ng),tt(ng,ng),tau(ng,ng)
double precision:: eps1(ng,ng),eps2(ng,ng)

double precision:: epi11(ng,ng),epi12(ng,ng)
double precision:: epi21(ng,ng),epi22(ng,ng)

double precision:: s1(ng,ng),s2(ng,ng)
double precision:: r1x(ng,ng),r1y(ng,ng)
double precision:: r2x(ng,ng),r2y(ng,ng)
double precision:: zpf1(ng,ng),zpf2(ng,ng)
double precision:: zx(ng,ng),zy(ng,ng)
double precision:: pnx(ng,ng),pny(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wkp(ng,ng)
double precision:: wkc(ng,ng),wkd(ng,ng)

double precision:: a1,r1,x1,y1,a2,r2,x2,y2,x,y,r
double precision:: perr

real:: tr4

integer:: ix,iy

!---------------------------------------------------------
! Initialise inversion constants and arrays:
call init_spectral

! Open output data files:
open(41,file='epn1.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(42,file='epn2.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)

open(51,file='apn1.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)
open(52,file='apn2.r4',form='unformatted',access='direct', &
                 & status='replace',recl=nbytes)

!======================================================================
! Define analytical non-hydrostatic pressure in order to set up source
! terms for determining the numerical non-hydrostatic pressure:
write(*,*) ' We consider a circular distribution of P_n in each layer'
write(*,*) ' having the form (A/2)*(1+cos(pi*r/R)) where A & R are'
write(*,*) ' constant, and r^2 = (x-X)^2+(y-Y)^2'
write(*,*)
write(*,*) ' Enter A, R, X & Y for layer 1:'
read(*,*) a1,r1,x1,y1
write(*,*) ' Enter A, R, X & Y for layer 2:'
read(*,*) a2,r2,x2,y2
a1=a1/two
a2=a2/two

! Check to make sure the circles are fully contained within the domain:
if (r1+abs(x1) .gt. pi) then
  write(*,*) ' Warning, circle 1 overlaps one of the x edges!  Stopping!'
  stop
endif

if (r1+abs(y1) .gt. pi) then
  write(*,*) ' Warning, circle 1 overlaps one of the y edges!  Stopping!'
  stop
endif

if (r2+abs(x2) .gt. pi) then
  write(*,*) ' Warning, circle 2 overlaps one of the x edges!  Stopping!'
  stop
endif

if (r2+abs(y2) .gt. pi) then
  write(*,*) ' Warning, circle 2 overlaps one of the y edges!  Stopping!'
  stop
endif

! Form arrays of NH pressure:
do ix=1,ng
  x=gl*dble(ix-1)-pi
  do iy=1,ng
    y=gl*dble(iy-1)-pi
    r=sqrt((x-x1)**2+(y-y1)**2)
    if (r .lt. r1) then
      epn1(iy,ix)=a1*(one+cos(pi*r/r1))
    else 
      epn1(iy,ix)=zero
    endif
    r=sqrt((x-x2)**2+(y-y2)**2)
    if (r .lt. r2) then
      epn2(iy,ix)=a2*(one+cos(pi*r/r2))
    else 
      epn2(iy,ix)=zero
    endif
  enddo
enddo

!======================================================================
! Define depth anomaly fields:
write(*,*) ' We also take tilde_h_1 = 0.3*cos(2x)*sin(3y) and'
write(*,*) '              tilde_h_2 = 0.2*sin(3x)*cos(2y).'

do ix=1,ng
  x=gl*dble(ix-1)-pi
  do iy=1,ng
    y=gl*dble(iy-1)-pi
    h1(iy,ix)=0.3d0*cos(two*x)*sin(three*y)
    h2(iy,ix)=0.2d0*sin(three*x)*cos(two*y)
  enddo
enddo

htot1=one+h1
htot2=one+h2

hinv1=one/htot1
call dealias(hinv1)
hinv2=one/htot2
call dealias(hinv2)

! mu & T arrays:
mu=mubar*htot2*hinv1
call dealias(mu)
tt=six/(four+three*mu)
call dealias(tt)

! tau array:
tau=tt-ttbar

! epsilon_1 & epsilon_2 arrays:
eps1=hinv1**2-one
call dealias(eps1)
eps2=hinv2**2-one
call dealias(eps2)

! Factors multiplying P_n1 & P_n2 in iteration below:
epi11=two*(tau+eps1*tt)
call dealias(epi11)
wka=eps2*tt
call dealias(wka) !wka = de-aliased eps2*tt
epi12=two*(two*tau+two*wka-three*eps2) !=two*(two*tau+eps2*(two*tt-three))
epi21=-f32*epi11
epi22=-six*(tau+wka-two*eps2)          !=-six*(tau+eps2*(tt-two))

! Z prefactors for pn1 & pn2:
zpf1=tt*f12*hinv1
call dealias(zpf1)
zpf2=-tt*(one+eps2)*htot1
call dealias(zpf2)

!======================================================================
! Determine `sources' h_1*tilde_gamma_1 & h_2*tilde_gamma_2 so that the
! pressure equations are satisfied by epn1 & epn2 above:
wka=epn1*(one+eps1)
wkb=epn2*(one+eps2)
s1=-tt*(two*wka-three*mu*wkb)
s2=-tt*(-three*wka+two*(one+three*mu)*wkb)

! Calculate H_1^2*grad(h_1)/h_1 and store in r1x & r1y:
wkp=h1*hbsq1
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,r1x,r1y)
r1x=r1x*hinv1
call dealias(r1x)
r1y=r1y*hinv1
call dealias(r1y)

! Define x & y components of H_1*Z:
wkc=zpf1*epn1+zpf2*epn2
call dealias(wkc)
zx=wkc*r1x
call dealias(zx)
zy=wkc*r1y
call dealias(zy)
! Note: H_1*Z = -H_1^2*T*(0.5*bar{P}_n1/h_1^2+bar{P}_n2/h_2^2)*grad{tilde{h_1}}

wkp=epn1
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
wkb=-hbsq1*rksq*wka ! Laplacian of H_1^2*P_n1
call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,pnx,pny)
! Compute divergence of mu*H_1*Z -> wkd:
wka=mu*zx
wkb=mu*zy
call divs(wka,wkb,wkc)
call spctop(ng,ng,wkc,wkd,xfactors,yfactors,xtrig,ytrig)
s1=s1+wkp-r1x*pnx-r1y*pny-htot1*wkd

! Calculate H_2^2*grad(h_2)/h_2 and store in r2x & r2y:
wkp=h2*hbsq2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,r2x,r2y)
r2x=r2x*hinv2
call dealias(r2x)
r2y=r2y*hinv2
call dealias(r2y)

wkp=epn2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
wkb=-hbsq2*rksq*wka ! Laplacian of H_2^2*P_n2
call spctop(ng,ng,wkb,wkp,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,pnx,pny)
! Compute divergence of Z*H_1 -> wkd:
call divs(zx,zy,wkc)
call spctop(ng,ng,wkc,wkd,xfactors,yfactors,xtrig,ytrig)
s2=s2+wkp-r2x*pnx-r2y*pny+htot2*hrat*wkd

!======================================================================
! Iterate to find approximate solution:

! Note, the fixed parts of the sources are in s1 & s2 above.
wkp=s1
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
! wka = s1 in spectral space

wkp=s2
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
! wkb = s2 in spectral space

! Obtain first approximation:
wkp=pi11*wka+pi12*wkb
call spctop(ng,ng,wkp,pn1,xfactors,yfactors,xtrig,ytrig)

wkp=pi21*wka+pi22*wkb
call spctop(ng,ng,wkp,pn2,xfactors,yfactors,xtrig,ytrig)

perr=one
do while (perr .gt. ptol)

  ! Define H_1*Z:
  wkc=zpf1*pn1+zpf2*pn2
  call dealias(wkc) 
  zx=wkc*r1x
  call dealias(zx)
  zy=wkc*r1y
  call dealias(zy)

  wkp=pn1
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  call gradient(wkb,pnx,pny)
  wka=mu*zx
  wkb=mu*zy
  call divs(wka,wkb,wkc)
  call spctop(ng,ng,wkc,wkd,xfactors,yfactors,xtrig,ytrig)
  wkp=s1+r1x*pnx+r1y*pny+epi11*pn1+epi12*pn2+htot1*wkd
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  ! wka = full S_1 in spectral space

  wkp=pn2
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  call gradient(wkb,pnx,pny)
  call divs(zx,zy,wkc)
  call spctop(ng,ng,wkc,wkd,xfactors,yfactors,xtrig,ytrig)
  wkp=s2+r2x*pnx+r2y*pny+epi21*pn1+epi22*pn2-htot2*hrat*wkd
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  ! wkb = full S_2 in spectral space

 !  Obtain next approximation (hold in wkc & wkd temporarily):
  wkp=pi11*wka+pi12*wkb
  call spctop(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)

  wkp=pi21*wka+pi22*wkb
  call spctop(ng,ng,wkp,wkd,xfactors,yfactors,xtrig,ytrig)

 !  Compute relative rms difference error:
  perr=sqrt(sum((wkc-pn1)**2+(wkd-pn2)**2)/sum(pn1**2+pn2**2))

 !  Copy wkc & wkd into pn1 & pn2:
  pn1=wkc
  pn2=wkd

  write(*,*) ' RMS NH pressure error = ',perr

enddo

!======================================================================
! Write data:
tr4=0.
write(41,rec=1) tr4,real(epn1)
write(42,rec=1) tr4,real(epn2)
write(51,rec=1) tr4,real(pn1)
write(52,rec=1) tr4,real(pn2)

! Close files:
close(41)
close(42)
close(51)
close(52)

write(*,*)
wka=(pn1-epn1)**2
wkb=epn1**2
write(*,*) ' R.m.s. error for P_n in layer 1 = ',sqrt(sum(wka)/sum(wkb))
wka=(pn2-epn2)**2
wkb=epn2**2
write(*,*) ' R.m.s. error for P_n in layer 2 = ',sqrt(sum(wka)/sum(wkb))
write(*,*)
write(*,*) ' See epn1.r4 & apn1.r4 for the exact and approximate forms'
write(*,*) ' of P_n in layer 1.'

write(*,*)
write(*,*) ' See epn2.r4 & apn2.r4 for the exact and approximate forms'
write(*,*) ' of P_n in layer 2.'


! End main program
end program ptest
!=======================================================================
