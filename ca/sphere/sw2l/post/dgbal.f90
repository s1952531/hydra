!######################################################################
!  Finds balanced fields given only the total PV in qq_init.r8 or
!  the PV anomaly in qAnom1.r4 & qAnom2.r4, using first-order delta-gamma
!  balance in which delta_t and gamma_t are set to zero to diagnose
!  all fields. Here, delta is the velocity divergence while gamma is
!  the acceleration divergence.

!  This routine can be used for both initialisation and post-processing.

!  Adapted to two layers on 12 November 2026 by D G Dritschel @ Kyoto
!######################################################################

program dgbal

 !Import spectral module, constants and parameters:
use spectral

implicit none

 !Global commons (to this program):

 !Physical fields (PV, delta, gamma, thickness and relative vorticity):
double precision:: qAnom1(numLatGridPts,numLongGridPts),qAnom2(numLatGridPts,numLongGridPts)
double precision:: d1(numLatGridPts,numLongGridPts),d2(numLatGridPts,numLongGridPts)
double precision:: g1(numLatGridPts,numLongGridPts),g2(numLatGridPts,numLongGridPts)
double precision:: h1(numLatGridPts,numLongGridPts),h2(numLatGridPts,numLongGridPts)
double precision:: z1(numLatGridPts,numLongGridPts),z2(numLatGridPts,numLongGridPts)

 !Physical velocity field:
double precision:: u1(numLatGridPts,numLongGridPts),u2(numLatGridPts,numLongGridPts)
double precision:: v1(numLatGridPts,numLongGridPts),v2(numLatGridPts,numLongGridPts)

 !Spectral fields (delta and gamma):
double precision:: deltaSpec1(numLatGridPts,numLongGridPts),deltaSpec2(numLatGridPts,numLongGridPts)
double precision:: gammaSpec1(numLatGridPts,numLongGridPts),gammaSpec2(numLatGridPts,numLongGridPts)

 !Bottom topography and associated spectral fields (if used):
double precision,allocatable,dimension(:,:):: botTopo,specDeltaTopoSrc,specGammaTopoSrc

 !Other quantities:
double precision:: ekin,epot,etot,angm,t
real:: qqr4(numLatGridPts,numLongGridPts),tr4
integer:: iopt,loop,iread,i

 !---------------------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !If fixed topography is present, allocate memory and read in its form:
if (topoFlag) then
  allocate(botTopo(numLatGridPts,numLongGridPts),specDeltaTopoSrc(numLatGridPts,numLongGridPts),specGammaTopoSrc(numLatGridPts,numLongGridPts))   
  open(11,file='bb_init.r8',form='unformatted',&
        access='direct',status='old',recl=2*nbytes)
  read(11,rec=1) t,botTopo
  close(11)
  !Store its contribution to the gamma and delta source terms:
  specGammaTopoSrc=csq1*botTopo
  call forfft(numLatGridPts,numLongGridPts,specGammaTopoSrc,trig,factors)
  call laplace(specGammaTopoSrc,specDeltaTopoSrc)
  call revfft(numLatGridPts,numLongGridPts,specDeltaTopoSrc,trig,factors)
endif

write(*,*) ' Choose one of the following options:'
write(*,*)
write(*,*) ' (0) initialise from data in qq_init.r8'
write(*,*) ' (1) post-process from data in qAnom1.r4 & qAnom2.r4'
write(*,*)
write(*,*) ' Option?'
read(*,*) iopt
write(*,*)

if (iopt .eq. 0) then
   !Only find balance at t = 0 to initialise; read total PV:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='old',recl=2*nbytes)
  read(11,rec=1) t,qAnom1
  read(11,rec=2) t,qAnom2
  close(11)
  
   !Find balanced fields:
  call balance

   !Save initial (total) PV:
  open(11,file='qq_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,qAnom1
  write(11,rec=2) zero,qAnom2
  close(11)
   !Note, the PV is modified by an overall constant offset due to the
   !mathematical requirement that the domain integral of the relative
   !vorticity be zero: see main_invert in spectral.f90.

   !Save initial velocity divergence:
  open(11,file='dd_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,d1
  write(11,rec=2) zero,d2
  close(11)

   !Save initial acceleration divergence:
  open(11,file='gg_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,g1
  write(11,rec=2) zero,g2
  close(11)

   !Save initial dimensionless height anomaly:
  open(11,file='hh_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,h1
  write(11,rec=2) zero,h2
  close(11)

   !Save initial relative vorticity:
  open(11,file='zz_init.r8',form='unformatted', &
        access='direct',status='replace',recl=2*nbytes)
  write(11,rec=1) zero,z1
  write(11,rec=2) zero,z2
  close(11)

  write(*,*)
  write(*,*) ' The balanced fields are available in xx_init.r8 where'
  write(*,*) ' xx = qq, dd, gg, hh or zz.'
  write(*,*)
  
else
   !Process all times using only the PV anomaly in qAnom1.r4 & qAnom2.r4:
  open(30,file='evolution/qAnom1.r4',form='unformatted', &
        access='direct',status='old',recl=nbytes)
  open(40,file='evolution/qAnom2.r4',form='unformatted', &
        access='direct',status='old',recl=nbytes)

  open(31,file='evolution/bq1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(41,file='evolution/bq2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(32,file='evolution/bd1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(42,file='evolution/bd2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(33,file='evolution/bh1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(43,file='evolution/bh2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(34,file='evolution/bg1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(44,file='evolution/bg2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(35,file='evolution/bz1.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)
  open(45,file='evolution/bz2.r4',form='unformatted', &
        access='direct',status='replace',recl=nbytes)

  open(51,file='evolution/becomp.asc',status='replace')

   !Read data at all available times and balance:
  loop=0
  do
    loop=loop+1
    iread=0
     !Read the PV anomaly from qAnom1.r4 & qAnom2.r4:
    read(30,rec=loop,iostat=iread) tr4,qqr4
    if (iread .ne. 0) exit
    qAnom1=dble(qqr4)
    t=dble(tr4)

    read(40,rec=loop,iostat=iread) tr4,qqr4
    if (iread .ne. 0) exit
    qAnom2=dble(qqr4)
    t=dble(tr4)

     !Define total PV by adding f:
    do i=1,numLongGridPts
      qAnom1(:,i)=qAnom1(:,i)+cof
      qAnom2(:,i)=qAnom2(:,i)+cof
    enddo

     !Balance using only the PV anomaly at time t = t:
    write(*,'(a,f9.2)') ' *** Processing t = ',t
    call balance

     !Redefine PV anomaly by subtracting f for output:
    do i=1,numLongGridPts
      qAnom1(:,i)=qAnom1(:,i)-cof
      qAnom2(:,i)=qAnom2(:,i)-cof
    enddo

     !Write data:
    write(31,rec=loop) tr4,real(qAnom1)
    write(41,rec=loop) tr4,real(qAnom2)

    write(32,rec=loop) tr4,real(d1)
    write(42,rec=loop) tr4,real(d2)

    write(33,rec=loop) tr4,real(h1)
    write(43,rec=loop) tr4,real(h2)

    write(34,rec=loop) tr4,real(g1)
    write(44,rec=loop) tr4,real(g2)

    write(35,rec=loop) tr4,real(z1)
    write(45,rec=loop) tr4,real(z2)

     !Write balanced energy components and angular momentum:
    write(51,'(f13.6,4(1x,f16.9))') t,ekin,epot,etot,angm
  enddo

  close(30)
  close(40)
  close(31)
  close(41)
  close(32)
  close(42)
  close(33)
  close(43)
  close(34)
  close(44)
  close(35)
  close(45)
  close(51)

  write(*,*)
  write(*,*) ' The balanced fields are available in evolution/bxx.r4'
  write(*,*) ' where xx = qAnom1, qAnom2, d1, d2, h1, h2, g1, g2, z1 & z2.'
  write(*,*)
  write(*,*) ' Also, becomp.asc contains the energy components.'
  write(*,*)

endif


 !Internal subroutine definitions (inherit global variables):

contains

!=======================================================================

subroutine balance
 !Given the *total* PV in qAnom1 & qAnom2, this routine finds all the other
 !fields by requiring delta_t = gamma_t = 0.
  
implicit none

 !Local parameters and variables:

double precision,parameter:: tole=1.d-10
 !tole: relative energy norm error in successive iterates when finding
 !      dj, gj, hj, uj & vj from qj(for j = 1 & 2). The energy norm is
 !      <du1^2+dv1^2+c_1^2*dh1^2> + (rho_2/rho_1)*(H_2/H_1)*
 !      <du2^2+dv2^2+c_2^2*dh2^2>
 !      where <:> means a horizontal domain average and (duj,dvj,dhj)
 !      is the difference from the previous guess.

 !Divergent part of the velocity field:
double precision:: ud1(numLatGridPts,numLongGridPts),ud2(numLatGridPts,numLongGridPts)
double precision:: vd1(numLatGridPts,numLongGridPts),vd2(numLatGridPts,numLongGridPts)

 !Previous guesses needed to compute the energy norm:
double precision:: h1pre(numLatGridPts,numLongGridPts),h2pre(numLatGridPts,numLongGridPts)
double precision:: u1pre(numLatGridPts,numLongGridPts),u2pre(numLatGridPts,numLongGridPts)
double precision:: v1pre(numLatGridPts,numLongGridPts),v2pre(numLatGridPts,numLongGridPts)

 !Used to store sources for the height inversion:
double precision:: sh1(numLatGridPts,numLongGridPts),sh2(numLatGridPts,numLongGridPts)

 !Spectral nonlinear delta dn gamma tendencies
double precision:: sds1(numLatGridPts,numLongGridPts),sds2(numLatGridPts,numLongGridPts)
double precision:: sgs1(numLatGridPts,numLongGridPts),sgs2(numLatGridPts,numLongGridPts)

 !Work arrays:
double precision:: wka(numLatGridPts,numLongGridPts),wkb(numLatGridPts,numLongGridPts),wkc(numLatGridPts,numLongGridPts),wkd(numLatGridPts,numLongGridPts)

 !Other quantities:
double precision:: enorm
integer:: i,m

 !----------------------------------------------------------------------
 !Start with zero divergence and acceleration divergence:
deltaSpec1=zero
deltaSpec2=zero
gammaSpec1=zero
gammaSpec2=zero

 !Initialise previous guesses for hj, uj & vj:
h1=zero
h2=zero
u1=zero
u2=zero
v1=zero
v2=zero

h1pre=zero
h2pre=zero
u1pre=zero
u2pre=zero
v1pre=zero
v2pre=zero

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 !Iteratively solve for hj, uj & vj (initial guess for dj = gj = 0):

 !Energy norm error (must be > tole to start):
enorm=one
do while (enorm .gt. tole)
   !Correct average PV by enforcing zero average vorticity:
  wkc=one+h1
  wka=qAnom1*wkc
  qAnom1=qAnom1-average(wka)
   !Recall qAnom1 = (z1 + cof)/(1 + h1)

  wkd=one+h2
  wkb=qAnom2*wkd
  qAnom2=qAnom2-average(wkb)
   !Recall qAnom2 = (z2 + cof)/(1 + h2)

   !---------------------------------------------------------------
   !Define initial "sources" sh1 & sh2 needed to invert h1 & h2:
  do i=1,numLongGridPts
    sh1(:,i)=beta*u1(:,i)-cof*(qAnom1(:,i)-cof)*wkc(:,i)
    sh2(:,i)=beta*u2(:,i)-cof*(qAnom2(:,i)-cof)*wkd(:,i)
  enddo
   !These use the fact that dsj = gsj = 0 in this first iterative loop.
   !Convert sources to semi-spectral space and de-alias:
  call dealias(sh1)
  call dealias(sh2)

   !Solve for h1 & h2 in mode space, where the equations are decoupled:
  wka=sh1*vm11+sh2*vm12
  call helm1(wka,sh1,wkb)
   !Solves f^2*h_1 - \tilde{c}_1^2*Lap(h_1) = wka, where h_1 is h for mode 1
   !The solution is stored in sh1.

  wka=sh1*vm21+sh2*vm22
  call helm2(wka,sh2,wkb)
   !Solves f^2*h_2 - \tilde{c}_2^2*Lap(h_2) = wka, where h_2 is h for mode 2
   !The solution is stored in sh2.

   !Return h1 & h2 to layers as semi-spectral arrays, then physical ones:
  h1=sh1*vt11+sh2*vt12
  call revfft(numLatGridPts,numLongGridPts,h1,trig,factors) 
  h2=sh1*vt21+sh2*vt22
  call revfft(numLatGridPts,numLongGridPts,h2,trig,factors) 

   !---------------------------------------------------------------
   !Correct average PV by enforcing zero average vorticity:
  wkc=one+h1
  wka=qAnom1*wkc
  qAnom1=qAnom1-average(wka)

  wkd=one+h2
  wkb=qAnom2*wkd
  qAnom2=qAnom2-average(wkb)

   !Compute relative vorticity (now guaranteed to have zero average):
  do i=1,numLongGridPts
    z1(:,i)=qAnom1(:,i)*wkc(:,i)-cof
    z2(:,i)=qAnom2(:,i)*wkd(:,i)-cof
  enddo

   !Compute non-divergent velocity field (as here ds = 0):

   !==>  Layer 1  <==!
  
   !Convert the relative vorticity to semi-spectral space as wkd and de-alias:
  wkd=z1
  call dealias(wkd)

   !Invert Laplace's operator on relative vorticity (wkd):
  call laplinv(wkd,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Compute d(wka)/dlon = wkc:
  call deriv(numLatGridPts,numLongGridPts,rk,wka,wkc)

   !Define velocity in semi-spectral space:
  do m=1,numLongGridPts
    u1(:,m)=-wkb(:,m)
    v1(:,m)=secantLat*wkc(:,m)
  enddo

   !Return to physical space:
  call revfft(numLatGridPts,numLongGridPts,u1,trig,factors)
  call revfft(numLatGridPts,numLongGridPts,v1,trig,factors)

   !==>  Layer 2  <==!
  
   !Convert the relative vorticity to semi-spectral space as wkd and de-alias:
  wkd=z2
  call dealias(wkd)

   !Invert Laplace's operator on relative vorticity (wkd):
  call laplinv(wkd,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Compute d(wka)/dlon = wkc:
  call deriv(numLatGridPts,numLongGridPts,rk,wka,wkc)

   !Define velocity in semi-spectral space:
  do m=1,numLongGridPts
    u2(:,m)=-wkb(:,m)
    v2(:,m)=+secantLat*wkc(:,m)
  enddo

   !Return to physical space:
  call revfft(numLatGridPts,numLongGridPts,u2,trig,factors)
  call revfft(numLatGridPts,numLongGridPts,v2,trig,factors)

   !Compute relative error (measured by the energy norm):
  wka=         (u1-u1pre)**2+(v1-v1pre)**2+csq1*(h1-h1pre)**2 &
       +mubar*((u2-u2pre)**2+(v2-v2pre)**2+csq2*(h2-h2pre)**2)
  wkb=u1**2+v1**2+csq1*h1**2+mubar*(u2**2+v2**2+csq2*h2**2)
  enorm=average(wka)/(average(wkb)+small)
   !Small = 1.d-12 is added to avoid dividing by zero for a rest state

   !Store current guess to recompute error in the next iteration:
  h1pre=h1
  h2pre=h2
  u1pre=u1
  u2pre=u2
  v1pre=v1
  v2pre=v2
enddo
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !Passing this, we have converged.

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !Iterate to find all the balanced fields, including delta & gamma:
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 !Energy norm error (must be > tole to start):
enorm=one
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
do while (enorm .gt. tole)

   !Compute nonlinear source terms for delta and gamma:
  call source(sds1,sds2,sgs1,sgs2)

   !Solve for delta in mode space after projection:
   !Mode 1:
  wka=sgs1*vm11+sgs2*vm12
  call helm1(wka,wkc,wkb) !wkb is not used
   !Solves f^2*d_1 - \tilde{c}_1^2*Lap(d_1) = sgs1,
   !where d_1 is delta for mode 1; the solution is stored in wkc.
   !ctsq1 = \tilde{c}_1^2 in parameters.f90

   !Mode 2:
  wka=sgs1*vm21+sgs2*vm22
  call helm2(wka,wkd,wkb) !wkb is not used
   !Solves f^2*d_2 - \tilde{c}_2^2*Lap(d_2) = sgs2,
   !where d_2 is delta for mode 2; the solution is stored in wkd.
   !ctsq2 = \tilde{c}_2^2 in parameters.f90

   !Return modes wkc & wkd to layers as the semi-spectral arrays deltaSpec1 & deltaSpec2:
  deltaSpec1=wkc*vt11+wkd*vt12
  deltaSpec2=wkc*vt21+wkd*vt22

   !gamma is just the negative of the divergence source term:
  gammaSpec1=-sds1
  gammaSpec2=-sds2

   !-----------------------------------------------------
   !Compute divergent part of the velocity field (ud,vd):

   !Invert Laplace's operator on the divergence (ds) in both layers:
  call laplinv(deltaSpec1,wka,vd1)
  call laplinv(deltaSpec2,wkb,vd2)
   !Here the divergence potential is wka or wkb, and its latitudinal
   !derivative is vd1 or vd2, the divergent meridional velocity
   !(all in semi-spectral space)

   !Compute divergent zonal velocity ud in both layers:
  call deriv(numLatGridPts,numLongGridPts,rk,wka,ud1)
  call deriv(numLatGridPts,numLongGridPts,rk,wkb,ud2)
  do m=1,numLongGridPts
    ud1(:,m)=secantLat*ud1(:,m)
    ud2(:,m)=secantLat*ud2(:,m)
  enddo
   !clati = 1/r, where r = cos(lat)

   !--------------------------------------------------------
   !Compute non-divergent velocity field and add to (ud,vd):

   !First correct average PV by enforcing zero average vorticity:
  wkc=one+h1
  wka=qAnom1*wkc
  qAnom1=qAnom1-average(wka)

  wkd=one+h2
  wka=qAnom2*wkd
  qAnom2=qAnom2-average(wka)

   !Compute relative vorticity (now guaranteed to have zero average):
  do i=1,numLongGridPts
    z1(:,i)=qAnom1(:,i)*wkc(:,i)-cof
    z2(:,i)=qAnom2(:,i)*wkd(:,i)-cof
  enddo

   !==>  Layer 1  <==!
  
   !Convert the relative vorticity to semi-spectral space as wkd and de-alias:
  wkd=z1
  call dealias(wkd)
   !Preserve wkd until it is used to define sh1 below

   !Invert Laplace's operator on relative vorticity:
  call laplinv(wkd,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Compute d(wka)/dlon = wkc:
  call deriv(numLatGridPts,numLongGridPts,rk,wka,wkc)

   !Define h1 in spectral space as sh1 for use in the loop below:
  sh1=h1
  call forfft(numLatGridPts,numLongGridPts,sh1,trig,factors)

  do m=1,numLongGridPts
    u1(:,m)=ud1(:,m)-wkb(:,m)
    v1(:,m)=vd1(:,m)+secantLat*wkc(:,m)
    sh1(:,m)=gammaSpec1(:,m)+beta*u1(:,m)-cof*(wkd(:,m)-cof*sh1(:,m))
  enddo
   !All the terms in this loop are semi-spectral. Note, we will solve
   !f^2*h_1 - Lap(c_1^2*h_1+(rho_2/rho_1)*c_2^2*h_2) = sh1, which follows
   !from the definition of gamma_1.

   !==>  Layer 2  <==!
  
   !Convert the relative vorticity to semi-spectral space as wkd and de-alias:
  wkd=z2
  call dealias(wkd)
   !Preserve wkd until it is used to define sh2 below
  
   !Invert Laplace's operator on relative vorticity:
  call laplinv(wkd,wka,wkb)
   !Here the streamfunction is wka while d(wka)/dlat = wkb.

   !Compute d(wka)/dlon = wkc:
  call deriv(numLatGridPts,numLongGridPts,rk,wka,wkc)

   !Define h2 in spectral space as sh2 for use in the loop below:
  sh2=h2
  call forfft(numLatGridPts,numLongGridPts,sh2,trig,factors)

  do m=1,numLongGridPts
    u2(:,m)=ud2(:,m)-wkb(:,m)
    v2(:,m)=vd2(:,m)+secantLat*wkc(:,m)
    sh2(:,m)=gammaSpec2(:,m)+beta*u2(:,m)-cof*(wkd(:,m)-cof*sh2(:,m))
  enddo
   !All the terms in this loop are semi-spectral. Note, we will solve
   !f^2*h_2 - Lap(c_1^2*h_1+c_2^2*h_2) = sh2, which follows from
   !the definition of gamma_2.

   !-------------------------------------------------------------------
   !Solve for h1 & h2 in mode space, where the equations are decoupled:
  wka=sh1*vm11+sh2*vm12
  call helm1(wka,sh1,wkb)
   !Solves f^2*h_1 - \tilde{c}_1^2*Lap(h_1) = wka, where h_1 is h for mode 1
   !The solution is stored in sh1.
   !ctsq1 = \tilde{c}_1^2 in parameters.f90

  wka=sh1*vm21+sh2*vm22
  call helm2(wka,sh2,wkb)
   !Solves f^2*h_2 - \tilde{c}_2^2*Lap(h_2) = wka, where h_2 is h for mode 2
   !The solution is stored in sh2.
   !ctsq2 = \tilde{c}_2^2 in parameters.f90

   !Return h1 & h2 to layers as semi-spectral arrays, then physical ones:
  h1=sh1*vt11+sh2*vt12
  call revfft(numLatGridPts,numLongGridPts,h1,trig,factors) 
  h2=sh1*vt21+sh2*vt22
  call revfft(numLatGridPts,numLongGridPts,h2,trig,factors) 

   !Return semi-spectral velocity field to physical space:
  call revfft(numLatGridPts,numLongGridPts,u1,trig,factors)
  call revfft(numLatGridPts,numLongGridPts,u2,trig,factors)

  call revfft(numLatGridPts,numLongGridPts,v1,trig,factors)
  call revfft(numLatGridPts,numLongGridPts,v2,trig,factors)

   !-----------------------------------------------------
   !Compute relative error (measured by the energy norm):
  wka=         (u1-u1pre)**2+(v1-v1pre)**2+csq1*(h1-h1pre)**2 &
       +mubar*((u2-u2pre)**2+(v2-v2pre)**2+csq2*(h2-h2pre)**2)
  wkb=u1**2+v1**2+csq1*h1**2+mubar*(u2**2+v2**2+csq2*h2**2)
  enorm=average(wka)/(average(wkb)+small)
   !Small = 1.d-12 is added to avoid dividing by zero for a rest state

   !Store current guess to recompute error in the next iteration:
  h1pre=h1
  h2pre=h2
  u1pre=u1
  u2pre=u2
  v1pre=v1
  v2pre=v2
enddo
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 !Passing this, we have converged.

 !Return delta and gamma to physical space:
d1=deltaSpec1
call revfft(numLatGridPts,numLongGridPts,d1,trig,factors)
d2=deltaSpec2
call revfft(numLatGridPts,numLongGridPts,d2,trig,factors)

g1=gammaSpec1
call revfft(numLatGridPts,numLongGridPts,g1,trig,factors)
g2=gammaSpec2
call revfft(numLatGridPts,numLongGridPts,g2,trig,factors)

 !Compute energy components, divided by rho_1*H where H = H_1+H_2:
wka=hAvgLayer1*(one+h1) !Physical thickness of layer 1
wkb=hAvgLayer2*(one+h2) !Physical thickness of layer 2
 !*** Do not re-use wka & wkb ***

 !Kinetic energy:
wkc=wka*(u1**2+v1**2)+densityRatio*wkb*(u2**2+v2**2) !kinetic energy density
ekin=twoPi*average(wkc)
 !This is the area integral of the kinetic energy density
 !densityRatio = rho_2/rho_1 (see constants.f90)

 !Potential energy:
if (topoFlag) then
   wkc=hAvgLayer1*botTopo+wka !Physical height of the top of layer 1
   wkc=wka*(two*wkc-wka)+densityRatio*wkb*(two*wkc+wkb)-dape
else
   wkc=wka**2+densityRatio*wkb*(two*wka+wkb)-dape
endif
epot=twoPi*gravity*average(wkc)
 !wkc: potential energy density
 !dape = hbar_1**2+(rho_2/rho_1)*hbar_2*(2*hbar_1+hbar_2): ensures that 
 !       a state of rest has zero potential energy (see constants.f90)
 !*** wka & wkb are now free to re-use ***

 !Total energy (sum of the above):
etot=ekin+epot

 !Compute angular momentum:
do i=1,numLongGridPts
  wka(:,i)=cosLat*((one+h1(:,i))*u1(:,i)+omega*cosLat*h1(:,i))
  wkb(:,i)=cosLat*((one+h2(:,i))*u2(:,i)+omega*cosLat*h2(:,i))
enddo
 !wka & wkb: angular momentum densities in layers 1 and 2
 !Omega = planetary rotation rate (= 2*pi, see constants.f90)
wkc=wka+densityRatio*wkb
angm=fourpi*average(wkc)

return
end subroutine balance

!=======================================================================

subroutine source(sds1,sds2,sgs1,sgs2)

! Gets the nonlinear parts of the source terms (sds,sgs) for divergence
! and acceleration divergence (ds,gs) --- all in semi-spectral space
! and for each *layer*. See further comments in caps/evolution.f90.

! The semi-spectral fields deltaSpec1, deltaSpec2, gammaSpec1 & gammaSpec2 are spectrally truncated.
! Note, h, u, v & z obtained by main_invert before calling this routine
! are also spectrally truncated. These variables are globally available.

implicit none

 !Passed variables:
double precision:: sds1(numLatGridPts,numLongGridPts),sds2(numLatGridPts,numLongGridPts),sgs1(numLatGridPts,numLongGridPts),sgs2(numLatGridPts,numLongGridPts)

 !Local variables:
double precision:: wka(numLatGridPts,numLongGridPts),wkb(numLatGridPts,numLongGridPts),wkc(numLatGridPts,numLongGridPts),wkd(numLatGridPts,numLongGridPts)
double precision:: wke(numLatGridPts,numLongGridPts),wkf(numLatGridPts,numLongGridPts),wkg(numLatGridPts,numLongGridPts),wkh(numLatGridPts,numLongGridPts)
double precision:: om1(numLatGridPts,numLongGridPts),om2(numLatGridPts,numLongGridPts)
double precision:: usq1(numLatGridPts,numLongGridPts),usq2(numLatGridPts,numLongGridPts)
double precision:: d1(numLatGridPts,numLongGridPts),d2(numLatGridPts,numLongGridPts),rhs(numLatGridPts),avgval
integer:: i

!--------------------------------------------------------------
 !First define Omega = u/r (where r = cos(lat)) for each layer:
do i=1,numLongGridPts
  om1(:,i)=secantLat*u1(:,i)
  om2(:,i)=secantLat*u2(:,i)
enddo
 !*** Do not re-use om1 & om2 below ***

!---------------------------------------------------------------
 !delta source (can be computed in each layer independently):

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain delta source for layer 1 (sds1):

 !Compute A = z*u/r + (1/r)*dv/d(lon) -> wka
 !    and B = z*v/r - (1/r)*du/d(lon) -> wkb:
 !where z = sin(lat) and r = cos(lat):

 !Convert v1 & om1 (= u_1/r) to semi-spectral space as wka & wkb:
wka=v1
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
wkb=om1
call forfft(numLatGridPts,numLongGridPts,wkb,trig,factors) 

 !Compute longitudinal derivatives of wka & wkb:
call deriv(numLatGridPts,numLongGridPts,rk,wka,wkc)
call deriv(numLatGridPts,numLongGridPts,rk,wkb,wkd)

 !Recover physical fields (inverse FFT):
call revfft(numLatGridPts,numLongGridPts,wkc,trig,factors) 
call revfft(numLatGridPts,numLongGridPts,wkd,trig,factors) 

 !Complete definition of A & B:
do i=1,numLongGridPts
  wka(:,i)=sinLat*om1(:,i)+secantLat*wkc(:,i)
  wkb(:,i)=tanLat*v1(:,i)-wkd(:,i)
enddo
 !tlat = tan(lat) above; *** Do not re-use wka or wkb ***

 !Get physical space velocity divergence -> d1:
d1=deltaSpec1
call revfft(numLatGridPts,numLongGridPts,d1,trig,factors) 
 !*** Do not re-use d1 ***

 !Get physical space derivatives of divergence:
call deriv(numLatGridPts,numLongGridPts,rk,deltaSpec1,wkg)
call revfft(numLatGridPts,numLongGridPts,wkg,trig,factors) ! wkg = d(delta)/d(lon)
call latder(d1,wkc)                 ! wkc = d(delta)/d(lat)

 !Compute the nonlinear part of the velocity divergence tendency:
usq1=u1**2+v1**2 ! *** Do not re-use usq1 ***
sds1=two*(wka*(z1-wka)-wkb*(wkb+d1))-usq1-d1**2-om1*wkg-v1*wkc
 !om1 = u_1/r on rhs above; *** wka & wkb now safe to re-use below ***

 !Add effect of fixed topography (if present) --- layer 1 only:
if (topoFlag) sds1=sds1-specDeltaTopoSrc

 !Convert to semi-spectral space and de-alias:
call dealias(sds1)

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain delta source for layer 2 (sds2):

 !Compute A = z*u/r + (1/r)*dv/d(lon) -> wka
 !    and B = z*v/r - (1/r)*du/d(lon) -> wkb:
 !where z = sin(lat) and r = cos(lat):

 !Convert v2 & om2 (= u_2/r) to semi-spectral space as wka & wkb:
wka=v2
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
wkb=om2
call forfft(numLatGridPts,numLongGridPts,wkb,trig,factors) 

 !Compute longitudinal derivatives of wka & wkb:
call deriv(numLatGridPts,numLongGridPts,rk,wka,wkc)
call deriv(numLatGridPts,numLongGridPts,rk,wkb,wkd)

 !Recover physical fields (inverse FFT):
call revfft(numLatGridPts,numLongGridPts,wkc,trig,factors) 
call revfft(numLatGridPts,numLongGridPts,wkd,trig,factors) 

 !Complete definition of A & B:
do i=1,numLongGridPts
  wka(:,i)=sinLat*om2(:,i)+secantLat*wkc(:,i)
  wkb(:,i)=tanLat*v2(:,i)-wkd(:,i)
enddo
 !tlat = tan(lat) above; *** Do not re-use wka or wkb ***

 !Get physical space velocity divergence -> d2:
d2=deltaSpec2
call revfft(numLatGridPts,numLongGridPts,d2,trig,factors) 
 !*** Do not re-use d2 ***

 !Get physical space derivatives of divergence:
call deriv(numLatGridPts,numLongGridPts,rk,deltaSpec2,wkg)
call revfft(numLatGridPts,numLongGridPts,wkg,trig,factors) ! wkg = d(delta)/d(lon)
call latder(d2,wkc)                 ! wkc = d(delta)/d(lat)

 !Compute the nonlinear part of the velocity divergence tendency:
usq2=u2**2+v2**2 ! *** Do not re-use usq2 ***
sds2=two*(wka*(z2-wka)-wkb*(wkb+d2))-usq2-d2**2-om2*wkg-v2*wkc
 !om2 = u_2/r on rhs above; *** wka & wkb now safe to re-use below ***

 !Convert to semi-spectral space and de-alias:
call dealias(sds2)

!---------------------------------------------------------------
 !gamma source (these are coupled through pressure effects)

 !First obtain c^2*Lap(div(h*(u,v))) in each layer as wkg & wkh;
 !both of these terms appear in each layer's gamma source term.

!::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Compute c_1^2*div((u_1,v_1)*h_1) -> wke:
wka=h1
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
call deriv(numLatGridPts,numLongGridPts,rk,wka,wkd)
call revfft(numLatGridPts,numLongGridPts,wkd,trig,factors) ! wkd = dh_1/d(lon)
call latder(h1,wkc)                 ! wkc = dh_1/d(lat)
wke=d1*h1+om1*wkd+v1*wkc
 !Note om1 = u_1/r on rhs above

 !If thermal damping is present, add its contribution:
if (thermal) wke=wke+thermDampRate*h1

 !Compute Laplacian of wke = div(h_1*(u_1,v_1)) after de-aliasing:
call dealias(wke)
call laplace(wke,wkg)
wkg=csq1*wkg
 !wkg = c_1^2*Lap(div(h_1*(u_1,v_1))) in semi-spectral space

!::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Compute c_2^2*div((u_2,v_2)*h_2) -> wke:
wka=h2
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors) 
call deriv(numLatGridPts,numLongGridPts,rk,wka,wkd)
call revfft(numLatGridPts,numLongGridPts,wkd,trig,factors) ! wkd = dh_2/d(lon)
call latder(h2,wkc)                 ! wkc = dh_2/d(lat)
wke=d2*h2+om2*wkd+v2*wkc
 !Note om2 = u_2/r on rhs above

 !If thermal damping is present, add its contribution:
if (thermal) wke=wke+thermDampRate*h2

 !Compute Laplacian of wke after de-aliasing:
call dealias(wke)
call laplace(wke,wkh)
wkh=csq2*wkh
 !wkh = c_2^2*Lap(div(h_2*(u_2,v_2))) in semi-spectral space

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain gamma source for layer 1 (sgs1):

 !Define wkf = f*zeta_1 & wkd = 2*f*beta*v_1:
do i=1,numLongGridPts
  wkf(:,i)=cof*z1(:,i)            !cof = f
  wkd(:,i)=dfb*v1(:,i)            !dfb = 2*f*beta
enddo

 !Compute div((u_1,v_1)*Z_1) where Z_1 = f*zeta_1:
wka=wkf
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors)
call deriv(numLatGridPts,numLongGridPts,rk,wka,wke)
call revfft(numLatGridPts,numLongGridPts,wke,trig,factors) ! wke = dZ/d(lon)
call latder(wkf,wkc)                ! wkc = dZ/d(lat)

 !Store div((u_1,v_1)*Z_1) + 2*f*beta*v_1 in sgs1 temporarily:
sgs1=d1*wkf+om1*wke+v1*wkc+wkd
 !De-alias and convert to semi-spectral space:
call dealias(sgs1)

 !Define B = c_1^2*h_1 + (rho_2/rho_1)*c_2^2*h_2 + (u_1^2 + v_1^2)/2 -> wkb:
wkb=csq1*h1+densityRatio*csq2*h2+f12*usq1
 !De-alias and convert to semi-spectral space:
call dealias(wkb)
 !Add effect of fixed topography (if present) --- layer 1 only:
if (topoFlag) wkb=wkb+specGammaTopoSrc
 !Calculate dB/d(lon) -> wkf (keep this semi-spectral):
call deriv(numLatGridPts,numLongGridPts,rk,wkb,wkf)

 !Complete calculation of (de-aliased) nonlinear part of gamma tendency:
sgs1=wkg+densityRatio*wkh+fpole*wkf-sgs1
 !wkg & wkh are defined above as the Laplacian of c^2*div(h*(u,v))
 !in layers 1 and 2, respectively.

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 !Obtain gamma source for layer 2 (sgs2):

 !Define wkf = f*zeta_2 & wkd = 2*f*beta*v_2:
do i=1,numLongGridPts
  wkf(:,i)=cof*z2(:,i)            !cof = f
  wkd(:,i)=dfb*v2(:,i)            !dfb = 2*f*beta
enddo

 !Compute div((u_2,v_2)*Z_2) where Z_2 = f*zeta_2:
wka=wkf
call forfft(numLatGridPts,numLongGridPts,wka,trig,factors)
call deriv(numLatGridPts,numLongGridPts,rk,wka,wke)
call revfft(numLatGridPts,numLongGridPts,wke,trig,factors) ! wke = dZ/d(lon)
call latder(wkf,wkc)                ! wkc = dZ/d(lat)

 !Store div((u_2,v_2)*Z_2) + 2*f*beta*v_2 in sgs2 temporarily:
sgs2=d2*wkf+om2*wke+v2*wkc+wkd
 !De-alias and convert to semi-spectral space:
call dealias(sgs2)

 !Define B = c_1^2*h_1 + c_2^2*h_2 + (u_2^2 + v_2^2)/2 -> wkb:
wkb=csq1*h1+csq2*h2+f12*usq2
 !De-alias and convert to semi-spectral space:
call dealias(wkb)
 !Calculate dB/d(lon) -> wkf (keep this semi-spectral):
call deriv(numLatGridPts,numLongGridPts,rk,wkb,wkf)

 !Complete calculation of (de-aliased) nonlinear part of gamma tendency:
sgs2=wkg+wkh+fpole*wkf-sgs2
 !wkg & wkh are defined above as the Laplacian of c^2*div(h*(u,v))
 !in layers 1 and 2, respectively.

!-----------------------------------------------------
 !Remove global mean values of sds & sgs:
rhs=sds1(:,1)*cosLat
avgval=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sds1(:,1)=sds1(:,1)-avgval

rhs=sds2(:,1)*cosLat
avgval=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sds2(:,1)=sds2(:,1)-avgval

rhs=sgs1(:,1)*cosLat
avgval=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sgs1(:,1)=sgs1(:,1)-avgval

rhs=sgs2(:,1)*cosLat
avgval=(f1112*(rhs(1)+rhs(numLatGridPts))+sum(rhs(2:ngm1)))*rsumi
sgs2(:,1)=sgs2(:,1)-avgval

return
end subroutine source

 ! End main program
end program dgbal
!=======================================================================
