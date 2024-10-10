program shape
!  -------------------------------------------------------------------------
!  |   Computes the shape parameter W = F/G^{3/2} defined by Wanming Qi,   |
!  |   together with the unsteadiness, defined as sqrt{<q_t^2>}/<q^2>,     |
!  |   and the energy from data in zz.r4.  Produces shape.asc.             |
!  |                                                                       |
!  |   For the shape parameter, the vorticity is projected on the          |
!  |   spherical harmonics of order 2:                                     |
!  |                                                                       |
!  |          w00*Y_20 + w21*Y_21 + c.c. + w22*Y_22 + c.c.  =              |
!  |                                                                       |
!  |   a0*Y_20 + a1r*Re(Y_21) - a1i*Im(Y_21) + a2r*Re(Y_22) - a2i*Im(Y_22) |
!  |                                                                       |
!  |   where a0 = w00, a1r = 2*Re(w21), a1i = 2*Im(w21), a2r = 2*Re(w22)   |
!  |   and a2i = 2*Im(w22).                                                |
!  |                                                                       |
!  |   Output is to the formatted file "shape.asc" which lists time vs     |
!  |   energy, log_10(sqrt{<q_t^2>}/<q^2>) & W.                            |
!  -------------------------------------------------------------------------
use spectral

implicit none

real:: tr4,zzr4(ng,nt)
double precision:: qq(ng,nt),zz(ng,nt),uu(ng,nt),vv(ng,nt),pp(ng,nt)
double precision:: wka(ng,nt),wkb(ng,nt),wkc(ng,nt)
double precision:: clon(nt),slon(nt)
double precision:: y20(ng)
double precision:: y21c(ng,nt),y21s(ng,nt)
double precision:: y22c(ng,nt),y22s(ng,nt)
double precision:: c20,c21,c22,w1,w2
double precision:: a0,a1r,a1i,a2r,a2i
double precision:: a1sq,a2sq,fpar,ff,gg,ww
double precision:: rlon,rlat,x,y
double precision:: ene,ens,uns
integer:: i,j,loop,iread

fpar=3.d0*sqrt(6.d0)/4.d0

!---------------------------------------------------------------
 !Define various fixed arrays:
do i=1,nt
  rlon=dble(i-1)*dl-pi
  clon(i)=cos(rlon)
  slon(i)=sin(rlon)
enddo

do j=1,ng
  rlat=(dble(j)-f12)*dl-hpi
  slat(j)=sin(rlat)
  clat(j)=cos(rlat)
enddo

 !Spherical harmonics of order 2 (real forms):
c20=f14*sqrt(five/pi)
do j=1,ng
  y20(j)=c20*(three*slat(j)**2-one)
enddo
c21=sqrt(15.d0/(8.d0*pi))
c22=f12*c21
do i=1,nt
  do j=1,ng
    x=clat(j)*clon(i)
    y=clat(j)*slon(i)
    y21c(j,i)=-c21*x*slat(j)
    y21s(j,i)=-c21*y*slat(j)
    y22c(j,i)= c22*(x**2-y**2)
    y22s(j,i)= c21*x*y
  enddo
enddo

 !Normalisation coefficients used below in projection:
c20=f1112*(y20(1)**2*clat(1)+y20(ng)**2*clat(ng))
do j=2,ngm1
  c20=c20+y20(j)**2*clat(j)
enddo
write(*,'(a,f11.7)') ' <Y_20>  = ',dl*c20*twopi
c20=one/(c20*dble(nt))

c21=zero
c22=zero
do i=1,nt
  c21=c21+f1112*(y21c(1,i)**2*clat(1)+y21c(ng,i)**2*clat(ng))
  c22=c22+f1112*(y22c(1,i)**2*clat(1)+y22c(ng,i)**2*clat(ng))
  do j=2,ngm1
    c21=c21+y21c(j,i)**2*clat(j)
    c22=c22+y22c(j,i)**2*clat(j)
  enddo
enddo
write(*,'(a,f11.7)') ' <Y_21c> = ',dl**2*c21
write(*,'(a,f11.7)') ' <Y_22c> = ',dl**2*c22
c21=one/c21
c22=one/c22

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

!---------------------------------------------------------------
 !Open input data file:
open(44,file='zz.r4',form='unformatted', & 
    & access='direct',status='old',recl=nbytes)

 !Open output file:
open(22,file='shape.asc',status='replace')

!---------------------------------------------------------------
 !Read data and process:
loop=0
do  
  loop=loop+1
  iread=0
  read(44,rec=loop,iostat=iread) tr4,zzr4
  if (iread .ne. 0) exit 

   !Convert to double precision:
  qq=dble(zzr4)

   !Get qq^2*cos(lat) -> wkc to compute <q^2>:
  do i=1,nt
    do j=1,ng
      zz(j,i)=qq(j,i)
      wkc(j,i)=clat(j)*qq(j,i)**2
    enddo
  enddo

   !Project zz onto l = 2 spherical harmonics:
  a0=zero
  a1r=zero
  a1i=zero
  a2r=zero
  a2i=zero
  do i=1,nt
    w1=zz(1,i)*clat(1)
    w2=zz(ng,i)*clat(ng)
    a0=a0+f1112*(y20(1)*w1+y20(ng)*w2)
    a1r=a1r+f1112*(y21c(1,i)*w1+y21c(ng,i)*w2)
    a1i=a1i-f1112*(y21s(1,i)*w1+y21s(ng,i)*w2)
    a2r=a2r+f1112*(y22c(1,i)*w1+y22c(ng,i)*w2)
    a2i=a2i-f1112*(y22s(1,i)*w1+y22s(ng,i)*w2)
    do j=2,ngm1
      w1=zz(j,i)*clat(j)
      a0=a0+y20(j)*w1
      a1r=a1r+y21c(j,i)*w1
      a1i=a1i-y21s(j,i)*w1
      a2r=a2r+y22c(j,i)*w1
      a2i=a2i-y22s(j,i)*w1
    enddo
  enddo
  a0=a0*c20
  a1r=a1r*c21
  a1i=a1i*c21
  a2r=a2r*c22
  a2i=a2i*c22
   !Compute shape parameter:
  a1sq=a1r**2+a1i**2
  a2sq=a2r**2+a2i**2
  gg=a0**2+f12*(a1sq+a2sq)
  ff=a0*(three*(a2sq-f12*a1sq)-two*a0**2) &
    &    -fpar*(a2r*(a1r**2-a1i**2)+two*a2i*a1r*a1i)
  ww=ff/(gg*sqrt(gg))

  ens=zero
  do i=1,nt
    ens=ens+f1112*(wkc(1,i)+wkc(ng,i))
    do j=2,ngm1
      ens=ens+wkc(j,i)
    enddo
  enddo
  ens=ens*dsumi

   !De-aliase and FFT vorticity in longitude (semi-spectral):
  call dealiase(qq)

   !Invert Laplace's operator on the vorticity (qq):
  call laplinv(qq,pp,uu)
   !Here the streamfunction psi is pp while uu = d(psi)/dlat.

   !Compute d(psi)/dlon = vv:
  call deriv(ng,nt,rk,pp,vv)

   !Get physical space velocity:
  call revfft(ng,nt,uu,trig,factors)
  call revfft(ng,nt,vv,trig,factors)  
  call revfft(ng,nt,pp,trig,factors)  
   !Here uu is minus the zonal velocity while vv is cos(lat) times
   !the meridional velocity

   !Calculate longitudinal derivative of vorticity spectrally:
  call deriv(ng,nt,rk,qq,wka)
   !Return to physical space:
  call revfft(ng,nt,wka,trig,factors) 
   !Calculate latitudinal derivative spectrally:
  call latder(qq,wkb)
   !wkb is returned in physical space here

   !Get zz*pp*cos(lat) -> wkc to compute energy:
  do i=1,nt
    do j=1,ng
      wkc(j,i)=clat(j)*zz(j,i)*pp(j,i)
    enddo
  enddo
  ene=zero
  do i=1,nt
    ene=ene+f1112*(wkc(1,i)+wkc(ng,i))
    do j=2,ngm1
      ene=ene+wkc(j,i)
    enddo
  enddo
  ene=-ene*dsumi

   !Get qq_t*cos(lat) -> wkc to compute <q_t^2>:
  do i=1,nt
    do j=1,ng
      wkc(j,i)=clati(j)*(uu(j,i)*wka(j,i)-vv(j,i)*wkb(j,i))**2
    enddo
  enddo

  uns=zero
  do i=1,nt
    uns=uns+f1112*(wkc(1,i)+wkc(ng,i))
    do j=2,ngm1
      uns=uns+wkc(j,i)
    enddo
  enddo
  uns=uns*dsumi

  uns=sqrt(uns)/ens
  write(22,'(f12.5,1x,f13.11,2(1x,f14.9))') tr4,ene,log10(uns),ww
  write( *,'(a,f12.5)') ' Processed t = ',tr4

enddo

 !Close all files:
close(44)
close(22)

write(*,*)
write(*,*) ' t vs energy, log_10(sqrt{<q_t^2>}/<q^2>) & W is listed in shape.asc;'
write(*,*) ' to view the results type'
write(*,*)
write(*,*) ' plotcol shape.asc'
write(*,*)

end program
