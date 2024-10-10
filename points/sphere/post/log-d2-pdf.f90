program d2pdf
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! Computes the PDF of squared separation distances between ALL 
! pairs of oppositely signed vortices, previously evolved with 
! sep-pvs.f90.  

! Uses logarithmically-spaced bins (in distance)

! It is assumed that positive vortices have odd indices (1,3,5,...)
! and negative vortices have even indices.
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 ! Import parameters and constants:
use parameters
use constants

implicit double precision (a-h,o-z)

integer,parameter:: nbins=50
 ! nbins: number of bins to form the PDF
double precision,parameter:: pi=3.1415926535897932385d0

double precision:: pdf(nbins)
double precision:: x(n), y(n), z(n)

character:: frame1*4,frame2*4

 !--------------------------------------------------------
 ! Select frame to process:
write(*,*) ' Beginning and ending time frames to examine (use 0 for t = 0)?'
read(*,*) kfr1,kfr2

write(*,*) ' Minimum distance, as a fraction of the mean intervortex distance?'
read(*,*) dminnd
davg=sqrt(4.d0*pi/dble(n))
dmin=dminnd*davg

write(*,*) ' Enter alpha, the initial dipole separation divided by L:'
read(*,*) alp
dini=alp*davg
write(*,'(a,f12.8)') ' log_10(initial dipole separation) = ',log10(dini)

 ! Initialise:
dbin=(log10(two)-log10(dmin))/dble(nbins-1)
do m=1,nbins
  pdf(m)=zero
enddo
aldmin=log10(dmin)

 ! Open data file, read and process:
open(11,file='points.dat',status='old')

if (kfr1 .gt. 0) then
   ! skip data:
  do k=1,(n+1)*kfr1
    read(11,*)
  enddo
endif

 ! Read selected frames:
do kfr=kfr1,kfr2
  read(11,*) t
  do k=1,n
    read(11,*) x(k),y(k),z(k)
  enddo

  write(*,'(a,f12.5)') ' *** Processing t = ',t

   ! Accumulate PDF:
  do j=1,n-1,2
    do i=2,n,2
      ald=half*log10(two*(one-x(i)*x(j)-y(i)*y(j)-z(i)*z(j)))-aldmin
      m=max(1,int(ald/dbin)+2)
      pdf(m)=pdf(m)+one
    enddo
  enddo
enddo

 ! Normalise PDF to give a probability density:
sumpdf=zero
do m=1,nbins
  sumpdf=sumpdf+pdf(m)
enddo
fac=one/(dmin*sumpdf)
pdf(1)=fac*pdf(1)
do m=2,nbins
  pdf(m)=fac*pdf(m)/(10.d0**(dbin*dble(m))-10.d0**(dbin*dble(m-1)))
enddo

 ! Write data:
write(frame1(1:4),'(i4.4)') kfr1
write(frame2(1:4),'(i4.4)') kfr2
open(80,file='lpdf'//frame1//'-'//frame2,status='unknown')
do m=1,nbins
  write(80,'(1x,f12.8,1x,f13.9)') aldmin+dbin*(dble(m)-half),pdf(m)
enddo

write(*,*) ' The PDF of log(d) is ready in lpdf'//frame1//'-'//frame2

 !End main program
end program
!=======================================================================
