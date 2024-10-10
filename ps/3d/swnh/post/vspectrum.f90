!#########################################################################
!  Computes the horizontally-averaged vertical spectrum of zeta, delta
!  and rho at a specified time.  Write zdr_vspecnnnn.asc where nnnn
!  is 0000, 0001 ... - the corresponding time frame.

!           Written 8/7/2019 by D G Dritschel @ St Andrews
!#########################################################################

program vspectrum

 !Import 1D FFT module:
use stafft

 !Import constants:
use constants

implicit none
integer,parameter:: ngsq=ng*ng

 !Original fields to process:
double precision:: z(ngsq,0:nz),d(ngsq,0:nz),r(ngsq,0:nz)

 !Other local variables:
double precision:: trig(2*nz)
integer:: factors(5)

double precision:: zspec(0:nz),dspec(0:nz),rspec(0:nz),fac
real:: tr4,qr4(ngsq,0:nz)
integer:: loop,i,m
character(len=17):: outfile

outfile='zdr_vspec0000.asc'

write(*,*) ' Enter the time you wish to analyse:'
read(*,*) tr4
loop=nint(tr4/tgsave)+1
write(outfile(10:13),'(i4.4)') loop-1

!----------------------------------------------------------------------
! Set up FFTs:
call initfft(nz,factors,trig)

 !Open input data files and read data:
open(31,file='3d/ql.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
read(31,rec=loop) tr4,qr4
z=dble(qr4)
close(31)

open(31,file='3d/r.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
read(31,rec=loop) tr4,qr4
r=dble(qr4)
close(31)
 !Define zeta from q_l (in z) and r:
z=z+cof*r

open(31,file= '3d/d.r4',form='unformatted',access='direct', &
                      status='old',recl=ntbytes)
read(31,rec=loop) tr4,qr4
d=dble(qr4)
close(31)

 !Perform cosine Fourier transforms (overwrites data by the coeffs):
call dct(ngsq,nz,z,trig,factors)
call dct(ngsq,nz,d,trig,factors)
call dct(ngsq,nz,r,trig,factors)

zspec=zero
dspec=zero
rspec=zero
do m=0,nz
  do i=0,ngsq
    zspec(m)=zspec(m)+z(i,m)**2
    dspec(m)=dspec(m)+d(i,m)**2
    rspec(m)=rspec(m)+r(i,m)**2
  enddo
enddo
fac=one/dble(ngsq)
zspec=fac*zspec
dspec=fac*dspec
rspec=fac*rspec
zspec(0)=f12*zspec(0)
dspec(0)=f12*dspec(0)
rspec(0)=f12*rspec(0)

open(60,file=outfile,status='replace')
do m=0,nz
  write(60,'(4(1x,f12.8))') log10(dble(m+1)),log10(zspec(m)+1.d-32), &
             log10(dspec(m)+1.d-32),log10(rspec(m)+1.d-32)
enddo
close(60)

write(*,*)
write(*,*) ' Spectra of zeta, delta & rho_tilde are now in '//outfile

 !End main program
end program vspectrum
!=======================================================================
