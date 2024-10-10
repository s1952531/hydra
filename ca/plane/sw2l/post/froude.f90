!#########################################################################
!     Calculates the Froude number for all data in h1.r4 ... z2.r4.

!         Written 10 December 2020 by D G Dritschel @ St Andrews
!#########################################################################

program froude

 !Import constants and parameters:
use constants
 !Import spectral module:
use spectral

implicit none
double precision:: h1(ng,ng),d1(ng,ng),u1(ng,ng),v1(ng,ng),z1(ng,ng)
double precision:: h2(ng,ng),d2(ng,ng),u2(ng,ng),v2(ng,ng),z2(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng),wkc(ng,ng),wkd(ng,ng)
double precision:: uds(ng,ng),vds(ng,ng)
double precision:: uio,vio,fr
real:: qr4(ng,ng),tr4
integer:: loop,iread

!---------------------------------------------------------
 !Initialise inversion constants and arrays:
call init_spectral

 !Open input data files:
open(31,file='evolution/h1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(41,file='evolution/h2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

open(32,file='evolution/d1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(42,file='evolution/d2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

open(33,file='evolution/z1.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)
open(43,file='evolution/z2.r4',form='unformatted', &
      access='direct',status='old',recl=nbytes)

 !Open output file:
open(51,file='evolution/froude.asc',status='replace')

 !Read data at all times and process:
loop=0
do
  loop=loop+1
  iread=0
  read(31,rec=loop,iostat=iread) tr4,qr4
  if (iread .ne. 0) exit
  h1=dble(qr4)
  read(41,rec=loop,iostat=iread) tr4,qr4
  h2=dble(qr4)

  read(32,rec=loop,iostat=iread) tr4,qr4
  d1=dble(qr4)
  read(42,rec=loop,iostat=iread) tr4,qr4
  d2=dble(qr4)

  read(33,rec=loop,iostat=iread) tr4,qr4
  z1=dble(qr4)
  read(43,rec=loop,iostat=iread) tr4,qr4
  z2=dble(qr4)

   !Invert to find velocity field:
  call ptospc(ng,ng,d1,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkc
  call xderiv(ng,ng,hrkx,wka,uds)
  call yderiv(ng,ng,hrky,wka,vds)

  call ptospc(ng,ng,z1,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkc
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)

  wkb=uds-wkb
  wkd=vds+wkd
  call spctop(ng,ng,wkb,u1,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,v1,xfactors,yfactors,xtrig,ytrig)

  call ptospc(ng,ng,d2,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkc
  call xderiv(ng,ng,hrkx,wka,uds)
  call yderiv(ng,ng,hrky,wka,vds)

  call ptospc(ng,ng,z2,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=rlap*wkc
  call xderiv(ng,ng,hrkx,wka,wkd)
  call yderiv(ng,ng,hrky,wka,wkb)

  wkb=uds-wkb
  wkd=vds+wkd
  call spctop(ng,ng,wkb,u2,xfactors,yfactors,xtrig,ytrig)
  call spctop(ng,ng,wkd,v2,xfactors,yfactors,xtrig,ytrig)

  uio=cio*sum(h1*u1+mubar*h2*u2)
  vio=cio*sum(h1*v1+mubar*h2*v2)
  u1=u1+uio
  u2=u2+uio
  v1=v1+vio
  v2=v2+vio

   !Calculate froude number:
  wka=((u2-u1)**2+(v2-v1)**2)/(hbar1*(one+h1)+hbar2*(one+h2))
  fr=sqrt(maxval(wka)/((one-alpha)*gravity))

  write(*,'(a,f9.2)') ' *** Processing t = ',tr4

   !Write data:
  write(51,'(1x,f9.2,1x,f13.9)') tr4,fr
enddo

close(31)
close(41)
close(32)
close(42)
close(33)
close(43)
close(51)

write(*,*)
write(*,*) ' Time vs Froude number is in froude.asc'
write(*,*)

end program froude
