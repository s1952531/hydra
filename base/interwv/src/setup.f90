module wavesetup

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! This module contains a routine to find the phase speeds c which 
! satisify the eigenvalue problem
! d^2(psi)/dy^2 + (N^2/c^2)psi = 0  ; psi(0) = psi(L_y) = 0
! for a three-layer configuration of piecewise-uniform N^2 
! bordered by error functions.

! This needs the lapack or an equivalent to carry out the matrix operations

! Adapted from older f77 source by S King, Sept 3 @ St Andrews
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

use constants

implicit double precision (a-h,o-z)

double precision:: bbarf(0:nyf),bfsf(0:nyf),dbfsfdy(0:nyf)
double precision:: psi(ny,nx)
double precision:: cp,cplong

contains

subroutine wave_init

implicit double precision (a-h,o-z)

 !Parameter declarations:
integer,parameter:: lwork=8*nym1
double precision,parameter:: bf1sq=bf1**2,bf2sq=bf2**2,bf3sq=bf3**2
double precision,parameter:: afac=f12*(bf1sq+bf3sq),bfac=f12*(bf2sq-bf1sq),cfac=f12*(bf3sq-bf2sq)
double precision,parameter:: wid=bfsmoond*gly

 !Array declarations:
double precision:: pbar(0:nyf),dpbardy(0:nyf)
double precision:: psik(0:ny),dpsikdy(0:ny)
double precision:: gam(nym1)
double precision:: smat(nym1,nym1)
double precision:: csq(nym1),dum(nym1)
double precision:: vecr(nym1,nym1),vecl(nym1,nym1)
double precision:: dgwork(lwork)
double precision:: lambda

!-------------------------------------------------------------------
 !Compute mean buoyancy and Bernoulli pressure:
cerf=two/(sqrt(pi)*wid)
do iy=0,nyf
  yg=glyf*dble(iy)
  arg1=(yg-y1)/wid
  arg2=(yg-y2)/wid
  bfsf(iy)=afac+bfac*erf(arg1)+cfac*erf(arg2)
  dbfsfdy(iy)=cerf*(bfac*exp(-arg1**2)+cfac*exp(-arg2**2))
  dpbardy(iy)=-yg*bfsf(iy)
enddo

bbarf(0)=zero
pbar(0)=zero
do iy=1,nyf
  bbarf(iy)=bbarf(iy-1)+hglyf*(bfsf(iy)+bfsf(iy-1))
  pbar(iy)=pbar(iy-1)+hglyf*(dpbardy(iy-1)+dpbardy(iy))
enddo
 !pbar and its derivative are used to calculate ape in any unsteady code run from 
 !here - they are not currently used but the calculation is retained in case needed
 !later - see old versions of clamb.F for how the calculation should be done.

!-------------------------------------------------------------------
 !Compute eigenmatrix (smat):
do m=1,nym1
  gam(m)=pi*dble(m)/elly
enddo

do m=1,nym1
  do n=1,nym1
    smat(n,m)=zero
  enddo
enddo

do iy=1,nym1
  yg=gly*dble(iy)
  iyf=mgf*iy
  do m=1,nym1
    do n=1,nym1    
      smat(n,m)=smat(n,m)+bfsf(iyf)*sin(gam(n)*yg)*sin(gam(m)*yg)
    enddo
  enddo
enddo

scfac=gly/hly
do m=1,nym1
  do n=1,nym1
    smat(n,m)=scfac*smat(n,m)/gam(m)**2
  enddo
enddo     

!-------------------------------------------------------------------
 !Solve eigenproblem:
ifail=0
call dgeev('N','V',nym1,smat,nym1,csq,dum,vecl,nym1,vecr,nym1,dgwork,lwork,ifail)

 !Select mode with largest csq:
csqmax=zero
do j=1,nym1
  if (csq(j) .gt. csqmax) then
    m=j
    csqmax=csq(j)
  endif
enddo

do iy=0,ny
  psik(iy)=zero
  dpsikdy(iy)=zero
enddo

do iy=0,ny
  yg=gly*dble(iy)
  do n=1,nym1
    psik(iy)=psik(iy)+vecr(n,m)*sin(gam(n)*yg)/gam(n)**2
    dpsikdy(iy)=dpsikdy(iy)+vecr(n,m)*cos(gam(n)*yg)/gam(n)
  enddo
enddo

 !Set psikmax for rescaling later:
psikmax=zero
do iy=1,ny
  if (abs(psik(iy)) .gt. abs(psikmax)) then
    psikmax=psik(iy)
  endif
enddo
 !Rescale psik(), and dpsikdy()	
do iy=0,ny
  psik(iy)=psik(iy)/psikmax
  dpsikdy(iy)=dpsikdy(iy)/psikmax
enddo
 
 !Write psik(y) etc to eigen.asc:
open(12,file='eigen.asc',status='unknown')
do iy=0,ny
  yg=gly*dble(iy)
  write(12,'(3(1x,f18.12))') yg,psik(iy),dpsikdy(iy)
enddo
close(12)

cp=sqrt(csqmax)
cplong=sqrt(csqmax)

write(*,'(a,2(1x,f16.12))') ' c, c_long = ',cp,cplong

!========================================================================= 
!----------------------------------------------------------------------
! Below are all the setup steps for a solitary wave based on the above
! vertical eigenfunction/value.
!----------------------------------------------------------------------

 !Set temp1=integral(dpsikdy**2),temp2=integral(dpsikdy**3)
temp1=f12*gly*(dpsikdy(0)**2.0d0+dpsikdy(ny)**2.0d0)
do iy=1,nym1
   temp1=temp1+gly*dpsikdy(iy)**2.0d0
enddo
temp2=f12*gly*(dpsikdy(0)**3.0d0+dpsikdy(ny)**3.0d0)
do iy=1,nym1
   temp2=temp2+gly*dpsikdy(iy)**3.0d0
enddo

 !Set alpha0:
alpha0=(3.0d0*cp*temp2)/(2.0d0*temp1)
 !Reset temp2=integral(psik**2)
temp2=f12*gly*((psik(0))**2.0d0+(psik(ny))**2.0d0)
do iy=1,nym1
   temp2=temp2+gly*(psik(iy))**2.0d0
enddo

 !Set beta:
beta=(cp*temp2)/(2.0d0*temp1)

 !Set lambda=sqrt(12.0d0*beta/(alpha0*amp))=hlx/6
lambda=f16*hlx
 !Set the amplitude to make this the correct choice for lambda:	
amp=12.0d0*beta/(alpha0*lambda**2)
 !Set the solitary wave speed and its inverse:
cp=cp+(amp*alpha0/3.0d0)

 !Then set the full psi values on the iy,ix grid:
do ix=1,nx
 xg=-hlx+dble(ix-1)*glx
   do iy=1,ny-1
     yg=dble(iy)*gly
     psi(iy,ix)=-cp*amp*((cosh((xg)/lambda))**(-2.0d0))*psik(iy)
   enddo
enddo      

return
end subroutine

end module
