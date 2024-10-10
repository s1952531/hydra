subroutine nhpsolve(d1,d2,s1,s2,wka,wkb,iopt)

! Finds the vertically-averaged non-hydrostatic pressure (pn1 & pn2) 
! in each layer:  pn1 = bar{P}_{n1}/h_1  and  pn2 = bar{P}_{n2}/h_2.
! On return, dj = divergence in layer j, and sj = -xi_j where
! x_j = D(delta_j)/Dt - delta_j^2, while (wka,wkb) = (J_1,J_2),
! the non-hydrostatic part of zeta_j tendency (j = 1 & 2).

! *** All fields are in physical space ***
  
implicit none

! Passed variables:
double precision:: d1(ng,ng),d2(ng,ng),s1(ng,ng),s2(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng)
integer:: iopt

! Local variables:
double precision,parameter:: ptol=1.d-7
! ptol: maximum relative rms NH pressure error

double precision:: htot1(ng,ng),htot2(ng,ng),hinv1(ng,ng),hinv2(ng,ng)
double precision:: h1x(ng,ng),h1y(ng,ng),tt(ng,ng),dt1c(ng,ng)
double precision:: r1x(ng,ng),r1y(ng,ng),t1x(ng,ng),t1y(ng,ng)
double precision:: r2x(ng,ng),r2y(ng,ng)
double precision:: ch1(ng,ng),ch2(ng,ng),ch3(ng,ng)
double precision:: ep11(ng,ng),ep12(ng,ng),ep21(ng,ng),ep22(ng,ng)
double precision:: xx11(ng,ng),xx12(ng,ng),xx21(ng,ng),xx22(ng,ng)
double precision:: yy11(ng,ng),yy12(ng,ng),yy21(ng,ng),yy22(ng,ng)
double precision:: pn1x(ng,ng),pn1y(ng,ng)
double precision:: pn2x(ng,ng),pn2y(ng,ng)
double precision:: wkc(ng,ng),wkd(ng,ng),wkp(ng,ng)
double precision:: perr

!---------------------------------------------------------------
! Get total layer thicknesses (h_1 & h_2 in the comments below):
htot1=hbar1*(one+h1)
htot2=hbar2*(one+h2)

! Get their inverses:
hinv1=one/htot1
call dealias(hinv1)
hinv2=one/htot2
call dealias(hinv2)

! Define T = 6/(4 + 3*mu) where mu = alpha*h_2/h_1:
wkp=alpha*htot2*hinv1
call dealias(wkp)
tt=six/(four+three*wkp)
call dealias(tt)

! Define c_hat_1,2,3:
wkp=three*tt*hinv1
call dealias(wkp)
ch1=wkp*hinv1-cona1
call dealias(ch1)
ch2=wkp*hinv2-cona2
call dealias(ch2)
wkp=two*tt*hinv2
call dealias(wkp)
ch3=wkp*(hinv2+alpha3*hinv1)-cona3   !alpha3 = 3*alpha
call dealias(ch3)

! Calculate (h1x,h1y) = grad{h_1} and (r1x,r1y) = grad{h_1}/h_1:
wkp=htot1
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,h1x,h1y)
r1x=hinv1*h1x
call dealias(r1x)
r1y=hinv1*h1y
call dealias(r1y)

! Calculate (t1x,t1y) = T*(r1x,r1y) = T*grad{h_1}/h_1 = tau in the notes:
t1x=tt*r1x
call dealias(t1x)
t1y=tt*r1y
call dealias(t1y)

! Calculate (r2x,r2y) = T*grad{h_2}:
wkp=htot2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,r2x,r2y)
r2x=hinv2*r2x
call dealias(r2x)
r2y=hinv2*r2y
call dealias(r2y)

! Define factors multiplying P_n1 & P_n2 and their derivatives
! needed in the iteration below:
xx11= -f23*t1x
yy11= -f23*t1y              ! (xx11,yy11) = -(2/3)*tau
xx12=alpha*t1x
yy12=alpha*t1y              ! (xx11,yy11) =  alpha*tau
xx21= -f12*t1x
yy21= -f12*t1y              ! (xx21,yy21) = -(1/2)*tau
call divs(t1x,t1y,wkc)
call spctop(ng,ng,wkc,wkd,xfactors,yfactors,xtrig,ytrig)
! Now wkd = div{tau}
wkp=wkd-ch1                 ! wkp  =        div{tau} - c_hat_1
ep11= -f23*wkp              ! ep11 = (2/3)*(c_hat_1 - div{tau})
ep12=alpha*wkp              ! ep12 = alpha*(div{tau} - c_hat_1)
ep21=-f12*wkd-ch2           ! ep21 = -(1/2)*div{tau} - c_hat_2

wkp=tt*hinv2
call dealias(wkp)           ! wkp = T/h_2 (temporarily)
xx22=-r2x-wkp*h1x
yy22=-r2y-wkp*h1y           ! (xx22,yy22) = -(grad{h_2}+T*grad{h_1})/h_2
call dealias(xx22)
call dealias(yy22)
call divs(xx22,yy22,wkc)
call spctop(ng,ng,wkc,wkd,xfactors,yfactors,xtrig,ytrig)
ep22=wkd+ch3                ! ep22 = c_hat_3 - div{(grad{h_2}+T*grad{h_1})/h_2}

!--------------------------------------------------------------------
! Form the fixed parts (s1,s2) of the sources needed to calculate the
! non-hydrostatic pressure:

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute gamma_tilde = gamma_l + 2*(J(u,v) - delta^2) in layer 1:
wka=ds1
call spctop(ng,ng,wka,d1,xfactors,yfactors,xtrig,ytrig)
! d1 contains delta_1 in physical space
call jacob(u1,v1,wkp)
! wkp contains J(u_1,v_1) in physical space
wkp=wkp-d1**2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
! Form eta_2 = f*(zeta_2 - f*h_2) -> wkb
wkb=cof*(z2-cof*h2)
call ptospc(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
wka=gs1+ee2*(gs2-wkc)+two*filt*wka
! Above gs1 + ee2*(gs2-wkc) = gamma_1l (spectral); ee2 = E_2 operator
call spctop(ng,ng,wka,s1,xfactors,yfactors,xtrig,ytrig)
! s1 now contains gamma_tilde_1 in physical space (de-aliased)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Compute gamma_tilde = gamma_l + 2*(J(u,v) - delta^2) in layer 2:
wka=ds2
call spctop(ng,ng,wka,d2,xfactors,yfactors,xtrig,ytrig)
! d2 contains delta_2 in physical space
call jacob(u2,v2,wkp)
! wkp contains J(u_2,v_2) in physical space
wkp=wkp-d2**2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
! Form eta_1 = f*(zeta_1 - f*h_1) -> wkb
wkb=cof*(z1-cof*h1)
call ptospc(ng,ng,wkb,wkc,xfactors,yfactors,xtrig,ytrig)
wka=gs2+ee1*(gs1-wkc)+two*filt*wka
! Above gs2 + ee1*(gs1-wkc) = gamma_2l (spectral); ee1 = E_1 operator
call spctop(ng,ng,wka,s2,xfactors,yfactors,xtrig,ytrig)
! s2 now contains gamma_tilde_2 in physical space (de-aliased)

!======================================================================
! Iterate to find approximate solution:

if (iopt .eq. 0) then
  ! this is done only when pn1 & pn2 are not yet defined:
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
endif
  
perr=one
do while (perr .gt. ptol)

  wkp=pn1
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  call gradient(wkb,pn1x,pn1y)

  wkp=pn2
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  call gradient(wkb,pn2x,pn2y)

  wkp=s1+ep11*pn1+ep12*pn2+xx11*pn1x+xx12*pn2x+yy11*pn1y+yy12*pn2y
  call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
  ! wka = full S_1 in spectral space

  wkp=s2+ep21*pn1+ep22*pn2+xx21*pn1x+xx22*pn2x+yy21*pn1y+yy22*pn2y
  call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
  ! wkb = full S_2 in spectral space

  ! Obtain next approximation (hold in wkc & wkd temporarily):
  wkp=pi11*wka+pi12*wkb
  call spctop(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)

  wkp=pi21*wka+pi22*wkb
  call spctop(ng,ng,wkp,wkd,xfactors,yfactors,xtrig,ytrig)

  ! Compute relative rms difference error:
  perr=sqrt(sum((wkc-pn1)**2+(wkd-pn2)**2)/sum(pn1**2+pn2**2))

  ! Copy wkc & wkd into pn1 & pn2:
  pn1=wkc
  pn2=wkd

enddo

!----------------------------------------------------------
! Compute -xi_1 = -[D(delta_1)/Dt - delta_1^2] -> s1
!     and -xi_2 = -[D(delta_2)/Dt - delta_2^2] -> s2:
wkp=f23*pn1-alpha*pn2
s1=(ch1+cona1)*wkp
s2=(ch3+cona3)*pn2-(ch2+cona2)*pn1

! Compute non-hydrostatic parts of zeta tendencies (J_1 & J_2):
wkp=tt*wkp
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call gradient(wkb,wkc,wkd)
wka=r1x*wkd-r1y*wkc                     ! wka = J_1

wkp=pn2
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
call gradient(wkb,pn2x,pn2y)

wkp=f12*pn1*hinv1+pn2*hinv2
call dealias(wkp)
wkp=tt*wkp                              ! wkp = kappa in the notes
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call gradient(wkb,wkc,wkd)              ! (wkc,wkd) = grad{kappa}
wkb=r2x*pn2y-r2y*pn2x+h1x*wkd-h1y*wkc   ! wkb = J_2

return
end subroutine nhpsolve
