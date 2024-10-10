subroutine nhpsolve(d1,d2,s1,s2,wka,wkb,iopt)

! Finds the scaled vertically-averaged non-hydrostatic pressure
! pnj = H_j^{-2}*bar{P}_{nj}/h_j, in layers j = 1 & 2.
! On return, dj = delta in layer j, and if ggen is true, 
! sj = -xi_j = -[D(delta_j)/Dt - delta_j^2], while 
! (wka,wkb) = (J_1,J_2) = non-hydrostatic part of zeta_j tendency.

! iopt = 0 on initialisation when a previous estimate for pnj
! is unavailable.  

! *** All passed variables are in physical space ***
  
implicit none

! Passed variables:
double precision:: d1(ng,ng),d2(ng,ng),s1(ng,ng),s2(ng,ng)
double precision:: wka(ng,ng),wkb(ng,ng)
integer:: iopt

! Local variables:
double precision,parameter:: ptol=1.d-7
! ptol: maximum relative rms NH pressure error

double precision:: htot1(ng,ng),htot2(ng,ng),hinv1(ng,ng),hinv2(ng,ng)
double precision:: h1x(ng,ng),h1y(ng,ng),tt(ng,ng)
double precision:: r1x(ng,ng),r1y(ng,ng),t1x(ng,ng),t1y(ng,ng)
double precision:: r2x(ng,ng),r2y(ng,ng),sgx(ng,ng),sgy(ng,ng)
double precision:: ht1x(ng,ng),ht1y(ng,ng)
double precision:: ch1(ng,ng),ch2(ng,ng),ch3(ng,ng)
double precision:: wkc(ng,ng),wkd(ng,ng),wkp(ng,ng)
double precision:: perr

!---------------------------------------------------------------
! Get total dimensionless layer thicknesses (1 + h_j):
htot1=one+h1
htot2=one+h2

! Get their inverses:
hinv1=one/htot1
call dealias(hinv1)
hinv2=one/htot2
call dealias(hinv2)

! Define T = 6/(4 + 3*mu) where mu = alpha*(H_2/H_1)*(1+h_2)/(1+h_1):
wkp=mubar*htot2*hinv1
call dealias(wkp)
tt=six/(four+three*wkp)
call dealias(tt)

! Define c_hat_1,2,3:
wkp=three*tt*hinv1
call dealias(wkp)
ch1=wkp*hinv1-cona1
call dealias(ch1)
ch2=hrati*wkp*hinv2-cona2
call dealias(ch2)
wkd=tt*hinv2          !wkd = T/(1+h_2): *** preserve to define sgx & sgy
call dealias(wkd)
ch3=two*wkd*(hinv2+mubar3*hinv1)-cona3   !mubar3 = 3*mubar
call dealias(ch3)

! Calculate (h1x,h1y) = grad{h_1} and (r1x,r1y) = H_1^2*grad{h_1}/(1+h_1):
wkp=h1
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,h1x,h1y)
r1x=hbsq1*hinv1*h1x
call dealias(r1x)
r1y=hbsq1*hinv1*h1y
call dealias(r1y)

! Calculate (t1x,t1y) = tau = T*(r1x,r1y) = H_1^2*T*grad{h_1}/(1+h_1):
t1x=tt*r1x
call dealias(t1x)
t1y=tt*r1y
call dealias(t1y)

! Define tau/2 for efficiency in pressure iteration below:
ht1x=f12*t1x
ht1y=f12*t1y

! Calculate (r2x,r2y) = H_2^2*grad{h_2}/(1+h_2):
wkp=hbsq2*h2
call ptospc(ng,ng,wkp,wka,xfactors,yfactors,xtrig,ytrig)
call gradient(wka,r2x,r2y)
r2x=hinv2*r2x
call dealias(r2x)
r2y=hinv2*r2y
call dealias(r2y)

! Calculate (sgx,sgy) = (r2x,r2y) + H_1*H_2*T*grad{h_1}/(1+h_2):
wka=wkd*h1x
call dealias(wka)
sgx=r2x+hb1hb2*wka
wka=wkd*h1y
call dealias(wka)
sgy=r2y+hb1hb2*wka
! Note: here wkd contains T/(1+h_2) from above; wkd can now be re-used

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
  ! This is done at t = 0 only when pn1 & pn2 are not yet defined:
  wkc=htot1*s1
  call dealias(wkc)
  wkd=htot2*s2
  call dealias(wkd)

  ! Use the exact solution valid in the long-wave limit (k*H_j)^2 -> 0:
  pn2=-htot2*(hhrati*wkc+f13*wkd)
  pn1=1.5d0*ahrsq*pn2-(f13*htot1+f14*mubar*htot2)*wkc
endif
  
perr=one
do while (perr .gt. ptol)

  ! Find full S_1 (rhs of 1st pressure equation) in spectral space (wka):
  wkp=f23*pn1-ahrsq*pn2     !P = (2/3)*pn1 - alpha*(H_2/H_1)^2*pn2
  wkc=wkp*t1x
  wkd=wkp*t1y
  call divs(wkc,wkd,wka)    !wka = div{P*tau} (spectral)
  wkp=s1+ch1*wkp            !wkp = gamma_tilde_1 + c_hat_1*P
  call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
  wka=wkc-wka               !wka = gamma_tilde_1 + c_hat_1*P - div{P*tau}

  ! Find full S_2 (rhs of 2nd pressure equation) in spectral space (wkb):
  wkc=pn1*ht1x+pn2*sgx
  wkd=pn1*ht1y+pn2*sgy
  call divs(wkc,wkd,wkb)    !wkb = div{(1/2)*pn1*tau+pn2*sigma} (spectral)
  wkp=s2-ch2*pn1+ch3*pn2    !wkp = gamma_tilde_2 - c_hat_2*pn1 + c_hat_3*pn2
  call ptospc(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)
  wkb=wkc-wkb               !wkb = gamma_tilde_2 - c_hat_2*pn1 + c_hat_3*pn2
                            !     -div{(1/2)*pn1*tau+pn2*sigma}

  ! Obtain next approximation for pn1 & pn2 (hold in wkc & wkd temporarily):
  wkp=pi11*wka+pi12*wkb
  call spctop(ng,ng,wkp,wkc,xfactors,yfactors,xtrig,ytrig)

  wkp=pi21*wka+pi22*wkb
  call spctop(ng,ng,wkp,wkd,xfactors,yfactors,xtrig,ytrig)

  ! Compute relative rms difference error:
  perr=sqrt(sum((wkc-pn1)**2+(wkd-pn2)**2)/sum(pn1**2+pn2**2))

  ! Copy wkc & wkd into pn1 & pn2 with relaxation (essential!):
  pn1=f12*(pn1+wkc)
  pn2=f12*(pn2+wkd)
  ! From tests, this 50/50 blend of old and new solutions leads to
  ! exponential convergence (50% error reduction each iteration).
enddo

! If this was called from savegrid, only pn1 & pn2 are needed:
if (.not. ggen) return

!----------------------------------------------------------
! Compute -xi_1 = -[D(delta_1)/Dt - delta_1^2] -> s1
!     and -xi_2 = -[D(delta_2)/Dt - delta_2^2] -> s2:
wkp=f23*pn1-ahrsq*pn2
s1=(ch1+cona1)*wkp
s2=(ch3+cona3)*pn2-(ch2+cona2)*pn1

! Compute non-hydrostatic parts of zeta_1 tendency (-> wka):
wkp=tt*wkp
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call gradient(wkb,wkc,wkd)
wka=r1x*wkd-r1y*wkc                     ! wka = J_1

! Compute non-hydrostatic parts of zeta_2 tendency (-> wkb):
wkp=pn2
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
call gradient(wkb,sgx,sgy)

wkp=pn2*hinv2+hhrati*pn1*hinv1          ! pn2/(1+h_2)+(H_1/2*H_2)*pn1/(1+h_1)
call dealias(wkp)
wkp=tt*wkp                              ! wkp = kappa in the notes
call ptospc(ng,ng,wkp,wkb,xfactors,yfactors,xtrig,ytrig)
wkb=filt*wkb
call gradient(wkb,wkc,wkd)              ! (wkc,wkd) = grad{kappa}
wkb=r2x*sgy-r2y*sgx+hb1hb2*(h1x*wkd-h1y*wkc)   ! wkb = J_2

!------------------------------------------------------------------
wkp=pn1*(one+h1)

return
end subroutine
