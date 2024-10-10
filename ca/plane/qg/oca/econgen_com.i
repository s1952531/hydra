      implicit real*8(a-h,o-z)

      include 'econgen_dim.i'
c      nz:   number of layers in the z direction
c      nh:   number of grid boxes in the x & y directions (inversion  grid)
c      nhf:  number of grid boxes in the x & y directions (ultra-fine grid)
c      nplm: maximum number of  nodes  in any given layer
c      nlm:  maximum number of contours "  "    "     "

      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0)
      parameter (half=one/two,third=one/three)
      parameter (four=4.d0,six=6.d0,quart=one/four,sixth=one/six)
      parameter (pi=3.1415926535897932385d0,twopi=two*pi)
c
      parameter (small=1.d-12,small3=small*small*small)
      parameter (hlxi=one/(pi+small),hlyi=hlxi)

c---------------------------------------------------------------
c      PV contour arrays:
      common /cont01/ xd(nplm),yd(nplm),dx(nplm),dy(nplm)
      common /cont02/  b(nplm), c(nplm), d(nplm)
      common /cont03/  u(nplm), v(nplm), q(nplm)
      common /cont04/ next(nplm)
      common /cont05/ i1d(nm),i2d(nm),npd(nm)

c      Grid -> Contour arrays:
      common /g2cc01/ qa(nhf,nhf)
      common /g2cc02/ xg(nhf+1)
      common /g2cc03/ qlev(2*nlevm)
      common /g2cc04/ ibx(nhf,0:1)

c      Basic parameters:
      common /base01/ amu,ell,elf,densf,dm,dmsq,dmi
      common /base02/ gl,gli,fnwhf
      common /base03/ dq,dqi,qoff,bfac
      common /base04/ na,nd,npta,nptd
      common /base05/ levbeg,levend

c      Control variables:
      logical corn(nplm),last
      common /ctrl01/ corn
      common /ctrl02/ last

