      implicit double precision(a-h,o-z)

      include 'g2c_dim.i'
c      nx, ny:   number of grid boxes in the x & y directions
c      nz:       number of layers in the z direction
      parameter (nxp1=nx+1,nyp1=ny+1)
      parameter (nxp2=nx+2,nyp2=ny+2)

      parameter (nm=npm/20)
c      npm:  maximum number of nodes on all contours
c      nm:   maximum number of contours
c      nplm: maximum number of nodes in any given layer

      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0)
      parameter (half=one/two,third=one/three)
      parameter (four=4.d0,six=6.d0,quart=one/four,sixth=one/six)

      parameter (small=1.d-12,small3=small*small*small)

c---------------------------------------------------------------
c      PV contour arrays:
      common /cont01/ xa(npm), ya(npm)
      common /cont02/ xd(nplm),yd(nplm),dx(nplm),dy(nplm)
      common /cont03/  a(nplm), b(nplm), c(nplm), d(nplm)
      common /cont04/  u(nplm), v(nplm), q(nplm)
      common /cont05/ next(nplm)
      common /cont06/ i1d(nm),i2d(nm),npd(nm),indd(nm)
      common /cont07/ i1a(nm),i2a(nm),npa(nm),inda(nm),laya(nm)
      common /cont08/ il1a(nz),il2a(nz),jl1a(nz),jl2a(nz)

c      Grid -> Contour arrays:
      common /g2cc01/ qa(0:nyp1,0:nxp1)
      common /g2cc02/ xg(0:nx),yg(0:ny)
      common /g2cc03/ qlev(2*nlevm)
      common /g2cc04/ ibx(0:nx,0:1),iby(0:ny,0:1)

c      Basic parameters:
      common /base01/ amu,ell,elf,densf,dm,dmsq,dmi
      common /base02/ glx,gly,dq,dqi,qoff
      common /base03/ xbeg,xend,ybeg,yend
      common /base04/ xwid,hlxi,xavg,ywid,hlyi,yavg
      common /base05/ ixbeg,iybeg
      common /base06/ na,nd,npta,nptd
      common /base07/ levbeg,levend

c      Control variables:
      logical corn(npm),last,perix,periy
      common /ctrl01/ corn
      common /ctrl02/ last,perix,periy
