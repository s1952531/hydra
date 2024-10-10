      program print
c--------------------------------------------------------------------
c      Displays masses on a surface of revolution S from a chosen
c      orthographic perspective
c--------------------------------------------------------------------
      parameter (n=2,npc=19,ng=12,nt=2*ng,nd=10)
c      n:   number of masses
c      npc: number of points used to plot circles over each mass
c      ng:  number of "co-latitude" divisions for plotting S
c      nt:  number of "longitude" divisions for plotting S
c      nd:  number of subdivisions in co-latitude or longitude when
c           plotting a circle of constant longitude or co-latitude;
c           ng*nd = ngf and nt*nd = ntf points are plotted, respectively
      parameter (ngf=ng*nd,ntf=nt*nd)
      parameter (wpage=7.68,rad=0.005,fdash=0.05)
c      wpage:  the width (in inches) of the plotted page that 
c              is assumed by the graphics package;
c      fdash:  the fractional length of the segment connecting
c              two points of constant phi or constant theta, 
c              for rendering the 3D shape of the surface
c      rad:    the "radius" given to the mass to make it visible
c              in the plot
      parameter (ald=0.5-0.5*fdash,bld=1.-ald)
      parameter (pi=3.14159265358979323846264)
      parameter (twopi=2.*pi,pif=pi/180.)

      real xv(n),yv(n),zv(n),t
      real u1x(n),u1y(n),u2x(n),u2y(n),u2z(n)
      real x(npc),y(npc),xp(npc),yp(npc)

c-------------------------------------------------------------
c      Compute various constants needed for rendering surface:
      dth=pi/float(ng)
      dthf=pi/float(ngf)
      dph=twopi/float(nt)
      dphf=twopi/float(ntf)

c---------------------------------------------------
      write(*,*) ' Co-latitude and longitude of the view (in degrees)?'
      read(*,*) thv,phv

c      Define rotation matrix:
      thv=pif*thv
      phv=pif*phv
      cthv=cos(thv)
      sthv=sin(thv)
      cphv=cos(phv)
      sphv=sin(phv)
      txx=cthv*cphv
      txy=cthv*sphv
      txz=-sthv
      tyx=-sphv
      tyy=cphv
      tzx=sthv*cphv
      tzy=sthv*sphv
      tzz=cthv

c      Work out min & max x & y of surface:
      xmin=100.
      ymin=100.
      xmax=-xmin
      ymax=-ymin
      do ith=1,ng-1
        th=dth*float(ith)
        rho=sin(th)
        zs=cos(th)
        do iph=1,ntf
          ph1=dphf*(float(iph)+ald)
          xs1=rho*cos(ph1)
          ys1=rho*sin(ph1)
          ph2=dphf*(float(iph)+bld)
          xs2=rho*cos(ph2)
          ys2=rho*sin(ph2)
c          Rotate points into chosen perspective:
          xh1=txx*xs1+txy*ys1+txz*zs
          yh1=tyx*xs1+tyy*ys1
          xh2=txx*xs2+txy*ys2+txz*zs
          yh2=tyx*xs2+tyy*ys2
          xmin=min(xmin,xh1,xh2)
          xmax=max(xmax,xh1,xh2)
          ymin=min(ymin,yh1,yh2)
          ymax=max(ymax,yh1,yh2)
        enddo
      enddo

      do iph=1,nt
        ph=dph*float(iph-1)
        cph=cos(ph)
        sph=sin(ph)
        do ith=1,ngf
          th1=dthf*(float(ith)+ald)
          rho=sin(th1)
          xs1=rho*cph
          ys1=rho*sph
          zs1=cos(th1)
          th2=dthf*(float(ith)+bld)
          rho=sin(th2)
          xs2=rho*cph
          ys2=rho*sph
          zs2=cos(th2)
c          Rotate points into chosen perspective:
          xh1=txx*xs1+txy*ys1+txz*zs1
          yh1=tyx*xs1+tyy*ys1
          xh2=txx*xs2+txy*ys2+txz*zs2
          yh2=tyx*xs2+tyy*ys2
          xmin=min(xmin,xh1,xh2)
          xmax=max(xmax,xh1,xh2)
          ymin=min(ymin,yh1,yh2)
          ymax=max(ymax,yh1,yh2)
        enddo
      enddo
      xdif=xmax-xmin
      xmax=xmax+0.005*xdif
      xmin=xmin-0.005*xdif
      ydif=ymax-ymin
      ymax=ymax+0.005*ydif
      ymin=ymin-0.005*ydif
      aspect=ydif/xdif

c      Create small circle to plot each mass:
      dthc=2.*pi/float(npc-1)
      do i=1,npc
        thc=dthc*float(i-1)
        xp(i)=rad*cos(thc)
        yp(i)=rad*sin(thc)
      enddo

c      Read energy file to work out how many time frames there are:
      open(44,file='energy.dat',status='old')
      nloops=0
3     read(44,*,end=4,err=4) t,ene
      nloops=nloops+1
      goto 3
4     close(44)

      if (n .gt. 1) then
        write(*,*) ' Which mass (use 0 for all)?'
        read(*,*) imass
      else
        imass=0
      endif

      if (imass .eq. 0) then
        ibeg=1
        iend=n
      else
        ibeg=imass
        iend=imass
      endif

      write(*,'(a,i5,a)') ' There are ',nloops,
     .                    ' time frames in the data.'
      write(*,*) ' Beginning and ending frames to use?'
      read(*,*) kfr1,kfr2

c      Open file containing point masses:
      nbytes=4*(3*n+1)
      open(10,file='points.r4',form='unformatted', 
     .       access='direct',status='old',recl=nbytes)

c      Initialize unix graphics:
      call gsetup

c      Start new image frame:
      call zframe
      call izoom
      call arzoom(aspect,aspect)
      call gscale(xmin,xmax,ymin,ymax)

c      Draw surface, first curves of constant theta, then constant phi:
      do ith=1,ng-1
        th=dth*float(ith)
        rho=sin(th)
        zs=cos(th)
        do iph=1,ntf
          ph1=dphf*(float(iph)+ald)
          xs1=rho*cos(ph1)
          ys1=rho*sin(ph1)
          ph2=dphf*(float(iph)+bld)
          xs2=rho*cos(ph2)
          ys2=rho*sin(ph2)
c          Rotate points into chosen perspective:
          xh1=txx*xs1+txy*ys1+txz*zs
          yh1=tyx*xs1+tyy*ys1
          xh2=txx*xs2+txy*ys2+txz*zs
          yh2=tyx*xs2+tyy*ys2
c          Draw segment rotated 90 degrees to come out upright in frame:
          call ggline(-xh1,-yh1,-xh2,-yh2,1)
        enddo
      enddo

      do iph=1,nt
        ph=dph*float(iph-1)
        cph=cos(ph)
        sph=sin(ph)
        do ith=1,ngf
          th1=dthf*(float(ith)+ald)
          rho=sin(th1)
          xs1=rho*cph
          ys1=rho*sph
          zs1=cos(th1)
          th2=dthf*(float(ith)+bld)
          rho=sin(th2)
          xs2=rho*cph
          ys2=rho*sph
          zs2=cos(th2)
c          Rotate points into chosen perspective:
          xh1=txx*xs1+txy*ys1+txz*zs1
          yh1=tyx*xs1+tyy*ys1
          xh2=txx*xs2+txy*ys2+txz*zs2
          yh2=tyx*xs2+tyy*ys2
c          Draw segment rotated 90 degrees to come out upright in frame:
          call ggline(-xh1,-yh1,-xh2,-yh2,1)
        enddo
      enddo

c      Draw tranjectory, as a series of circles with connecting lines:
      do loop=kfr1,kfr2

c        Read Cartesian coordinates of points (xv,yv,zv), then define local 
c        surface tangent unit vectors (u1x,u1y,0) and (u2x,u2y,u2z):
        read(10,rec=loop) t,xv,yv,zv
        write(*,'(a,f12.5)') ' t = ',t

        do i=ibeg,iend
          cth=zv(i)
          sth=sqrt(abs(1.-zv(i)**2))
          cph=xv(i)/sth
          sph=yv(i)/sth
          u1x(i)=-sph
          u1y(i)= cph
          u2x(i)= cph*cth
          u2y(i)= sph*cth
          u2z(i)=-sth
        enddo

c        Draw points as projected circles lying on the surface:
        do i=ibeg,iend
          do j=1,npc
c            Form circle (nearly) on surface:
            xs=xv(i)+xp(j)*u1x(i)+yp(j)*u2x(i)
            ys=yv(i)+xp(j)*u1y(i)+yp(j)*u2y(i)
            zs=zv(i)             +yp(j)*u2z(i)
c            Rotate point into chosen perspective:
            xh=txx*xs+txy*ys+txz*zs
            yh=tyx*xs+tyy*ys
c            Rotate by 90 degrees to come out upright in frame:
            x(j)=-xh
            y(j)=-yh
          enddo
c          Draw curve:
          call gcurv(npc,x,y,1)
        enddo

      enddo

c      terminate graphics:
      call gclear
      call gend
      end
