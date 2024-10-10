      program print
c--------------------------------------------------------------------
c      Displays masses on the 2D plane
c--------------------------------------------------------------------
      parameter (n=N_POINTS,npc=19)
c      n:   number of masses
c      npc: number of points used to plot circles over each mass
      parameter (wpage=7.68,rad=0.02)
c      wpage:  the width (in inches) of the plotted page that 
c              is assumed by the graphics package;
c      rad:    the "radius" given to the mass to make it visible
c              in the plot
      parameter (pi=3.14159265358979323846264)

      real xv(n),yv(n),t
      real x(npc),y(npc),xp(npc),yp(npc)

c-------------------------------------------------------------
      write(*,*) ' Plot scale (max{|x|,|y|})?'
      read(*,*) sc

      srad=sc*rad
c      Create small circle to plot each mass:
      dthc=2.*pi/float(npc-1)
      do i=1,npc
        thc=dthc*float(i-1)
        xp(i)=srad*cos(thc)
        yp(i)=srad*sin(thc)
      enddo

c      Read energy file to work out how many time frames there are:
      open(44,file='energy.dat',status='old')
      nloops=0
3     read(44,*,end=4,err=4) t,ene
      nloops=nloops+1
      goto 3
4     close(44)

c      Open file containing point masses:
      nbytes=4*(2*n+1)
      open(10,file='points.r4',form='unformatted', 
     .       access='direct',status='old',recl=nbytes)

c      Initialize unix graphics:
      call gsetup

c      Start repeated loop over data (infinite loop):

      do loop=1,nloops

c      Read Cartesian coordinates of points (xv,yv):
      read(10,rec=loop) t,xv,yv
      write(*,'(a,f12.5)') ' t = ',t

c      Start new image frame:
      call zframe
      call izoom
      call arzoom(1.,1.)
      call gscale(-sc,sc,-sc,sc)

c      Draw points as circles:
      do i=1,n
        do j=1,npc
c          Form circle (nearly) on surface:
          xs=xv(i)+xp(j)
          ys=yv(i)+yp(j)
c          Rotate to come out upright in frame:
          x(j)=-ys
          y(j)= xs
        enddo
c        Draw curve:
        call gcurv(npc,x,y,1)
      enddo

c      Clear image:
      call gclear

      enddo

c      terminate graphics:
      call gend
      end
