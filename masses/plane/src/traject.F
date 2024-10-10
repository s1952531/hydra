      program print
c--------------------------------------------------------------------
c      Displays masses on a surface of revolution S from a chosen
c      orthographic perspective
c--------------------------------------------------------------------
      parameter (n=N_POINTS,npc=19)
c      n:   number of masses
c      npc: number of points used to plot circles over each mass
      parameter (wpage=7.68,rad=0.005)
c      wpage:  the width (in inches) of the plotted page that is
c              assumed by the graphics package;
c      rad:    the "radius" given to the mass to make it visible
c              in the plot
      parameter (pi=3.14159265358979323846264)

      real xv(n),yv(n),t
      real x(npc),y(npc),xp(npc),yp(npc)

c--------------------------------------------------------------------
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

      write(*,'(a,i5,a)') ' There are ',nloops,
     .                    ' time frames in the data.'
      write(*,*) ' Beginning and ending frames to use?'
      read(*,*) kfr1,kfr2

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

c      Open file containing point masses:
      nbytes=4*(2*n+1)
      open(10,file='points.r4',form='unformatted', 
     .       access='direct',status='old',recl=nbytes)

c      Initialize unix graphics:
      call gsetup

c      Start new image frame:
      call zframe
      call izoom
      call arzoom(1.,1.)
      call gscale(-sc,sc,-sc,sc)

c      Draw tranjectory, as a series of circles:
      do loop=kfr1,kfr2

c        Read Cartesian coordinates of points (xv,yv):
        read(10,rec=loop) t,xv,yv
        write(*,'(a,f12.5)') ' t = ',t

c        Draw points as circles:
        do i=ibeg,iend
          do j=1,npc
            xs=xv(i)+xp(j)
            ys=yv(i)+yp(j)
c            Rotate to come out upright in frame:
            x(j)=-xs
            y(j)=-ys
          enddo
c          Draw curve:
          call gcurv(npc,x,y,1)
        enddo

      enddo

c      terminate graphics:
      call gclear
      call gend
      end
