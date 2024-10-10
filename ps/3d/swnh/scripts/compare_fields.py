#!/usr/bin/env python

# This script plots a selected field from the SW, GN and 3D directories,
# at three times specified below and for a selected resolution.  

# NOTE: The differences SW-3D and GN-3D are plotted.

#========== Perform the generic imports =========
import warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)

# set tick label size:
label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
# set x tick width and size:
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
# set y tick width and size:
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
# set axes width:
mpl.rcParams['axes.linewidth'] = 2

#====================== Function definitions =======================
def contint(fmin,fmax):
    #Determines a nice contour interval (giving 10-20 divisions with
    #interval 1, 2 or 5x10^m for some m) given the minimum & maximum
    #values of the field data (fmin & fmax).

    fmax=0.9999999*fmax
    fmin=0.9999999*fmin
    #The 0.99... factor avoids having a superfluous tick interval
    #in cases where fmax-fmin is 10^m or 2x10^m

    mpow=0
    rmult=fmax-fmin
    while rmult < 10.0:
       mpow+=1
       rmult=rmult*10.0

    while rmult >= 100.0:
       mpow-=1
       rmult=rmult/10.0

    emag=10.0**(float(-mpow))

    kmult=int(rmult/10.0)

    if kmult < 1:
       ci=emag
    elif kmult < 2:
       ci=2.0*emag
    elif kmult < 4:
       ci=4.0*emag
    elif kmult < 8:
       ci=10.0*emag
    else:
       ci=20.0*emag

    return ci

#=================================================================
# Select data to compare:
field_list=['h','zeta','delta','w']
field_2d=['hh','zz','dd','ww']
field_3d=['h','zeta','d','w']

print
t_list=[5.0,15.0,25.0]
print ' Showing times ',t_list[0],', ',t_list[1],' and ',t_list[2]
t_list=np.array(t_list)

print
ng=int(raw_input(' Resolution (default 256)? ') or 256)
print

# Define grid:
#xg=np.linspace(-np.pi,np.pi,ng+1)
#yg=xg
#X,Y=np.meshgrid(xg,yg)
N=ng*ng

#=================================================================
# List of directories to compare (need the final /):
direc=['../plane/sw/swps/vstrip/ng'+str(ng)+'ld0.5r001/','../plane/gn/gnps/vstrip/ng'+str(ng)+'ld0.5r001/','../3d/swnh/ng'+str(ng)+'ld0.5r001/2d/']

# Open ecomp.asc file in one directory to get time between frames (dt):
in_file=open(direc[0]+'ecomp.asc','r')
time, e1, e2, e3, e4, e5 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

dt=time[1]-time[0]

# Frames to display:
frame_list=[]
for t in t_list:
   frame_list.append(int((t+1.e-6)/dt))

#=================================================================
# Loop over fields (h, zeta, delta and w):
print ' To view the images, type'
for k in range(len(field_list)):
# Set up figure:
   fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(figsize=[20,13], nrows=2, ncols=3)

   # Customise tick values:
   ax1.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax2.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax3.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax4.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax5.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax6.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax2.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax3.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax4.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax5.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax6.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

   ax1.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax2.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax3.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax4.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax5.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax6.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax1.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax2.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax3.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax4.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax5.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax6.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

   # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
   plt.setp(ax1.get_xticklabels(), visible=False)
   plt.setp(ax2.get_xticklabels(), visible=False)
   plt.setp(ax3.get_xticklabels(), visible=False)
   plt.setp(ax2.get_yticklabels(), visible=False)
   plt.setp(ax3.get_yticklabels(), visible=False)
   plt.setp(ax5.get_yticklabels(), visible=False)
   plt.setp(ax6.get_yticklabels(), visible=False)

   # Add information in the left panels:
   ax1.text( 0.08, 0.86, 'SW - 3D', transform=ax1.transAxes, fontsize=25)
   ax4.text( 0.08, 0.86, 'GN - 3D', transform=ax4.transAxes, fontsize=25)

   # Add times above top panels:
   ax1.set_title('$t = {x:.0f}$'.format(x=t_list[0]), fontsize=30)
   ax2.set_title('$t = {x:.0f}$'.format(x=t_list[1]), fontsize=30)
   ax3.set_title('$t = {x:.0f}$'.format(x=t_list[2]), fontsize=30)

   #fig.subplots_adjust(wspace=0.8, hspace=0.5, top=0.6, bottom=0.05)
   fig.subplots_adjust(wspace=0.8, hspace=-0.1)

   #Select field:
   field=field_list[k]
   acro2=field_2d[k]
   acro3=field_3d[k]

   #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   # First time:
   frame=frame_list[0]

   # Read data into arrays for plotting:
   Z1=np.empty([ng+1,ng+1])
   Z2=np.empty([ng+1,ng+1])
   Z3=np.empty([ng+1,ng+1])

   # Read full fields in all directories:
   in_file=open(direc[0]+acro2+'.r4','r')
   raw_array=np.fromfile(in_file,dtype=np.float32)
   in_file.close()
   Z1[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

   in_file=open(direc[1]+acro2+'.r4','r')
   raw_array=np.fromfile(in_file,dtype=np.float32)
   in_file.close()
   Z2[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

   in_file=open(direc[2]+acro3+'.r4','r')
   raw_array=np.fromfile(in_file,dtype=np.float32)
   in_file.close()
   Z3[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

   # Add periodic edges:
   Z1[ng,0:ng]=Z1[0,0:ng]
   Z1[0:ng+1,ng]=Z1[0:ng+1,0]
   Z2[ng,0:ng]=Z2[0,0:ng]
   Z2[0:ng+1,ng]=Z2[0:ng+1,0]
   Z3[ng,0:ng]=Z3[0,0:ng]
   Z3[0:ng+1,ng]=Z3[0:ng+1,0]

   # Compute differences from the vertically-averaged 3D field:
   Z1=Z1-Z3
   Z2=Z2-Z3

   # Work out the overall min/max values:
   zmin1=np.amin(Z1)
   zmin2=np.amin(Z2)
   zmax1=np.amax(Z1)
   zmax2=np.amax(Z2)

   # Make all limits the same and centre the limits on 0:
   zmag=max(abs(zmin1),abs(zmin2),abs(zmax1),abs(zmax2))
   zmin1=-zmag
   zmax1=zmag
   zmin2=-zmag
   zmax2=zmag

   # Obtain contour levels for plotting the colorbars:
   dz=contint(zmin1,zmax1)
   jmin=-int(-zmin1/dz)
   jmax= int( zmax1/dz)
   clevels1=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

   dz=contint(zmin2,zmax2)
   jmin=-int(-zmin2/dz)
   jmax= int( zmax2/dz)
   clevels2=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

   # Plot the images in an array with individual colourbars:
   im1=ax1.imshow(Z1,cmap=cm.seismic,vmin=zmin1,vmax=zmax1,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   divider = make_axes_locatable(ax1)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im1, cax=cax, ticks=clevels1)
   setp(cbar.ax.yaxis.set_ticklabels(clevels1), fontsize=18)

   im2=ax4.imshow(Z2,cmap=cm.seismic,vmin=zmin2,vmax=zmax2,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   divider = make_axes_locatable(ax4)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im2, cax=cax, ticks=clevels2)
   setp(cbar.ax.yaxis.set_ticklabels(clevels2), fontsize=18)

   #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   # Second time:
   frame=frame_list[1]

   # Read data into arrays for plotting:
   Z1=np.empty([ng+1,ng+1])
   Z2=np.empty([ng+1,ng+1])
   Z3=np.empty([ng+1,ng+1])

   # Read full fields in all directories:
   in_file=open(direc[0]+acro2+'.r4','r')
   raw_array=np.fromfile(in_file,dtype=np.float32)
   in_file.close()
   Z1[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

   in_file=open(direc[1]+acro2+'.r4','r')
   raw_array=np.fromfile(in_file,dtype=np.float32)
   in_file.close()
   Z2[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

   in_file=open(direc[2]+acro3+'.r4','r')
   raw_array=np.fromfile(in_file,dtype=np.float32)
   in_file.close()
   Z3[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

   # Add periodic edges:
   Z1[ng,0:ng]=Z1[0,0:ng]
   Z1[0:ng+1,ng]=Z1[0:ng+1,0]
   Z2[ng,0:ng]=Z2[0,0:ng]
   Z2[0:ng+1,ng]=Z2[0:ng+1,0]
   Z3[ng,0:ng]=Z3[0,0:ng]
   Z3[0:ng+1,ng]=Z3[0:ng+1,0]

   # Compute differences from the vertically-averaged 3D field:
   Z1=Z1-Z3
   Z2=Z2-Z3

   # Work out the overall min/max values:
   zmin1=np.amin(Z1)
   zmin2=np.amin(Z2)
   zmax1=np.amax(Z1)
   zmax2=np.amax(Z2)

   # Make all limits the same and centre the limits on 0:
   zmag=max(abs(zmin1),abs(zmin2),abs(zmax1),abs(zmax2))
   zmin1=-zmag
   zmax1=zmag
   zmin2=-zmag
   zmax2=zmag

   # Obtain contour levels for plotting the colorbars:
   dz=contint(zmin1,zmax1)
   jmin=-int(-zmin1/dz)
   jmax= int( zmax1/dz)
   clevels1=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

   dz=contint(zmin2,zmax2)
   jmin=-int(-zmin2/dz)
   jmax= int( zmax2/dz)
   clevels2=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

   # Plot the images in an array with individual colourbars:
   im1=ax2.imshow(Z1,cmap=cm.seismic,vmin=zmin1,vmax=zmax1,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   divider = make_axes_locatable(ax2)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im1, cax=cax, ticks=clevels1)
   setp(cbar.ax.yaxis.set_ticklabels(clevels1), fontsize=18)

   im2=ax5.imshow(Z2,cmap=cm.seismic,vmin=zmin2,vmax=zmax2,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   divider = make_axes_locatable(ax5)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im2, cax=cax, ticks=clevels2)
   setp(cbar.ax.yaxis.set_ticklabels(clevels2), fontsize=18)

   #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   # Last time:
   frame=frame_list[2]

   # Read data into arrays for plotting:
   Z1=np.empty([ng+1,ng+1])
   Z2=np.empty([ng+1,ng+1])
   Z3=np.empty([ng+1,ng+1])

   # Read full fields in all directories:
   in_file=open(direc[0]+acro2+'.r4','r')
   raw_array=np.fromfile(in_file,dtype=np.float32)
   in_file.close()
   Z1[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

   in_file=open(direc[1]+acro2+'.r4','r')
   raw_array=np.fromfile(in_file,dtype=np.float32)
   in_file.close()
   Z2[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

   in_file=open(direc[2]+acro3+'.r4','r')
   raw_array=np.fromfile(in_file,dtype=np.float32)
   in_file.close()
   Z3[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

   # Add periodic edges:
   Z1[ng,0:ng]=Z1[0,0:ng]
   Z1[0:ng+1,ng]=Z1[0:ng+1,0]
   Z2[ng,0:ng]=Z2[0,0:ng]
   Z2[0:ng+1,ng]=Z2[0:ng+1,0]
   Z3[ng,0:ng]=Z3[0,0:ng]
   Z3[0:ng+1,ng]=Z3[0:ng+1,0]

   # Compute differences from the vertically-averaged 3D field:
   Z1=Z1-Z3
   Z2=Z2-Z3

   # Work out the overall min/max values:
   zmin1=np.amin(Z1)
   zmin2=np.amin(Z2)
   zmax1=np.amax(Z1)
   zmax2=np.amax(Z2)

   # Make all limits the same and centre the limits on 0:
   zmag=max(abs(zmin1),abs(zmin2),abs(zmax1),abs(zmax2))
   zmin1=-zmag
   zmax1=zmag
   zmin2=-zmag
   zmax2=zmag

   # Obtain contour levels for plotting the colorbars:
   dz=contint(zmin1,zmax1)
   jmin=-int(-zmin1/dz)
   jmax= int( zmax1/dz)
   clevels1=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

   dz=contint(zmin2,zmax2)
   jmin=-int(-zmin2/dz)
   jmax= int( zmax2/dz)
   clevels2=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

   # Plot the images in an array with individual colourbars:
   im1=ax3.imshow(Z1,cmap=cm.seismic,vmin=zmin1,vmax=zmax1,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   divider = make_axes_locatable(ax3)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im1, cax=cax, ticks=clevels1)
   setp(cbar.ax.yaxis.set_ticklabels(clevels1), fontsize=18)

   im2=ax6.imshow(Z2,cmap=cm.seismic,vmin=zmin2,vmax=zmax2,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
   divider = make_axes_locatable(ax6)
   cax = divider.append_axes("right", size="7%", pad=0.1)
   cbar=fig.colorbar(im2, cax=cax, ticks=clevels2)
   setp(cbar.ax.yaxis.set_ticklabels(clevels2), fontsize=18)

   outfile=field+str(ng)+'diff.eps'
   # Save image:
   fig.savefig(outfile, format='eps', dpi=300)

   print ' gv',outfile
