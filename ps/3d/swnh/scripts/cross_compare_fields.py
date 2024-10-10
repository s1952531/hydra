#!/usr/bin/env python

# This script plots the difference of a selected field from the 3D data
# for the SW, GN and VA models.

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

    emag=10.0**(-mpow)

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
# Specify directories:
dir3d='../3d/swnh/hbar0.4ng256ld0.5/2d/'
dirsw='../plane/sw/ng256/'
dirgn='../plane/gn/hbar0.4ng256ld0.5/'
dirva='../plane/va/hbar0.4ng256ld0.5/'

# Select data to compare:
print ' Choose one of the following fields to plot:'
print
print ' (1) h;  (2) zeta;  (3) delta;  (4) gamma-tilde;  (5) P_n'
print
iopt=int(raw_input(' Choice (default 3): ') or 3)
print

t=float(raw_input(' Time to compare (default 25): ') or 25)
print

ng=int(raw_input(' Resolution (default 256)? ') or 256)
N=ng*ng
print

hbar=float(raw_input('Mean depth H (default 0.4)? ') or 0.4)
print

# Open input files:
in_file3d=open(dir3d+'h.r4','r')
in_filesw=open(dirsw+'hh.r4','r')
in_filegn=open(dirgn+'hh.r4','r')
in_fileva=open(dirva+'hh.r4','r')

if iopt == 1:
   title_label='$\\tilde{h}$'
   field='h'
elif iopt == 2:
   title_label='$\\zeta$'
   field='zeta'
elif iopt == 3:
   title_label='$\\delta$'
   field='delta'
elif iopt == 4:
   title_label='$\\tilde{\\gamma}$'
   field='gamma-tilde'
else:
   title_label='$\\bar{P}_n$'
   field='pn'

title_label=title_label+',$\ \ \ \ \ t = {x:.0f}$'.format(x=t)+',$\ \ \ \ \ H = {x:.1f}$'.format(x=hbar)+',$\ \ \ \ \ n = {x:.0f}$'.format(x=ng)
   
# Open ecomp.asc file in one directory to get time between frames (dt):
in_file=open(dirsw+'ecomp.asc','r')
time, e1, e2, e3, e4, e5 = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

dt=time[1]-time[0]

# Frame to display:
frame=int(t/dt+0.001)

#=================================================================
# Set up figure:
fig, (ax1, ax2, ax3) = plt.subplots(figsize=[20,7], nrows=1, ncols=3)

# Customise tick values:
ax1.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax2.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax3.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax2.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax3.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

ax1.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax2.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax3.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax1.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax2.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax3.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

# Fine-tune figure; hide y ticks for right plots
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)

# Add titles:
ax1.set_title('SW - 3D', fontsize=30, fontname='Times New Roman')
ax2.set_title('GN - 3D', fontsize=30, fontname='Times New Roman')
ax3.set_title('VA - 3D', fontsize=30, fontname='Times New Roman')

# Add an overall title:
fig.suptitle(title_label, fontsize=32)

#fig.subplots_adjust(wspace=0.8, hspace=-0.1)

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Read data into arrays for plotting:

# 3D data:
Z0=np.empty([ng+1,ng+1])
raw_array=np.fromfile(in_file3d,dtype=np.float32)
in_file3d.close()
Z0[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
# Add periodic edges:
Z0[ng,0:ng]=Z0[0,0:ng]
Z0[0:ng+1,ng]=Z0[0:ng+1,0]

# SW data:
Z1=np.empty([ng+1,ng+1])
raw_array=np.fromfile(in_filesw,dtype=np.float32)
in_filesw.close()
Z1[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z1[ng,0:ng]=Z1[0,0:ng]
Z1[0:ng+1,ng]=Z1[0:ng+1,0]
# Compute differences from the vertically-averaged 3D field:
Z1=Z1-Z0
# Work out the overall min/max values for 3D and 2D fields:
zmin1=np.amin(Z1)
zmax1=np.amax(Z1)
# Make all limits the same and centre the limits on 0:
zmag=max(abs(zmin1),abs(zmax1))
zmin1=-zmag
zmax1=zmag

# GN data:
Z2=np.empty([ng+1,ng+1])
raw_array=np.fromfile(in_filegn,dtype=np.float32)
in_filegn.close()
Z2[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z2[ng,0:ng]=Z2[0,0:ng]
Z2[0:ng+1,ng]=Z2[0:ng+1,0]
# Compute differences from the vertically-averaged 3D field:
Z2=Z2-Z0
# Work out the overall min/max values for 3D and 2D fields:
zmin2=np.amin(Z2)
zmax2=np.amax(Z2)
# Make all limits the same and centre the limits on 0:
zmag=max(abs(zmin2),abs(zmax2))
zmin2=-zmag
zmax2=zmag

# VA data:
Z3=np.empty([ng+1,ng+1])
raw_array=np.fromfile(in_fileva,dtype=np.float32)
in_fileva.close()
Z3[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T
Z3[ng,0:ng]=Z3[0,0:ng]
Z3[0:ng+1,ng]=Z3[0:ng+1,0]
# Compute differences from the vertically-averaged 3D field:
Z3=Z3-Z0
# Work out the overall min/max values for 3D and 2D fields:
zmin3=np.amin(Z3)
zmax3=np.amax(Z3)
# Make all limits the same and centre the limits on 0:
zmag=max(abs(zmin3),abs(zmax3))
zmin3=-zmag
zmax3=zmag

# Obtain contour levels for plotting the colorbars:
dz=contint(zmin1,zmax1)
jmin=-int(-zmin1/dz)
jmax= int( zmax1/dz)
clevels1=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

dz=contint(zmin2,zmax2)
jmin=-int(-zmin2/dz)
jmax= int( zmax2/dz)
clevels2=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

dz=contint(zmin3,zmax3)
jmin=-int(-zmin3/dz)
jmax= int( zmax3/dz)
clevels3=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

# Plot the images in an array with individual colourbars:
im1=ax1.imshow(Z1,cmap=cm.seismic,vmin=zmin1,vmax=zmax1,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im1, cax=cax)
#setp(cbar.ax.yaxis.set_ticklabels(clevels1), fontsize=18)
#cbar=fig.colorbar(im1, cax=cax, ticks=clevels1)
#setp(cbar.ax.yaxis.set_ticklabels(clevels1), fontsize=18)

im2=ax2.imshow(Z2,cmap=cm.seismic,vmin=zmin2,vmax=zmax2,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax2)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im2, cax=cax)
#cbar=fig.colorbar(im2, cax=cax, ticks=clevels2)
#setp(cbar.ax.yaxis.set_ticklabels(clevels2), fontsize=18)

im3=ax3.imshow(Z3,cmap=cm.seismic,vmin=zmin3,vmax=zmax3,extent=(-np.pi,np.pi,-np.pi,np.pi),origin='lower',interpolation='bilinear')
divider = make_axes_locatable(ax3)
cax = divider.append_axes("right", size="7%", pad=0.1)
cbar=fig.colorbar(im3, cax=cax)
#cbar=fig.colorbar(im3, cax=cax, ticks=clevels3)
#setp(cbar.ax.yaxis.set_ticklabels(clevels3), fontsize=18)

#=========================================================================
# Save figure:
outfile='cross_'+field+'_H{x:.1f}'.format(x=hbar)+'n'+str(ng)+'t{x:.0f}'.format(x=t)+'.eps'
fig.savefig(outfile, format='eps', dpi=1200)

print
print ' To view the image, type'
print
print ' gv',outfile,'&'
print
