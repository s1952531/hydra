#!/usr/bin/env python3

# This script plots the PV field at a selected time or frame
# from fine-grid data created using genfg.

#========== Perform the generic imports =========
import sys,os,warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib.colors as clrs
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

#===========================================================================
print
nx=int(raw_input(' Grid resolution (default 256)? ') or 256)
ng=nx
N=ng*ng
#dpi=int(raw_input(' DPI to save plot with (default 600)? ') or 600)
dpi=300

field=str(raw_input(' Field to show (a, d, g, h, g, q or zeta; default d)? ') or 'd')

frame=int(raw_input(' Frame to show (1, 2, ...; default 51)? ') or 51)
yam=int(raw_input(' Max |y| to show (default 1)? ') or 1.0)
cbopt=str(raw_input(' Add a colourbar (default y)? ') or 'y')

# Read data:
in_file=open('../3d/swnh/ng'+str(ng)+'ld0.5r001/2d/'+field+'.r4','r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()

Z=np.empty([ng+1,ng+1])
Z[0:ng,0:ng]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(ng,ng).T

# Add periodic edges:
Z[ng,0:ng]=Z[0,0:ng]
Z[0:ng+1,ng]=Z[0:ng+1,0]

# Determine pixel range in y:
iy1=int(0.5*float(ng)*(1.0-yam/np.pi))
iy2=ng+1-iy1

# Get min/max field values:
zmin=np.amin(Z)
zmax=np.amax(Z)

# Get aspect ratio of the plot:
aspect=yam/np.pi

# Image width:
width=8.0

# Set up figure:
fig = plt.figure(1,figsize=[width,width])
ax = fig.add_subplot(111)

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

if field == 'q':
   cmap='terrain'
elif field == 'a':
   cmap='Greys'
   zmin=0.0
else:
   cmap='seismic'
   zmax=max(abs(zmin),abs(zmax))
   zmin=-zmax

# Plot the image (optionally with a colourbar):
if cbopt == 'n':
   ax.imshow(Z[iy1:iy2,:],cmap=cmap,vmin=zmin,vmax=zmax,extent=(-np.pi,np.pi,-yam,yam),origin='lower',interpolation='bilinear')
else:
   im=ax.imshow(Z[iy1:iy2,:],cmap=cmap,vmin=zmin,vmax=zmax,extent=(-np.pi,np.pi,-yam,yam),origin='lower',interpolation='bilinear')
   # Obtain contour levels for plotting the colorbar:
   dz=contint(zmin,zmax)
   jmin=-int(-zmin/dz)
   jmax= int( zmax/dz)
   clevels=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)
   # Add a colourbar with nice intervals:
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("bottom", size="7%", pad=0.1)
   cbar=fig.colorbar(im, cax=cax, ticks=clevels, orientation='horizontal')
   setp(cbar.ax.xaxis.set_ticklabels(clevels), fontsize=20)

#=========================================================================
# Save image:
fig.savefig(field+str(frame)+'.eps', format='eps', dpi=dpi, bbox_inches = 'tight')

print
print ' To view the image, type'
print
print ' gv '+field+str(frame)+'.eps'
print
