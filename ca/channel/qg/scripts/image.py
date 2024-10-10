#!/usr/bin/env python

# This script plots data in qq.r4, qa.r4, qe.r4 or qo.r4,
# the latter three obtained after running pvsep.f90.

#  @@@@   Run from the current job directory   @@@@

#========== Perform the generic imports =========
import sys,os,warnings
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
mpl.rcParams['axes.linewidth'] = 1

#====================== Function definitions =======================
def contint(fmin,fmax):
    #Determines a nice contour interval (giving 10-20 divisions with
    #interval 1, 2 or 5x10^m for some m) given the minimum & maximum
    #values of the field data (fmin & fmax).

    fmax=0.9999999*fmax
    fmin=0.9999999*fmin
    #The 0.99... factor avoids having a superfluous tick interval
    #in cases where fmax-fmin is 10^m or 2x10^m

    emag=1.0
    rmult=fmax-fmin
    while rmult < 10:
       emag=emag/10
       rmult=rmult*10

    while rmult >= 100:
       emag=emag*10
       rmult=rmult/10

    kmult=int(rmult/10)

    if kmult < 1:
       ci=emag
    elif kmult < 2:
       ci=2*emag
    elif kmult < 4:
       ci=4*emag
    elif kmult < 8:
       ci=10*emag
    else:
       ci=20*emag

    return ci

#=================================================================
# Set dimensions:
nx=2048
nyf=513
nyh=(nyf-1)/2+1

ellx=20.0*np.pi
elly=5.0*np.pi

hlx=ellx/2.0
hly=elly/2.0

# Time between frames:
dtsave=5.0

# Total number of frames:
nfmax=120

# Select data to compare:
field_list=['qq','qa','qe','qo']
field_acro=['q','q_a','q_e','q_o']
ny_list=[nyf,nyf,nyh,nyh]

print()
print(' Which field do you wish to image')
print()
print(' (1) q')
print(' (2) q_a = q - beta*y')
print(' (3) q_e (even part of q_a), or')
print(' (4) q_o (odd  part of q_a)')
print()
op_in = input(' Option (default 3)? ')
option = int(op_in or 3)
k=option-1

ny=int(ny_list[k])
N=nx*ny

fname=field_list[k]+'.r4'
figname=field_acro[k]+'.eps'

nf_in = input(' Number of frames to show (default 3)? ')
nf = int(nf_in or 3)

ncol = 1
if nf%2 == 0:
    ncol_in = input(' Number of columns in figure layout (default 2)? ')
    ncol = int(ncol_in or 2)
nrow = int(nf/ncol+0.1)

frame_list=np.empty(nf)
print(' Enter the frames (note: 0 corresponds to t = 0)')
dframe=int(nfmax/nf+0.1)
frame_def=0
for j in range(nf):
    frame_def+=dframe
    frame_in = input(' Frame '+str(j)+' (default '+str(frame_def)+')? ')
    frame_list[j] = int(frame_in or frame_def)
print()

f_in = input(' Starting fraction of x domain to plot (default 0.25)? ')
fx1 = float(f_in or 0.25)
f_in = input(' Ending fraction of x domain to plot (default 0.515625)? ')
fx2 = float(f_in or 0.515625)

if k > 1:
    f_in = input(' Starting fraction of y domain to plot (default 0)? ')
    fy1 = float(f_in or 0.0)
    f_in = input(' Ending fraction of y domain to plot (default 0.25)? ')
    fy2 = float(f_in or 0.25)
else:
    f_in = input(' Starting fraction of y domain to plot (default 0.25)? ')
    fy1 = float(f_in or 0.25)
    f_def = 1.0-fy1
    f_in = input(' Ending fraction of y domain to plot (default '+str(f_def)+')? ')
    fy2 = float(f_in or f_def)

ix1=int(fx1*float(nx)+0.5)
ix2=min(int(fx2*float(nx)+0.5),nx-1)

iy1=int(fy1*float(ny)+0.5)
iy2=min(int(fy2*float(ny)+0.5),ny-1)

print()
print(' Range of x grid points:',ix1,'to',ix2)
print(' Range of y grid points:',iy1,'to',iy2)
print()

xmin=fx1*ellx-hlx
xmax=fx2*ellx-hlx
if k < 2:
    ymin=fy1*elly-hly
    ymax=fy2*elly-hly
else:
    ymin=fy1*hly
    ymax=fy2*hly

# Formatting for displayed time:
time_template = '$t =$ %.0f'

#==============================================================================
# Set up figure:
aspect=float(iy2-iy1-1)/float(ix2-ix1-1)

fig, ax = plt.subplots(figsize=[10,float(nrow)*(7.8*aspect/float(ncol)+0.5)], nrows=nrow, ncols=ncol)
ax = ax.flatten()

# Read data into array for plotting:
in_file=open(fname,'r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()

for j,frame in enumerate(frame_list):
    t=float(frame)*dtsave
    print(' Time t =',t)
    ax1=ax[j]

    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])

    if ncol == 1:
        if j == nf-1:
            ax1.set_xlabel('$x$', fontsize=20)
        else:
            plt.setp(ax1.get_xticklabels(), visible=False)
            ax1.set_ylabel('$y$', fontsize=20)
    else:
        if j == nf-1 or j == nf-2:
            ax1.set_xlabel('$x$', fontsize=20)
        else:
            plt.setp(ax1.get_xticklabels(), visible=False)
        if j%2 == 0:
            ax1.set_ylabel('$y$', fontsize=20)

    Z=np.empty([nx,ny])
    Z=raw_array[int(frame)*(N+1)+1:(int(frame)+1)*(N+1)].reshape(nx,ny)

    # Work out the overall min/max values:
    zmin=np.amin(Z)
    zmax=np.amax(Z)
    print()
    print(' Min & max field levels = '+str(zmin)+', '+str(zmax))

    zmag=max(abs(zmin),abs(zmax))
    zmin=-zmag
    zmax= zmag

    zm_in = input(' Min level to show (default '+str(zmin)+')? ')
    zmin = float(zm_in or zmin)
    zm_in = input(' Max level to show (default '+str(zmax)+')? ')
    zmax = float(zm_in or zmax)

    # Obtain contour levels for plotting the colorbars:
    dz=2.0*contint(zmin,zmax)
    jmin=-int(-zmin/dz)
    jmax= int( zmax/dz)
    clevels=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

    # Plot the image in an array with an optional colourbar:
    im1=ax1.imshow(Z[ix1:ix2,iy1:iy2].T,cmap=cm.seismic,vmin=zmin,vmax=zmax,extent=(xmin,xmax,ymin,ymax),origin='lower',interpolation='bilinear')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="4%", pad=0.1)
    cbar=fig.colorbar(im1, cax=cax, ticks=clevels)

    ax1.set_title('$t = {x:.0f}$'.format(x=t), fontsize=25)

fig.subplots_adjust(wspace=0.8, hspace=0.0)
    
#=========================================================================
# Save image:
fig.savefig(figname, format='eps', bbox_inches='tight', dpi=300)

print()
print(' To view the image, type')
print()
print(' gv '+figname+' &')
print()
