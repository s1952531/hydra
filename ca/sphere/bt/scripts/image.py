#!/usr/bin/env python3

# This script plots the vorticity field as a function of longitude and
# latitude at selected times.

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
print()
ncol = 3
q_in = input(' Number of columns in the figure (default '+str(ncol)+')? ')
ncol = int(q_in or ncol)
nrow = 3
q_in = input(' Number of rows in the figure (default '+str(nrow)+')? ')
nrow = int(q_in or nrow)

nim = ncol*nrow

#---------------------------------------------------------------------------
# Work out grid resolution and data save interval by reading parameters.f90:
with open('src/parameters.f90','r') as in_file:
    fread = in_file.readlines()
    for line in fread:
        if ':: ng=' in line:
            ng = int(line.split("=")[1].split(",")[0])
        if ',tsave=' in line:
            tsave = float(line.split("e=")[1].split("d0,tsim")[0])

nt = 2*ng
N = ng*nt

# Read all data:
in_file = open('zz.r4','r')
raw_array = np.fromfile(in_file,dtype = np.float32)
in_file.close()

d = np.empty((nim,nt+1,ng))
time = np.empty(nim)

print()
tp = 0.0
ts = 0.0
dts = 0.0
for i in range(nim):
    ts = ts+dts
    q_in = input(' Time to show for image '+str(i+1)+' (default '+str(ts)+')? ')
    ts = float(q_in or ts)
    frame = int(ts/tsave+0.5)
    time[i] = raw_array[frame*(N+1)]
    d[i,0:nt,0:ng] = raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(nt,ng)
    if i == 0:
        dts = tsave
    else:
        dts = ts-tp
    tp = ts
    # Add periodic edge:
    d[i,nt,:] = d[i,0,:]
    print(' Found time '+str(time[i])+' in the data.')
    
# Work out the overall max abs value:
dmax = abs(d.max())

# Obtain contour levels for plotting the colorbar:
dd = contint(-dmax,dmax)
jmax = int(dmax/dd)
clevels = np.linspace(-dd*float(jmax),dd*float(jmax),2*jmax+1)

#==============================================================================
# Set up figure:
fig, ax = plt.subplots(figsize = [20,10*float(nrow)/float(ncol)+1.5], nrows = nrow, ncols = ncol)
ax = ax.flatten()

# Plot all images:
for i in range(nim):
   ax0 = ax[i]
   # Customise tick values:
   ax0.set_xticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
   ax0.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=24)
   ax0.set_yticks([-np.pi/2.0,0.0,np.pi/2.0])
   ax0.set_yticklabels([r'$-\pi/2$',r'$0$',r'$\pi/2$'],fontsize=24)

   # Plot image:
   im = ax0.imshow(d[i].T,cmap=cm.seismic,vmin=-dmax,vmax=dmax,extent=(-np.pi,np.pi,-np.pi/2,np.pi/2),origin='lower',interpolation='bilinear')

   # Add title:
   ax0.set_title('$t = {x:.1f}$'.format(x=time[i]), fontsize=24)

# Fine-tune figure; hide x ticks for top plots and y ticks for right plots:
if nrow > 1:
    for i in range(nim-ncol):
        ax0 = ax[i]
        plt.setp(ax0.get_xticklabels(), visible=False)

if ncol > 1:
    for j in range(nrow):
        for i in range(j*ncol+1,(j+1)*ncol):
            ax0 = ax[i]
            plt.setp(ax0.get_yticklabels(), visible=False)

# Add colourbar underneath:
#ax_cbar = fig.add_axes([0.125, -0.05, 0.75, 0.03])
#cbar = fig.colorbar(im, cax=ax_cbar, ticks=clevels, orientation='horizontal')
#setp(cbar.ax.xaxis.set_ticklabels(clevels), fontsize=24)

#=========================================================================
# Save image:

print()
print(' To view the image, type')
print()
fig.savefig('vorticity.pdf', format='pdf', dpi=300)
print(' ev vorticity.pdf &')
print()
