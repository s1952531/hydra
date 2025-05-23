#!/usr/bin/env python3

# This script plots the relative rms difference norms for h, zeta and delta
# versus time, considering the differences SW - 3D and GN - 3D.

# The data directories are specified below

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

#=================================================================
# Specify data directories:
dir3d='../3d/swnh/hbar0.4ng256ld0.5/2d/'
dirsw='../plane/sw/ng256/'
dirgn='../plane/gn/hbar0.4ng256ld0.5/'

ng=int(raw_input(' Resolution (default 256)? ') or 256)
N=ng*ng
print

hbar=float(raw_input('Mean depth H (default 0.4)? ') or 0.4)
print

# Read field data and compute rms differences:
Zsw=np.empty(N)
Zgn=np.empty(N)
Z3d=np.empty(N)

# Height anomaly:
sw_file=open(dirsw+'hh.r4','r')
sw_array=np.fromfile(sw_file,dtype=np.float32)
sw_file.close()

gn_file=open(dirgn+'hh.r4','r')
gn_array=np.fromfile(gn_file,dtype=np.float32)
gn_file.close()

three_file=open(dir3d+'h.r4','r')
three_array=np.fromfile(three_file,dtype=np.float32)
three_file.close()

nframes = int(len(sw_array)/(N+1))
time=np.empty(nframes)
hsw=np.empty(nframes)
hgn=np.empty(nframes)

for frame in range(nframes):
   time[frame]=sw_array[frame*(N+1)]
   Zsw=sw_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Zgn=gn_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Z3d=three_array[frame*(N+1)+1:(frame+1)*(N+1)]
   hsw[frame]=np.sqrt(sum((Zsw-Z3d)**2)/sum(Z3d**2))
   hgn[frame]=np.sqrt(sum((Zgn-Z3d)**2)/sum(Z3d**2))

# Relative vorticity:
sw_file=open(dirsw+'zz.r4','r')
sw_array=np.fromfile(sw_file,dtype=np.float32)
sw_file.close()

gn_file=open(dirgn+'zz.r4','r')
gn_array=np.fromfile(gn_file,dtype=np.float32)
gn_file.close()

three_file=open(dir3d+'zeta.r4','r')
three_array=np.fromfile(three_file,dtype=np.float32)
three_file.close()

nframes = int(len(sw_array)/(N+1))
time=np.empty(nframes)
zsw=np.empty(nframes)
zgn=np.empty(nframes)

for frame in range(nframes):
   time[frame]=sw_array[frame*(N+1)]
   Zsw=sw_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Zgn=gn_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Z3d=three_array[frame*(N+1)+1:(frame+1)*(N+1)]
   zsw[frame]=np.sqrt(sum((Zsw-Z3d)**2)/sum(Z3d**2))
   zgn[frame]=np.sqrt(sum((Zgn-Z3d)**2)/sum(Z3d**2))

# Divergence:
sw_file=open(dirsw+'dd.r4','r')
sw_array=np.fromfile(sw_file,dtype=np.float32)
sw_file.close()

gn_file=open(dirgn+'dd.r4','r')
gn_array=np.fromfile(gn_file,dtype=np.float32)
gn_file.close()

three_file=open(dir3d+'d.r4','r')
three_array=np.fromfile(three_file,dtype=np.float32)
three_file.close()

nframes = int(len(sw_array)/(N+1))
time=np.empty(nframes)
dsw=np.empty(nframes)
dgn=np.empty(nframes)

for frame in range(nframes):
   time[frame]=sw_array[frame*(N+1)]
   Zsw=sw_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Zgn=gn_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Z3d=three_array[frame*(N+1)+1:(frame+1)*(N+1)]
   dsw[frame]=np.sqrt(sum((Zsw-Z3d)**2)/sum(Z3d**2))
   dgn[frame]=np.sqrt(sum((Zgn-Z3d)**2)/sum(Z3d**2))

#=================================================================
# Set up figure:
fig, (ax1, ax2, ax3) = plt.subplots(figsize=[20,6], nrows=1, ncols=3)

tmax=time[-1]
ax1.set_xlim(0.0,tmax)
ax2.set_xlim(0.0,tmax)
ax3.set_xlim(0.0,tmax)

ax1.set_xlabel('$t$', fontsize=30)
ax2.set_xlabel('$t$', fontsize=30)
ax3.set_xlabel('$t$', fontsize=30)

ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')

yrange=1.5e2
ymult=1.0+0.05*np.log10(yrange)
ymax=ymult*max(np.amax(hsw),np.amax(hgn))
ymin=ymax/yrange
ax1.set_ylim(ymin,ymax)

yrange=3.6e3
ymult=1.0+0.05*np.log10(yrange)
ymax=ymult*max(np.amax(zsw),np.amax(zgn))
ymin=ymax/yrange
ax2.set_ylim(ymin,ymax)

yrange=15.0
ymult=1.0+0.05*np.log10(yrange)
ymax=ymult*max(np.amax(dsw),np.amax(dgn))
ymin=ymax/yrange
ax3.set_ylim(ymin,ymax)

ax1.set_ylabel('$\\langle{(\\Delta\\tilde{h})^2}\\rangle^{1/2}/\\langle{\\tilde{h}^2}\\rangle^{1/2}$', fontsize=30)
ax2.set_ylabel('$\\langle{(\\Delta\\zeta)^2}\\rangle^{1/2}/\\langle{\\zeta^2}\\rangle^{1/2}$', fontsize=30)
ax3.set_ylabel('$\\langle{(\\Delta\\delta)^2}\\rangle^{1/2}/\\langle{\\delta^2}\\rangle^{1/2}$', fontsize=30)

ax1.plot(time,hsw,'b-',lw=3, label='SW - 3D')
ax1.plot(time,hgn,'r-',lw=3, label='GN - 3D')

ax2.plot(time,zsw,'b-',lw=3)
ax2.plot(time,zgn,'r-',lw=3)

ax3.plot(time,dsw,'b-',lw=3)
ax3.plot(time,dgn,'r-',lw=3)

ax1.legend(loc='lower right',prop={'size':20}, shadow=True)

# Add information about the resolution in the first panel:
ax1.text( 0.06, 0.92, '$n = {x:.0f}$'.format(x=ng), transform=ax1.transAxes, fontsize=25)
ax1.text( 0.06, 0.84, '$H = {x:.1f}$'.format(x=hbar), transform=ax1.transAxes, fontsize=25)

#=========================================================================
# Save figure:
outfile='hzd_rmsdiff_H{x:.1f}'.format(x=hbar)+'n'+str(ng)+'.eps'
fig.savefig(outfile, format='eps', dpi=1200)

print
print ' To view the image, type'
print
print ' gv',outfile,'&'
print
