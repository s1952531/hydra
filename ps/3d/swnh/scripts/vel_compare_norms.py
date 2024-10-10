#!/usr/bin/env python

# This script plots the rms differences between SW and 3D, and between
# GN and 3D, for the fields of delta, zeta, u & v.

#========== Perform the generic imports =========
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)

# set tick label size:
label_size = 24
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
# set x tick width and size:
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 5
mpl.rcParams['xtick.minor.width'] = 1
# set y tick width and size:
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 5
mpl.rcParams['ytick.minor.width'] = 1
# set axes width:
mpl.rcParams['axes.linewidth'] = 3

#=================================================================
# Select resolution:
print
ng=int(raw_input(' Resolution (default 256)? ') or 256)

#=================================================================
# Set up figure:
tmax=25.0
# tmax: maximum time to show

# Adjust height to get square images:
height=14.6

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(figsize=[16,height], nrows=2, ncols=2)

ax1.set_xlim(0.0,tmax)
ax2.set_xlim(0.0,tmax)
ax3.set_xlim(0.0,tmax)
ax4.set_xlim(0.0,tmax)

ax3.set_xlabel('$t$', fontsize=30)
ax4.set_xlabel('$t$', fontsize=30)

#=================================================================
outfile='vel_norms_n'+str(ng)+'.eps'
ax1.set_ylabel(r'$\langle{\Delta{u}^2}\rangle^{1/2}$', fontsize=30)
ax2.set_ylabel(r'$\langle{\Delta{v}^2}\rangle^{1/2}$', fontsize=30)
ax3.set_ylabel(r'$\langle{\Delta\gamma^2}\rangle^{1/2}$', fontsize=30)
ax4.set_ylabel(r'$\langle{\Delta{\tilde{h}}^2}\rangle^{1/2}$', fontsize=30)

#=================================================================
# List of directories to compare (need the final /):
direc=['../plane/sw/swps/vstrip/ng'+str(ng)+'ld0.5r001/','../plane/gn/gnps/vstrip/ng'+str(ng)+'ld0.5r001/','../3d/swnh/ng'+str(ng)+'ld0.5r001/']
# Corresponding labels on the curves plotted:
label_list=['SW -- 3D','GN -- 3D']
# Corresponding colours (or shades of grey):
colorlist=[(0.0,0.0,0.0),(0.4,0.4,0.4)]

#=================================================================
# Read field data and compute rms differences:
N=ng*ng
Zsw=np.empty(N)
Zgn=np.empty(N)
Z3d=np.empty(N)
sumfac=1.0/float(N)

# x velocity:
sw_file=open(direc[0]+'uu.r4','r')
sw_array=np.fromfile(sw_file,dtype=np.float32)
sw_file.close()

gn_file=open(direc[1]+'uu.r4','r')
gn_array=np.fromfile(gn_file,dtype=np.float32)
gn_file.close()

three_file=open(direc[2]+'2d/u.r4','r')
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
   dsw[frame]=np.sqrt(sumfac*sum((Zsw-Z3d)**2))
   dgn[frame]=np.sqrt(sumfac*sum((Zgn-Z3d)**2))

# y velocity:
sw_file=open(direc[0]+'vv.r4','r')
sw_array=np.fromfile(sw_file,dtype=np.float32)
sw_file.close()

gn_file=open(direc[1]+'vv.r4','r')
gn_array=np.fromfile(gn_file,dtype=np.float32)
gn_file.close()

three_file=open(direc[2]+'2d/v.r4','r')
three_array=np.fromfile(three_file,dtype=np.float32)
three_file.close()

nframes = int(len(sw_array)/(N+1))
zsw=np.empty(nframes)
zgn=np.empty(nframes)

for frame in range(nframes):
   Zsw=sw_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Zgn=gn_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Z3d=three_array[frame*(N+1)+1:(frame+1)*(N+1)]
   zsw[frame]=np.sqrt(sumfac*sum((Zsw-Z3d)**2))
   zgn[frame]=np.sqrt(sumfac*sum((Zgn-Z3d)**2))

# Gamma:
sw_file=open(direc[0]+'gg.r4','r')
sw_array=np.fromfile(sw_file,dtype=np.float32)
sw_file.close()

gn_file=open(direc[1]+'gg.r4','r')
gn_array=np.fromfile(gn_file,dtype=np.float32)
gn_file.close()

three_file=open(direc[2]+'2d/g.r4','r')
three_array=np.fromfile(three_file,dtype=np.float32)
three_file.close()

nframes = int(len(sw_array)/(N+1))
gsw=np.empty(nframes)
ggn=np.empty(nframes)

for frame in range(nframes):
   Zsw=sw_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Zgn=gn_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Z3d=three_array[frame*(N+1)+1:(frame+1)*(N+1)]
   gsw[frame]=np.sqrt(sumfac*sum((Zsw-Z3d)**2))
   ggn[frame]=np.sqrt(sumfac*sum((Zgn-Z3d)**2))

# Dimensionless depth anomaly:
sw_file=open(direc[0]+'hh.r4','r')
sw_array=np.fromfile(sw_file,dtype=np.float32)
sw_file.close()

gn_file=open(direc[1]+'hh.r4','r')
gn_array=np.fromfile(gn_file,dtype=np.float32)
gn_file.close()

three_file=open(direc[2]+'2d/h.r4','r')
three_array=np.fromfile(three_file,dtype=np.float32)
three_file.close()

nframes = int(len(sw_array)/(N+1))
wsw=np.empty(nframes)
wgn=np.empty(nframes)

for frame in range(nframes):
   Zsw=sw_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Zgn=gn_array[frame*(N+1)+1:(frame+1)*(N+1)]
   Z3d=three_array[frame*(N+1)+1:(frame+1)*(N+1)]
   wsw[frame]=np.sqrt(sumfac*sum((Zsw-Z3d)**2))
   wgn[frame]=np.sqrt(sumfac*sum((Zgn-Z3d)**2))

#=========================================================================
# Plot results in the appropriate panel:
ax1.plot(time,dsw,c=colorlist[0],lw=3)
ax1.plot(time,dgn,c=colorlist[1],lw=3)

ax2.plot(time,zsw,c=colorlist[0],lw=3,label=label_list[0])
ax2.plot(time,zgn,c=colorlist[1],lw=3,label=label_list[1])

ax2.legend(loc='upper left',prop={'size':20}, shadow=True)

ax3.plot(time,gsw,c=colorlist[0],lw=3)
ax3.plot(time,ggn,c=colorlist[1],lw=3)

ax4.plot(time,wsw,c=colorlist[0],lw=3)
ax4.plot(time,wgn,c=colorlist[1],lw=3)

# Add information about the resolution in the first panel:
ax2.text( 0.10, 0.75, r'\textit{n} = '+str(ng), transform=ax2.transAxes, fontsize=25)

# Share x axes in top and bottom panels:
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)

# Add spacing between panels:
plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=1.0)

# Save figure:
fig.savefig(outfile, format='eps', dpi=1200)

print
print ' To view the image, type'
print
print ' gv ',outfile
print
