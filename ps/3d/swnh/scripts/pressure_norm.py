#!/usr/bin/env python3

# This script plots the rms non-hydrostatic pressure divided by the
# rms total pressure (using g(h-H) as the hydrostatic part) for two
# resolutions, provided below.

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
# Set up figure:
fig = plt.figure(1,figsize=[8,7.5])
ax = fig.add_subplot(111)
ax.set_xlabel('$t$', fontsize=30)

#=================================================================
outfile='nh-pressure-rms.eps'
ax.set_ylabel('$\\langle{P_n}\\rangle^{1/2}/\\langle{P}\\rangle^{1/2}$', fontsize=30)

#=================================================================
# Resolutions to compare:
resol=np.array([256,512])
# Corresponding labels on the curves plotted:
label_list=['$n_g=256$','$n_g=512$']
# Corresponding colours (or shades of grey):
#colorlist=[(0.0,0.0,0.0),(0.4,0.4,0.4)]
# Corresponding line styles:
dashlist=[(8,3),(1,0.0001)]

#=================================================================
# Read data and plot for each resolution:
for m,ng in enumerate(resol):
   direc='../3d/swnh/ng'+str(ng)+'ld0.5r001/'

   inputfile = open(direc+'pn_rms.asc','r')
   t, pnrms = np.loadtxt(inputfile,dtype=float,unpack=True)
   inputfile.close()

#   ax.plot(t,pnrms,c=colorlist[m],label=label_list[m])
   ax.plot(t,pnrms,dashes=dashlist[m],c='k',label=label_list[m])
   ax.legend(loc='upper left',prop={'size':20}, shadow=True)

ax.set_xlim(0.0,t[-1])

# Save figure:
fig.savefig(outfile, format='eps', dpi=1200)

print
print ' To view the image, type'
print
print ' gv ',outfile
print
