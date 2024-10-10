#!/usr/bin/env python

# This script plots the rms non-hydrostatic acceleration divided by the
# rms acceleration for two resolutions and for both the 3D and GN results.

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

# Adjust height to get square images:
height=7.9

fig, (ax1, ax2) = plt.subplots(figsize=[16,height], nrows=1, ncols=2)
ax1.set_xlabel('$t$', fontsize=30)
ax2.set_xlabel('$t$', fontsize=30)

#=================================================================
outfile='nh-accel-rms.eps'
ax1.set_ylabel('$\\langle{\\|\\mathbf{a}_n\\|}\\rangle^{1/2}/\\langle{\\|\\mathbf{a}\\|}\\rangle^{1/2}$', fontsize=30)
#ax2.set_ylabel('$\\langle{\\|\\mathbf{a}_n\\|}\\rangle^{1/2}/\\langle{\\|\\mathbf{a}\\|}\\rangle^{1/2}$', fontsize=30)

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

   inputfile = open(direc+'an_rms.asc','r')
   t, anrms = np.loadtxt(inputfile,dtype=float,unpack=True)
   inputfile.close()

#   ax1.plot(t,anrms,c=colorlist[m],label=label_list[m])
   ax1.plot(t,anrms,dashes=dashlist[m],c='k',label=label_list[m])

   ax1.legend(loc='upper left',prop={'size':20}, shadow=True)

   direc='../plane/gn/gnps/vstrip/ng'+str(ng)+'ld0.5r001/'

   inputfile = open(direc+'an_rms.asc','r')
   t, anrms = np.loadtxt(inputfile,dtype=float,unpack=True)
   inputfile.close()

#   ax2.plot(t,anrms,c=colorlist[m],label=label_list[m])
   ax2.plot(t,anrms,dashes=dashlist[m],c='k',label=label_list[m])
   ax2.legend(loc='upper left',prop={'size':20}, shadow=True)

ax1.set_xlim(0.0,t[-1])
ax2.set_xlim(0.0,t[-1])

# Add information about the model in each panel:
ax1.text( 0.14, 0.75, '3D', transform=ax1.transAxes, fontsize=25)
ax2.text( 0.14, 0.75, 'GN', transform=ax2.transAxes, fontsize=25)

# Save figure:
fig.savefig(outfile, format='eps', dpi=1200)

print
print ' To view the image, type'
print
print ' gv ',outfile
print
