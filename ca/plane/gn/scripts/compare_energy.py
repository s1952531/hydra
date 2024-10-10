#!/usr/bin/env python

# This script plots the kinetic, potential and total energy from data
# in e.asc in 4 separate directories specified below.

#========== Perform the generic imports =========
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
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
# Set up figures:
tmax=25.0
# tmax: maximum time to show

# Set up figure:
fig, (ax1, ax2, ax3) = plt.subplots(figsize=[24.0,7.53], nrows=1, ncols=3)

ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('${\\mathcal{K}}$', fontsize=30)
ax1.set_xlim(0.0,tmax)
ax1.yaxis.set(major_locator=MultipleLocator(0.01),
                 major_formatter=FormatStrFormatter('%1.2f'))

ax2.set_xlabel('$t$', fontsize=30)
ax2.set_ylabel('${\\mathcal{P}}$', fontsize=30)
ax2.set_xlim(0.0,tmax)
ax2.yaxis.set(major_locator=MultipleLocator(0.01),
                 major_formatter=FormatStrFormatter('%1.2f'))

ax3.set_xlabel('$t$', fontsize=30)
ax3.set_ylabel('${\\mathcal{E}}$', fontsize=30)
ax3.set_xlim(0.0,tmax)
ax3.yaxis.set(major_locator=MultipleLocator(0.001),
                 major_formatter=FormatStrFormatter('%1.3f'))

#=================================================================
ng=int(raw_input(' Resolution (default 256)? ') or 256)

# List of directories to compare (need the final /):
dirend=str(ng)+'/'
dir_list=['sw/ng'+dirend,'sw/bal_ng'+dirend,'gn/ng'+dirend,'gn/bal_ng'+dirend]
# Corresponding labels on the curves plotted:
label_list=['SW','SW-bal','GN','GN-bal']
# Corresponding line styles:
dashlist=[(1,0.0001),(8,3),(1,0.0001),(8,3)]
# Corresponding colours (or shades of grey):
colorlist=[(0.4,0.4,0.4),(0.4,0.4,0.4),(0.0,0.0,0.0),(0.0,0.0,0.0)]
# Datafile name:
datafile='e.asc'

#=================================================================
# Loop over directories and plot results:
for m,dir in enumerate(dir_list):
   # Open input file and read data:
   in_file=open(dir+datafile,'r')
   time, dum1, dum2, ekin, epot, etot = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   ax1.plot(time,ekin,dashes=dashlist[m],c=colorlist[m],lw=3,label=label_list[m])
   ax1.legend(loc='upper right',prop={'size':20}, shadow=True)
   ax2.plot(time,epot,dashes=dashlist[m],c=colorlist[m],lw=3)
   ax3.plot(time,etot,dashes=dashlist[m],c=colorlist[m],lw=3)

#=========================================================================
# Add information about the resolution in the second panel:
ax1.text( 0.74, 0.64, r'\textit{n} = '+str(ng), transform=ax1.transAxes, fontsize=25)

# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0)

fig.savefig('energy'+str(ng)+'.eps',  format='eps', dpi=1200)


