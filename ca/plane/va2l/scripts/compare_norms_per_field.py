#!/usr/bin/env python

# This script plots the full, balanced and imbalanced field values 
# for h, zeta, delta or gamma from data in ?norm.asc in 4 separate 
# directories specified below.

#========== Perform the generic imports =========
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

## global settings

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
ax1.set_xlim(0.0,tmax)
#ax1.yaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

ax2.set_xlabel('$t$', fontsize=30)
ax2.set_xlim(0.0,tmax)
#ax2.yaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

ax3.set_xlabel('$t$', fontsize=30)
ax3.set_xlim(0.0,tmax)
#ax3.yaxis.set(major_locator=MultipleLocator(0.001),
#                 major_formatter=FormatStrFormatter('%1.3f'))

#=================================================================
# Select data to compare:
print
print ' The following fields may be compared:'
print
print ' (1) h;'
print ' (2) zeta;'
print ' (3) delta;'
print ' (4) gamma.'
print
k=int(raw_input(' Option (default 1)? ') or 1)

ng=int(raw_input(' Resolution (default 256)? ') or 256)

if k==1:
   datafile='hnorms.asc'
   outfile='h_n'+str(ng)+'.eps'
   ax1.set_ylabel('$\\langle{h}\\rangle$', fontsize=30)
   ax2.set_ylabel('$\\langle{h_b}\\rangle$', fontsize=30)
   ax3.set_ylabel('$\\langle{h_i}\\rangle$', fontsize=30)
   legloc='lower right'
elif k==2:
   datafile='znorms.asc'
   outfile='z_n'+str(ng)+'.eps'
   ax1.set_ylabel('$\\langle{\\zeta}\\rangle$', fontsize=30)
   ax2.set_ylabel('$\\langle{\\zeta_b}\\rangle$', fontsize=30)
   ax3.set_ylabel('$\\langle{\\zeta_i}\\rangle$', fontsize=30)
   legloc='upper right'
elif k==3:
   datafile='dnorms.asc'
   outfile='d_n'+str(ng)+'.eps'
   ax1.set_ylabel('$\\langle{\\delta}\\rangle$', fontsize=30)
   ax2.set_ylabel('$\\langle{\\delta_b}\\rangle$', fontsize=30)
   ax3.set_ylabel('$\\langle{\\delta_i}\\rangle$', fontsize=30)
   legloc='lower right'
else:
   datafile='gnorms.asc'
   outfile='g_n'+str(ng)+'.eps'
   ax1.set_ylabel('$\\langle{\\gamma}\\rangle$', fontsize=30)
   ax2.set_ylabel('$\\langle{\\gamma_b}\\rangle$', fontsize=30)
   ax3.set_ylabel('$\\langle{\\gamma_i}\\rangle$', fontsize=30)
   legloc='lower right'

#=================================================================
# List of directories to compare (need the final /):
dirend=str(ng)+'/'
dir_list=['sw/ng'+dirend,'sw/bal_ng'+dirend,'gn/ng'+dirend,'gn/bal_ng'+dirend]
# Corresponding labels on the curves plotted:
label_list=['SW','SW-bal','GN','GN-bal']
# Corresponding line styles:
dashlist=[(1,0.0001),(8,3),(1,0.0001),(8,3)]
# Corresponding colours (or shades of grey):
colorlist=[(0.4,0.4,0.4),(0.4,0.4,0.4),(0.0,0.0,0.0),(0.0,0.0,0.0)]

#=================================================================
# Loop over directories and plot results:
for m,dir in enumerate(dir_list):
   # Open input file and read data:
   in_file=open(dir+datafile,'r')
   time, rms, brms, irms = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   ax1.plot(time,rms,dashes=dashlist[m],c=colorlist[m],lw=3)

   ax2.plot(time,brms,dashes=dashlist[m],c=colorlist[m],lw=3,label=label_list[m])
   ax2.legend(loc=legloc,prop={'size':20}, shadow=True)

   ax3.plot(time,irms,dashes=dashlist[m],c=colorlist[m],lw=3)

#=========================================================================
# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0)

fig.savefig(outfile, format='eps', dpi=1200)

print
print ' To view the image, type'
print
print ' gv ',outfile
print


