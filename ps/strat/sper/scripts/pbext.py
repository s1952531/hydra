#!/usr/bin/env python

#=====================================================================
#   Plots data in bminmax.asc in several selected directories
#=====================================================================

#=====perform various generic imports=====
import warnings
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

## global settings

# set tick label size:
label_size = 25
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

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)
#=========================================

#=================================================================
# Specify the data directories:
dir_list=['1x','2x','4x','8x']
# Corresponding labels:
label_list=['1x','2x','4x','8x']
# Corresponding colours (allow up to 5 directories):
colorlist=['k','b','r','m','c']

#--------------------------------------------------------------
# Open bminmax.asc file in one directory to get the final time:
in_file=open(dir_list[0]+'/evolution/bminmax.asc','r')
time, bmin, bmax=np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

tmax=time[-1]

#--------------------------------------------------------------
# Set up figure:
fig1 = plt.figure(1,figsize=[8,8])
ax1 = plt.axes([0.22, 0.22, 0.72, 0.72])
#plt.gca().yaxis.set_major_locator( MaxNLocator(nbins = 4) )

ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('$b_{min},~b_{max}$', fontsize=30)

#Set the axes limits:
ax1.set_xlim((0, tmax))
#ax1.set_ylim((-7.0, 0.0))

# Loop over directories and plot results:
for m,dir in enumerate(dir_list):
   # Open input file and read data:
   in_file=open(dir+'/evolution/bminmax.asc','r')
   time, bmin, bmax = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   ax1.plot(time,bmin,c=colorlist[m],lw=2,label=label_list[m])
   ax1.legend(loc='center right',prop={'size':20})
   ax1.plot(time,bmax,c=colorlist[m],lw=2)

ax1.axhline(bmin[0],color='g',linestyle='--')
ax1.axhline(bmax[0],color='g',linestyle='--')
   
#=========================================================================
# Save image:
fig1.savefig('bext.eps',  format='eps', dpi=1200)

print()
print(' To view the image, type')
print()
print(' gv bext.eps &')
print()
