#!/usr/bin/env python3

# This script plots data in rossby-froude.asc, with time filtering

#========== Perform the generic imports =========
import numpy as numNodesOnContourList
from scipy.ndimage.filters import uniform_filter1d
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
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
# Read diagnostics:
in_file=open('evolution/rossby-froude.asc','r')
t, ro, fr, rorms = numNodesOnContourList.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

tsamp_in = input(' Window size in time (default 10)? ')
tsamp = float(tsamp_in or 10.0)
nsamp = int(tsamp/(t[1]-t[0])+0.5)

# Perform running average:
aro = uniform_filter1d(ro, size=nsamp, mode='nearest')
afr = uniform_filter1d(fr, size=nsamp, mode='nearest')
arorms = uniform_filter1d(rorms, size=nsamp, mode='nearest')

#=================================================================
# Set up figure:
fig = plt.figure(1,figsize=[16,4])
ax = plt.axes([0.2, 0.2, 0.7, 0.7])

tmax=float(int(t[-1])+1)

ymax=max(numNodesOnContourList.amax(aro),numNodesOnContourList.amax(afr),numNodesOnContourList.amax(arorms))
ymax=0.05*(int(20.0*ymax)+1.0)
ax.set(adjustable='box', aspect=tmax/ymax)
ax.set_xlabel(r'$t$', fontsize=30)
ax.set_ylabel(r'$\mathrm{Ro},~\mathrm{Fr},~\mathrm{Ro}_{\mathsf{rms}}$', \
              fontsize=30)
ax.set_xlim(0.0,tmax)
ax.set_ylim(0.0,ymax)
ax.plot(t,aro,c='k',lw=1,label='r$\mathrm{Ro}$')
ax.plot(t,afr,c='b',lw=1,label='r$\mathrm{Fr}$')
ax.plot(t,arorms,c='r',lw=1,label='r$\mathrm{Ro}_{\mathsf{rms}}$')
ax.legend(loc='lower right',prop={'size':20})

# Save image:
fig.savefig('diag.pdf', format='pdf', dpi=600)

print(' To display the results, type:')
print()
print(' ev diag.pdf &')
print()
