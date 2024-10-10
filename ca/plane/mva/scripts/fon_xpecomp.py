#!/usr/bin/env python

# This script plots the kinetic, potential, magnetic and total energy
# versus time, with the possibility of comparing data in different
# directories as specified below.

#========== Perform the generic imports =========
import warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib import rc
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

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
# Specify data directories here:d
dir_list=['fon0.25/','fon1.0/','fon4.0/']
ndir=len(dir_list)

# Corresponding line styles:
dashlist=[(1,0.0001),(8,3),(1,3,3,1)]

# Corresponding Rossby numbers for plotting dimensionless time:
rossby=np.array([0.25,0.25,0.25])

# Corresponding values of f/N:
fon=np.array([0.25,1,4])

# Corresponding labels for the plots:
lab_list=['$f/N=1/4$','$f/N=1$','$f/N=4$']

# Set up figure:
fig, (ax1, ax2, ax3, ax4) = plt.subplots(figsize=[24.3,6], nrows=1, ncols=4)

# For working out maximum energy and time:
emax=0.0
emin=1.e20
ebmax=0.0
ebmin=0.0
ekmax=0.0
ekmin=1.e20
ehmax=0.0
ehmin=1.e20
tmax=0.0

# Read data and plot:
for m,direc in enumerate(dir_list):
   in_file=open(direc+'ecomp.asc','r')
   t, ekin, epot, emag, etot = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()
   t=rossby[m]*t
   epot=epot-epot[0]
   etot=ekin+epot+emag

   emin=min(emin,np.amin(etot))
   emax=max(emax,np.amax(etot))
   ebmax=max(ebmax,np.amax(emag))
   ekmin=min(ekmin,np.amin(ekin))
   ekmax=max(ekmax,np.amax(ekin))
   ehmin=min(ehmin,np.amin(epot))
   ehmax=max(ehmax,np.amax(epot))
   tmax=max(tmax,t[-1])

   ax1.plot(t,etot,c='k',lw=2,dashes=dashlist[m],label=lab_list[m])
   ax2.plot(t,epot,c='k',lw=2,dashes=dashlist[m])
   ax3.plot(t,ekin,c='k',lw=2,dashes=dashlist[m])
   ax4.plot(t,emag,c='k',lw=2,dashes=dashlist[m])

ax1.set_xlim(0.0,tmax)
ax2.set_xlim(0.0,tmax)
ax3.set_xlim(0.0,tmax)
ax4.set_xlim(0.0,tmax)
sfac=0.02
de=sfac*(emax-emin)
emin=emin-de
emax=emax+de
de=sfac*(ebmax-ebmin)
ebmax=ebmax+de
de=sfac*(ekmax-ekmin)
ekmin=ekmin-de
ekmax=ekmax+de
de=sfac*(ehmax-ehmin)
ehmin=ehmin-de
ehmax=ehmax+de
ax1.set_ylim(emin,emax)
ax2.set_ylim(ehmin,ehmax)
ax3.set_ylim(ekmin,ekmax)
ax4.set_ylim(ebmin,ebmax)

ax1.set(adjustable='box', aspect=tmax/(emax-emin))
ax2.set(adjustable='box', aspect=tmax/(ehmax-ehmin))
ax3.set(adjustable='box', aspect=tmax/(ekmax-ekmin))
ax4.set(adjustable='box', aspect=tmax/(ebmax-ebmin))
ax1.legend(loc='best',prop={'size':20})

ax1.set_xlabel('$\\varepsilon t$', fontsize=30)
ax1.set_ylabel('$\mathcal{E}$', fontsize=30)
ax2.set_xlabel('$\\varepsilon t$', fontsize=30)
ax2.set_ylabel('$\mathcal{E}_h$', fontsize=30)
ax3.set_xlabel('$\\varepsilon t$', fontsize=30)
ax3.set_ylabel('$\mathcal{E}_u$', fontsize=30)
ax4.set_xlabel('$\\varepsilon t$', fontsize=30)
ax4.set_ylabel('$\mathcal{E}_b$', fontsize=30)

# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0)

fig.savefig('fon_ecomp.eps',  format='eps', dpi=1200)

print(' To display the results, type:')
print()
print(' gv fon_ecomp.eps &')
print()
