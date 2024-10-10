#!/usr/bin/env python

# This script plots the balanced profiles of d, h, u & zeta for all data
# in a chosen directory (see default parameters below)

#========== Perform the generic imports =========
import warnings,os,sys
import numpy as np
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
# Set up figure:
fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(figsize=[24.0,15.2], nrows=2, ncols=3)

ax1.set_xlabel('$\phi-a$', fontsize=30)
ax1.set_ylabel('$a$', fontsize=30)
ax1.set_ylim(-np.pi/2.0,np.pi/2.0)
#ax1.xaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

ax2.set_xlabel('$h/H$', fontsize=30)
ax2.set_ylabel('$a$', fontsize=30)
ax2.set_ylim(-np.pi/2.0,np.pi/2.0)
#ax2.xaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

ax3.set_xlabel('$u/\Omega$', fontsize=30)
ax3.set_ylabel('$a$', fontsize=30)
ax3.set_ylim(-np.pi/2.0,np.pi/2.0)
#ax3.xaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

ax4.set_xlabel('$\zeta\sqrt{1-A}$', fontsize=30)
ax4.set_ylabel('$\phi$', fontsize=30)
ax4.set_ylim(-np.pi/2.0,np.pi/2.0)
#ax4.xaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

ax5.set_xlabel('$h/H$', fontsize=30)
ax5.set_ylabel('$\phi$', fontsize=30)
ax5.set_ylim(-np.pi/2.0,np.pi/2.0)
#ax5.xaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

ax6.set_xlabel('$u/\Omega$', fontsize=30)
ax6.set_ylabel('$\phi$', fontsize=30)
ax6.set_ylim(-np.pi/2.0,np.pi/2.0)
#ax6.xaxis.set(major_locator=MultipleLocator(0.01),
#                 major_formatter=FormatStrFormatter('%1.2f'))

#=================================================================
print()
gamma_in = input(' Enter gamma = 1/L_d (default: 2.0)? ')
gamma = float(gamma_in or 2.0)
print()
print( ' Processing ...')
print()

if gamma==1:
   gamma=1.0
elif gamma==2:
   gamma=2.0
elif gamma==5:
   gamma=5.0
elif gamma==10:
   gamma=10.0
elif gamma==20:
   gamma=20.0
elif gamma==50:
   gamma=50.0
elif gamma==100:
   gamma=100.0
elif gamma==200:
   gamma=200.0

# Default width:
wid=0.01

# Values of A:
amp_list=['0.2','0.4','0.6','0.8','0.999']

# Corresponding labels:
label_list=['$A=0.2$','$A=0.4$','$A=0.6$','$A=0.8$','$A=0.999$']

# Corresponding colours:
col=['r','b','c','m','k','g']

# Open input files and read data:
for m,amp in enumerate(amp_list):
   data='results/g'+str(gamma)+'/g'+str(gamma)+'w'+str(wid)+'A'+str(amp)+'.asc'
   in_file=open(data,'r')
   y, x1, x2, x3, x4 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   ax1.plot(x1,y,c=col[m],lw=2)
   ax2.plot(x3,y,c=col[m],lw=2)
   ax3.plot(x2,y,c=col[m],lw=2)

   y=y+x1
   x4=np.array(x4)*np.sqrt(1.0-float(amp))
   
   ax4.plot(x4,y,c=col[m],lw=2)
   ax5.plot(x3,y,c=col[m],lw=2,label=label_list[m])
   ax6.plot(x2,y,c=col[m],lw=2)
   
   ax5.legend(loc='lower left',prop={'size':25})

# Save image:
fig.subplots_adjust(wspace=0.1, hspace=0.1)

ofile='g'+str(gamma)+'w'+str(wid)+'.eps'
fig.savefig(ofile,  format='eps', dpi=300)

print(' Type:')
print()
print(' gv',ofile)
print()
