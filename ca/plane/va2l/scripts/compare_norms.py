#!/usr/bin/env python3

# This script plots the full, balanced or imbalanced rms field values 
# for h, zeta, delta and gamma from data in ?norms.asc in 4 separate 
# directories specified below.

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
# Select full, balanced or imbalanced results:
print ' The following three options are available:'
print
print ' (1) Compare full fields;'
print ' (2) Compare balanced fields;'
print ' (3) Compare imbalanced fields;'
print
option=int(raw_input('Option (default 1)? ') or 1)
print
ng=int(raw_input(' Resolution (default 256)? ') or 256)

#=================================================================
# Set up figure:
tmax=25.0
# tmax: maximum time to show

if option==1:
   height=14.6
elif option==2:
   height=14.6
else:
   height=14.2

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(figsize=[16,height], nrows=2, ncols=2)

ax1.set_xlim(0.0,tmax)
ax2.set_xlim(0.0,tmax)
ax3.set_xlim(0.0,tmax)
ax4.set_xlim(0.0,tmax)

ax3.set_xlabel('$t$', fontsize=30)
ax4.set_xlabel('$t$', fontsize=30)

#=================================================================
if option==1:
   outfile='norms_n'+str(ng)+'.eps'
   ax1.set_ylabel('$\\langle{h^2}\\rangle^{1/2}$', fontsize=30)
   ax2.set_ylabel('$\\langle{\\zeta^2}\\rangle^{1/2}$', fontsize=30)
   ax3.set_ylabel('$\\langle{\\delta^2}\\rangle^{1/2}$', fontsize=30)
   ax4.set_ylabel('$\\langle{\\gamma^2}\\rangle^{1/2}$', fontsize=30)

elif option==2:
   outfile='norms_bal_n'+str(ng)+'.eps'
   ax1.set_ylabel('$\\langle{h_b^2}\\rangle^{1/2}$', fontsize=30)
   ax2.set_ylabel('$\\langle{\\zeta_b^2}\\rangle^{1/2}$', fontsize=30)
   ax3.set_ylabel('$\\langle{\\delta_b^2}\\rangle^{1/2}$', fontsize=30)
   ax4.set_ylabel('$\\langle{\\gamma_b^2}\\rangle^{1/2}$', fontsize=30)

else:
   outfile='norms_imb_n'+str(ng)+'.eps'
   ax1.set_ylabel('$\\langle{h_i^2}\\rangle^{1/2}$', fontsize=30)
   ax2.set_ylabel('$\\langle{\\zeta_i^2}\\rangle^{1/2}$', fontsize=30)
   ax3.set_ylabel('$\\langle{\\delta_i^2}\\rangle^{1/2}$', fontsize=30)
   ax4.set_ylabel('$\\langle{\\gamma_i^2}\\rangle^{1/2}$', fontsize=30)

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
   # Open input data files and read data:

   in_file=open(dir+'hnorms.asc','r')
   time, dum1, dum2, dum3 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()
   if option==1:
      h=dum1
   elif option==2:
      h=dum2
   else:
      h=dum3

   in_file=open(dir+'znorms.asc','r')
   time, dum1, dum2, dum3 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()
   if option==1:
      z=dum1
   elif option==2:
      z=dum2
   else:
      z=dum3

   in_file=open(dir+'dnorms.asc','r')
   time, dum1, dum2, dum3 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()
   if option==1:
      d=dum1
   elif option==2:
      d=dum2
   else:
      d=dum3

   in_file=open(dir+'gnorms.asc','r')
   time, dum1, dum2, dum3 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()
   if option==1:
      g=dum1
   elif option==2:
      g=dum2
   else:
      g=dum3

   # Plot results in the appropriate panel:
   ax1.plot(time,h,dashes=dashlist[m],c=colorlist[m],lw=3)
   ax2.plot(time,z,dashes=dashlist[m],c=colorlist[m],lw=3,label=label_list[m])
   ax2.legend(loc='upper right',prop={'size':20}, shadow=True)
   ax3.plot(time,d,dashes=dashlist[m],c=colorlist[m],lw=3)
   ax4.plot(time,g,dashes=dashlist[m],c=colorlist[m],lw=3)

#=========================================================================
# Add information about the resolution in the first panel:
ax2.text( 0.74, 0.62, r'\textit{n} = '+str(ng), transform=ax2.transAxes, fontsize=25)

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
