#!/usr/bin/env python3

# This script plots spectra for the height anomaly from data in 
# 4 separate directories specified below.

#========== Perform the generic imports =========
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
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
# Set up figure:
fig1 = plt.figure(1,figsize=[10,10])
ax1 = plt.axes([0.2, 0.2, 0.7, 0.7])
ax1.set_xlabel('$\log_{10}k$', fontsize=30)
ax1.set_xlim(0.0,2.25)

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
t=float(raw_input('Time to show (default 25)? ') or 25.0)
print
hbar=float(raw_input('Value of H used in GN simulations (default 0.2)? ') or 0.2)
print
ng=int(raw_input(' Resolution (default 256)? ') or 256)
kgn=np.sqrt(3.0)/hbar
xgn=np.log10(kgn)

#=================================================================
if option==1:
   suffix='_n'+str(ng)+'_t'+str(int(t+0.01))+'.png'
   datafile='alt-spectra.asc'
   ax1.set_ylabel('$\log_{10}k|\\hat{h}|^2$', fontsize=30)
   print
   # Range in log_10 spectra to show:
   yrange=12.0
   ymax=float(raw_input('Maximum value in log_10   h   spectrum to show (default -3)? ') or -3.0)
   ax1.set_ylim(ymax-yrange,ymax)

elif option==2:
   suffix='_n'+str(ng)+'_bal_t'+str(int(t+0.01))+'.png'
   datafile='alt-bspectra.asc'
   ax1.set_ylabel('$\log_{10}k|\\hat{h}_b|^2$', fontsize=30)
   print
   # Range in log_10 spectra to show:
   yrange=12.0
   ymax=float(raw_input('Maximum value in log_10   h   spectrum to show (default -3)? ') or -3.0)
   ax1.set_ylim(ymax-12.0,ymax)

else:
   suffix='_n'+str(ng)+'_imb_t'+str(int(t+0.01))+'.png'
   datafile='alt-ispectra.asc'
   ax1.set_ylabel('$\log_{10}k|\\hat{h}_i|^2$', fontsize=30)
   print
   # Range in log_10 spectra to show:
   yrange=12.0
   ymax=float(raw_input('Maximum value in log_10   h   spectrum to show (default -8)? ') or -8.0)
   ax1.set_ylim(ymax-yrange,ymax)

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
   # Open input file:
   in_file=open(dir+datafile,'r')
   # Read the first header line to get kmax:
   first_line=in_file.readline()
   kmax=int(first_line.split()[-1])
   in_file.seek(0)

   nx=kmax

   # Read in the full data to a 1d array and close input file:
   raw_data = np.fromfile(file=in_file,dtype=float,sep='\n')
   in_file.close()

   # Determine the number of frames:
   nframes = int(len(raw_data)/(4*kmax+2))  
   print 'Number of frames found %i' %nframes

   # Shape the data array into a useful shape for plotting:
   frames=range(0,nframes)
   time=[raw_data[i*(4*kmax+2)] for i in frames]
   dt=time[1]-time[0]
   tim_eles = [i*(4*kmax+2)+j for i in frames for j in range(2)]
   shaped_data = np.delete(raw_data,tim_eles)[0:(4*kmax+2)*nframes].reshape((nframes,kmax,4))
   k=np.zeros((nframes,nx))
   h=np.zeros((nframes,nx))
   for i in frames:
      k[i,:]=shaped_data[i].transpose()[0][0:nx]
      h[i,:]=shaped_data[i].transpose()[1][0:nx]

   # Select frame:
   ic=int(t/dt+0.01)
   print
   print ' Showing results at t = ',time[ic]

   ax1.plot(k[ic],h[ic],dashes=dashlist[m],c=colorlist[m],lw=3,label=label_list[m])
   ax1.legend(loc='lower left',prop={'size':20}, shadow=True)
   ax1.axvline(xgn,color='k',ls='--',lw=1)

#=========================================================================
# Save image:
fig1.savefig('h'+suffix,  bbox_inches='tight', pad_inches = 0.1, dpi=600)

plt.show()
