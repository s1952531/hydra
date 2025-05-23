#!/usr/bin/env python3

#=====================================================================
#                Animates data created by pam.f90
#=====================================================================

#=====perform various generic imports=====
import sys,os,warnings
import numpy as np

from matplotlib import pyplot as plt
import matplotlib.animation as anim

import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

# set tick label size:
label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
# set x tick width and size:
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
# set y tick width and size:
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 2
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
# set axes width:
mpl.rcParams['axes.linewidth'] = 2
#=========================================

opt_in = input(' Animate (1) h,  (2) u,  or (3) v (default: 1)? ')
opt=int(opt_in or 1)

if opt == 1:
   field='h'
   xlabel='$h$'
elif opt == 2:
   field='u'
   xlabel='$u$'
elif opt == 3:
   field='v'
   xlabel='$v$'
else:
   print (' Not a valid option.')
   exit

# Get number of phi intervals from length of init.asc:
with open('init.asc') as init_file:
   n = 0
   for line in init_file:
      n += 1

# Read data:
datafile=field+'.r8'
file_bytes = os.path.getsize(datafile)
nframes = int(file_bytes/((n+1)*8))
print (' Number of frames generated: %d' %nframes)
frames = range(nframes)
in_file = open(datafile,'r')
raw_array = np.fromfile(in_file,dtype=np.float32)
x = np.empty((nframes,n))
for i in frames:
   x[i,:] = raw_array[i*(n+1)+1:(i+1)*(n+1)]
in_file.close()    
xmax=np.amax(abs(x))

datafile='phi.r8'
in_file = open(datafile,'r')
raw_array = np.fromfile(in_file,dtype=np.float32)
y = np.empty((nframes,n))
for i in frames:
   y[i,:] = raw_array[i*(n+1)+1:(i+1)*(n+1)]
in_file.close()    

time=[raw_array[i*(n+1)] for i in frames]

# Create figure:
fig = plt.figure(figsize=[8.0,8.0])

# Add axes to figure
ax = fig.add_subplot(111,autoscale_on=False)

# Fix domain of view:
ax.set_xlim([-xmax,xmax])
ax.set_ylim([-np.pi/2.0,np.pi/2.0])

ax.set_xlabel(xlabel, fontsize=30)
ax.set_ylabel('$\\phi$', fontsize=30)

# Formatting for displayed time:
time_template = 't = %.3f'

# Set up list of images for animation:
ims=[]
for i,t in enumerate(time):
   # On each characteristic starting from x = x0, u remains u0:
   x0=x[i,:]
   y0=y[i,:]

   cims = ax.plot(x0,y0,c='k',lw=2)

   cims.append(ax.text(0.72, 0.96, time_template%(t), transform=ax.transAxes, weight='bold', size=20))

   ims.append(cims)

# Run animation:
ani = anim.ArtistAnimation(fig, ims, interval=40, repeat_delay=400, blit=False)

#print ()
#print (' Creating movie (evo.mp4) ... (this can take some time)')
#ani.save('evo.mp4', fps=10, extra_args=['-vcodec', 'libx264'])
#
#print ()
#print (' Display the movie by typing')
#print ()
#print (' mplayer evo.mp4 -loop 0')
#print ()

plt.show()
