#!/usr/bin/env python3

# This script plots spectra for either r, (u,v), zeta or delta,
# and compares the 3D decomposed spectra with the 2D SW and GN spectra.

# Specify the 3D, SW and GN directories below

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
# Specify directories:
dir3d='../3d/swnh/hbar0.4ng256ld0.5/'
dirsw='../plane/sw/ng256/'
dirgn='../plane/gn/hbar0.4ng256ld0.5/'

print ' Plot spectra for one of the following fields:'
print
print '   (1) r;    (2) zeta;    or (3) delta.'
print

iopt=int(raw_input(' Choice (default 3): ') or 3)
print

hbar=float(raw_input('Value of H used in the 3D & GN simulations (default 0.4)? ') or 0.4)
print
kgn=2.0*np.pi/hbar
xgn=np.log10(kgn)

# Open input files:
if iopt == 1:
   in_file3d=open(dir3d+'rspec.asc','r')
   in_filesw=open(dirsw+'alt-spectra.asc','r')
   in_filegn=open(dirgn+'alt-spectra.asc','r')
   col2d=1
   ylimits=np.array([-15.0,-3.0])
   field='rho'
elif iopt == 2:
   in_file3d=open(dir3d+'zspec.asc','r')
   in_filesw=open(dirsw+'spectra.asc','r')
   in_filegn=open(dirgn+'spectra.asc','r')
   col2d=1
   ylimits=np.array([-12.0,0.0])
   field='zeta'
else:
   in_file3d=open(dir3d+'dspec.asc','r')
   in_filesw=open(dirsw+'spectra.asc','r')
   in_filegn=open(dirgn+'spectra.asc','r')
   col2d=2
   ylimits=np.array([-16.0,-4.0])
   field='delta'

t=float(raw_input('Time to show (default 25)? ') or 25.0)
print

# Read the first header line to get kmax (number of rows):
first_line=in_filesw.readline()
kmax=int(first_line.split()[-1])
in_filesw.seek(0)

ng=kmax
xmax=np.log10(2.0*float(ng)/3.0)
dx=0.1
xmax=dx*(float(int(xmax/dx)+1))
xlimits=np.array([0.0,xmax])

#=================================================================
# Set up figure:
fig = plt.figure(1,figsize=[8,7.5])
ax = plt.axes([0.2, 0.2, 0.7, 0.7])
ax.set_xlim(xlimits)
ax.set_ylim(ylimits)
ax.set_xlabel('$\\log_{10}k$', fontsize=30)
if iopt == 1:
   ax.set_ylabel('$\\log_{10}S_{\\tilde{\\rho}_\\theta}$', fontsize=30)
elif iopt == 2:
   ax.set_ylabel('$\\log_{10}S_{\\zeta}$', fontsize=30)
else:
   ax.set_ylabel('$\\log_{10}S_{\\delta}$', fontsize=30)

#=================================================================
# Read 3D data:
raw_data = np.fromfile(file=in_file3d,dtype=float,sep='\n')
in_file3d.close()

# Determine the number of frames:
nframes = int(len(raw_data)/(4*kmax+2))  

# Shape the data array into a useful shape for plotting:
frames=range(0,nframes)
time=[raw_data[i*(4*kmax+2)] for i in frames]
dt=time[1]-time[0]
tim_eles = [i*(4*kmax+2)+j for i in frames for j in range(2)]
shaped_data = np.delete(raw_data,tim_eles)[0:(4*kmax+2)*nframes].reshape((nframes,kmax,4))
k=np.zeros((nframes,ng))
z3d1=np.zeros((nframes,ng))
z3d2=np.zeros((nframes,ng))
z3d3=np.zeros((nframes,ng))
for i in frames:
   k[i,:]=shaped_data[i].transpose()[0][0:ng]
   z3d1[i,:]=shaped_data[i].transpose()[1][0:ng]
   z3d2[i,:]=shaped_data[i].transpose()[2][0:ng]
   z3d3[i,:]=shaped_data[i].transpose()[3][0:ng]

#=================================================================
# Read SW data:
raw_data = np.fromfile(file=in_filesw,dtype=float,sep='\n')
in_filesw.close()

# Shape the data array into a useful shape for plotting:
shaped_data = np.delete(raw_data,tim_eles)[0:(4*kmax+2)*nframes].reshape((nframes,kmax,4))
zsw=np.zeros((nframes,ng))
for i in frames:
   zsw[i,:]=shaped_data[i].transpose()[col2d][0:ng]

#=================================================================
# Read GN data:
raw_data = np.fromfile(file=in_filegn,dtype=float,sep='\n')
in_filegn.close()

# Shape the data array into a useful shape for plotting:
shaped_data = np.delete(raw_data,tim_eles)[0:(4*kmax+2)*nframes].reshape((nframes,kmax,4))
zgn=np.zeros((nframes,ng))
for i in frames:
   zgn[i,:]=shaped_data[i].transpose()[col2d][0:ng]

#=================================================================
# Plot chosen time:
ic=int(t/dt+0.01)

if iopt == 1:
   ax.plot(k[ic],z3d1[ic],'b-',lw=2,label='$\\tilde{h}_{3D}$')
   ax.plot(k[ic],z3d2[ic],'r-',lw=2,label='$\\tilde{\\rho}_{\\theta}$')
   ax.plot(k[ic],z3d3[ic],'m-',lw=2,label='${\\rho}^\\prime_{\\theta}$')
   ax.plot(k[ic], zsw[ic],'k-',lw=2,label='$\\tilde{h}_{SW}$')
   ax.plot(k[ic], zgn[ic],'g-',lw=2,label='$\\tilde{h}_{GN}$')
elif iopt == 2:
   ax.plot(k[ic],z3d1[ic],'b-',lw=2,label='$\\bar{\\zeta}_{3D}$')
   ax.plot(k[ic],z3d2[ic],'r-',lw=2,label='$\\zeta$')
   ax.plot(k[ic],z3d3[ic],'m-',lw=2,label='$\\zeta^\\prime$')
   ax.plot(k[ic], zsw[ic],'k-',lw=2,label='$\\bar{\\zeta}_{SW}$')
   ax.plot(k[ic], zgn[ic],'g-',lw=2,label='$\\bar{\\zeta}_{GN}$')
else:
   ax.plot(k[ic],z3d1[ic],'b-',lw=2,label='$\\bar{\\delta}_{3D}$')
   ax.plot(k[ic],z3d2[ic],'r-',lw=2,label='$\\delta$')
   ax.plot(k[ic],z3d3[ic],'m-',lw=2,label='$\\delta^\\prime$')
   ax.plot(k[ic], zsw[ic],'k-',lw=2,label='$\\bar{\\delta}_{SW}$')
   ax.plot(k[ic], zgn[ic],'g-',lw=2,label='$\\bar{\\delta}_{GN}$')

ax.legend(loc='lower left',prop={'size':20}, shadow=True)
ax.axvline(xgn,color='k',ls='--',lw=1)

ax.set_title('$H = {x:.1f}$'.format(x=hbar)+'$\ \ \ \ \ n = {x:.0f}$'.format(x=ng)+'$\ \ \ \ \ t = {x:.0f}$'.format(x=t), fontsize=30)

#=========================================================================
# Save figure:
outfile='zeta_hole_t{x:.0f}_zoom.eps'.format(x=t)

outfile=field+'_spec_H{x:.1f}'.format(x=hbar)+'n'+str(ng)+'t{x:.0f}'.format(x=t)+'.eps'
fig.savefig(outfile, format='eps', dpi=1200)

print
print ' To view the image, type'
print
print ' gv',outfile,'&'
print
