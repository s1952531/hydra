#!/usr/bin/python

# This script plots the PV, streamfunction and vorticity in
# evolution/qq.r4, pp.r4 and zz.r4, for all layers or for all
# modes --- either at a single time or as a time average.

#  @@@@   Run from the current job directory   @@@@

#========== Perform the generic imports =========
import sys,os,warnings
import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
import matplotlib.cm as cm
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

# Ensure latex fonts throughout:
rc('font',**{'family': 'Times New Roman'})
rc('text',usetex=True)

# set tick label size:
label_size=20
mpl.rcParams['xtick.labelsize']=label_size
mpl.rcParams['ytick.labelsize']=label_size
# set x tick width and size:
mpl.rcParams['xtick.major.size']=5
mpl.rcParams['xtick.major.width']=2
mpl.rcParams['xtick.minor.size']=3
mpl.rcParams['xtick.minor.width']=1
# set y tick width and size:
mpl.rcParams['ytick.major.size']=5
mpl.rcParams['ytick.major.width']=2
mpl.rcParams['ytick.minor.size']=3
mpl.rcParams['ytick.minor.width']=1
# set axes width:
mpl.rcParams['axes.linewidth']=1

#====================== Function definitions =======================
def contint(fmin,fmax):
    #Determines a nice contour interval (giving 10-20 divisions with
    #interval 1, 2 or 5x10^m for some m) given the minimum & maximum
    #values of the field data (fmin & fmax).

    fmax=0.9999999*fmax
    fmin=0.9999999*fmin
    #The 0.99... factor avoids having a superfluous tick interval
    #in cases where fmax-fmin is 10^m or 2x10^m

    emag=1.0
    rmult=max(1.0E-12,fmax-fmin)
    while rmult < 10:
        emag=emag/10
        rmult=rmult*10

    while rmult >= 100:
        emag=emag*10
        rmult=rmult/10

    kmult=int(rmult/10)

    if kmult < 1:
        ci=emag
    elif kmult < 2:
        ci=2*emag
    elif kmult < 4:
        ci=4*emag
    elif kmult < 8:
        ci=10*emag
    else:
        ci=20*emag

    return ci

#=================================================================
# Work out grid resolution (nx, ny & nz) by reading parameters.f90:
in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: nx=' in line:
      pline=line

line=pline.split("=")[1]
nx=int(line.split(",")[0])
in_file.close()

in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: ny=' in line:
      pline=line

line=pline.split("=")[1]
ny=int(line.split(",")[0])
in_file.close()

in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: nz=' in line:
      pline=line

line=pline.split("=")[1]
nz=int(line.split(",")[0])
in_file.close()

# Increase nx & ny by 1 to include boundary points:
nx=nx+1
ny=ny+1

#-------------------------------------------------
# Work out x & y limits by reading parameters.f90:
in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: xmin=' in line:
      pline=line

line=pline.split("=")[1]
xmin=float(line.split(",")[0])
in_file.close()

in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: xmax=' in line:
      pline=line

line=pline.split("=")[1]
xmax=float(line.split(",")[0])
in_file.close()

in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: ymin=' in line:
      pline=line

line=pline.split("=")[1]
ymin=float(line.split(",")[0])
in_file.close()

in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: ymax=' in line:
      pline=line

line=pline.split("=")[1]
ymax=float(line.split(",")[0])
in_file.close()

#-----------------------------------------------------
# Work out data save interval:
in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: tgsave=' in line:
      pline=line

line=pline.split("=")[1]
dtsave=float(line.split(",")[0])
in_file.close()

#-----------------------------------------------------
# Read energy data to get final time in data:
in_file=open('evolution/energy.asc','r')
time, ekin, epot, etot = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()
tsim = time[-1]

#=================================================================
# Select viewing options:
print()
op_in = input(' View layers or modes (l/m) (default l)? ')
option = str(op_in or "l")

# Read in vertical mode matrix if required:
if option == "m":
    vec=np.empty((nz,nz))
    in_file=open('modes.asc','r')
    for m in range(nz):
        line=in_file.readline()
    for m in range(nz):
        for iz in range(nz):
            line=in_file.readline()
            string=line.split()
            vec[iz,m]=string[0]
    in_file.close()

q_in = input(' Starting time (default '+str(tsim)+')? ')
t1 = float(q_in or tsim)
frame1 = int(t1/dtsave+0.5)

q_in = input('   Ending time (default '+str(t1)+')? ')
t2 = float(q_in or t1)
frame2 = int(t2/dtsave+0.5)

if frame2 > frame1:
    div=1.0/float(frame2-frame1+1)

#=================================================================
# Set up figure:

# Layout is 3 rows and nz columns:
nrow=3
ncol=nz
nim=nrow*ncol

# Ensure square images:
aspect=(ymax-ymin)/(xmax-xmin)

fig, ax = plt.subplots(figsize=[15,float(nz)*(12.5*aspect/float(nrow)+0.5)], nrows=nrow, ncols=ncol)
ax = ax.flatten()

# Set up an array to store all images to be plotted:
d=np.empty((nim,nx,ny))

# Number of horizontal grid points:
NH=nx*ny

# Total number of data elements per time frame:
N=nz*NH+1
# The "+1" includes the time element

# Field titles over each column:
field=['$q$','$\\psi$','$\\zeta$']

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read PV data into array for plotting:
in_file=open('evolution/qq.r4','r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()

# Loop over layers:
for iz in range(nz):
    # Frame where image will appear (0 1 2 first row, 3 4 5 second, ...):
    k=3*iz
    if frame2 > frame1:
        # Time average:
        for frame in range(frame1,frame2):
            offset=frame*N
            d[k,:,:]+=raw_array[offset+iz*NH+1:offset+(iz+1)*NH+1].reshape(nx,ny)
        d[k,:,:]=div*d[k,:,:]

    else:
        # Single time:
        offset=frame1*N
        d[k,:,:]=raw_array[offset+iz*NH+1:offset+(iz+1)*NH+1].reshape(nx,ny)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read streamfunction data into array for plotting:
in_file=open('evolution/pp.r4','r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()

# Loop over layers:
for iz in range(nz):
    # Frame where image will appear (0 1 2 first row, 3 4 5 second, ...):
    k=3*iz+1
    if frame2 > frame1:
        # Time average:
        for frame in range(frame1,frame2):
            offset=frame*N
            d[k,:,:]+=raw_array[offset+iz*NH+1:offset+(iz+1)*NH+1].reshape(nx,ny)
        d[k,:,:]=div*d[k,:,:]

    else:
        # Single time:
        offset=frame1*N
        d[k,:,:]=raw_array[offset+iz*NH+1:offset+(iz+1)*NH+1].reshape(nx,ny)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read vertical vorticity data into array for plotting:
in_file=open('evolution/zz.r4','r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()

# Loop over layers:
for iz in range(nz):
    # Frame where image will appear (0 1 2 first row, 3 4 5 second, ...):
    k=3*iz+2
    if frame2 > frame1:
        # Time average:
        for frame in range(frame1,frame2):
            offset=frame*N
            d[k,:,:]+=raw_array[offset+iz*NH+1:offset+(iz+1)*NH+1].reshape(nx,ny)
        d[k,:,:]=div*d[k,:,:]

    else:
        # Single time:
        offset=frame1*N
        d[k,:,:]=raw_array[offset+iz*NH+1:offset+(iz+1)*NH+1].reshape(nx,ny)

#=================================================================
# Plot each image with its own colourbar:
for j in range(nim):
    ax1=ax[j]

    ax1.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])

    # Label x axis only for images in bottom row:
    row=int(j/3)
    if row < nz-1:
        plt.setp(ax1.get_xticklabels(), visible=False)
    else:
        ax1.set_xlabel('$x$', fontsize=20)

    # Label y axis only for images in leftmost column:
    col=j-3*row
    if col > 0:
        plt.setp(ax1.get_yticklabels(), visible=False)
    else:
        ax1.set_ylabel('$y$', fontsize=20)

    if option == "m":
        # Project data onto vertical modes:
        m=row
        Z=np.zeros((nx,ny))
        for iz in range(nz):
            k=3*iz+col
            Z+=vec[iz,m]*d[k]
    else:
        # Use existing layer data:
        Z=d[j]

    # Work out the overall min/max values:
    zmin=np.amin(Z)
    zmax=np.amax(Z)
    #print()
    #print(' Min & max field levels = '+str(zmin)+', '+str(zmax))

    zmag=max(abs(zmin),abs(zmax))
    zmin=-zmag
    zmax= zmag

    #zm_in = input(' Min level to show (default '+str(zmin)+')? ')
    #zmin = float(zm_in or zmin)
    #zm_in = input(' Max level to show (default '+str(zmax)+')? ')
    #zmax = float(zm_in or zmax)

    # Obtain contour levels for plotting the colorbars:
    dz=2.0*contint(zmin,zmax)
    jmin=-int(-zmin/dz)
    jmax= int( zmax/dz)
    clevels=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

    # Plot the image in an array with an optional colourbar:
    im1=ax1.imshow(Z.T,cmap=cm.seismic,vmin=zmin,vmax=zmax,extent=(xmin,xmax,ymin,ymax),origin='lower',interpolation='bilinear')
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="4%", pad=0.1)
    cbar=fig.colorbar(im1, cax=cax, ticks=clevels)

    if row == 0:
        # Add title:
        ax1.set_title(field[col], fontsize=36)

fig.subplots_adjust(wspace=0.8, hspace=0.8)
    
#=========================================================================
# Save image:
if frame2 > frame1:
    figname=option+'_t'+str(t1)+'-'+str(t2)+'.png'
else:
    figname=option+'_t'+str(t1)+'.png'

fig.savefig(figname, bbox_inches='tight', pad_inches=0.025, dpi=300)

print()
print(' To view the image, type')
print()
print(' eom '+figname+' &')
print()
