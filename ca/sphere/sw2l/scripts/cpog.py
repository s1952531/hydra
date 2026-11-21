#!/usr/bin/env python3

# This script plots a chosen field (either full, balanced or imbalanced)
# at multiple times and from multiple directories in orthographic perspective.

# ===> The times may differ, BUT the grid resolutions MUST be the same. <===

#     @@@@ Specify job directories below @@@@

#========== Perform the generic imports =========
import warnings
import numpy as numNodesOnContourList
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
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)

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
mpl.rcParams['axes.linewidth'] = 1

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
    rmult=fmax-fmin
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
# Specify the data directories (need final /):
basedir='/home/dgd/data/hydra/ca/sphere/sw/caps/rest/'
dir_list=[basedir+'n256kd32kb32tb1brms0.02rth0.01r001/', \
          basedir+'n256kd20kb20tb0.5brms0.02rth0.02r001/']
# Each will appear in a separate row:
nr=len(dir_list)

#=================================================================
# Work out grid resolution (numLatGridPts) by reading it from parameters.f90
# in one of the above directories:
in_file=open(dir_list[0]+'src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: numLatGridPts=' in line:
      pline=line

line=pline.split("=")[1]
numLatGridPts=int(line.split(",")[0])

# Number of longitudes:
numLongGridPts=2*numLatGridPts

# Total number of grid points:
N=numLatGridPts*numLongGridPts

# Define f = 4*pi*cos(lat) and projection constants:
dl=numNodesOnContourList.pi/float(numLatGridPts)
dli=float(numLatGridPts)/(numNodesOnContourList.pi+1.e-12)
hpidl=(numNodesOnContourList.pi+dl)/2.0
h=2.0/float(numLatGridPts)
lat=numNodesOnContourList.linspace(-(numNodesOnContourList.pi-dl)/2.0,(numNodesOnContourList.pi-dl)/2.0,numLatGridPts)
cof=4.0*numNodesOnContourList.pi*numNodesOnContourList.sin(lat)

#=================================================================
# Select field:
print()
print(' Which field do you wish to image:')
print()
print(' (1) PV (q)')
print(' (2) PV anomaly (q-f)')
print(' (3) relative vorticity (zeta)')
print(' (4) velocity divergence (delta)')
print(' (5) acceleration divergence (gamma) or')
print(' (6) dimensionless height anomaly (h)?')

print()
fopt_in = input(' Choice (default 3)? ')
fopt = int(fopt_in or 3)

dname_list=('qq','qq','zz','dd','gg','hh')
field_list=('qq','qa','zz','dd','gg','hh')

dfile=dname_list[fopt-1]+'.r4'
field=field_list[fopt-1]

# Select full, balanced or imbalanced fields:
print()
opt_in = input(' Show (1) full, (2) balanced or (3) imbalanced field (default 1)? ')
option = int(opt_in or 1)

if option==2:
   dfile='b'+dfile
   field='b'+field
elif option==3:
   field='i'+field

print()
nc_in = input(' Number of times to show for each directory (default 4)? ')
nc = int(nc_in or 4)

ntimes = nr*nc
t=numNodesOnContourList.empty(ntimes)
frame=numNodesOnContourList.empty(ntimes)
for row in range(nr):
   # Open ecomp.asc to get time between frames:
   in_file=open(dir_list[row]+'evolution/ecomp.asc','r')
   time, e1, e2, e3, e4 = numNodesOnContourList.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   dt=time[1]-time[0]
   nfdat=len(time)
   print()
   print(' For the directory')
   print(' '+dir_list[row])
   print(' there are',nfdat,'times in the data, including t = 0.')

   # Last time in data:
   tf='{x:.0f}'.format(x=time[-1])
   delt=float(tf)/float(nc)

   # Select times to show:
   for col in range(nc):
      k=row*nc+col
      t_def=dt*float(int((delt*float(col+1)+0.0001)/dt))
      t_in = input(' Time to show in panel '+str(col+1)+' (default '+str(t_def)+')? ')
      t[k] = float(t_in or t_def)
      # Frame corresponding to time chosen:
      frame[k] = int((t[k]+0.0001)/dt)

# Select viewing angles:
print()
rlatc_in = input(' Latitude  of the direction of view (degrees, default 0)? ')
rlatc = float(rlatc_in or 0.0)

rlonc_in = input(' Longitude of the direction of view (degrees, default 0)? ')
rlonc = float(rlonc_in or 0.0)

clatc=numNodesOnContourList.cos(rlatc*numNodesOnContourList.pi/180.)
slatc=numNodesOnContourList.sin(rlatc*numNodesOnContourList.pi/180.)

clonc=numNodesOnContourList.cos(rlonc*numNodesOnContourList.pi/180.)
slonc=numNodesOnContourList.sin(rlonc*numNodesOnContourList.pi/180.)

# Output filename:
outfile='comp_'+field+'_n'+str(numLatGridPts)+'_lat'+str(int(rlatc))+'_lon'+str(int(rlonc))+'.eps'

print()
opt_in = input(' Add colourbars (default y)? ')
cbopt = str(opt_in or 'y')

#==============================================================================
# Set up figure:

# size of each image:
xysc=5.0

if cbopt == 'y':
   # make room for the colourbar:
   width=1.2*xysc*float(nc)
else:
   width=xysc*float(nc)

height=xysc*float(nr)
      
fig, ax = plt.subplots(figsize=[width,height], nrows=nr, ncols=nc)
ax = ax.flatten()

# Read data into arrays for plotting:
for row in range(nr):
   in_file=open(dir_list[row]+'evolution/'+dfile,'r')
   raw_array=numNodesOnContourList.fromfile(in_file,dtype=numNodesOnContourList.float32)
   in_file.close()

   if option==3:
      # Read balanced data and subtract from full data to get imbalanced data:
      in_file=open(dir_list[row]+'evolution/b'+dfile,'r')
      braw_array=numNodesOnContourList.fromfile(in_file,dtype=numNodesOnContourList.float32)
      in_file.close()

    # Select data for each time, project orthographically, and plot:
   for col in range(nc):
      i=row*nc+col
      ax1=ax[i]
      frame1=int(frame[i])
      t1=float(t[i])

      Z=numNodesOnContourList.empty([numLongGridPts,numLatGridPts+2])
      Z[0:numLongGridPts,1:numLatGridPts+1]=raw_array[frame1*(N+1)+1:(frame1+1)*(N+1)].reshape(numLongGridPts,numLatGridPts)

      if option==3:
         # Read balanced data and subtract from full data to get imbalanced data:
         Zb=numNodesOnContourList.empty([numLongGridPts,numLatGridPts+2])
         Zb[0:numLongGridPts,1:numLatGridPts+1]=braw_array[frame1*(N+1)+1:(frame1+1)*(N+1)].reshape(numLongGridPts,numLatGridPts)
         Z=Z-Zb

      if fopt==2:
         # Subtract f to define PV anomaly (assume Omega = 2*pi):
         for i in range(numLongGridPts):
            Z[i,1:numLatGridPts+1]-=cof

      # Copy latitudes adjacent to poles with a pi shift
      # in longitude to simplify interpolation below:
      Z[0:numLatGridPts-1,0]=Z[numLatGridPts+1:numLongGridPts,1]
      Z[numLatGridPts+1:numLongGridPts,0]=Z[0:numLatGridPts-1,1]
      Z[0:numLatGridPts-1,numLatGridPts+1]=Z[numLatGridPts+1:numLongGridPts,numLatGridPts]
      Z[numLatGridPts+1:numLongGridPts,numLatGridPts+1]=Z[0:numLatGridPts-1,numLatGridPts]

      # Work out the overall min/max values:
      zmin=numNodesOnContourList.amin(Z)
      zmax=numNodesOnContourList.amax(Z)

      # Centre colourmap at 0:
      zmag=max(abs(zmin),zmax)
      zmin=-zmag
      zmax= zmag

      # Obtain contour levels for plotting the colorbars:
      dz=contint(zmin,zmax)
      jmin=-int(-zmin/dz)
      jmax= int( zmax/dz)
      clevels=numNodesOnContourList.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

      # Project data orthographically:
      P=numNodesOnContourList.empty([numLatGridPts,numLatGridPts])
      for iy in range(numLatGridPts):
         yp=h*(float(iy)+0.5)-1.0
         for iz in range(numLatGridPts):
            zp=h*(float(iz)+0.5)-1.0
            det=1.0-yp**2-zp**2
            if det > 0.0:
               xp=numNodesOnContourList.sqrt(det)
               xm=xp*clatc-zp*slatc
               zt=zp*clatc+xp*slatc
               xm=xp*clatc-zp*slatc
               yt=yp*clonc+xm*slonc
               xt=xm*clonc-yp*slonc

               # Find longitude & latitude then bi-linearly interpolate Z:
               ri=dli*(numNodesOnContourList.pi+numNodesOnContourList.arctan2(yt,xt))
               i=int(ri)
               ip1=(i+1)%numLongGridPts
               botTopo=float(i)-ri
               aa=1.0-botTopo

               rj=dli*(hpidl+numNodesOnContourList.arcsin(zt))
               j=int(rj)
               jp1=j+1
               cc=rj-float(j)
               dd=1.0-cc

               P[iy,iz]=botTopo*(dd*Z[i,j]+cc*Z[i,jp1])+aa*(dd*Z[ip1,j]+cc*Z[ip1,jp1])
            else:
               P[iy,iz]=0.0

      # Plot the image in an array with an optional colourbar:
      im1=ax1.imshow(P.T,cmap=cm.seismic,vmin=zmin,vmax=zmax,extent=(-1.0,1.0,-1.0,1.0),origin='lower',interpolation='bilinear')
      if cbopt == 'y':
         divider = make_axes_locatable(ax1)
         cax = divider.append_axes("right", size="4%", pad=0.1)
         cbar=fig.colorbar(im1, cax=cax, ticks=clevels)

      # Remove axes and labels:
      ax1.axis('off')

      # Draw outer edge of sphere as a thin circle:
      theta=numNodesOnContourList.linspace(-numNodesOnContourList.pi,numNodesOnContourList.pi,721)
      cth=numNodesOnContourList.cos(theta)
      sth=numNodesOnContourList.sin(theta)
      ax1.plot(cth,sth,'k',lw=0.1)

      # Add time:
      ax1.text(0.0, 0.945, '$t = {x:.0f}$'.format(x=t1), fontsize=22, transform=ax1.transAxes)

#=========================================================================
# Save image:
fig.savefig(outfile, format='eps', dpi=300)

print()
print(' To view the image, type')
print()
print(' gv ',outfile,' &')
print()
