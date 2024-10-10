#!/usr/bin/env python

# This script plots zeta, delta and h_tilde at the final time
# from multiple directories in orthographic perspective.

# ===> The times may differ, BUT the grid resolutions MUST be the same. <===

#     @@@@ Specify job directories below @@@@

#========== Perform the generic imports =========
import warnings
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
basedir='/local_raid1/dgd1/hydra/ca/sphere/sw/caps/rest/'
dir_list=[basedir+'n256kd5n1-72tb50tau100brms0.005/', \
          basedir+'n256kd25n1-72tb50tau100brms0.005/']
# Each will appear in a separate row:
nr=len(dir_list)

#=================================================================
# Work out grid resolution (ng) by reading it from parameters.f90
# in one of the above directories:
in_file=open(dir_list[0]+'src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: ng=' in line:
      pline=line

line=pline.split("=")[1]
ng=int(line.split(",")[0])
# *** ng is assumed to be equal for all directories!!!

# Number of longitudes:
nt=2*ng

# Total number of grid points:
N=ng*nt

# Define f = 4*pi*cos(lat) and projection constants:
dl=np.pi/float(ng)
dli=float(ng)/(np.pi+1.e-12)
hpidl=(np.pi+dl)/2.0
h=2.0/float(ng)
lat=np.linspace(-(np.pi-dl)/2.0,(np.pi-dl)/2.0,ng)
cof=4.0*np.pi*np.sin(lat)

#=================================================================
field=['zz','dd','hh']
nc=len(field)
ntimes = nr*nc
t=np.empty(nr)
frame=np.empty(nr)
for row in range(nr):
   # Open ecomp.asc to get final time and corresponding frame:
   in_file=open(dir_list[row]+'evolution/ecomp.asc','r')
   time, e1, e2, e3, e4 = np.loadtxt(in_file,dtype=float,unpack=True)
   in_file.close()

   dt=time[1]-time[0]
   t[row]='{x:.0f}'.format(x=time[-1])
   frame[row] = int((t[row]+0.0001)/dt)
   print(' For the directory')
   print(' '+dir_list[row])
   print(' showing t =',t[row])

# Select viewing angles:
print()
rlatc_in = input(' Latitude  of the direction of view (degrees, default 30)? ')
rlatc = float(rlatc_in or 30.0)

rlonc_in = input(' Longitude of the direction of view (degrees, default 0)? ')
rlonc = float(rlonc_in or 0.0)

clatc=np.cos(rlatc*np.pi/180.)
slatc=np.sin(rlatc*np.pi/180.)

clonc=np.cos(rlonc*np.pi/180.)
slonc=np.sin(rlonc*np.pi/180.)

# Output filename:
outfile='zdh_n'+str(ng)+'_lat'+str(int(rlatc))+'_lon'+str(int(rlonc))+'.eps'

#print()
#opt_in = input(' Add colourbars (default y)? ')
#cbopt = str(opt_in or 'y')
cbopt='y'

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
   frame1=int(frame[row])
   t1=float(t[row])
    # Select data for each field, project orthographically, and plot:
   for col in range(nc):
      in_file=open(dir_list[row]+'evolution/'+field[col]+'.r4','r')
      raw_array=np.fromfile(in_file,dtype=np.float32)
      in_file.close()
      i=row*nc+col
      ax1=ax[i]

      Z=np.empty([nt,ng+2])
      Z[0:nt,1:ng+1]=raw_array[frame1*(N+1)+1:(frame1+1)*(N+1)].reshape(nt,ng)

      # Copy latitudes adjacent to poles with a pi shift
      # in longitude to simplify interpolation below:
      Z[0:ng-1,0]=Z[ng+1:nt,1]
      Z[ng+1:nt,0]=Z[0:ng-1,1]
      Z[0:ng-1,ng+1]=Z[ng+1:nt,ng]
      Z[ng+1:nt,ng+1]=Z[0:ng-1,ng]

      # Work out the overall min/max values:
      zmin=np.amin(Z)
      zmax=np.amax(Z)

      # Centre colourmap at 0:
      zmag=max(abs(zmin),zmax)
      zmin=-zmag
      zmax= zmag

      # Obtain contour levels for plotting the colorbars:
      dz=contint(zmin,zmax)
      jmin=-int(-zmin/dz)
      jmax= int( zmax/dz)
      clevels=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

      # Project data orthographically:
      P=np.empty([ng,ng])
      for iy in range(ng):
         yp=h*(float(iy)+0.5)-1.0
         for iz in range(ng):
            zp=h*(float(iz)+0.5)-1.0
            det=1.0-yp**2-zp**2
            if det > 0.0:
               xp=np.sqrt(det)
               xm=xp*clatc-zp*slatc
               zt=zp*clatc+xp*slatc
               xm=xp*clatc-zp*slatc
               yt=yp*clonc+xm*slonc
               xt=xm*clonc-yp*slonc

               # Find longitude & latitude then bi-linearly interpolate Z:
               ri=dli*(np.pi+np.arctan2(yt,xt))
               i=int(ri)
               ip1=(i+1)%nt
               bb=float(i)-ri
               aa=1.0-bb

               rj=dli*(hpidl+np.arcsin(zt))
               j=int(rj)
               jp1=j+1
               cc=rj-float(j)
               dd=1.0-cc

               P[iy,iz]=bb*(dd*Z[i,j]+cc*Z[i,jp1])+aa*(dd*Z[ip1,j]+cc*Z[ip1,jp1])
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
      theta=np.linspace(-np.pi,np.pi,721)
      cth=np.cos(theta)
      sth=np.sin(theta)
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
