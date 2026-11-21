#!/usr/bin/env python3

#================================================================
#  Contours a given field at a given time
#================================================================

#=====perform various generic imports=====
import os,warnings
import numpy as numNodesOnContourList
import matplotlib.patches as mpatches
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
# Set default plot save resolution to a large value:
mpl.rcParams['savefig.dpi'] = 200
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.artist import setp
#import matplotlib.cm as cm
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
mpl.rcParams['axes.linewidth'] = 2

# --------------------- Function definitions -----------------------
def contint(fmin,fmax):
    #Determines a nice contour interval (giving 10-20 divisions with
    #interval 1, 2 or 5x10^m for some m) given the minimum & maximum
    #values of the field data (fmin & fmax).

    fmax=0.9999999*fmax
    fmin=0.9999999*fmin
    #The 0.99... factor avoids having a superfluous tick interval
    #in cases where fmax-fmin is 10^m or 2x10^m

    mpow=0
    rmult=fmax-fmin
    while rmult < 10.0:
       mpow+=1
       rmult=rmult*10.0

    while rmult >= 100.0:
       mpow-=1
       rmult=rmult/10.0

    emag=10.0**(float(-mpow))

    kmult=int(rmult/10.0)

    if kmult < 1:
       ci=emag
    elif kmult < 2:
       ci=2.0*emag
    elif kmult < 4:
       ci=4.0*emag
    elif kmult < 8:
       ci=10.0*emag
    else:
       ci=20.0*emag

    return ci

#-------------------------------------------------------------------------
#Get image dimensions (numLongGridPts x numLatGridPts):
file_bytes=os.path.getsize('qq_init.r8')-8
numLatGridPts=int(numNodesOnContourList.sqrt(file_bytes/16))
numLongGridPts=2*numLatGridPts

dl=numNodesOnContourList.pi/float(numLatGridPts)
lat=numNodesOnContourList.linspace((dl-numNodesOnContourList.pi)/2.0,(numNodesOnContourList.pi-dl)/2.0,numLatGridPts)
lat=numNodesOnContourList.insert(lat,0,-numNodesOnContourList.pi/2.0)
lat=numNodesOnContourList.append(lat,numNodesOnContourList.pi/2.0)
lon=numNodesOnContourList.linspace(-numNodesOnContourList.pi,numNodesOnContourList.pi,numLongGridPts+1)

X, Y = numNodesOnContourList.meshgrid(lon, lat)

print()
q_in = input(' Data file prefix (h1, h2, d1, d2, ...; default h1)? ')
prefix = str(q_in or 'h1')

filename='evolution/'+prefix+'.r4'

#Read in data and find min/max values:
in_file = open(filename,'r')
raw_array = numNodesOnContourList.fromfile(in_file,dtype=numNodesOnContourList.float32)
in_file.close()

#Open ecomp.asc to get final time and corresponding frame:
in_file=open('evolution/ecomp.asc','r')
time, e1, e2, e3, e4 = numNodesOnContourList.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

frame_def = int((round(time[-1],0)+0.0001)/(time[1]-time[0]))
q_in = input(' Frame to read (0 for first; default '+str(frame_def)+')? ')
frame = int(q_in or frame_def)

N=numLatGridPts*numLongGridPts
Z=numNodesOnContourList.zeros([numLatGridPts+2,numLongGridPts+1])
#Grab frame:
Z[1:numLatGridPts+1,0:numLongGridPts]=raw_array[frame*(N+1)+1:(frame+1)*(N+1)].reshape(numLongGridPts,numLatGridPts).T
#Add polar values by extrapolation:
Z[0,0:numLongGridPts]=1.5*Z[1,0:numLongGridPts]-0.5*Z[2,0:numLongGridPts]
Z[numLatGridPts+1,0:numLongGridPts]=1.5*Z[numLatGridPts,0:numLongGridPts]-0.5*Z[numLatGridPts-1,0:numLongGridPts]
#Add periodic edge in longitude:
Z[0:numLatGridPts+2,numLongGridPts]=Z[0:numLatGridPts+2,0]

q_in = input(' Add axis labels and ticks (y/n; default y)? ')
opt = str(q_in or 'y')

#-----------------------------------------------------------------------------
#Next plot results:
Zmin=numNodesOnContourList.amin(Z)
Zmax=numNodesOnContourList.amax(Z)
print(' Min field value = ',round(Zmin,5),'  Max field value = ',round(Zmax,5))

Zmin=-Zmax
dZ_def=contint(Zmin,Zmax)/2.0
q_in = input(' Contour interval (default '+str(dZ_def)+')? ')
dZ = float(q_in or dZ_def)

jmax=int(Zmax/dZ+0.5)
Zmax=dZ*(float(jmax)-0.5)
clevels=numNodesOnContourList.linspace(-Zmax,Zmax,2*jmax)

fig1 = plt.figure(1,figsize=[12,6])
ax1 = fig1.add_subplot(111)

ax1.xaxis.set_ticks([-numNodesOnContourList.pi,-numNodesOnContourList.pi/2.0,0.0,numNodesOnContourList.pi/2.0,numNodesOnContourList.pi])

ax1.yaxis.set_ticks([-numNodesOnContourList.pi/2.0,0.0,numNodesOnContourList.pi/2.0])

if opt=='y':
   ax1.set_xlabel('$\\lambda$', fontsize=24)
   ax1.set_ylabel('$\\phi$', fontsize=24)
   ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
   ax1.set_yticklabels([r'$-\pi/2$',r'$0$',r'$\pi/2$'],fontsize=20)
else:
   plt.setp(ax1.get_xticklabels(), visible=False)
   plt.setp(ax1.get_yticklabels(), visible=False)

# Contour:
ax1.contour(X,Y,Z,clevels,colors='k',origin='lower',extent=(-numNodesOnContourList.pi,numNodesOnContourList.pi,-numNodesOnContourList.pi/2,numNodesOnContourList.pi/2),linewidths=1)

# Save figure after cropping:
figname=prefix+'-'+str(frame)+'.pdf'
plt.savefig(figname, format='pdf', dpi=600)

print(' To display the results, type:')
print()
print(' ev',figname,'&')
print()
