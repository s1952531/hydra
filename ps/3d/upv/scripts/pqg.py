#!/usr/bin/env python3

# This script plots xz cross sections at y = 0 of b, zeta, static
# stability and u from data in 3d/qgb.r4, qgz.r4, qgg.r4 and qgu.r4.

#     @@@@   Run from the current job directory   @@@@

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
# Work out grid resolution (ng & nz) by reading it from parameters.f90:
in_file=open('src/parameters.f90','r')
fread=in_file.readlines()
for line in fread:
   if ':: ng=' in line:
      pline=line

line=pline.split("=")[1]
ng=int(line.split(",")[0])
line=pline.split("=")[2]
nz=int(line.split(",")[0])+1

# To show more contour lines in main image:
mult=2

#---------------------------------------------------
# Read buoyancy field (cross section at constant x):
in_file=open('3d/qgb.r4','r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()
A=np.empty([nz,ng,ng])
A=raw_array[1:nz*ng*ng+1].reshape(nz,ng,ng)
ix=int(ng/2+0.5)
b=A[:,ix,:]

# Work out the overall min/max values:
bmin=np.amin(b)
bmax=np.amax(b)

# Obtain contour levels for plotting the colorbars:
db=contint(bmin,bmax)
jmin=-int(-bmin/db)
jmax= int( bmax/db)
clevels1=np.linspace(db*float(jmin),db*float(jmax),jmax-jmin+1)
db=db/mult
jmin=-int(-bmin/db)
jmax= int( bmax/db)
clevels1f=np.linspace(db*float(jmin),db*float(jmax),jmax-jmin+1)

#---------------------------------------------------
# Read vorticity field (cross section at constant x):
in_file=open('3d/qgz.r4','r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()
A=np.empty([nz,ng,ng])
A=raw_array[1:nz*ng*ng+1].reshape(nz,ng,ng)
ix=int(ng/2+0.5)
z=A[:,ix,:]

# Work out the overall min/max values:
zmin=np.amin(z)
zmax=np.amax(z)
zmag=max(-zmin,zmax)
zmin=-zmag
zmax=zmag

# Obtain contour levels for plotting the colorbars:
dz=contint(zmin,zmax)
jmin=-int(-zmin/dz)
jmax= int( zmax/dz)
clevels2=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)
dz=dz/mult
jmin=-int(-zmin/dz)
jmax= int( zmax/dz)
clevels2f=np.linspace(dz*float(jmin),dz*float(jmax),jmax-jmin+1)

#-----------------------------------------------------------
# Read static stability field (cross section at constant x):
in_file=open('3d/qgg.r4','r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()
A=np.empty([nz,ng,ng])
A=raw_array[1:nz*ng*ng+1].reshape(nz,ng,ng)
ix=int(ng/2+0.5)
g=A[:,ix,:]

# Work out the overall min/max values:
gmin=np.amin(g)
gmax=np.amax(g)

# Obtain contour levels for plotting the colorbars:
dg=contint(gmin,gmax)
jmin=-int(-gmin/dg)
jmax= int( gmax/dg)
clevels3=np.linspace(dg*float(jmin),dg*float(jmax),jmax-jmin+1)
dg=dg/mult
jmin=-int(-gmin/dg)
jmax= int( gmax/dg)
clevels3f=np.linspace(dg*float(jmin),dg*float(jmax),jmax-jmin+1)

#-----------------------------------------------------
# Read y velocity field (cross section at constant x):
in_file=open('3d/qgu.r4','r')
raw_array=np.fromfile(in_file,dtype=np.float32)
in_file.close()
A=np.empty([nz,ng,ng])
A=raw_array[1:nz*ng*ng+1].reshape(nz,ng,ng)
ix=int(ng/2+0.5)
u=A[:,ix,:]

# Work out the overall min/max values:
umin=np.amin(u)
umax=np.amax(u)
umag=max(-umin,umax)
umin=-umag
umax=umag

# Obtain contour levels for plotting the colorbars:
du=contint(umin,umax)
jmin=-int(-umin/du)
jmax= int( umax/du)
clevels4=np.linspace(du*float(jmin),du*float(jmax),jmax-jmin+1)
du=du/mult
jmin=-int(-umin/du)
jmax= int( umax/du)
clevels4f=np.linspace(du*float(jmin),du*float(jmax),jmax-jmin+1)

#==============================================================================
# Set up figure:

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(figsize=[14,7], nrows=2, ncols=2)

ax1.set_xlim([-np.pi,np.pi])
ax1.set_ylim([-np.pi,0.0])
ax2.set_xlim([-np.pi,np.pi])
ax2.set_ylim([-np.pi,0.0])
ax3.set_xlim([-np.pi,np.pi])
ax3.set_ylim([-np.pi,0.0])
ax4.set_xlim([-np.pi,np.pi])
ax4.set_ylim([-np.pi,0.0])

# Customise tick values:
ax1.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax2.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax3.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax4.xaxis.set_ticks([-np.pi,-np.pi/2.0,0.0,np.pi/2.0,np.pi])
ax1.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax2.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax3.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)
ax4.set_xticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$',r'$\pi/2$',r'$\pi$'],fontsize=20)

ax1.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0])
ax2.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0])
ax3.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0])
ax4.yaxis.set_ticks([-np.pi,-np.pi/2.0,0.0])
ax1.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$'],fontsize=20)
ax2.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$'],fontsize=20)
ax3.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$'],fontsize=20)
ax4.set_yticklabels([r'$-\pi$',r'$-\pi/2$',r'$0$'],fontsize=20)

# Fine-tune figure; hide x ticks for upper plots and y ticks for right plots:
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
plt.setp(ax4.get_yticklabels(), visible=False)

ax1.set_title('$b(0,y,z)/N^2 H$', fontsize=20)
ax2.set_title('$\\zeta(0,y,z)/f$', fontsize=20)
ax3.set_title('$\\Gamma(0,y,z)$', fontsize=20)
ax4.set_title('$u(0,y,z)/f$', fontsize=20)

ax3.set_xlabel('$y$', fontsize=20)
ax4.set_xlabel('$y$', fontsize=20)
ax1.set_ylabel('$Nz/f$', fontsize=20)
ax3.set_ylabel('$Nz/f$', fontsize=20)

extent=(-np.pi,np.pi,-np.pi,0.0)

# Plot the image in an array with an optional colourbar:
im1=ax1.imshow(b,cmap=cm.seismic,vmin=bmin,vmax=bmax,extent=extent,origin='lower',interpolation='bilinear')
cs1=ax1.contour(b, clevels1f, colors='k', extent=extent, linewidths=1)
divider = make_axes_locatable(ax1)
cax1 = divider.append_axes("right", size="4%", pad=0.1)
cbar=fig.colorbar(im1, cax=cax1, ticks=clevels1)
cbar.add_lines(cs1)

im2=ax2.imshow(z,cmap=cm.seismic,vmin=zmin,vmax=zmax,extent=extent,origin='lower',interpolation='bilinear')
cs2=ax2.contour(z, clevels2f, colors='k', extent=extent, linewidths=1)
divider = make_axes_locatable(ax2)
cax2 = divider.append_axes("right", size="4%", pad=0.1)
cbar=fig.colorbar(im2, cax=cax2, ticks=clevels2)
cbar.add_lines(cs2)

im3=ax3.imshow(g,cmap=cm.seismic,vmin=gmin,vmax=gmax,extent=extent,origin='lower',interpolation='bilinear')
cs3=ax3.contour(g, clevels3f, colors='k', extent=extent, linewidths=1)
divider = make_axes_locatable(ax3)
cax3 = divider.append_axes("right", size="4%", pad=0.1)
cbar=fig.colorbar(im3, cax=cax3, ticks=clevels3)
cbar.add_lines(cs3)

im4=ax4.imshow(u,cmap=cm.seismic,vmin=umin,vmax=umax,extent=extent,origin='lower',interpolation='bilinear')
cs4=ax4.contour(u, clevels4f, colors='k', extent=extent, linewidths=1)
divider = make_axes_locatable(ax4)
cax4 = divider.append_axes("right", size="4%", pad=0.1)
cbar=fig.colorbar(im4, cax=cax4, ticks=clevels4)
cbar.add_lines(cs4)

#=========================================================================
# Save image:
fig.savefig('qg.eps', format='eps', dpi=600)

print()
print(' To view the image, type')
print()
print(' gv qg.eps &')
print()
