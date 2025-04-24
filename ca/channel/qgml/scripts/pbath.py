#!/usr/bin/python

# This script plots the bathymetry bath.r8.

#  @@@@   Run from the current job directory   @@@@

#========== Perform the generic imports =========
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#-------------------------------------------------
# Work out x & y limits and grid resolution (nx, ny) by reading parameters.f90:
with open('src/parameters.f90','r') as in_file:
    fread=in_file.readlines()
    for line in fread:
        if ':: ellx=' in line:
            ellx=float(line.split("=")[1].split(",")[0])
        if ':: elly=' in line:
            elly=float(line.split("=")[1].split(",")[0])
        if ':: nx=' in line:
            nx=int(line.split("=")[1].split(",")[0])
        if ':: ny=' in line:
            ny=int(line.split("=")[1].split(",")[0])

xmax=0.5*ellx
ymax=0.5*elly
xmin=-xmax
ymin=-ymax

# Increase ny by 1 to include boundary points:
ny=ny+1

with open('bath.r8', 'rb') as in_file:
    zero=np.fromfile(in_file,dtype=np.float64,count=1)[0]
    bb_array=np.fromfile(in_file,dtype=np.float64)

# bb_array now contains ny*nx values, matching Fortran's (ny+1) x nx grid
bb_array=bb_array.reshape((ny, nx), order='F')  # Fortran order

fig,ax=plt.subplots(figsize=(8,6))
bbi=ax.imshow(bb_array.T,cmap=cm.jet,vmin=bb_array.min(),vmax=bb_array.max(),
              extent=(xmin,xmax,ymin,ymax),origin='lower',interpolation='bilinear')
cbar=fig.colorbar(bbi,ax=ax)
ax.set_xlabel('$x$',fontsize=20)
ax.set_ylabel('$y$',fontsize=20)
ax.set_title("Bathymetry",fontsize=36)
plt.savefig("bath.png",dpi=150)
