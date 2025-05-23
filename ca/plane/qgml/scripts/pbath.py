#!/usr/bin/python

# This script plots the bathymetry bath.r8.

#  @@@@   Run from the current job directory   @@@@

#========== Perform the generic imports =========
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

#-------------------------------------------------
# Work out x & y limits, grid resolution (ng & nz),
# and data save interval by reading parameters.f90:
with open('src/parameters.f90','r') as in_file:
    fread=in_file.readlines()
    for line in fread:
        if ':: ng=' in line:
            ng=int(line.split("=")[1].split(",")[0])
        if ':: nz=' in line:
            nz=int(line.split("=")[1].split(",")[0])
        if ':: tgsave=' in line:
            dtsave=float(line.split("=")[1].split(",")[0])

xmax=np.pi
ymax=np.pi
xmin=-xmax
ymin=-ymax

# For including periodic edges at x = pi & y = pi:
ngp1=ng+1

# Set up the image array:
bb_array=np.empty((ngp1,ngp1))

with open('bath.r8', 'rb') as in_file:
    # Assuming the first value is 'zero' and the rest is the bb array
    zero = np.fromfile(in_file, dtype=np.float64, count=1)
    bb_array[0:ng,0:ng] = np.fromfile(in_file, dtype=np.float64).reshape(ng,ng)

# Add periodic edges:
bb_array[ng,0:ng]=bb_array[0,0:ng]
bb_array[0:ng+1,ng]=bb_array[0:ng+1,0]

fig, ax = plt.subplots(figsize=(8, 6))
bbi = ax.imshow( bb_array.T, cmap=cm.seismic, vmin=bb_array.min(), vmax=bb_array.max(),
                 extent=(xmin, xmax, ymin, ymax), origin='lower', interpolation='bilinear' )
cbar = fig.colorbar(bbi, ax=ax)
ax.set_xlabel('$x$', fontsize=20)
ax.set_ylabel('$y$', fontsize=20)
ax.set_title("Bathymetry", fontsize=36)
plt.savefig("bath.png", dpi=150)
