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
        if ':: xmin=' in line:
            xmin=float(line.split("=")[1].split(",")[0])
        if ':: xmax=' in line:
            xmax=float(line.split("=")[1].split(",")[0])
        if ':: ymin=' in line:
            ymin=float(line.split("=")[1].split(",")[0])
        if ':: ymax=' in line:
            ymax=float(line.split("=")[1].split(",")[0])
        if ':: nx=' in line:
            nx=int(line.split("=")[1].split(",")[0])
        if ':: ny=' in line:
            ny=int(line.split("=")[1].split(",")[0])

# Increase nx & ny by 1 to include boundary points:
nx+=1
ny+=1

with open('bath.r8', 'rb') as in_file:
    # Assuming the first value is 'zero' and the rest is the bb array
    zero = np.fromfile(in_file, dtype=np.float64, count=1)
    bb_array = np.fromfile(in_file, dtype=np.float64)

bb_array = bb_array[0:nx*ny]
bb_array = bb_array.reshape((ny, nx))

fig, ax = plt.subplots(figsize=(8, 6))
bbi = ax.imshow( bb_array.T, cmap=cm.seismic, vmin=bb_array.min(), vmax=bb_array.max(),
                 extent=(xmin, xmax, ymin, ymax), origin='lower', interpolation='bilinear' )
cbar = fig.colorbar(bbi, ax=ax)
ax.set_xlabel('$x$', fontsize=20)
ax.set_ylabel('$y$', fontsize=20)
ax.set_title("Bathymetry", fontsize=36)
plt.savefig("bath.png", dpi=150)
