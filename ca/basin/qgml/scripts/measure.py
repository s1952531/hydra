#!/usr/bin/python

#--------------------------------------------------------------------
# This script calculates the PV measure in each layer, defined to be
# the fractional area of the domain with PV above a given threshold.
# We use ncontq discrete thresholds spanning the range [q_min,q_max],
# where ncontq is read from parameters.f90.

# Outputs the results (for all times saved in qq.r4) in the file
# evolution/far.asc. This is an ascii file containing the following:
# the time,
# t
# then,
# for j in range(nz-1):
#    farea[:,j]
# where the inner dimension of farea is [0,ncontq].
#--------------------------------------------------------------------

#==================================================
#  @@@@   Run from the current job directory   @@@@
#==================================================

import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt

#-------------------------------------------------------
# Extract nx, ny, nz and ncontq from src/parameters.f90:
with open('src/parameters.f90','r') as in_file:
    fread=in_file.readlines()
    for line in fread:
        if ':: nx=' in line:
            nx=int(line.split("=")[1].split(",")[0])
        if ':: ny=' in line:
            ny=int(line.split("=")[1].split(",")[0])
        if ':: nz=' in line:
            nz=int(line.split("=")[1].split(",")[0])
        if ':: ncontq=' in line:
            nq=int(line.split("=")[1].split(",")[0])

# Increase nx & ny by 1 to include boundary points:
nx=nx+1
ny=ny+1

# Number of horizontal grid points:
NH=nx*ny
da=1.0/float(NH)

# Total number of data elements per time frame:
N=nz*NH+1
# The "+1" includes the time element.

#=================================================================
# Read energy data to get the data save times:
in_file=open('evolution/energy.asc','r')
time, ekin, epot, etot = np.loadtxt(in_file,dtype=float,unpack=True)
in_file.close()

# Number of times in the data:
nt=len(time)

# Read all PV data into an array for processing:
with open('evolution/qq.r4','rb') as in_file:
    qq_array=np.fromfile(in_file,dtype=np.float32)

print(' Determining min/max values for each layer.')

global_farea_min = float('inf')
global_farea_max = float('-inf')

global_qmin = [float('inf')] * nz
global_qmax = [float('-inf')] * nz

for frame in range(nt):
    offset=frame*N
    for iz in range(nz):
        q = qq_array[offset + iz * NH + 1:offset + (iz + 1) * NH + 1]
        qmin = np.min(q)
        qmax = np.max(q)
        global_qmin[iz] = min(global_qmin[iz], qmin)
        global_qmax[iz] = max(global_qmax[iz], qmax)

        dq = (qmax - qmin) / float(nq)
        dqi = 1.0 / dq

        qmin = qmin - dq / 2.0
        qmax = qmax + dq / 2.0

        farea = np.zeros([nq + 1])
        for i in range(NH):
            k = int(dqi * (q[i] - qmin))
            farea[k] += da

        for k in range(1, nq + 1):
            farea[k] += farea[k - 1]

        global_farea_min = min(global_farea_min, np.min(farea))
        global_farea_max = max(global_farea_max, np.max(farea))

qmax_values = [[] for _ in range(nz)]
qmin_values = [[] for _ in range(nz)]

# Open output file:
diag_file = open('evolution/far.asc','w+')

# Loop over all times in the data:
for frame in range(nt):
    print(time[frame], file=diag_file)
    print(' Processing t =',time[frame])

    # Create a new figure for this time frame:
    plt.figure(figsize=(10,5*nz))

    # Process each layer in turn:
    offset=frame*N
    for iz in range(nz):
        # Extract PV field in each layer:
        q=qq_array[offset+iz*NH+1:offset+(iz+1)*NH+1]

        # Compute min/max values of PV:
        qmin=np.min(q)
        qmax=np.max(q)

        qmin_values[iz].append(qmin)
        qmax_values[iz].append(qmax)

        dq=(qmax-qmin)/float(nq)
        dqi=1.0/dq

        # Increase range to catch min/max values correctly in PV measure:
        qmin=qmin-dq/2.0
        qmax=qmax+dq/2.0

        # Count number of grid points having PV in each interval
        # [qmin+k*dq,qmin+(k+1)*dq):
        farea=np.zeros([nq+1])
        for i in range(NH):
            k=int(dqi*(q[i]-qmin))
            farea[k]+=da # da = 1/NH

        # Form cumulative PDF:
        for k in range(nq+1):
            farea[k]+=farea[k-1]

        # Write data for this iz:
        print(farea, file=diag_file)

        # Plot farea for this layer and time frame:
        plt.subplot(nz,2,iz*2+1)
        plt.plot(range(nq + 1), farea)
        plt.xlabel('Threshold Index')
        plt.ylabel('Fractional Area')
        plt.title(f'farea at t={time[frame]:.2f}, Layer {iz + 1}')
        plt.xlim(0, nq)
        plt.ylim(global_farea_min, global_farea_max)
        plt.grid()

        # Plot qmax and qmin evolution for all layers:
        plt.subplot(nz,2,iz*2+2)
        plt.plot(time[:frame + 1], qmax_values[iz])
        plt.plot(time[:frame + 1], qmin_values[iz])
        plt.xlabel('Time')
        plt.ylabel('qmax / qmin')
        plt.title('Evolution of qmax and qmin')
        plt.xlim(time[0], time[-1])
        plt.ylim(global_qmin[iz], global_qmax[iz])

    # Save the plot for this time frame:
    plt.tight_layout()
    plt.savefig(f'{frame:04d}.png')
    plt.close()

diag_file.close()

print()
print(' All done. Results are in evolution/far.asc')
