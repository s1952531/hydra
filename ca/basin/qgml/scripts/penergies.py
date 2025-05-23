#!/usr/bin/python

# This script uses the streamfunction in evolution/qq.r4
# to obtain and plot the energies for all layers/interfaces.

#=====perform various generic imports=====
import warnings
import numpy as np

from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)

## global settings

# set tick label size:
label_size = 18
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

# Ensure latex fonts throughout:
rc('font', **{'family': 'Times New Roman'})
rc('text', usetex=True)
#=========================================

#=================================================================
# Re-arrange the data for all frames
def arrange_data(field_data):
    data_array=np.zeros((nx,ny,nz,nt))

    for frame in range(nt):
        for iz in range(nz):
            # Use existing layer data:
            offset=frame*N+iz*NH+1
            data_array[:,:,iz,frame]=field_data[offset:offset+NH].reshape(nx,ny)

    return data_array

#-------------------------------------------------
# Work out x & y limits, grid resolution (nx, ny & nz),
# and data save interval by reading parameters.f90:
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
        if ':: nz=' in line:
            nz=int(line.split("=")[1].split(",")[0])
        if ':: tgsave=' in line:
            dtsave=float(line.split("=")[1].split(",")[0])

dsumi=1/(nx*ny)
danorm=np.full((ny+1,nx+1),dsumi)
danorm[1:ny,0   ]= 0.5*dsumi
danorm[1:ny,nx  ]= 0.5*dsumi
danorm[0   ,1:nx]= 0.5*dsumi
danorm[ny  ,1:nx]= 0.5*dsumi
danorm[0   ,0   ]=0.25*dsumi
danorm[0   ,nx  ]=0.25*dsumi
danorm[ny  ,0   ]=0.25*dsumi
danorm[ny  ,nx  ]=0.25*dsumi

# Increase nx & ny by 1 to include boundary points:
nx=nx+1
ny=ny+1

# Read energy data to get final time in data:
with open('evolution/energy.asc','r') as in_file:
    time,ekin,epot,etot=np.loadtxt(in_file,dtype=float,unpack=True)
tsim=time[-1]

#=================================================================
# Read vertical structure file:
hhat = np.zeros(nz)
kdsq = np.zeros(nz)

with open('vertical.asc', 'r') as in_file:
    for iz in range(nz):
        line = in_file.readline()
        hhat[iz], kdsq[iz] = map(float, line.split())

# Number of horizontal grid points:
NH=nx*ny

# Total number of data elements per time frame:
N=nz*NH+1
# The "+1" includes the time element

#=================================================================
# Read streamfunction data into array for plotting:
with open('evolution/pp.r4','rb') as in_file:
    pp_array=np.fromfile(in_file,dtype=np.float32)

nt=len(pp_array) // N
pp_array=arrange_data(pp_array)

# Compute the wave numbers
kx = np.fft.fftfreq(nx, d=(xmax - xmin) / (nx - 1)) * 2 * np.pi
ky = np.fft.fftfreq(ny, d=(ymax - ymin) / (ny - 1)) * 2 * np.pi
kx, ky = np.meshgrid(kx, ky, indexing='ij')

ke=np.zeros((nt,nz))
pe=np.zeros((nt,nz-1))
if nz==2:
    bt=np.zeros(nt)
    bc=np.zeros(nt)

t=0.0
while t<=tsim:
    frame=int(t/dtsave+0.5)

    for iz in range(nz):
        pp_fft = np.fft.fft2(pp_array[:, :, iz, frame])
        uu = np.fft.ifft2(-1j * ky * pp_fft).real
        vv = np.fft.ifft2( 1j * kx * pp_fft).real

        ke[frame, iz] = 0.5 * hhat[iz] * np.sum((uu**2 + vv**2) * danorm)
        if iz < nz-1:
            pe[frame,iz] = 0.5 * kdsq[iz] * np.sum((pp_array[:,:,iz+1,frame] - pp_array[:,:,iz,frame])**2 * danorm)

    if nz==2:
        psiB = hhat[0]*pp_array[:,:,0,frame] + hhat[1]*pp_array[:,:,1,frame]
        psiB_fft = np.fft.fft2(psiB)
        psiB_x = np.fft.ifft2(1j * kx * psiB_fft).real
        psiB_y = np.fft.ifft2(1j * ky * psiB_fft).real
        bt[frame] = 0.5 * np.sum((psiB_x**2 + psiB_y**2) * danorm)

        psiT = pp_array[:,:,0,frame] - pp_array[:,:,1,frame]
        psiT_fft = np.fft.fft2(psiT)
        psiT_x = np.fft.ifft2(1j * kx * psiT_fft).real
        psiT_y = np.fft.ifft2(1j * ky * psiT_fft).real
        bc[frame] = 0.5 * np.sum((psiT_x**2 + psiT_y**2 + kdsq[0]*psiT**2) * danorm)

    t+=dtsave

#=================================================================
# Save the energies in a single file:
with open('evolution/energies.asc', 'w') as out_file:
    for i in range(nt):
        line = f"{time[i]}"
        for iz in range(nz):
            line += f"\t{ke[i, iz]}"
        for iz in range(nz-1):
            line += f"\t{pe[i, iz]}"
        line += f"\t{etot[i]}\n"
        out_file.write(line)

#-------------------------------------------------------------------------
# Set up figure:
fig = plt.figure(1,figsize=[8,5])
ax1 = fig.add_subplot(111)

ax1.set_xlabel('$t$', fontsize=30)
ax1.set_ylabel('Energy', fontsize=30)
styles = ['--', '-.', ':', '-']

for iz in range(nz):
    ax1.plot(time, ke[:, iz], color='b', linestyle=styles[iz % len(styles)], lw=2, label=f'KE$_{iz+1}$')
for iz in range(nz-1):
    ax1.plot(time, pe[:, iz], color='r', linestyle=styles[iz % len(styles)], lw=2, label=f'PE$_{iz+1}$')

# Calculate the total energy
te = np.zeros(nt)
for i in range(nt):
    te[i] = np.sum(ke[i, :]) + np.sum(pe[i, :])

ax1.plot(time, te, color='k', linestyle='-', lw=2, label='TE')

if nz==2:
    ax1.plot(time, bt, color='g', linestyle=':', lw=2, label='BT')
    ax1.plot(time, bc, color='g', linestyle='-', lw=2, label='BC')
    print("Maximum relative error in difference between layerwise and modal total energies:", max(abs((bt+bc-te)/te)))

ax1.legend(loc='best',prop={'size':18})
plt.savefig('energy.png', bbox_inches='tight', pad_inches = 0.025, dpi=300)
