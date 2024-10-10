#!/usr/bin/python
#Program to call a matlab function from python and save the data 23 May 2019:

#WARNING! Uses Python3.0 or above. 

import matlab.engine
import numpy as np
import os 
from os.path import expanduser

# Keep track of original directory in cwd:
cwd = os.getcwd()

# Change directory to location of MATLAB routines:
home = expanduser("~")
os.chdir(home+"/math/conformal/MATLAB/sc-toolbox-master")

# Ask for domain shape parameters:
print (' The domain has a width of L_x and a height of L_y, and along the')
print (' bottom there is a ramp of width w and height h.')

lx = float(input(" Enter L_x (default 6.4): ") or ("6.4"))
ly = float(input(" Enter L_y (default 0.3): ") or ("0.3"))
w = float(input(" Enter w (default 1.5): ") or ("1.5"))
h = float(input(" Enter h (default 0.1): ") or ("0.1"))

ny = int(input(" Resolution in y (in the conformal domain, default 128): ") or ("128"))
nxd=16*ny
nx = int(input(" Resolution in x (in the conformal domain, default "+str(nxd)+"): ") or (nxd))

# Set up prevertices:
z1=str(ly)+"i"
z2="0.0"
z3=str(lx-w)
z4=str(lx)+"+"+str(h)+"i"
z5=str(lx)+"+"+str(ly)+"i"
prevert=[z1,z2,z3,z4,z5]

print ()
print (' Computing the conformal map using MATLAB...')

eng = matlab.engine.start_matlab()
[X,Xxt,Yxt,xt,yt]=eng.hydra_slope(nx,ny,prevert,nargout=5)

# Change back to original directory:
os.chdir(cwd)

print ()
print (' Writing data files...')

# Open files to read in data:
out_file1 = open("coords.r8", "bw")
out_file2 = open("derivs.r8", "bw")
out_file3 = open("domdim.asc", "w")
out_file4 = open("resolution.asc", "w")

X=np.array(X,dtype=np.complex128)
Xxt=np.array(Xxt,dtype=np.float64)
Yxt=np.array(Yxt,dtype=np.float64)
xt=np.array(xt)
yt=np.array(yt)

# Get min/max domain values:
xtmin=np.amin(xt)
xtmax=np.amax(xt)
ytmin=np.amin(yt)
ytmax=np.amax(yt)

# Write data to files:
X.real.astype('float64').tofile(out_file1)
X.imag.astype('float64').tofile(out_file1)
Yxt.astype('float64').tofile(out_file2)
Xxt.astype('float64').tofile(out_file2)
out_file3.write(str(xtmin)+' '+str(xtmax)+'\n'+str(ytmin)+' '+str(ytmax)+'\n')
out_file4.write(str(nx)+' '+str(ny)+'\n')

out_file1.close()
out_file2.close()
out_file3.close()
out_file4.close()

print ()
print (' Note that the conformal domain width (in x) is ',xtmax-xtmin)
print ('                       and the height (in y) is ',ytmax-ytmin)
