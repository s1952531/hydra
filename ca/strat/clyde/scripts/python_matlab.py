#!/usr/bin/python
#Program to call a matlab function from python and save the data 10 jan 2019:

#WARNING! Uses Python3.0 or above. 

import matlab.engine
import numpy as np
import os 
from os.path import expanduser

exists = os.path.isfile('resolution.asc')
if exists:
    # Read dimensions in this file
    f = open("resolution.asc", "r")
    nx=int(f.readline())
    ny=int(f.readline())
    f.close()
else:
    # Keep presets
    print ('Using default resolution')

# Keep track of original directory in cwd:
cwd = os.getcwd()

# Change directory to location of MATLAB routines:
home = expanduser("~")
os.chdir(home+"/math/conformal/MATLAB/sc-toolbox-master")

print ('Computing the conformal map using MATLAB...')

eng = matlab.engine.start_matlab()
if exists:
    [X,Xxt,Yxt,xt,yt]=eng.heidi_originalgrid(nx,ny,nargout=5)
else:
    [X,Xxt,Yxt,xt,yt]=eng.heidi_originalgrid(nargout=5)

# Change back to original directory:
os.chdir(cwd)

print ('Writing data files...')

# Open files to read in data:
out_file1 = open("coords.r8", "bw")
out_file2 = open("derivs.r8", "bw")
out_file3 = open("domdim.asc", "w")

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


