#!/usr/bin/python
#Program to call a matlab function from python and save the data 10 jan 2019:

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
print (' The lower left corner of the domain is at x = y = 0, and the bottom')
print (' has the shape y = -a*x^2.')
print ()
xrt = float(input(" Total length of domain in x (default 8): ") or ("8.0"))
xw = float(input(" x coordinate of weir (default 4): ") or ("4.0"))
yw = float(input(" y coordinate of bottom of the weir (default 0.25): ") or ("0.25"))
ww = float(input(" Width of the weir (default 0.1): ") or ("0.1"))
ylt = float(input(" y coordinate of top left corner of the domain (default 1): ") or ("1.0"))
yrt = float(input(" y coordinate of top right corner of the domain (default 2): ") or ("2.0"))
yrb = float(input(" y coordinate of bottom right corner of the domain (default -1): ") or ("-1.0"))
print ()
nx = int(input(" Resolution in x (in the conformal domain, default 1000): ") or ("1000"))
ny = int(input(" Resolution in y (in the conformal domain, default 400): ") or ("400"))

# Set up prevertices:
xwl=xw-ww/2.0
xwr=xw+ww/2.0
ymin=yw+ww/2.0
z1=str(xwl)+"+"+str(ylt)+"i"
z2=str(ylt)+"i"
z3=str(0)
z4=str(xrt)+"+"+str(yrb)+"i"
z5=str(xrt)+"+"+str(yrt)+"i"
z6=str(xwr)+"+"+str(yrt)+"i"
prevert=[z1,z2,z3,z4,z5,z6]

print ()
print (' Computing the conformal map using MATLAB...')

eng = matlab.engine.start_matlab()
[X,Xxt,Yxt,xt,yt]=eng.heidi_originalgrid(nx,ny,prevert,ymin,nargout=5)

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
