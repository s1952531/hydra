#!/usr/bin/python3
#Program to call a matlab function from python and save the data 10 jan 2019:

#WARNING! Uses Python3.0 or above. 

import matlab.engine
import numpy as np
import os 
from os.path import expanduser

# Check if resolution is supplied in a file:
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

# ------------------------------------------------------------------------
# Initialise arrays to generate several conformal maps about inital point:
kx=nx+1
ky=ny+1
X=np.zeros((9,kx,ky))
Y=np.zeros((9,kx,ky))
Yx=np.zeros((9,kx,ky))
Yy=np.zeros((9,kx,ky))
Lx=np.zeros(9)
Ly=np.zeros(9)

Yso=1.0
Yro=1.5
dYs=0.01
dYr=0.01

k=0
for js in range(-1,2):
   Ys=Yso+dYs*float(js)
   cYs=str(Ys)+'i'
   for jr in range(-1,2):
      Yr=Yro+dYr*float(jr)
      #define prevert:
      cYr=str(Yr)+'i'
      prevert=['0.95+'+cYs, cYs, '0', '2-0.5i', '2+'+cYr, '1+'+cYr]
      print(' Generating map',k+1)
      [Z,dYdy,dYdx,xt,yt]=eng.heidi_originalgrid(nx,ny,prevert,nargout=5)
      Z=np.array(Z,dtype=complex)
      X[k]=Z.real
      Y[k]=Z.imag
      Yx[k]=np.array(dYdx)
      Yy[k]=np.array(dYdy)
      xt=np.array(xt)
      yt=np.array(yt)
      Lx[k]=np.amax(xt)-np.amin(xt)
      Ly[k]=np.amax(yt)-np.amin(yt)
      k+=1

# ------------------------------------------------------------------------
# Change back to original directory and write data:
os.chdir(cwd)
print()

# ------------------------------------------------------------------------
# Write conformal domain length and its Yr & Ys derivatives to domdim.asc:
print (' Writing domdim.asc')
out_file = open("domdim.asc", "w")
out_file.write(str(Lx[4])+' '+str(Ly[4])+'\n')
out_file.write(str((Lx[5]-Lx[3])/(2.0*dYr))+'\n')
out_file.write(str((Lx[7]-Lx[1])/(2.0*dYs))+'\n')
out_file.write(str((Lx[5]-2.0*Lx[4]+Lx[3])/dYr**2)+'\n')
out_file.write(str((Lx[7]-2.0*Lx[4]+Lx[1])/dYs**2)+'\n')
out_file.write(str((Lx[8]-Lx[6]+Lx[0]-Lx[2])/(4.0*dYr*dYs))+'\n')
out_file.close()

# -------------------------------------------------------
# Write X & Y and their Yr & Ys derivatives to coords.r8:
print (' Writing coords.r8')
out_file = open("coords.r8", "bw")
A=np.array(X[4],dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array(Y[4],dtype=np.float64)
A.astype('float64').tofile(out_file)

# First derivatives:
A=np.array((X[5]-X[3])/(2.0*dYr),dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Y[5]-Y[3])/(2.0*dYr),dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((X[7]-X[1])/(2.0*dYs),dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Y[7]-Y[1])/(2.0*dYs),dtype=np.float64)
A.astype('float64').tofile(out_file)

# Second derivatives:
A=np.array((X[5]-2.0*X[4]+X[3])/dYr**2,dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Y[5]-2.0*Y[4]+Y[3])/dYr**2,dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((X[7]-2.0*X[4]+X[1])/dYr**2,dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Y[7]-2.0*Y[4]+Y[1])/dYr**2,dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((X[8]-X[6]+X[0]-X[2])/(4.0*dYr*dYs),dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Y[8]-Y[6]+Y[0]-Y[2])/(4.0*dYr*dYs),dtype=np.float64)
A.astype('float64').tofile(out_file)
out_file.close()

# -----------------------------------------------------------
# Write Y_x & Y_y and their Yr & Ys derivatives to derivs.r8:
print (' Writing derivs.r8')
out_file = open("derivs.r8", "bw")
A=np.array(Yx[4],dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array(Yy[4],dtype=np.float64)
A.astype('float64').tofile(out_file)

# First derivatives:
A=np.array((Yx[5]-Yx[3])/(2.0*dYr),dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Yy[5]-Yy[3])/(2.0*dYr),dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Yx[7]-Yx[1])/(2.0*dYs),dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Yy[7]-Yy[1])/(2.0*dYs),dtype=np.float64)
A.astype('float64').tofile(out_file)

# Second derivatives:
A=np.array((Yx[5]-2.0*Yx[4]+Yx[3])/dYr**2,dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Yy[5]-2.0*Yy[4]+Yy[3])/dYr**2,dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Yx[7]-2.0*Yx[4]+Yx[1])/dYr**2,dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Yy[7]-2.0*Yy[4]+Yy[1])/dYr**2,dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Yx[8]-Yx[6]+Yx[0]-Yx[2])/(4.0*dYr*dYs),dtype=np.float64)
A.astype('float64').tofile(out_file)
A=np.array((Yy[8]-Yy[6]+Yy[0]-Yy[2])/(4.0*dYr*dYs),dtype=np.float64)
A.astype('float64').tofile(out_file)
out_file.close()

print()
print (' All done!')
