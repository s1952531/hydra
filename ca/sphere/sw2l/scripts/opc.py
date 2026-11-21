#!/usr/bin/env python

# This plots contours from evolution/contours.dat, here just the
# first frame, but both layers. Adapted from cs/new_sphere/scripts/mpc.py

#=====perform the various imports========
import warnings
import numpy as numNodesOnContourList
import matplotlib as mpl
import matplotlib.pyplot as plt
# Set default plot save resolution to a large value:
mpl.rcParams['savefig.dpi'] = 800
warnings.simplefilter("ignore",DeprecationWarning)

# Fading value of alpha for parts of contours on the far side:
fade=0.5

# Approximate buffer width around edge of plot:
buf=0.1

## global settings                       
# set tick label size:                   
label_size = 0.001                       
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
# set x tick width and size:                
mpl.rcParams['xtick.major.size'] = 0.001    
mpl.rcParams['xtick.major.width'] = 0.001   
mpl.rcParams['xtick.minor.size'] = 0.001    
mpl.rcParams['xtick.minor.width'] = 0.001   
# set y tick width and size:                
mpl.rcParams['ytick.major.size'] = 0.001    
mpl.rcParams['ytick.major.width'] = 0.001   
mpl.rcParams['ytick.minor.size'] = 0.001    
mpl.rcParams['ytick.minor.width'] = 0.001
mpl.rcParams['axes.linewidth'] = 0.5  
#========================================

hw=buf*float(1+int(1.0/buf+0.5))

purecol = numNodesOnContourList.array([[0,0,0],[255,0,0],[0,0,255],[0,255,255],[255,0,255],[255,140,0],[0,255,0],[255,105,30],[155,20,147],[46,139,87]])/255.0
# RGB values for black, red, blue, cyan, magenta, orange, lime, chocolate, deep pink, sea green

#-----------------------------------------------------------------
# For plotting the edge of the visible sphere:
ange=numNodesOnContourList.linspace(0.0,2.0*numNodesOnContourList.pi,181)
xe=numNodesOnContourList.cos(ange)
ye=numNodesOnContourList.sin(ange)

#-----------------------------------------------------------------
# Set viewing direction:
q_in = input(' Enter co-latitude (deg) of viewing direction (default: 30) ')
thv = float(q_in or 30.0)
q_in = input(' Enter  longitude  (deg) of viewing direction (default: 60) ')
phv = float(q_in or 0.0)

thv=numNodesOnContourList.pi*thv/180.0
phv=numNodesOnContourList.pi*phv/180.0

#-----------------------------------------------------------------
# Define rotation matrix:
cthv=numNodesOnContourList.cos(thv)
sthv=numNodesOnContourList.sin(thv)
cphv=numNodesOnContourList.cos(phv)
sphv=numNodesOnContourList.sin(phv)
txx=-cthv*cphv
txy=-cthv*sphv
txz=sthv
tyx=-sphv
tyy=cphv
tzx=sthv*cphv
tzy=sthv*sphv
tzz=cthv

#-----------------------------------------------------------------
#Define lines of constant latitude:
nphi = 90
phi = numNodesOnContourList.linspace(-numNodesOnContourList.pi, numNodesOnContourList.pi, nphi)
cphi = numNodesOnContourList.cos(phi)
sphi = numNodesOnContourList.sin(phi)
nlat = 11
dthetag = numNodesOnContourList.pi/float(nlat+1)
thetag = numNodesOnContourList.linspace(dthetag, numNodesOnContourList.pi-dthetag, nlat)

xa=numNodesOnContourList.zeros((nlat,nphi))
ya=numNodesOnContourList.zeros((nlat,nphi))
za=numNodesOnContourList.zeros((nlat,nphi))

for j in range(nphi):
   for i,th in enumerate(thetag):
      rg=numNodesOnContourList.sin(th)
      xg=rg*cphi[j]
      yg=rg*sphi[j]
      zg=numNodesOnContourList.cos(th)
      # Apply rotation matrix:
      xa[i,j]=tyx*xg+tyy*yg
      ya[i,j]=txx*xg+txy*yg+txz*zg
      za[i,j]=tzx*xg+tzy*yg+tzz*zg

#-----------------------------------------------------------------
#Define lines of constant longitude:
ntheta = 45
theta = numNodesOnContourList.linspace(0.0, numNodesOnContourList.pi, ntheta)
cth = numNodesOnContourList.cos(theta)
sth = numNodesOnContourList.sin(theta)
nlon = 24
dphig = 2.0*numNodesOnContourList.pi/float(nlon)
phig = numNodesOnContourList.linspace(0.0, 2.0*numNodesOnContourList.pi-dphig, nlon)

xb=numNodesOnContourList.zeros((nlon,ntheta))
yb=numNodesOnContourList.zeros((nlon,ntheta))
zb=numNodesOnContourList.zeros((nlon,ntheta))

for j in range(ntheta):
   for i,ph in enumerate(phig):
      xg=sth[j]*numNodesOnContourList.cos(ph)
      yg=sth[j]*numNodesOnContourList.sin(ph)
      zg=cth[j]
      # Apply rotation matrix:
      xb[i,j]=tyx*xg+tyy*yg
      yb[i,j]=txx*xg+txy*yg+txz*zg
      zb[i,j]=tzx*xg+tzy*yg+tzz*zg

#-----------------------------------------------------------------
# Create figure:
fig = plt.figure(figsize=[10.0,20.0])

##########################################################################
# Loop over rows (corresponds to layers here):
for row in range(2):
   # Add axes to figure
   ax = fig.add_subplot(1,2,row+1)

   ax.set(adjustable='box', aspect=1.0)

   # Plot edge of visible sphere:
   ax.plot(xe,ye,c='g',lw=0.5,alpha=fade)

   # Plot surface:
   for i in range(nlat):
      ax.plot(xa[i,:],ya[i,:],c='k',lw=0.5,alpha=0.2)

   for i in range(nlon):
      ax.plot(xb[i,:],yb[i,:],c='k',lw=0.5,alpha=0.2)
      
   # Open input file:
   contour_file=open('evolution/contours.asc','r')

   if row > 0:
      # Skip frames:
      for k in range(row):
         record=contour_file.readline()
         nc=int(record.split()[0])
         numNodesTot=int(record.split()[1])
         for i in range(nc+numNodesTot):
            record=contour_file.readline()

   # Read the chosen frame and plot:
   record=contour_file.readline()
   nc=int(record.split()[0])
   numNodesTot=int(record.split()[1])
   t=float(record.split()[2])

   print(' iz =',row,'  nc =',nc,'  numNodesTot =',numNodesTot)
   print()

   # Get number of nodes on each contour (numNodesOnContourList) and vorticity levels (fieldLvlIdx):
   npc=[]
   fieldLvlIdx=[]
   for j in range(nc):
      record=contour_file.readline()
      npc.append(int(record.split()[0]))
      fieldLvlIdx.append(int(record.split()[2]))

   for j in range(nc):
      x=[]
      y=[]
      z=[]
      # Read and assign the x, y & z coordinates of each contour:
      record=contour_file.readline()
      x0=float(record.split()[0])
      y0=float(record.split()[1])
      z0=float(record.split()[2])
      x.append(x0)
      y.append(y0)
      z.append(z0)

      for i in range(1,npc[j]):
         record=contour_file.readline()
         x.append(float(record.split()[0]))
         y.append(float(record.split()[1]))
         z.append(float(record.split()[2]))

      x.append(x0)
      y.append(y0)
      z.append(z0)

      x=numNodesOnContourList.array(x)
      y=numNodesOnContourList.array(y)
      z=numNodesOnContourList.array(z)

      # Apply rotation matrix:
      xh=tyx*x+tyy*y
      yh=txx*x+txy*y+txz*z
      zh=tzx*x+tzy*y+tzz*z

      # Plot contours at this time:
      c=purecol[fieldLvlIdx[j]%10]
      #ax.plot(xh,yh,color=c,lw=0.5)
      for i in range(len(zh)-1):
         if zh[i]<0.0 and zh[i+1]<0.0:
            ax.plot([xh[i],xh[i+1]],[yh[i],yh[i+1]],color=c,lw=0.5,alpha=fade)
         else:
            ax.plot([xh[i],xh[i+1]],[yh[i],yh[i+1]],color=c,lw=1)

   # Fix domain of view:
   ax.set_xlim([-hw,hw])
   ax.set_ylim([-hw,hw])

   # Hide axes:
   plt.axis('off')

#######################################################################
# Finalise:
fig.subplots_adjust(wspace=0.0, hspace=0.0)

fig.savefig('contours.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.025, dpi=600)

print(' To view the results, type')
print()
print(' ev contours.pdf &')
print()
