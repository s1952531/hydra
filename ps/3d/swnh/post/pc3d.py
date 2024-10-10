#!/usr/bin/env python

# This plots contours from congen.asc

#=====perform the various imports========
import os,warnings
import numpy as npy
import matplotlib as mpl
import matplotlib.pyplot as plt
# Set default plot save resolution to a large value:
mpl.rcParams['savefig.dpi'] = 200
warnings.simplefilter("ignore",DeprecationWarning)

# Maximum colour saturation level (1 = black):
sat=0.7

#========================================
# Open input files:
infile0=open('3d/vstruct.asc','r')
record0=infile0.readline()
nl=int(record0.split()[0])
depth=float(record0.split()[1])
dz=depth/float(nl-1)
z=npy.linspace(-depth/2.0,depth/2.0,nl)

# Get list of files for imaging:
filelist=[]
for file in os.listdir("3d"):
   if file.endswith("contours.asc"):
      filelist.append(file)

print ' Choose one of the following files to image:'
print
for i,file in enumerate(filelist):
   print ' ('+str(i+1)+')',file

print
i=int(raw_input(' Choice (default 1): ') or 1)
file=filelist[i-1]
print
print ' Imaging data in',file,'...'
print

in_file=open('3d/'+file,'r')

thh=float(raw_input('Horizontal rotation angle (default 45) ') or 45.0)
thv=float(raw_input('  Vertical rotation angle (default 45) ') or 45.0)
#ski=int(raw_input('Contours to skip for front ')) 
ski = 1

thh = npy.pi*thh/180.0
thv = npy.pi*thv/180.0
coh = npy.cos(thh)
cov = npy.cos(thv)
sih = npy.sin(thh)
siv = npy.sin(thv)

# Read the chosen frame and plot:
record=in_file.readline()
nc=int(record.split()[0])
npt=int(record.split()[1])
t=float(record.split()[2])

print
print ' nc = ',nc,'   npt = ',npt,'   t = ',t

print
xmax=float(raw_input('Maximum |x| to show (default: pi)? ') or npy.pi)
ymax=float(raw_input('Maximum |y| to show (default: pi)? ') or npy.pi)
print

#========================================================
# Create figure:
fig = plt.figure(figsize=[5.0,5.0*ymax/xmax])

# Add axes to figure
ax = fig.add_subplot(111,autoscale_on=False)

# Fix domain of view:
ax.set_xlim([-xmax,xmax])
ax.set_ylim([-ymax,ymax])

# Loop over contours and plot:
np=[]
ind=[]
nlr=[]
# Get minimum and maximum "levels" for choosing contour colours:
indmin=1000000
indmax=-1000000
for j in range(0,nc):
   record=in_file.readline()
   np.append(int(record.split()[0]))
   lev=int(record.split()[2])
   ind.append(lev)
   indmin=min(indmin,lev)
   indmax=max(indmax,lev)
   nll=int(record.split()[3])
   nlr.append(nll)

for j in range(0,nc):
   # Choose contour colour depending on sign and magnitude of ind:
   z0=z[nlr[j]]
   zfac=(z0+depth/2.0)/depth
   if ind[j] > 0:
     color=plt.cm.Reds(1.0-sat*zfac)
   elif ind[j] < 0:
     color=plt.cm.Blues(1.0-sat*zfac)

   zi = siv*z0

   x=[]
   y=[]

   record=in_file.readline()
   x0=float(record.split()[0])
   y0=float(record.split()[1])
   x.append(x0)
   y.append(y0)

   for i in range(1,np[j]):
      record=in_file.readline()
      x.append(float(record.split()[0]))
      y.append(float(record.split()[1]))

   x.append(x0)
   y.append(y0)
     
   i=0
   while (i < np[j]):
     xp=[]
     yp=[]
     iflagx=0
     iflagy=0
     while ( (iflagx==0) and (iflagy==0) and (i<np[j]) ):
       ip = ((i+1) % np[j])
       dx = abs(x[ip]-x[i])
       dy = abs(y[ip]-y[i])
       if (dx > npy.pi):
         iflagx = 1

       if (dy > npy.pi):
         iflagy = 1

       if ((iflagx==0) and (iflagy==0) ):
         xx = x[i]*coh-y[i]*sih
         yy = (y[i]*coh+x[i]*sih)*cov+zi
         xp.append(xx)
         yp.append(yy)

       i = i+1

     if ((iflagx == 0) and (iflagy==0)):
       xx = x[0]*coh-y[0]*sih
       yy = (y[0]*coh+x[0]*sih)*cov+zi
       xp.append(xx)
       yp.append(yy)

     if ( (ind[j]<2) or ( (nlr[j] % ski) == 0 ) ):   
        ax.plot(xp,yp,color=color,lw=0.1)


ax.set_xticks([],[])
ax.set_yticks([],[])

basename = file.split("contours")[0]
fig.savefig('3d/'+basename+'.eps', format='eps', dpi=300)

print
print ' To view the image, type'
print
print ' gv 3d/'+basename+'.eps'
print
