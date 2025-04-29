#!/usr/bin/python

# This script plots the PV, streamfunction and vorticity in
# evolution/qq.r4, pp.r4 and zz.r4, for all layers.

#  @@@@   Run from the current job directory   @@@@

#========== Perform the generic imports =========
import sys,os,warnings
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
from matplotlib import rcParams
from matplotlib import rc
rcParams.update({'figure.autolayout': True})
warnings.simplefilter("ignore",DeprecationWarning)
sys.path.append(os.path.expandvars('${HOME}/hydra/lib/'))
from wbgyr import cmap_wbgyr

# Ensure latex fonts throughout:
rc('font',**{'family': 'Times New Roman'})
rc('text',usetex=True)

# set tick label size:
label_size=20
mpl.rcParams['xtick.labelsize']=label_size
mpl.rcParams['ytick.labelsize']=label_size
# set x tick width and size:
mpl.rcParams['xtick.major.size']=5
mpl.rcParams['xtick.major.width']=2
mpl.rcParams['xtick.minor.size']=3
mpl.rcParams['xtick.minor.width']=1
# set y tick width and size:
mpl.rcParams['ytick.major.size']=5
mpl.rcParams['ytick.major.width']=2
mpl.rcParams['ytick.minor.size']=3
mpl.rcParams['ytick.minor.width']=1
# set axes width:
mpl.rcParams['axes.linewidth']=1

#====================== Function definitions =======================
def contint(fmin,fmax):
    #Determines a nice contour interval (giving 10-20 divisions with
    #interval 1, 2 or 5x10^m for some m) given the minimum & maximum
    #values of the field data (fmin & fmax).

    fmax=0.9999999*fmax
    fmin=0.9999999*fmin
    #The 0.99... factor avoids having a superfluous tick interval
    #in cases where fmax-fmin is 10^m or 2x10^m

    emag=1.0
    rmult=max(1.0E-12,fmax-fmin)
    while rmult < 10:
        emag=emag/10
        rmult=rmult*10

    while rmult >= 100:
        emag=emag*10
        rmult=rmult/10

    kmult=int(rmult/10)

    if kmult < 1:
        ci=emag
    elif kmult < 2:
        ci=2*emag
    elif kmult < 4:
        ci=4*emag
    elif kmult < 8:
        ci=10*emag
    else:
        ci=20*emag

    return ci

#=================================================================
# Function to calculate global min and max for a field across all frames
def get_layer_min_max(field_data,is_pv=False):
    layer_min=np.full(nz, np.inf)
    layer_max=np.full(nz,-np.inf)
    data_array=np.zeros((nx,ny,nz,nt))

    for frame in range(nt):
        for iz in range(nz):
            if option=="m":
                # Project data onto vertical modes:
                for jz in range(nz):
                    offset=frame*N+jz*NH+1
                    data_array[:,:,iz,frame]+=vec[jz,iz]*field_data[offset:offset+NH].reshape(nx,ny)
                if is_pv and pv_vis=="a" and iz==0:
                    data_array[:,:,0,frame]=data_array[:,:,0,frame]-bety

            else:
                # Use existing layer data:
                offset=frame*N+iz*NH+1
                data_array[:,:,iz,frame]=field_data[offset:offset+NH].reshape(nx,ny)
                if is_pv and pv_vis=="a":
                    data_array[:,:,iz,frame]=data_array[:,:,iz,frame]-bety

            layer_min[iz]=min(layer_min[iz],np.min(data_array[:,:,iz,frame]))
            layer_max[iz]=max(layer_max[iz],np.max(data_array[:,:,iz,frame]))

    return layer_min,layer_max,data_array

# Ensure the script is called with the correct number of arguments
if len(sys.argv) != 3:
    print("Usage: python p_temporal_avg.py <PV option: a/t> <View option: l/m>")
    sys.exit(1)

# Read the arguments
pv_vis = sys.argv[1].lower()  # First argument: PV option ('a' or 't')
option = sys.argv[2].lower()  # Second argument: View option ('l' or 'm')

# Validate the arguments
if pv_vis not in ['a', 't']:
    print("Error: PV option must be 'a' (anomaly) or 't' (total).")
    sys.exit(1)

if option not in ['l', 'm']:
    print("Error: View option must be 'l' (layers) or 'm' (modes).")
    sys.exit(1)

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

# Increase nx & ny by 1 to include boundary points:
nx=nx+1
ny=ny+1

if pv_vis == "a":
    with open('src/parameters.f90','r') as in_file:
        fread=in_file.readlines()
        for line in fread:
            if ':: beta=' in line:
                beta=float(line.split("=")[1].split(",")[0])
                bety=np.tile(beta*np.linspace(ymin,ymax,ny),(nx,1))

# Read energy data to get final time in data:
with open('evolution/energy.asc','r') as in_file:
    time,ekin,epot,etot=np.loadtxt(in_file,dtype=float,unpack=True)
tsim=time[-1]

# Read in vertical mode matrix if required:
vec=np.empty((nz,nz))
if option == "m":
    in_file=open('modes.asc','r')
    for m in range(nz):
        line=in_file.readline()
    for m in range(nz):
        for iz in range(nz):
            line=in_file.readline()
            string=line.split()
            vec[iz,m]=string[0]
    in_file.close()

#=================================================================
# Set up figure:
# Layout is nz rows and 3 columns:
nrow=nz
ncol=3
nim=nrow*ncol

# The following figure might need to be adjusted depending on the problem
fig,ax=plt.subplots(figsize=[12+1*(nz>2),1+3*nrow],nrows=nrow,ncols=ncol)
ax=ax.flatten()

# Set up an array to store all images to be plotted:
d=np.empty((nim,nx,ny))

# Number of horizontal grid points:
NH=nx*ny

# Total number of data elements per time frame:
N=nz*NH+1
# The "+1" includes the time element

# Field titles over each column:
field=['$\\psi$','$q$','$\\zeta$']

min_vals=np.full((nz,ncol), np.inf)
max_vals=np.full((nz,ncol),-np.inf)

#=================================================================
# Read streamfunction data into array for plotting:
with open('evolution/pp.r4','rb') as in_file:
    pp_array=np.fromfile(in_file,dtype=np.float32)

# Read PV data into array for plotting:
with open('evolution/qq.r4','rb') as in_file:
    qq_array=np.fromfile(in_file,dtype=np.float32)

# Read relative vorticity data into array for plotting:
with open('evolution/zz.r4','rb') as in_file:
    zz_array=np.fromfile(in_file,dtype=np.float32)

nt=len(pp_array) // N

min_vals[:,0],max_vals[:,0],pp_array=get_layer_min_max(pp_array)
min_vals[:,1],max_vals[:,1],qq_array=get_layer_min_max(qq_array,True)
min_vals[:,2],max_vals[:,2],zz_array=get_layer_min_max(zz_array)

# compute temporal averages for each layer
for iz in range(nz):
    d[iz*ncol  ,:,:]=np.mean(pp_array[:,:,iz,:], axis=-1)
    d[iz*ncol+1,:,:]=np.mean(qq_array[:,:,iz,:], axis=-1)
    d[iz*ncol+2,:,:]=np.mean(zz_array[:,:,iz,:], axis=-1)

im=[] # To store the images
cbar=[] # To store the colorbars
cax=[0]*nim
zmin_arr=[0]*nim
zmax_arr=[0]*nim
clevels_arr=[[] for _ in range(nim)]

for j in range(nim):
    ax[j].set_aspect('equal')
    row=int(j/ncol)
    col=j-ncol*row

    # Work out the overall min/max values:
    zmin=min_vals[row,col]
    zmax=max_vals[row,col]

    zmag=max(abs(zmin),abs(zmax))
    zmin_arr[j]=-zmag
    zmax_arr[j]=zmag

    dz=2.0*contint(zmin,zmax)
    jj=int(zmag/dz)
    z=dz*float(jj)
    clevels_arr[j]=np.linspace(-z,z,max(3,jj+1))

    ax[j].set_xlim([xmin,xmax])
    ax[j].set_ylim([ymin,ymax])

    # Label x axis only for images in the bottom row:
    if row < nz - 1:
        plt.setp(ax[j].get_xticklabels(),visible=False)
    else:
        ax[j].set_xlabel('$x$',fontsize=20)

    # Label y axis only for images in the leftmost column:
    if col > 0:
        plt.setp(ax[j].get_yticklabels(),visible=False)
    else:
        ax[j].set_ylabel('$y$',fontsize=20)

    if row == 0:
        ax[j].set_title(field[col],fontsize=36)

    # Initialize imshow with placeholder data
    img=ax[j].imshow( d[j,:,:].T,cmap=cmap_wbgyr(),vmin=zmin_arr[j],vmax=zmax_arr[j],
                      extent=(xmin,xmax,ymin,ymax),origin='lower',interpolation='bilinear' )

    # Create colorbar for each subplot (only once)
    divider=make_axes_locatable(ax[j])
    cax[j]=divider.append_axes("right",size="4%",pad=0.1)
    cbar_j=fig.colorbar(img,cax=cax[j],ticks=clevels_arr[j])

    im.append(img)
    cbar.append(cbar_j)

#=========================================================================
# Adjust figure size to make it divisible by 2
# Get the current size of the figure in pixels
fig.canvas.draw()
width,height=fig.canvas.get_width_height()

# Check if width and height are divisible by 2
new_width=width if width % 2 == 0 else width-1
new_height=height if height % 2 == 0 else height-1

# Resize the figure if necessary
if new_width != width or new_height != height:
    fig.set_size_inches(new_width/fig.dpi,new_height/fig.dpi)
#=========================================================================

fig.subplots_adjust(wspace=0.0,hspace=0.0)
fig.savefig(f"temporal_avg.png",dpi=150)
