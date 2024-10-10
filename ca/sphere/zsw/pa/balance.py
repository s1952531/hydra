#!/usr/bin/env python

#=================================================================
#   Searches for a balanced zonal spherical shallow water flow
#=================================================================

#=====perform various generic imports=====
import numpy as np
#=========================================

#-------------------------------------------------------------------------
print()
print(' This routine finds a balanced zonal shallow-water flow on a sphere.')
print()

gamma_in = input(' Enter gamma = 1/L_d (default: 20)? ')
gamma = float(gamma_in or 20.0)
gam2=(gamma/2.0)**2

Omega=2.0*np.pi
fpole=2.0*Omega
cgw=fpole/gamma

print()
print( ' We consider h_bar(a) = 1 - A*tanh(z_bar/w) where z_bar = sin(a).')

amp_in = input(' Enter A (default: 0.1)? ')
amp = float(amp_in or 0.1)

wid_in = input(' Enter w (default: 0.1)? ')
wid = float(wid_in or 0.1)
wi=1.0/wid

n_def = max(500,int(25.0*gamma+0.5),int(50.0*wi+0.5))
n_in = input(' Enter the resolution, n (default: '+str(n_def)+')? ')
n = int(n_in or n_def)
print()

nm2=n-2
nm1=n-1
np1=n+1

a=np.linspace(-np.pi/2.0,np.pi/2.0,np1)
da=a[1]-a[0]

ddai=0.5/da
da2i=1.0/da**2
tda2i=2.0*da2i

rb=np.cos(a)
rb4=(rb**2)**2
zb=np.sin(a)
hb=1.0-amp*np.tanh(wi*zb)
# r_bar*h_bar:
rbhb=rb*hb
# d(r_bar*h_bar)/da:
drbhbda=-zb*hb-rb*amp*wi/np.cosh(wi*zb)**2

#---------------------------------------------------------------------
# Iterate to find displacement d:

# Maximum error between successive guesses:
toler=abs(amp)*1.e-11

# Initial guess:
phi=a
phipre=a
dpda=np.ones(np1)
d2pda2=np.zeros(np1)

a0=np.zeros(n)
ap=np.zeros(n)
am=np.zeros(n)
htd=np.zeros(n)
etd=np.zeros(n)

relax=0.5
ddmax=1.0

kit=1
while ddmax > toler:
   r=np.cos(phi)
   r[0]=1.0
   r[n]=1.0
   ri=1.0/r
   r[0]=0.0
   r[n]=0.0
   r2=r*r
   r2i=ri*ri
   z=np.sin(phi)
   zdr=z*ri
   
   dpi=1.0/dpda
   dpi2=dpi**2

   tt=r2-rb4*r2i
   dd22=d2pda2*dpi2
   F=gam2*z*tt+rbhb*dpi*(dd22-zdr)-drbhbda*dpi2
   F[0]=0.0
   F[n]=0.0
   # dF/dphi:
   G0=gam2*r*(tt-2.0*z*z*(1.0+rb4*r2i**2))-rbhb*dpi*r2i
   # dF/dphi':
   G1=dpi2*(rbhb*(zdr-3.0*dd22)+2.0*drbhbda*dpi)
   # dF/dphi'':
   G2=rbhb*dpi*dpi2

   a0[1:n]=G0[1:n]-tda2i*G2[1:n]
   ap[1:nm1]=da2i*G2[1:nm1]+ddai*G1[1:nm1]
   am[2:n]=da2i*G2[2:n]-ddai*G1[2:n]

   htd[1]=1.0/a0[1]
   etd[1]=-ap[1]*htd[1]

   for j in range(2,nm1):
      htd[j]=1.0/(a0[j]+am[j]*etd[j-1])
      etd[j]=-ap[j]*htd[j]

   htd[nm1]=1.0/(a0[nm1]+am[nm1]*etd[nm2])

   dphi=np.zeros(np1)
   dphi[1]=-F[1]*htd[1]

   for j in range(2,n):
      dphi[j]=-(F[j]+am[j]*dphi[j-1])*htd[j]

   for j in range(nm2,0,-1):
      dphi[j]=etd[j]*dphi[j+1]+dphi[j]

   phi=phi+relax*dphi
   dpda[1:n]=ddai*(phi[2:np1]-phi[0:nm1])
   d2pda2[1:n]=da2i*(phi[2:np1]-2.0*phi[1:n]+phi[0:nm1])

   ddmax=np.amax(abs(phi-phipre))
   print (' Iter ',kit,'  Error = ',ddmax)

   if ddmax > 1.0:
      break
   phipre=phi
   kit+=1

#-----------------------------------------------------------
if ddmax > 1.0:
   print()
   print( ' *** Not converging!  Stopping!')
   print()

else:
   # Write out various results:

   # Displacement:
   dphi=phi-a
   # Zonal velocity u / Omega:
   u=rb*rb*ri-r
   # Dimensionless height:
   h=rb*hb*ri
   h[1:n]=h[1:n]/dpda[1:n]
   h[0]=(4.0*h[1]-h[2])/3.0
   h[n]=(4.0*h[nm1]-h[nm2])/3.0
   # Vorticity (h*q-f) / (2*Omega):
   zeta=h*zb/hb-z
   ofile='g'+str(gamma)+'w'+str(wid)+'A'+str(amp)+'.asc'
   out_file = open(ofile,'w')
   for j in range(np1):
      out_file.write("%15.12f %15.12f %15.12f %15.12f %15.12f \n" % \
                     (a[j], dphi[j], u[j], h[j], zeta[j]))
   out_file.close()

   # Pointwise error in equilibrium equation:
   out_file = open('F.asc','w')
   for j in range(np1):
      out_file.write("%15.12f %15.12f \n" % (F[j], a[j]))
   out_file.close()

   # Prepare initialisation file for pam.f90 to test equilibrium:
   out_file = open('init.asc','w')
   fac=amp/wi
   mass=zb[1:np1]-zb[0:n]-fac*np.log(np.cosh(zb[1:np1]*wi)/np.cosh(zb[0:n]*wi))
   for j in range(1,np1):
      qm=Omega*(zb[j]**2-zb[j-1]**2)
      out_file.write("%15.12f %15.12f %15.12f %15.12f \n" % \
                     (phi[j], mass[j-1], qm, 0.0))
   out_file.close()

   # Compute energy components and total:
   ekin=0.25*Omega**2*sum(mass*(u[0:n]**2+u[1:np1]**2))
   epot=0.25*cgw**2*sum(mass*(h[0:n]+h[1:np1]-2.0))
   etot=ekin+epot
   print()
   print('   Kinetic energy = ',ekin)
   print(' Potential energy = ',epot)
   print('     Total energy = ',etot)
   print()
   print(' Use init.asc to initialise pam.f90')
   print()
   print(' Latitude displacement phi-a, zonal velocity u and scaled height')
   print(' h/H versus a are listed in ',ofile)
   print()
