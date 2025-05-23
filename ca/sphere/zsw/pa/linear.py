#!/usr/bin/env python3

#=====perform various generic imports=====
import numpy as np
from numpy import linalg as LA
#=========================================

gamma_in = input(' Enter gamma = 1/L_d (default: 20)? ')
gamma = float(gamma_in or 20.0)
gsq=gamma**2

fpole=4.0*np.pi
cgw=fpole/gamma

n_in = input(' Enter n (default: 500)? ')
n = int(n_in or 500)

n1=n-1
n2=n-2

a=np.linspace(-np.pi/2.0,np.pi/2.0,n+1)
r=np.cos(a)
rsq=r*r

z=np.sin(a)
zsq=z*z

da=a[1]-a[0]
dasq=da*da

ah=np.linspace(-(np.pi-da)/2.0,(np.pi-da)/2.0,n)
rh=np.cos(ah)

Amat = np.zeros([n1,n1])
for j in range(n1):
   Amat[j,j]=(rh[j]+rh[j+1])/(dasq*r[j+1])+1.0/rsq[j+1]+gsq*zsq[j+1]

for j in range(1,n1):
   Amat[j-1,j]=-rh[j]/(dasq*r[j])
   Amat[j,j-1]=-rh[j]/(dasq*r[j+1])

#Compute eigenvalues (sig) and eigenvectors (vec):
sig, vec = LA.eig(Amat)

sfac=cgw/(2.0*np.pi)
sig = sfac*np.sqrt(sig)

freq_file = open('freq.asc','w')

print('')
print(' Minimum frequencies, omega/(2*pi):')

nlow=10
for m in range(nlow):

   i=np.argmin(sig)
   print (sig[i])
   freq_file.write("%f \n" % sig[i])

   if sum(r[1:n]*vec[:,i]) < 0.0:
      vec[:,i]=-vec[:,i]

   vnorm=np.sqrt(da*sum(r[1:n]*vec[:,i]**2))
   vec[:,i]=vec[:,i]/vnorm
      
   vect_file = open('d'+str(m+1)+'.asc','w')
   vect_file.write("%f %f \n" % (0.0, a[0]))
   for j in range(n1):
      vect_file.write("%f %f \n" % (vec[j,i], a[j+1]))
   vect_file.write("%f %f \n" % (0.0, a[n]))
   vect_file.close()

   vect_file = open('u'+str(m+1)+'.asc','w')
   vect_file.write("%f %f \n" % (0.0, a[0]))
   for j in range(n1):
      vect_file.write("%f %f \n" % (z[j+1]*vec[j,i], a[j+1]))
   vect_file.write("%f %f \n" % (0.0, a[n]))
   vect_file.close()

   vect_file = open('h'+str(m+1)+'.asc','w')
   vect_file.write("%f %f \n" % (-r[1]*vec[0,i]/(da*rh[0]), ah[0]))
   for j in range(1,n1):
      vect_file.write("%f %f \n" % ((r[j]*vec[j-1,i]-r[j+1]*vec[j,i])/(da*rh[j]), a[j]))
   vect_file.write("%f %f \n" % (r[n1]*vec[n2,i]/(da*rh[n1]), ah[n1]))
   vect_file.close()

   sig[i]=1.e9

freq_file.close()
