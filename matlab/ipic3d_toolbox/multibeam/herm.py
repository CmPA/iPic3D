import numpy as np
import math as mt
import matplotlib.pyplot as plt

Np =100000
Nh=80
vp=np.random.randn(Np,1)
vect=np.zeros(Nh)
m=np.arange(0,Nh,1)
for i in range(0,Nh):
	print(i)
	c=np.zeros(i)
	c=np.hstack((c,[1]))
	hv=np.polynomial.hermite.hermval(vp/np.sqrt(2),c)
	vect[i]=np.mean(hv)/mt.factorial(i+1)
	print(vect[i])

vt=np.sqrt(2)
k=90
nu=2*k-1
sig=np.sqrt(k*vt*vt/nu)
vpk=np.sqrt(nu*(pow(np.random.rand(Np,1),-2/nu) -1.0));
vpk=sig*vpk*np.cos(2*np.pi*np.random.rand(Np,1));
vectk=np.zeros(Nh)
for i in range(0,Nh):
	print(i)
	c=np.zeros(i)
	c=np.hstack((c,[1]))
	hv=np.polynomial.hermite.hermval(vpk/np.sqrt(2),c)
	vectk[i]=np.mean(hv)/mt.factorial(i+1)
	print(vectk[i])
    
ax= plt.subplot(2,2,1)
n, bins, patches = plt.hist(vp, 50, log=True, color='r')

plt.subplot(2,2,2)
n, bins, patches = plt.hist(vpk, 50, log=True, color='b')

plt.subplot(2,2,3)
n, bins, patches = plt.hist(vp, 50, log=True, color='r')
n, bins, patches = plt.hist(vpk, bins=bins, log=True, color='b',fill=False)
plt.subplot(2,2,4)
plt.semilogy(m,np.abs(vect),'r',m,np.abs(vectk),'b')
plt.ylim((pow(10,-10),pow(10,0)))
plt.show()

