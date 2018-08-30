from scipy.io import loadmat
import numpy as np
from sklearn.cluster import KMeans
x = loadmat('bimodal_ions_FAC.mat')

fcutoff = x['fcutoff']

Nx = fcutoff.shape[0]
Ny = fcutoff.shape[1]
Nz = fcutoff.shape[2]

print(Nx,Ny,Nz)
eg=np.linspace(-560,560,Nx)

xg, yg, zg = np.meshgrid(eg,eg,eg, indexing='ij')

print('x',xg[0:4,0,0])

print('y',yg[0,0:4,0])

print('z',zg[0,0,0:4])

Uxdata=np.mean(xg*fcutoff)/np.mean(fcutoff)
Uydata=np.mean(yg*fcutoff)/np.mean(fcutoff)
Uzdata=np.mean(zg*fcutoff)/np.mean(fcutoff)


print(Uxdata,Uydata,Uzdata)

a=10*fcutoff/np.max(fcutoff);

X=np.zeros((1000000,3));
ip=0
for i in range (0,Nx):
    for j in range(0,Ny):
        for k in range (0,Nz):
            Nv=int(np.floor(a[i,j,k]))
            r1= np.random.rand(1)
            if a[i,j,k]-Nv>r1:
                Nv=Nv+1   
            #print(a[i,j,k]-Nv,r1, Nv)
            #if Nv > 0:
            #    print('Nv',Nv)
            for ip2 in range (0,Nv):
            #    b=np.array([i, j, k])
            #    X=np.vstack((X, b))
            #X=np.vstack((X, np.tile(b,(Nv,1))));
                X[ip,0]=i
                X[ip,1]=j
                X[ip,2]=k
                ip = ip +1
    print( i,ip) #,  X.shape)
X=X[:ip-1,:]/(Nx-1)*560*2-560
kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
print(kmeans.cluster_centers_)