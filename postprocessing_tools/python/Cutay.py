import numpy as np

# Uncommnt to make plots only on file
import matplotlib
matplotlib.use('Agg')

import matplotlib.cm as cm
import matplotlib.mlab as mlab # for setting up the data
import matplotlib.pyplot as plt
from scipy import ndimage
import h5py


""" core routines
"""

#==============================================================================
#  Read Variable
#==============================================================================

def read_variables(f):
    global p, b, gb, db, gyro_radius, vth, e, ex, ey, ez, bx, by, bz, az, epar, rho, j
    Bx = np.array(f['/Step#0/Block/Bx/0'])
    By = np.array(f['/Step#0/Block/By/0'])
    Bz = np.array(f['/Step#0/Block/Bz/0'])
    Nz=Bx.shape[0]
    Ny=Bx.shape[1]
    Nx=Bx.shape[2]
    
    dx = Lx/(Nx-1.)
    dy = Ly/(Ny-1.)

    
    print(Nx, Ny, Nz)
    
    Ncut= int(Nz/2)
    print(Ncut)
    Ncut=0
    
    Nxstart=0 
    Nxend= Nx 
    
    Nystart=0
    Nyend= Ny
    
    bx = Bx[Ncut][:][:]
    by = By[Ncut][:][:]
    bz = Bz[Ncut][:][:]

    b = np.sqrt(pow(bx,2) + pow(by,2) + pow(bz,2))

    #dBx = np.array(f['/Step#0/Block/Bx/0'])
    #dBy = np.array(f['/Step#0/Block/By/0'])
    #dBz = np.array(f['/Step#0/Block/Bz/0'])
    #dbx = dBx[Ncut][:][:]
    #dby = dBy[Ncut][:][:]
    #dbz = dBz[Ncut][:][:]
    #db = np.sqrt(pow(dbx,2) + pow(dby,2) + pow(dbz,2))
    
    az = vecpot(dx,dy,bx,by)

    #gbx,gby = np.gradient(db)
    
    #gb = np.sqrt(pow(gbx,2) + pow(gby,2))/b

    Ex = np.array(f['/Step#0/Block/Ex/0'])
    Ey = np.array(f['/Step#0/Block/Ey/0'])
    Ez = np.array(f['/Step#0/Block/Ez/0'])
    ex = Ex[Ncut][:][:]
    ey = Ey[Ncut][:][:]
    ez = Ez[Ncut][:][:]
    e = np.sqrt(pow(ex,2) + pow(ey,2) + pow(ez,2))
    
    epar= ex*bx+ey*by+ez*bz
    
    Jx = np.array(f['/Step#0/Block/Jx_0/0']) + np.array(f['/Step#0/Block/Jx_2/0'])
    Jy = np.array(f['/Step#0/Block/Jy_0/0']) + np.array(f['/Step#0/Block/Jy_2/0'])
    Jz = np.array(f['/Step#0/Block/Jz_0/0']) + np.array(f['/Step#0/Block/Jz_2/0'])
    jex = Jx[Ncut][:][:]
    jey = Jy[Ncut][:][:]
    jez = Jz[Ncut][:][:]
    
    j = np.sqrt(pow(jex,2) + pow(jey,2) + pow(jez,2))
    
    Pxx = np.array(f['/Step#0/Block/Pxx_1/0'])+np.array(f['/Step#0/Block/Pxx_3/0'])
    Pyy = np.array(f['/Step#0/Block/Pyy_1/0'])+np.array(f['/Step#0/Block/Pyy_3/0'])
    Pzz = np.array(f['/Step#0/Block/Pzz_1/0'])+np.array(f['/Step#0/Block/Pzz_3/0'])

    Rho = np.array(f['/Step#0/Block/rho_1/0'])+np.array(f['/Step#0/Block/rho_3/0'])
    
    
    pxx = Pxx[Ncut][:][:]
    pyy = Pyy[Ncut][:][:]
    pzz = Pzz[Ncut][:][:]
    rho = Rho[Ncut,:,:]*4*np.pi
    
    p = pxx + pyy + pzz
    
    qom=1
    c=1
    wc=qom*b/c
#    vth=np.sqrt(np.abs(p/(rho+1e-10)*qom))
#    gyro_radius=vth/wc

    return True

def vecpot(dx,dy,bx,by):
    nx = bx.shape[0]
    ny = bx.shape[1]
    nymezzo = int(np.ceil(ny/2))
    print('vecpot',nx,ny,nymezzo)
    az=np.zeros((nx,ny))
    for i in range(1,nx):
      az[nymezzo,i] = az[nymezzo,i-1]+ (by[nymezzo,i-1]+by[nymezzo,i])*dx/2.0

    for ind in range(nymezzo+1,ny):
        for j in range(0,nx):
            az[ind,j]=az[ind-1,j]- (bx[ind,j] + bx[ind-1,j])*dy/2.0
	
    for ind in range(nymezzo-1,-1,-1):
        for j in range(0,nx):
            az[ind,j]=az[ind+1,j]+ (bx[ind+1,j] + bx[ind,j])*dy/2.0

    return az
  
import math
def round_to_n(x, n):
    if not x: return 0
    power = -int(math.floor(math.log10(abs(x)))) + (n - 1)
    factor = (10 ** power)
    return round(x * factor) / factor
	
#==============================================================================
#  Plot Variable
#==============================================================================
def plot_variables(time_label, Lx, Ly, variable, az, name_variable, minval, maxval):
	levels = 20

	variablesm=ndimage.filters.gaussian_filter(variable,2, mode='nearest')
	if minval==0 and maxval==0:
		minval = variable.min()
		maxval = variable.max()
		
	print(minval,maxval)
	if minval==0 and maxval==0:
		return False
# plot the filled contour
# using a colormap (jet)
#    CF = plt.contourf(gyro_radius, levels=np.linspace(1,100, levels),
#    CF = plt.contourf(vth, levels=np.linspace(0,.005, levels),
	#CF = plt.contourf(variable, levels=np.linspace(minval,maxval, levels),
        #          extent=(0,Lx,0,Ly),cmap=cm.jet)
	#CF = plt.pcolor(variable, cmap='RdBu', vmin=minval, vmax=maxval)
	
	CF = plt.imshow(variable, cmap=cm.jet, vmin=minval, vmax=maxval,
	   extent=[0, Lx, 0, Ly],
	   interpolation='nearest', origin='lower')
# plot the contour lines
# using gray scale
	levels = 60
	CL = plt.contour(az, levels,
                 linewidths=0.1,
                 extent=(0,Lx,0,Ly),cmap=cm.gray)

# plot color bars for both contours (filled and lines)
	CB = plt.colorbar(CL, extend='both')
	mn=minval      # colorbar min value
	mx=maxval       # colorbar max value
	md=round_to_n((mx+mn)/2,2)                     # colorbar midpoint value
	print('round',mn,md,mx)
	CBI = plt.colorbar(CF, ticks=[mn, md,mx],orientation='horizontal')
	CBI.set_ticklabels([mn,md,mx])

# Plotting the second colorbar makes
# the original colorbar look a bit out of place,
# so let's improve its position.

	l,b,w,h = plt.gca().get_position().bounds
	ll,bb,ww,hh = CB.ax.get_position().bounds
	CB.ax.set_position([ll, b, ww, h])

	plt.axes().set_aspect('equal')
	print('plotted')
	plt.show()
	fname=name_variable+time_label+'.png' # .svg also works very well
	plt.savefig(fname, dpi=300) #1200 for very good rez
	plt.clf()
	return True


#==============================================================================
#  Main program
#==============================================================================
#0 6000 ogni 100

#directory='/discover/nobackup/jpark35/ecsim/build/wb7_0'
directory='cyl2_5/'
directory='data/'
ini=800
nt=ini
dt=10
Lx=40
Ly=30

for i in range(ini,nt+dt,dt):
    print(i)
    time_label = str(i).zfill(6)
    filename = directory+'GEM-Fields_' + time_label +'.h5'
    #print filename
    f = h5py.File(filename, 'r')
    read_variables(f)

    plot_variables(time_label,Lx, Ly, rho, az, 'rho',0,0)#-3.5e-7,3.5e-7)
    
print('done')
