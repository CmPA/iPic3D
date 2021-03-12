import numpy as np

# Uncommnt to make plots only on file
import matplotlib
matplotlib.use('Agg')

import matplotlib.cm as cm
import matplotlib.mlab as mlab # for setting up the data
import matplotlib.pyplot as plt
from scipy import ndimage
import h5py
import os
import sys


""" core routines
"""

#==============================================================================
#  Read Variable
#==============================================================================

def read_variables(f,fB):
    global dx, dy
    global rho,p, b, gb, db, gyro_radius, vth, e, ex, ey, ez, bx, by, bz, az, jx, jy, jz, rhoc, rho0, rho1, rho2, rho3, bx_ext, by_ext, bz_ext, lamb, divB
    global rhoc, divE
    global spec, uth, vth, wth, small
    global Pxx, Pyy, Pzz, Pxy, Pxz, Pyz, Jx, Jy, Jz, Rho, qom, rhoc_avg, divE_avg
    
    Bx = np.array(f['/Step#0/Block/Bx/0'])
    #Bx_ext = np.array(f['/Step#0/Block/Bx_ext/0'])
    By = np.array(f['/Step#0/Block/By/0'])
    #By_ext = np.array(f['/Step#0/Block/By_ext/0'])
    Bz = np.array(f['/Step#0/Block/Bz/0'])
    #Bz_ext = np.array(f['/Step#0/Block/Bz_ext/0'])
    Nz=Bx.shape[0]
    Ny=Bx.shape[1]
    Nx=Bx.shape[2]
    
    dx = Lx/(Nx-1.)
    dy = Ly/(Ny-1.)

    
    #print(Nx, Ny, Nz)
    
    Ncut= int(Nz/2)
    #print(Ncut)
    Ncut=0
    
    Nxstart=0 
    Nxend= Nx 
    
    Nystart=0
    Nyend= Ny
    
    bx = Bx[Ncut][:][:]
    by = By[Ncut][:][:]
    bz = Bz[Ncut][:][:]
    
    #bx_ext = Bx_ext[Ncut][:][:]
    #by_ext = By_ext[Ncut][:][:]
    #bz_ext = Bz_ext[Ncut][:][:]

    #print(bx.shape[0],bx.shape[1])
    b = np.sqrt(pow(bx,2) + pow(by,2) + pow(bz,2))+small


    
    az = vecpot(dx,dy,bx,by)

    gbx,gby = np.gradient(b)
    
    #gb = np.sqrt(pow(gbx,2) + pow(gby,2))/b

    Ex = np.array(f['/Step#0/Block/Ex/0'])
    Ey = np.array(f['/Step#0/Block/Ey/0'])
    Ez = np.array(f['/Step#0/Block/Ez/0'])
    ex = Ex[Ncut][:][:]
    ey = Ey[Ncut][:][:]
    ez = Ez[Ncut][:][:]
    e = np.sqrt(pow(ex,2) + pow(ey,2) + pow(ez,2))

    Pxx = np.array(f['/Step#0/Block/Pxx_'+spec+'/0']) 
    Pyy = np.array(f['/Step#0/Block/Pyy_'+spec+'/0']) 
    Pzz = np.array(f['/Step#0/Block/Pzz_'+spec+'/0']) 
    Pxy = np.array(f['/Step#0/Block/Pxy_'+spec+'/0']) 
    Pxz = np.array(f['/Step#0/Block/Pxz_'+spec+'/0']) 
    Pyz = np.array(f['/Step#0/Block/Pyz_'+spec+'/0']) 

    Jx = np.array(f['/Step#0/Block/Jx_'+spec+'/0'])
    Jy = np.array(f['/Step#0/Block/Jy_'+spec+'/0'])
    Jz = np.array(f['/Step#0/Block/Jz_'+spec+'/0'])
    
    Rho = np.array(f['/Step#0/Block/rho_'+spec+'/0'])
    Rho0 = np.array(f['/Step#0/Block/rho_'+'0'+'/0'])
    if ns>1:
        Rho1 = np.array(f['/Step#0/Block/rho_'+'1'+'/0'])
    #Rho2 = np.array(f['/Step#0/Block/rho_'+'2'+'/0'])
    #Rho3 = np.array(f['/Step#0/Block/rho_'+'3'+'/0'])
    
    #Lambda = np.array(f['/Step#0/Block/Lambda/0'])
    
    compute_pressure()
    

    pxx = Pxx[Ncut][:][:]
    pyy = Pyy[Ncut][:][:]
    pzz = Pzz[Ncut][:][:]
    rho = Rho[Ncut][:][:]+small
    
    rho0 = Rho0[Ncut][:][:]
    if ns>1:
        rho1 = Rho1[Ncut][:][:]
        rhoc = Rho0[Ncut][:][:] + Rho1[Ncut][:][:]
    #rho2 = Rho2[Ncut][:][:]
    #rho3 = Rho3[Ncut][:][:]
    uth = np.sqrt(abs(qom*pxx/rho))
    vth = np.sqrt(abs(qom*pyy/rho))
    wth = np.sqrt(abs(qom*pzz/rho))
 #   p = pxx + pyy + pzz
    

    jx = Jx[Ncut][:][:]
    jy = Jy[Ncut][:][:]
    jz = Jz[Ncut][:][:]
    
    #lamb = Lambda[Ncut][:][:]
    
    #DivB = np.array(f['/Step#0/Block/divB/0'])
    #divB = DivB[Ncut][:][:]
    
    #Rhoc = np.array(fB['/Step#0/Block/rhoc_avg/0'])
    #DivE = np.array(fB['/Step#0/Block/divE_avg/0'])
    #rhoc = Rhoc[Ncut][:][:]
    #divE = DivE[Ncut][:][:]
        
    c=1
    wc=qom*b/c
#    vth=np.sqrt(np.abs(p/(rho+1e-10)*qom))
#    gyro_radius=vth/wc

#Rhoc_avg = np.array(f['/Step#0/Block/rho_avg/0'])
#DivE_avg = np.array(f['/Step#0/Block/divE/0'])
#rhoc_avg = Rhoc_avg[Ncut][:][:]
#divE_avg = DivE_avg[Ncut][:][:]

    return True
    
def compute_pressure():
	global Pxx, Pyy, Pzz, Pxy, Pxz, Pyz, Jx, Jy, Jz, Rho, qom, small
	Pxx = (Pxx - Jx*Jx / (Rho+small) ) /qom;
	Pyy = (Pyy - Jy*Jy / (Rho+small) ) /qom;
	Pzz = (Pzz - Jz*Jz / (Rho+small) ) /qom;
	Pxy = (Pxy - Jx*Jy / (Rho+small) ) /qom;
	Pxz = (Pxz - Jx*Jz / (Rho+small) ) /qom;
	Pyz = (Pyz - Jy*Jz / (Rho+small) ) /qom;
	return True
	
def vecpot(dx,dy,bx,by):
    nx = bx.shape[0]
    ny = bx.shape[1]
    nymezzo = int(np.ceil(ny/2))
    #print('vecpot',nx,ny)
    az=np.zeros((nx,ny))
    for i in range(1,nx):
      az[i][nymezzo] = az[i-1][nymezzo]- (bx[i-1][nymezzo]+bx[i][nymezzo])*dy/2.0

    for ind in range(nymezzo+1,ny):
        for j in range(0,nx):
            az[j,ind]=az[j][ind-1]+ (by[j][ind-1] + by[j][ind-1])*dx/2.0
	
    for ind in range(nymezzo-1,-1,-1):
        for j in range(0,nx):
            az[j,ind]=az[j][ind+1]- (by[j][ind+1] + by[j][ind])*dx/2.0
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
	global directory
	levels = 20
	#variable=ndimage.filters.gaussian_filter(variable, 2, mode='nearest')
	if minval==0 and maxval==0:
		minval = round_to_n(variable.min(),4)
		maxval = round_to_n(variable.max(),4)
		
	print(name_variable, minval,maxval)
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
	plt.title(name_variable+'   Cycle='+time_label)

# plot color bars for both contours (filled and lines)
#CB = plt.colorbar(CL, extend='both')
	mn=minval      # colorbar min value
	mx=maxval       # colorbar max value
	md=round_to_n((mx+mn)/2,2)                     # colorbar midpoint value
	#print('round',mn,md,mx)
	CBI = plt.colorbar(CF, ticks=[mn, md,mx],orientation='horizontal')
	CBI.set_ticklabels([mn,md,mx])

# Plotting the second colorbar makes
# the original colorbar look a bit out of place,
# so let's improve its position.

#l,b,w,h = plt.gca().get_position().bounds
#ll,bb,ww,hh = CB.ax.get_position().bounds
#CB.ax.set_position([ll, b, ww, h])

	#plt.axes().set_aspect('equal')
	#print('plotted')
	#plt.show()
	fname=dirout+name_variable+time_label+'.png' # .svg also works very well
	#plt.savefig(fname, dpi=1200)
	plt.savefig(fname, dpi=300,bbox_inches='tight')
	plt.clf()
	return True


#==============================================================================
#  Main program
#==============================================================================
#0 6000 ogni 100

directory='data'+sys.argv[2]+'/'
dirout='png'+sys.argv[2]+'/'
returned_value = os.system('mkdir '+dirout)
nt=int(sys.argv[1])
ini=nt
dt=100

B0 = 0.0097
vthe = 0.045
n0 = 1

ns = 1

Lx=10 #99.58
Ly=10 #200
# label of species
spec='0' # 0=electrons, 1=ions
# enter correct qom here
qom=-256
small=1e-50*qom;

de = np.sqrt(n0/-qom)
rhoe = vthe / (-B0*qom)



for i in range(ini,nt+dt,dt):
    #print(i)
    time_label = str(i).zfill(6)
    filename = directory+'DoubleHarris-Fields_' + time_label +'.h5'
    filenameB = directory+'DoubleHarris-fieldB_' + time_label +'.h5'
    print(filename)
    f = h5py.File(filename, 'r')
    fB = f
    #fB = h5py.File(filenameB, 'r')
    spec='0' # 0=electrons, 1=ions
    read_variables(f,fB)
    
    print('Scales de=', de, '  rhoe=',rhoe,'   dx=',dx)


    plot_variables(time_label,Lx, Ly, az, az, 'az',0,0)
    #plot_variables(time_label,Lx, Ly, lamb, az, 'lambda',0,0)
    #plot_variables(time_label,Lx, Ly, divB, az, 'divB',0,0)
    
    #plot_variables(time_label,Lx, Ly, divE, divE, 'divE',0,0)
    #plot_variables(time_label,Lx, Ly, rhoc, rhoc*4*np.pi, 'rhoc',0,0)
        
    plot_variables(time_label,Lx, Ly, bx, az, 'bx',-0.01,0.01)
    plot_variables(time_label,Lx, Ly, by, az, 'by',0,0)
    plot_variables(time_label,Lx, Ly, bz, az, 'bz',0,0)
    #plot_variables(time_label,Lx, Ly, bx_ext, az, 'bx_ext',0,0)
    #plot_variables(time_label,Lx, Ly, by_ext, az, 'by_ext',0,0)
    #plot_variables(time_label,Lx, Ly, bz_ext, az, 'bz_ext',0,0)
    
    plot_variables(time_label,Lx, Ly, ex, az, 'ex',0,0)
    plot_variables(time_label,Lx, Ly, ey, az, 'ey',0,0)
    plot_variables(time_label,Lx, Ly, ez, az, 'ez',0,0)

    for ispec in range(0,ns):
        spec=str(ispec).zfill(1) # 0=electrons, 1=ions

        read_variables(f,fB)
    
        plot_variables(time_label,Lx, Ly, uth, az, 'uth'+spec,0,0)
        plot_variables(time_label,Lx, Ly, vth, az, 'vth'+spec,0,0)
        plot_variables(time_label,Lx, Ly, wth, az, 'wth'+spec,0,0)


        plot_variables(time_label,Lx, Ly, uth, az, 'uth'+spec,0,0)
        plot_variables(time_label,Lx, Ly, vth, az, 'vth'+spec,0,0)
        plot_variables(time_label,Lx, Ly, wth, az, 'wth'+spec,0,0)

    
        plot_variables(time_label,Lx, Ly, jx, az, 'jx'+spec,0,0)
        plot_variables(time_label,Lx, Ly, jy, az, 'jy'+spec,0,0)
        plot_variables(time_label,Lx, Ly, jz, az, 'jz'+spec,0,0)
        plot_variables(time_label,Lx, Ly, jz/rho, az, 'Vz'+spec,0,0)
        plot_variables(time_label,Lx, Ly, rho, az, 'rho'+spec,0,0)
        #rhosc=(rho[1:,1:]+rho[1:,:-1]+rho[-1:,1:]+rho[-1:,-1:])/4
        #print('Sum rhoc =', np.sum(rhoc),'Sum rho =', np.sum(rho))
        #plot_variables(time_label,Lx, Ly, rhoc, rhoc/rhosc, 'rhocOrho',0,0)

print('done')
montage = False
if(montage):
    os.system("montage data/bx*.png -density 300 -tile 2x3  -geometry 1000 tile_bx.pdf")
    os.system("montage data/by*.png -density 300 -tile 2x3  -geometry 1000 tile_by.pdf")
    os.system("montage data/bz*.png -density 300 -tile 2x3  -geometry 1000 tile_bz.pdf")
    os.system("montage data/uth0*.png -density 300 -tile 2x3  -geometry 1000 tile_uth0.pdf")
    os.system("montage data/uth1*.png -density 300 -tile 2x3  -geometry 1000 tile_uth1.pdf")
    os.system("montage data/vth0*.png -density 300 -tile 2x3  -geometry 1000 tile_vth0.pdf")
    os.system("montage data/vth1*.png -density 300 -tile 2x3  -geometry 1000 tile_vth1.pdf")
    os.system("montage data/wth0*.png -density 300 -tile 2x3  -geometry 1000 tile_wth0.pdf")
    os.system("montage data/wth1*.png -density 300 -tile 2x3  -geometry 1000 tile_wth1.pdf")
    os.system("montage data/uth2*.png -density 300 -tile 2x3  -geometry 1000 tile_uth2.pdf")
    os.system("montage data/uth3*.png -density 300 -tile 2x3  -geometry 1000 tile_uth3.pdf")
    os.system("montage data/vth2*.png -density 300 -tile 2x3  -geometry 1000 tile_vth2.pdf")
    os.system("montage data/vth3*.png -density 300 -tile 2x3  -geometry 1000 tile_vth3.pdf")
    os.system("montage data/wth2*.png -density 300 -tile 2x3  -geometry 1000 tile_wth2.pdf")
    os.system("montage data/wth3*.png -density 300 -tile 2x3  -geometry 1000 tile_wth3.pdf")

