################################################################################
################################################################################
####                                                                        ####
####  Compute the omega-k spectrum from the field files (H5hut)             ####
####                                                                        ####
####  Diego Gonzalez                                                        ####
####  February 2018                                                         ####
####                                                                        ####
####------------------------------------------------------------------------####
####  USAGE:                                                                ####
####    Set the name  an the path of the results                            ####
####    Set the initial and the last cycle and the step between files       ####
####    Set the time step, the size of the domain and other parameters      ####
####    used in the simulation.                                             ####
####    Select a field component and a direction to do the analysis.        ####
####    Execute the script                                                  ####
####                                                                        ####
################################################################################
################################################################################

from readhdf import DataSet, DataSetH5hut
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys

################################################################################
#---------- INPUT DATA --------------------------------------------------------

h5hut = True
path = "data_peri"
# Name of the files (without -Fields_######.h5)
name = "SolarWind"
#name = "/home/diego/Codes/ECsim/Results/Whistler/Paper/Case01/data03/Whistler"
# Inital cycle, last cycle, and step between files
cStart   = 0    
cEnd     = 800  
cStep    = 5
# Time step used in the simulation
dt  = 0.1
# Size of the domain used in the simulation
Lx  = 20.0
Ly  = 20.0
Lz  = 1.00

# Other simulation parameters
qom_ion = 1.0
qom_ele = 1836
d_ion = 1.0
d_ele = d_ion/np.sqrt(qom_ele)
wp_ion = 1./d_ion
wp_ele = 1./d_ele

# Field to use for the analysis
Key = "Bz"

# Direction for the K spectrum
direction = "Y"

# Scale of the plot:
# k: di, de
k_scale = "di" 
# w: wpi, wpe
w_scale = "wpi" 

# Range of the visualization: [max, min]
range_k = []
range_w = []

# Save the plot as png
savepng = True




################################################################################
##    	MAIN
################################################################################

if (h5hut):
  dataSet = DataSetH5hut(path, name)
  Ex = dataSet.readField("Ex", cStart)
else:
  dataSet = DataSet(path)
  Ex = dataSet.readField("/fields/Ex", cStart)

n_tp = int((cEnd - cStart)/cStep + 1)
n_sp0 = Ex.shape[0]
n_sp1 = Ex.shape[1]
n_sp2 = Ex.shape[2]
n_sp = 0
data = np.empty((n_tp, n_sp0, n_sp1, n_sp2))

axFT = [] # Axes to do the FT
axSU = [] # Axes to sum over
Dt = dt*cStep
dx = 0
label_x = ""
if (direction == "X"):
  n_sp = n_sp2
  dx = Lx/n_sp
  axFT = (0,3)
  axSU = (1,2)
  label_x = r'$k_x$'
elif (direction == "Y"):
  n_sp = n_sp1
  dx = Ly/n_sp
  axFT = (0,2)
  axSU = (1,3)
  label_x = r'$k_y$'
elif (direction == "Z"):
  n_sp = n_sp0
  dx = Lz/n_sp
  axFT = (0,1)
  axSU = (2,3)
  label_x = r'$k_z$'

K_max = np.pi/dx
W_max = np.pi/Dt

print("Kmax = " + str(K_max))
print("Wmax = " + str(W_max))

# Read the data
it = 0
for c in range(cStart,cEnd+1,cStep):
  print ("Cycle " +str(c))
  if (h5hut):
    values = dataSet.readField(Key,c)
  else:
    values = dataSet.readField("/fields/"+Key,c) 
  data[it,:,:,:] = values[:,:,:]
  it += 1

# Fourier transform
ft = np.fft.fft2(data, axes=axFT)
fts = np.fft.fftshift(ft)
spectrum = np.sum(np.abs(fts), axis=axSU)

# Plot
f2 = plt.figure()
plt.title(r'$k-\omega$ spectrum. (' + Key + ') ' + direction + '-direction') 

if (k_scale == "de"):
  K_max = K_max*d_ele
  plt.xlabel(label_x + r' $d_{ele}$')
elif (k_scale == "di"):
  K_max = K_max*d_ion
  plt.xlabel(label_x + r' $d_{ion}$')

if (w_scale == "wpe"):
  W_max = W_max/wp_ele
  plt.ylabel(r'$\omega / \omega_{pe}$')
elif (w_scale == "wpi"):
  W_max = W_max/wp_ion
  plt.ylabel(r'$\omega / \omega_{pi}$')

cax = plt.imshow(np.log(spectrum), cmap='gist_ncar', interpolation='nearest', extent = [-K_max, K_max,-W_max, W_max], aspect="auto")
#plt.yscale('log')
#plt.xscale('log')
if (len(range_k) > 0):
  plt.xlim(range_k)
if (len(range_w) > 0):
  plt.ylim(range_w)

f2.colorbar(cax, orientation="vertical")

plt.tight_layout()

if (savepng):
  plt.savefig( name + "_spectrum-"+Key+"-" + direction + ".png", format='png', dpi=300)
plt.show()

