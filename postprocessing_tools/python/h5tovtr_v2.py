#!/usr/bin/python

from evtk.vtk import VtkFile, VtkRectilinearGrid

import numpy as np
import h5py
import sys

# Set the basename, range and final filename

h5filename  = str(sys.argv[1])
istep       = int(sys.argv[2])
istart      = int(sys.argv[3])
iend        = int(sys.argv[4])
vtkfilename = str(sys.argv[5])

# Loop over the simulation cycles

for i in np.arange(istart, iend+istep, istep):

  # Open the HDF5 file

  print "Processing file", i

  cycle = "%05d" % i
  file  = h5py.File(h5filename+"_"+cycle+".h5",'r')

  # Dimensions 
  
  ncellxyz = file['/Parameters/ncell']
  nx, ny, nz = ncellxyz[...]
  
  LxLyLz = file['/Parameters/LxLyLz']
  lx, ly, lz = LxLyLz[...]
  dx, dy, dz = lx/nx, ly/ny, lz/nz 
  
  ncells  = nx * ny * nz 
  npoints = (nx + 1) * (ny + 1) * (nz + 1) 
  
  # Coordinates 
  
  xc = np.arange(0, lx + 0.1*dx, dx, dtype='float64') 
  yc = np.arange(0, ly + 0.1*dy, dy, dtype='float64') 
  zc = np.arange(0, lz + 0.1*dz, dz, dtype='float64') 
  
  # Read the datasets to a 'data' vector
  
  data1  = file['/Fields/Bx'][...]
  data2  = file['/Fields/By'][...]
  data3  = file['/Fields/Bz'][...]
  data4  = file['/Fields/Ex'][...]
  data5  = file['/Fields/Ey'][...]
  data6  = file['/Fields/Ez'][...]
  data7  = file['/Fields/Rho_0'][...]
  data8  = file['/Fields/Rho_1'][...]
  data9  = file['/Fields/Jx_0'][...]
  data10 = file['/Fields/Jy_0'][...]
  data11 = file['/Fields/Jz_0'][...]
  data12 = file['/Fields/Jx_1'][...]
  data13 = file['/Fields/Jy_1'][...]
  data14 = file['/Fields/Jz_1'][...]
  
  # Write the VTK file

  start, end  = (0, 0, 0), (nx, ny, nz)
  w = VtkFile(vtkfilename+"_"+cycle, VtkRectilinearGrid)
  w.openGrid (start = start, end = end)
  w.openPiece(start = start, end = end)

  # Cell data
  w.openData("Cell", scalars = "Rho_0, Rho_1", vectors = "B, E, J_0, J_1")
  w.addData ("Rho_0", data7)
  w.addData ("Rho_1", data8)
  w.addData ("B",   (data1,  data2,  data3 ))
  w.addData ("E",   (data4,  data5,  data6 ))
  w.addData ("J_0", (data9,  data10, data11))
  w.addData ("J_1", (data12, data13, data14))
  w.closeData("Cell")

  # Coordinates of cell vertices
  w.openElement("Coordinates")
  w.addData("x_coordinates", xc)
  w.addData("y_coordinates", yc)
  w.addData("z_coordinates", zc)
  w.closeElement("Coordinates")

  w.closePiece()
  w.closeGrid()

  # Apped the data to the file
  w.appendData(data = data7)
  w.appendData(data = data8)
  w.appendData(data = (data1,  data2,  data3 ))
  w.appendData(data = (data4,  data5,  data6 ))
  w.appendData(data = (data9,  data10, data11))
  w.appendData(data = (data12, data13, data14))
  w.appendData(xc).appendData(yc).appendData(zc)
  w.save()
