################################################################################
################################################################################
####                                                                        ####
#### Create a Xdmf file from iPic3D adn ECsim result in HDF5 format         ####
#### (written with the "default" option, no with H5hut)                     ####
####                                                                        ####
#### Diego Gonzalez                                                         ####
#### November 2017                                                          ####
####                                                                        ####
#### ---------------------------------------------------------------------- ####
####                                                                        ####
#### To generate the Xdmf file just execute this script in the same folder  ####
#### where the results are (file "seetings.hdf" is required)                ####
#### It will generate a file called "Results.xmf". This file can be opened  ####
#### directly with paraview.                                                ####
####                                                                        ####
#### ---------------------------------------------------------------------- ####
####                                                                        ####
#### WARNING: Xdmf assumes that the data in HDF5 files is saved with        ####
#### slowest varying dimension first (KJI). In the current version of       ####
#### iPic3D this is wrong.                                                  ####
#### In consequenece in the visualization the Z axis means X axis in the    ####
#### simulation.                                                            ####
####                                                                        ####
################################################################################
################################################################################
import numpy as np
import h5py

def exploreTree(listKeys, group, name="", s=""):
	for g in group:
		if (g.find("cycle") != -1):
			D = [group.name, name+s]
			listKeys.append(D)
			return True
		else:
			if (g.find("species_") != -1):
				txt, txtS = g.split("_")
				s = "_"+txtS
			subgroup = group.get(g)
			exploreTree(listKeys, subgroup,g,s)

################################################################################
####                              MAIN                                      ####
################################################################################
# Read the general settings
filename = "settings.hdf"
fileIn  = h5py.File(filename)
dx = fileIn.get("collective/Dx")[0]
dy = fileIn.get("collective/Dy")[0]
dz = fileIn.get("collective/Dz")[0]
Nxc = fileIn.get("collective/Nxc")[0]
Nyc = fileIn.get("collective/Nyc")[0]
Nzc = fileIn.get("collective/Nzc")[0]
Ns  = fileIn.get("collective/Ns")[0]
nprocs = fileIn.get("topology/Nprocs")[0]
XLEN = fileIn.get("topology/XLEN")[0]
YLEN = fileIn.get("topology/YLEN")[0]
ZLEN = fileIn.get("topology/ZLEN")[0]
dt = fileIn.get("collective/Dt")[0]
nxc = Nxc/XLEN
nyc = Nyc/YLEN
nzc = Nzc/ZLEN
nxn = nxc+1
nyn = nyc+1
nzn = nzc+1

# Check the last cycle and the distance between datasets
fileIn   = h5py.File("proc0.hdf")
groupBx  = fileIn.get("fields/Bx")
nsets = len(groupBx)
kk = groupBx.get("cycle_0")
cycle_max = 0
for g in groupBx:
	txtS, txtC = g.split('_')
	cycle = int(txtC);
	if (cycle > cycle_max):
		cycle_max = cycle

cycle_del = cycle_max
cycle_min = cycle_max
for g in groupBx:
	txtS, txtC = g.split('_')
	cycle = int(txtC);
	if (cycle != cycle_max and cycle_max - cycle < cycle_del):
		cycle_del = cycle_max-cycle
	if (cycle < cycle_min):
		cycle_min = cycle

# Explore the whole tree to find all the fields:
listKeys = []
group  = fileIn.get("fields")
exploreTree(listKeys, group)
group  = fileIn.get("moments")
exploreTree(listKeys, group)


# WRITE THE XDMF FILE	
fileOut = open("Results.xmf","w")

fileOut.write("<?xml version=\"1.0\" ?>\n")
fileOut.write("<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\">\n")
fileOut.write("<Xdmf Version=\"2.0\">\n")

fileOut.write("  <Domain>\n")
fileOut.write("    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n")
fileOut.write("      <Time TimeType=\"List\">\n")
fileOut.write("        <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"1\">\n")
fileOut.write("        ")
for c in range(nsets):
	cycle = cycle_min + c*cycle_del
	fileOut.write(" " +  str(dt*cycle))
fileOut.write("\n        </DataItem>\n")
fileOut.write("      </Time>\n")
for c in range(nsets):
	cycle = cycle_min + c*cycle_del
	fileOut.write("      <Grid Name=\"" + str(c) + "\" GridType=\"Collection\" CollectionType=\"Spatial\">\n")

	for p in range(nprocs):
		localfile = "proc"+str(p)+".hdf"
		fileIn   = h5py.File(localfile)
		cartesian_coord = fileIn.get("topology/cartesian_coord")
		x0 = cartesian_coord[0]*nxc*dx
		y0 = cartesian_coord[1]*nyc*dy
		z0 = cartesian_coord[2]*nzc*dz
		fileOut.write("        <Grid Name=\"Structured mesh\" GridType=\"Uniform\"> \n")
		fileOut.write("          <Topology TopologyType=\"3DRectMesh\" Dimensions=\"" + str(nxn) + " " + str(nyn) + " " + str(nzn) + "\"/> \n")
		fileOut.write("          <Geometry GeometryType=\"ORIGIN_DXDYDZ\"> \n")
		fileOut.write("            <DataItem Format=\"XML\" Dimensions=\"3\" NumberType=\"Float\"> \n")
		fileOut.write("             " + str(z0) + " " + str(y0) + " " + str(x0) + "\n")
		fileOut.write("            </DataItem> \n")
		fileOut.write("            <DataItem Format=\"XML\" Dimensions=\"3\" NumberType=\"Float\"> \n")
		fileOut.write("             " + str(dz) + " " + str(dy) + " " + str(dx) + "\n")
		fileOut.write("            </DataItem> \n")
		fileOut.write("          </Geometry> \n")

		for k in listKeys:
			fileOut.write("          <Attribute Name=\"" + k[1] + "\" Center=\"Node\">\n") 
			fileOut.write("          <DataItem Format=\"HDF\" Dimensions=\"" + str(nxn) + " " + str(nyn) + " " + str(nzn) + "\" NumberType=\"Float\"> \n")
			fileOut.write("           " + localfile + ":"+ k[0] +"/cycle_"+str(cycle)+"\n")
			fileOut.write("          </DataItem> \n")
			fileOut.write("          </Attribute> \n")

		fileOut.write("        </Grid>\n")
	fileOut.write("      </Grid>\n") # Spatial series
fileOut.write("    </Grid>\n")  # Time series
fileOut.write("  </Domain>\n")
fileOut.write("</Xdmf>\n")

fileOut.close()
