""" Create multiple vtr files that needs to be reassemble by a vtr file
    My spin-off of Arnaud's read_proc_para.py routine.

    This one reads densities and currents separately for different species.
    Also, full pressure tensor will be implemented.

    NOTE: if you, e.g., work with currents, where the sum: Jx = Jx0+Jx1.
    This sum would be different when computed here comparing to what you'll get in ParaView!
    So the result could be different from what you'd get from read_proc_para.py!

    USAGE:
    See the code below, you have to specify a few parameters.
    Modify, run (change the number of cores): 
    mpirun -np 4 python read_proc_spec.py

    (c) V. Olshevsky: sya@mao.kiev.ua

"""
from mpi4py import MPI
import xml.dom.minidom
import scipy
import tables
import vtk
import sys

comm = MPI.COMM_WORLD
rankproc = comm.Get_rank() # Proc rank
numproc = comm.Get_size() # Total number of procs
######### INPUT THE FOLLOWING PARAMETERS ############################################
first_cycle = 0
last_cycle = 1001
step = 100
directory = "../runs/gem0/data/" #Directory of the simulation data
output_pvd_fields_name = directory+"fields.pvd"
output_pvd_rho_name = directory+"rho.pvd"
output_pvd_currents_name = directory+"current.pvd"
#shape = (9,9,9) # Shape of one field array for 1 proc*.hdf file. TEMPORARY. Need to read the shape of each proc domain.
shape = (33, 33, 2) # for GEM
ns = 2 # number of species, no more than 6 (due to the declaration of Fields_local)
####################################################################################

## Some necessary constants
str_moments = 'moments'
str_fields = 'fields'
str_species = 'species_'
str_rho = 'rho'
current_comp_names = ['Jx', 'Jy', 'Jz']
field_comp_names = ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez']
field_comp_count = len(field_comp_names)
field_comp = xrange(field_comp_count)
species = xrange(ns)
list_cycle = scipy.arange(first_cycle,last_cycle+1,step).tolist()

## Gathering information and initialization
setting_filename=directory + "settings.hdf" # Setting file name
h5setting_file = tables.openFile(setting_filename, mode = "r", title = "Setting_file")
collective_group = h5setting_file.root.collective
topology_group = h5setting_file.root.topology
Nprocs = topology_group._f_getChild("Nprocs").read()[0]
dx = collective_group._f_getChild("Dx").read()[0]
dy = collective_group._f_getChild("Dy").read()[0]
dz = collective_group._f_getChild("Dz").read()[0]
Nxc = collective_group._f_getChild("Nxc").read()[0]
Nyc = collective_group._f_getChild("Nyc").read()[0]
Nzc = collective_group._f_getChild("Nzc").read()[0]
XLEN = topology_group._f_getChild("XLEN").read()[0]
YLEN = topology_group._f_getChild("YLEN").read()[0]
ZLEN = topology_group._f_getChild("ZLEN").read()[0]
h5setting_file.close()

if rankproc == 0:
    print "Reading directory " + directory

Zbyproc = ZLEN/numproc #Minimum Nomber of Z you have to read
reste = ZLEN - Zbyproc*numproc
if rankproc < reste :
    Zrange = scipy.arange(rankproc*(Zbyproc+1),(rankproc+1)*(Zbyproc+1))
else:
    Zrange = scipy.arange((reste*(Zbyproc+1)+(rankproc-reste)*Zbyproc),(reste*(Zbyproc+1)+(rankproc-reste+1)*Zbyproc))
list_read_proc = []
for Z in Zrange:
    list_read_proc.extend(range(Z,Nprocs,ZLEN))

list_read_proc = scipy.array(list_read_proc)
zread_proc = len(list_read_proc)/(YLEN*XLEN) #Number of proc file you want to read in the z direction

Nzlocal = zread_proc*(shape[2] - 1) + 1
Fields_local = scipy.zeros((field_comp_count,Nxc+1,Nyc+1,Nzlocal), dtype="float32")
numproc = min(numproc,ZLEN) # To account for specific case where more processors are used here than in the Z direction of the simulation
if len(list_read_proc) > 0:
    numpoints = (Nxc+1)*(Nyc+1)*Nzlocal
    Xcoordinates = scipy.arange(Nxc+1).flatten().astype('float32')*dx
    Ycoordinates = scipy.arange(Nyc+1).flatten().astype('float32')*dy
    Zcoordinates = (scipy.arange(Nzlocal).flatten().astype('float32') + list_read_proc[0]*(shape[2]-1))*dz
    vtkXcoordinates = vtk.vtkFloatArray()
    vtkXcoordinates.SetNumberOfComponents(1)
    vtkXcoordinates.SetNumberOfTuples(Nxc+1)
    vtkXcoordinates.SetVoidArray(Xcoordinates,Nxc+1,1)
    vtkYcoordinates = vtk.vtkFloatArray()
    vtkYcoordinates.SetNumberOfComponents(1)
    vtkYcoordinates.SetNumberOfTuples(Nyc+1)
    vtkYcoordinates.SetVoidArray(Ycoordinates,Nyc+1,1)
    vtkZcoordinates = vtk.vtkFloatArray()
    vtkZcoordinates.SetNumberOfComponents(1)
    vtkZcoordinates.SetNumberOfTuples(Nzlocal)
    vtkZcoordinates.SetVoidArray(Zcoordinates,Nzlocal,1)
    ###########################################################################################
    ##                       FIELDS                                                         ##
    ###########################################################################################
    for cycle in list_cycle:    
        for read_proc in list_read_proc:
            proc_filename = directory + "proc"+repr(read_proc)+".hdf" # Proc file name
            
            h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")
            list_group = h5proc_file.root.fields._f_listNodes("Group")
            
            p_coordinate=h5proc_file.root.topology.cartesian_coord.read() 
            p_coordinate[2]=p_coordinate[2]-list_read_proc[0] 
            x_start=p_coordinate[0]*(shape[0]-1) #-1 to account for the 1 node overlapping
            x_end=x_start+shape[0]
            y_start=p_coordinate[1]*(shape[1]-1)
            y_end=y_start+shape[1]
            z_start=p_coordinate[2]*(shape[2]-1)
            z_end=z_start+shape[2]
            ## Cycle over field components
            for c in field_comp :                                                                                                    
                for node in h5proc_file.walkNodes('/' + str_fields + '/'+ field_comp_names[c]) :
                    if type(node) == tables.array.Array :
                        try :
                            num = int(node.name[6:])
                        except :
                            print 'Node name conversion error for node: ', node
                        if num == cycle :
                            Fields_local[c, x_start:x_end, y_start:y_end, z_start:z_end] = node.read()               
                        #end if                                                                                                    
                    #end if                                                                                                        
                #end for                                                                                                           
            #end for
            h5proc_file.close()
        #end for
        #Write the .vtr file
        rtg = vtk.vtkRectilinearGrid()        
        rtg.SetDimensions(Nxc+1,Nyc+1,Nzlocal)
        rtg.SetExtent(0,Nxc,0,Nyc,rankproc*Nzc/numproc,(rankproc+1)*Nzc/numproc)
        rtg.SetXCoordinates(vtkXcoordinates)
        rtg.SetYCoordinates(vtkYcoordinates)
        rtg.SetZCoordinates(vtkZcoordinates)
        ## Lists of arrays for writing
        vtk_data = [] 
        array_list = []
        for c in field_comp :
            vtk_data.append(vtk.vtkFloatArray())
            vtk_data[c].SetNumberOfTuples(numpoints)                              
            vtk_data[c].SetNumberOfComponents(1)
            array_list.append(Fields_local[c,:,:,:].swapaxes(0,2).flatten().astype("float32"))
            vtk_data[c].SetVoidArray(array_list[c], numpoints, 1)    
            vtk_data[c].SetName(field_comp_names[c])
            if c == 0 :
                rtg.GetPointData().SetScalars(vtk_data[c])
            else :
                rtg.GetPointData().AddArray(vtk_data[c])
        #end for       
        writer =vtk.vtkXMLRectilinearGridWriter()
        writer.SetFileName(directory+"fields_cycle"+repr(cycle).rjust(8,"0")+"-"+repr(rankproc)+".vtr")
        writer.SetInput(rtg)
        writer.Write()
        ## Dispose the lists; garbage collector SHOULD also clear memory from arrays
        array_list = None
        vtk_data = None
        ## Write the pvtr file
        if rankproc == 0:
            output_pvtr_name=directory+"fields_cycle"+repr(cycle).rjust(8,"0")+".pvtr"
            liste_proc=range(numproc)
            output_pvtr_id = open(output_pvtr_name, 'w')
            
            pvtr = xml.dom.minidom.Document()
            pvtr_root = pvtr.createElementNS("VTK", "VTKFile")
            pvtr_root.setAttribute("type", "PRectilinearGrid")
            pvtr_root.setAttribute("version", "0.1")
            pvtr_root.setAttribute("byte_order", "LittleEndian")
            pvtr.appendChild(pvtr_root)
            
            grid = pvtr.createElementNS("VTK", "PRectilinearGrid")
            grid.setAttribute("WholeExtent", "0 "+ repr(Nxc)+" 0 "+repr(Nyc)+" 0 "+repr(Nzc))
            grid.setAttribute("GhostLevel","0")
            pvtr_root.appendChild(grid)
            
            coordinate = pvtr.createElementNS("VTK", "PCoordinates")
            for k in range(3):
                data_array = pvtr.createElementNS("VTK", "PDataArray")
                data_array.setAttribute("type","Float32")
                data_array.setAttribute("format","Appended")
                coordinate.appendChild(data_array)
            #end for
            grid.appendChild(coordinate)
            
            p_data = pvtr.createElementNS("VTK","PPointData")
            p_data.setAttribute("Scalars", field_comp_names[0])
            for c in field_comp :
                data_array = pvtr.createElementNS("VTK", "PDataArray")
                data_array.setAttribute("type","Float32")
                data_array.setAttribute("Name", field_comp_names[c])
                p_data.appendChild(data_array)
            #end for
            grid.appendChild(p_data)
            
            for i in liste_proc:
                num = repr(i) 
                filename = "fields_cycle"+repr(cycle).rjust(8,"0")+"-"+repr(i)+".vtr"
                piece = pvtr.createElementNS("VTK", "Piece")
                piece.setAttribute("Extent", "0 "+ repr(Nxc)+" 0 "+repr(Nyc)+" "+repr(i*Nzc/numproc)+" "+repr((i+1)*Nzc/numproc))
                piece.setAttribute("Source",filename)
                grid.appendChild(piece)
            #end for
            pvtr.writexml(output_pvtr_id, newl='\n')
            output_pvtr_id.close() 
        #end if
    #end for
    ## Write the pvd file
    if rankproc == 0:
        output_pvd_id = open(output_pvd_fields_name, 'w')
        
        pvd = xml.dom.minidom.Document()
        pvd_root = pvd.createElementNS("VTK", "VTKFile")
        pvd_root.setAttribute("type", "Collection")
        pvd_root.setAttribute("version", "0.1")
        pvd_root.setAttribute("byte_order", "LittleEndian")
        pvd.appendChild(pvd_root)
        
        collection = pvd.createElementNS("VTK", "Collection")
        pvd_root.appendChild(collection)
        
        for i in list_cycle:
                num=repr(i).rjust(8,"0") #Must match the rjust from read_proc
                filename="fields_cycle"+num+".pvtr"
                dataSet = pvd.createElementNS("VTK", "DataSet")
                dataSet.setAttribute("timestep", str(i))
                dataSet.setAttribute("group", "")
                dataSet.setAttribute("part", "0")
                dataSet.setAttribute("file", filename)
                collection.appendChild(dataSet)
        #end for
        pvd.writexml(output_pvd_id, newl='\n')
        output_pvd_id.close()
        print 'Fields written'
    #end if    
    ###########################################################################################
    ##                       RHO                                                             ##
    ###########################################################################################
    for cycle in list_cycle:    
        Fields_local = scipy.zeros((ns+1,Nxc+1,Nyc+1,Nzlocal), dtype="float32")
        ## Cycle over processors; for each of them read data for the current cycle
        for read_proc in list_read_proc:
            proc_filename = directory + "proc"+repr(read_proc)+".hdf" # Proc file name            
            h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")
            #list_group = h5proc_file.root.moments._f_listNodes("Group")
            p_coordinate = h5proc_file.root.topology.cartesian_coord.read() 
            p_coordinate[2] = p_coordinate[2] - list_read_proc[0] 
            x_start = p_coordinate[0]*(shape[0]-1) #-1 to account for the 1 node overlapping
            x_end = x_start+shape[0]
            y_start = p_coordinate[1]*(shape[1]-1)
            y_end = y_start+shape[1]
            z_start = p_coordinate[2]*(shape[2]-1)
            z_end = z_start+shape[2]
            # Read the total density first
            for node in h5proc_file.walkNodes('/' + str_moments + '/' + str_rho) :
                if type(node) == tables.array.Array :
                    try : 
                        num = int(node.name[6:])
                    except :
                        print 'Node name conversion error for node: '
                        print node
                    #end try
                    if num == cycle :                    
                        Fields_local[0, x_start:x_end, y_start:y_end, z_start:z_end] = node.read()
                    #end if
                #end if
            #end for                
            ## Cycle over species for each particle type density
            for specie in species :
                for node in h5proc_file.walkNodes('/' + str_moments + '/'+ str_species + str(specie) + '/' + str_rho) :
                    if type(node) == tables.array.Array :
                        try : 
                            num = int(node.name[6:])
                        except :
                            print 'Node name conversion error for node: ', node
                        if num == cycle :
                            Fields_local[specie + 1, x_start:x_end, y_start:y_end, z_start:z_end] = node.read()
                        #end if
                    #end if
                #end for
            #end for            
            h5proc_file.close()
        #end for
        ## Write the .vtr file
        rtg = vtk.vtkRectilinearGrid()        
        rtg.SetDimensions(Nxc+1,Nyc+1,Nzlocal)
        rtg.SetExtent(0,Nxc,0,Nyc,rankproc*Nzc/numproc,(rankproc+1)*Nzc/numproc)        
        rtg.SetXCoordinates(vtkXcoordinates)
        rtg.SetYCoordinates(vtkYcoordinates)
        rtg.SetZCoordinates(vtkZcoordinates)
        ## Lists of arrays for writing
        vtk_data = [] 
        array_list = []
        vtk_data.append(vtk.vtkFloatArray())
        vtk_data[0].SetNumberOfTuples(numpoints)                              
        vtk_data[0].SetNumberOfComponents(1)
        array_list.append(Fields_local[0,:,:,:].swapaxes(0,2).flatten().astype("float32"))
        vtk_data[0].SetVoidArray(array_list[0], numpoints, 1)    
        vtk_data[0].SetName('Rho')
        rtg.GetPointData().SetScalars(vtk_data[0])        
        ## Cycle over species, and add the components of rho to the vtk grid
        for specie in species :
            vtk_data.append(vtk.vtkFloatArray())
            vtk_data[specie+1].SetNumberOfTuples(numpoints)
            vtk_data[specie+1].SetNumberOfComponents(1)
            array_list.append(Fields_local[specie + 1,:,:,:].swapaxes(0,2).flatten().astype("float32"))
            vtk_data[specie+1].SetVoidArray(array_list[specie + 1], numpoints, 1)	
            vtk_data[specie+1].SetName('Rho' + repr(specie))
            rtg.GetPointData().AddArray(vtk_data[specie+1])
        #end for
        writer = vtk.vtkXMLRectilinearGridWriter()
        writer.SetFileName(directory+"rho_cycle"+repr(cycle).rjust(8,"0")+"-"+repr(rankproc)+".vtr")
        writer.SetInput(rtg)
        writer.Write()
        ## Dispose the lists; garbage collector SHOULD also clear memory from arrays
        array_list = None
        vtk_data = None
        ## Write the pvtr file
        if rankproc == 0:
            output_pvtr_name = directory+"rho_cycle"+repr(cycle).rjust(8,"0")+".pvtr"
            liste_proc = range(numproc)
            output_pvtr_id = open(output_pvtr_name, 'w')
            
            pvtr = xml.dom.minidom.Document()
            pvtr_root = pvtr.createElementNS("VTK", "VTKFile")
            pvtr_root.setAttribute("type", "PRectilinearGrid")
            pvtr_root.setAttribute("version", "0.1")
            pvtr_root.setAttribute("byte_order", "LittleEndian")
            pvtr.appendChild(pvtr_root)
            
            grid = pvtr.createElementNS("VTK", "PRectilinearGrid")
            grid.setAttribute("WholeExtent", "0 "+ repr(Nxc)+" 0 "+repr(Nyc)+" 0 "+repr(Nzc))
            grid.setAttribute("GhostLevel","0")
            pvtr_root.appendChild(grid)
            
            coordinate = pvtr.createElementNS("VTK", "PCoordinates")
            for k in range(3):
                data_array = pvtr.createElementNS("VTK", "PDataArray")
                data_array.setAttribute("type","Float32")
                data_array.setAttribute("format","Appended")
                coordinate.appendChild(data_array)
            #end for
            grid.appendChild(coordinate)
            
            p_data = pvtr.createElementNS("VTK","PPointData")
            p_data.setAttribute("Scalars","Rho")
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","Rho")
            p_data.appendChild(data_array)

            for specie in species :
                data_array = pvtr.createElementNS("VTK", "PDataArray")
                data_array.setAttribute("type","Float32")
                data_array.setAttribute("Name","Rho" + repr(specie))
                p_data.appendChild(data_array)
            #end for
            grid.appendChild(p_data)
            
            for i in liste_proc:
                num=repr(i) 
                filename="rho_cycle"+repr(cycle).rjust(8,"0")+"-"+repr(i)+".vtr"
                piece = pvtr.createElementNS("VTK", "Piece")
                piece.setAttribute("Extent", "0 "+ repr(Nxc)+" 0 "+repr(Nyc)+" "+repr(i*Nzc/numproc)+" "+repr((i+1)*Nzc/numproc))
                piece.setAttribute("Source",filename)
                grid.appendChild(piece)
            #end for
            pvtr.writexml(output_pvtr_id, newl='\n')
            output_pvtr_id.close() 
        #end if
    #end for
    ## Write the pvd file
    if rankproc == 0:
        output_pvd_id = open(output_pvd_rho_name, 'w')
        
        pvd = xml.dom.minidom.Document()
        pvd_root = pvd.createElementNS("VTK", "VTKFile")
        pvd_root.setAttribute("type", "Collection")
        pvd_root.setAttribute("version", "0.1")
        pvd_root.setAttribute("byte_order", "LittleEndian")
        pvd.appendChild(pvd_root)
        
        collection = pvd.createElementNS("VTK", "Collection")
        pvd_root.appendChild(collection)
        
        for i in list_cycle:
                num = repr(i).rjust(8,"0") #Must match the rjust from read_proc
                filename = "rho_cycle" + num + ".pvtr"
                dataSet = pvd.createElementNS("VTK", "DataSet")
                dataSet.setAttribute("timestep", str(i))
                dataSet.setAttribute("group", "")
                dataSet.setAttribute("part", "0")
                dataSet.setAttribute("file", filename)
                collection.appendChild(dataSet)
        #end for
        pvd.writexml(output_pvd_id, newl='\n')
        output_pvd_id.close()
        print 'Rho written'
    #end if
    ##########################################################################################
    ##                       CURRENTS                                                        ##
    ###########################################################################################    
    for cycle in list_cycle:    
        Fields_local = scipy.zeros((ns*3,Nxc+1,Nyc+1,Nzlocal), dtype="float32")
        for read_proc in list_read_proc:
            proc_filename = directory + "proc"+repr(read_proc)+".hdf" # Proc file name            
            h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")            
            p_coordinate = h5proc_file.root.topology.cartesian_coord.read() 
            p_coordinate[2] = p_coordinate[2]-list_read_proc[0] 
            x_start = p_coordinate[0]*(shape[0]-1) #-1 to account for the 1 node overlapping
            x_end = x_start+shape[0]
            y_start = p_coordinate[1]*(shape[1]-1)
            y_end = y_start+shape[1]
            z_start = p_coordinate[2]*(shape[2]-1)
            z_end = z_start+shape[2]
            ## Cycle over species to get current for each specie
            for specie in species :
                for c in xrange(3) :
                    for node in h5proc_file.walkNodes('/' + str_moments + '/'+ str_species + repr(specie) + '/' + current_comp_names[c]) :
                        if type(node) == tables.array.Array :
                            try : 
                                num = int(node.name[6:])
                            except :
                                print 'Node name conversion error for node: ', node
                            if num == cycle :
                                #Fields_local[c, x_start:x_end, y_start:y_end, z_start:z_end] += node.read()
                                Fields_local[specie*3 + c, x_start:x_end, y_start:y_end, z_start:z_end] = node.read()                                
                            #end if
                        #end if
                    #end for
                #end for
            #end for            
            h5proc_file.close()
        #end for
        ## Write the .vtr file        
        rtg = vtk.vtkRectilinearGrid()        
        rtg.SetDimensions(Nxc+1,Nyc+1,Nzlocal)
        rtg.SetExtent(0,Nxc,0,Nyc,rankproc*Nzc/numproc,(rankproc+1)*Nzc/numproc)        
        rtg.SetXCoordinates(vtkXcoordinates)        
        rtg.SetYCoordinates(vtkYcoordinates)        
        rtg.SetZCoordinates(vtkZcoordinates)
        ## Lists of arrays for writing
        vtk_data = [] 
        array_list = []
        for specie in species :
            for c in xrange(3) :
                i = 3*specie+c
                vtk_data.append(vtk.vtkFloatArray())
                vtk_data[i].SetNumberOfTuples(numpoints)                              
                vtk_data[i].SetNumberOfComponents(1)
                array_list.append(Fields_local[i,:,:,:].swapaxes(0,2).flatten().astype("float32"))
                vtk_data[i].SetVoidArray(array_list[i], numpoints, 1)    
                vtk_data[i].SetName(current_comp_names[c] + repr(specie))
                if i == 0 :
                    rtg.GetPointData().SetScalars(vtk_data[i])
                else :
                    rtg.GetPointData().AddArray(vtk_data[i])
            #end for
        #end for       
        writer =vtk.vtkXMLRectilinearGridWriter()
        writer.SetFileName(directory+"currents_cycle"+repr(cycle).rjust(8,"0")+"-"+repr(rankproc)+".vtr")
        writer.SetInput(rtg)
        writer.Write()
        ## Dispose the lists; garbage collector SHOULD also clear memory from arrays
        vtk_data = None
        array_list = None
        ## Write the pvtr file
        if rankproc == 0:
            output_pvtr_name=directory+"currents_cycle"+repr(cycle).rjust(8,"0")+".pvtr"
            liste_proc=range(numproc)
            output_pvtr_id = open(output_pvtr_name, 'w')
            
            pvtr = xml.dom.minidom.Document()
            pvtr_root = pvtr.createElementNS("VTK", "VTKFile")
            pvtr_root.setAttribute("type", "PRectilinearGrid")
            pvtr_root.setAttribute("version", "0.1")
            pvtr_root.setAttribute("byte_order", "LittleEndian")
            pvtr.appendChild(pvtr_root)
            
            grid = pvtr.createElementNS("VTK", "PRectilinearGrid")
            grid.setAttribute("WholeExtent", "0 "+ repr(Nxc)+" 0 "+repr(Nyc)+" 0 "+repr(Nzc))
            grid.setAttribute("GhostLevel","0")
            pvtr_root.appendChild(grid)
            
            coordinate = pvtr.createElementNS("VTK", "PCoordinates")
            for k in range(3):
                data_array = pvtr.createElementNS("VTK", "PDataArray")
                data_array.setAttribute("type","Float32")
                data_array.setAttribute("format","Appended")
                coordinate.appendChild(data_array)
            #end for
            grid.appendChild(coordinate)
            
            p_data = pvtr.createElementNS("VTK","PPointData")
            p_data.setAttribute("Scalars", current_comp_names[0] + "0")

            for specie in species :
                for c in xrange(3) :
                    data_array = pvtr.createElementNS("VTK", "PDataArray")
                    data_array.setAttribute("type","Float32")
                    data_array.setAttribute("Name", current_comp_names[c] + repr(specie))
                    p_data.appendChild(data_array)
                #end for
            #end for    
            grid.appendChild(p_data)
            
            for i in liste_proc:
                num = repr(i) 
                filename = "currents_cycle"+repr(cycle).rjust(8,"0")+"-"+repr(i)+".vtr"
                piece = pvtr.createElementNS("VTK", "Piece")
                piece.setAttribute("Extent", "0 "+ repr(Nxc)+" 0 "+repr(Nyc)+" "+repr(i*Nzc/numproc)+" "+repr((i+1)*Nzc/numproc))
                piece.setAttribute("Source",filename)
                grid.appendChild(piece)
            #end for
            
            pvtr.writexml(output_pvtr_id, newl='\n')
            output_pvtr_id.close() 
        #end if
    #end for
    ## Write the pvd file
    if rankproc == 0:
        output_pvd_id = open(output_pvd_currents_name, 'w')
        
        pvd = xml.dom.minidom.Document()
        pvd_root = pvd.createElementNS("VTK", "VTKFile")
        pvd_root.setAttribute("type", "Collection")
        pvd_root.setAttribute("version", "0.1")
        pvd_root.setAttribute("byte_order", "LittleEndian")
        pvd.appendChild(pvd_root)
        
        collection = pvd.createElementNS("VTK", "Collection")
        pvd_root.appendChild(collection)
        
        for i in list_cycle:
                num=repr(i).rjust(8,"0") #Must match the rjust from read_proc
                filename="currents_cycle"+num+".pvtr"
                dataSet = pvd.createElementNS("VTK", "DataSet")
                dataSet.setAttribute("timestep", str(i))
                dataSet.setAttribute("group", "")
                dataSet.setAttribute("part", "0")
                dataSet.setAttribute("file", filename)
                collection.appendChild(dataSet)
        #end for
        pvd.writexml(output_pvd_id, newl='\n')
        output_pvd_id.close()
        print 'Currents written'
    #end if
#end if