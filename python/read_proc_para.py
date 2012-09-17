#Create multiple vtr files that needs to be reassemble by a vtr file

from mpi4py import MPI
import xml.dom.minidom
import scipy
import tables
import vtk
import sys
import time

comm = MPI.COMM_WORLD
rankproc = comm.Get_rank() # Proc rank
numproc = comm.Get_size() # Total number of procs

if rankproc == 0:
    start_time = time.clock()
######### INPUT THE FOLLOWING PARAMETERS ############################################
first_cycle=0
last_cycle=2001
step=100
directory = "../runs/gem0/data/" #Directory of the simulation data
output_pvd_fields_name = directory+"fields.pvd"
output_pvd_rho_name = directory+"rho.pvd"
output_pvd_currents_name = directory+"current.pvd"
#shape=(9,9,9) # Shape of one field array for 1 proc*.hdf file. TEMPORARY. Need to read the shape of each proc domain.
shape = (33, 33, 2) # for GEM
####################################################################################
list_cycle=scipy.arange(first_cycle,last_cycle+1,step).tolist()

###Gathering in formation and initialization
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

Zbyproc=ZLEN/numproc #Minimum Nomber of Z you have to read
reste=ZLEN-Zbyproc*numproc
if (rankproc<reste):
    Zrange=scipy.arange(rankproc*(Zbyproc+1),(rankproc+1)*(Zbyproc+1))
else:
    Zrange=scipy.arange((reste*(Zbyproc+1)+(rankproc-reste)*Zbyproc),(reste*(Zbyproc+1)+(rankproc-reste+1)*Zbyproc))
list_read_proc=[]
for Z in Zrange:
    list_read_proc.extend(range(Z,Nprocs,ZLEN))

list_read_proc=scipy.array(list_read_proc)
zread_proc=len(list_read_proc)/(YLEN*XLEN) #Number of proc file you want to read in the z direction


Nzlocal=zread_proc*(shape[2]-1)+1
Fields_local = scipy.zeros((6,Nxc+1,Nyc+1,Nzlocal),dtype="float32")
# 6 for 6 components of the field
numproc=min(numproc,ZLEN) # To account for specific case where more processors are used here than in the Z direction of the simulation
if len(list_read_proc) > 0:
    numpoints = (Nxc+1)*(Nyc+1)*Nzlocal
    Xcoordinates = scipy.arange(Nxc+1).flatten().astype('float32')*dx
    Ycoordinates = scipy.arange(Nyc+1).flatten().astype('float32')*dy
    Zcoordinates = (scipy.arange(Nzlocal).flatten().astype('float32')+list_read_proc[0]*(shape[2]-1))*dz
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
        ### Start reading #######
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
            k=0 #Counts the 6 components
            for group in list_group:
                list_array=group._f_listNodes("Array") #List of all recorded arrays in "group"
                for array in list_array:
                    num=scipy.int64(repr(array).split()[0][17:])
                    if num==cycle:
                        num_final=num
                        Fields_local[k,x_start:x_end,y_start:y_end,z_start:z_end]=array.read()
                k=k+1
            h5proc_file.close()
    #Write the .vtr file
        rtg = vtk.vtkRectilinearGrid()
        
        rtg.SetDimensions(Nxc+1,Nyc+1,Nzlocal)
        rtg.SetExtent(0,Nxc,0,Nyc,rankproc*Nzc/numproc,(rankproc+1)*Nzc/numproc)
        rtg.SetXCoordinates(vtkXcoordinates)
        rtg.SetYCoordinates(vtkYcoordinates)
        rtg.SetZCoordinates(vtkZcoordinates)
        
        Bx=Fields_local[0,:,:,:].swapaxes(0,2).flatten().astype("float32")
        By=Fields_local[1,:,:,:].swapaxes(0,2).flatten().astype("float32")
        Bz=Fields_local[2,:,:,:].swapaxes(0,2).flatten().astype("float32")
        Ex=Fields_local[3,:,:,:].swapaxes(0,2).flatten().astype("float32")
        Ey=Fields_local[4,:,:,:].swapaxes(0,2).flatten().astype("float32")
        Ez=Fields_local[5,:,:,:].swapaxes(0,2).flatten().astype("float32")
        
        vtk_dataBx = vtk.vtkFloatArray()
        vtk_dataBx.SetNumberOfTuples(numpoints)
        vtk_dataBx.SetNumberOfComponents(1)
        vtk_dataBx.SetVoidArray(Bx, numpoints, 1)	
        vtk_dataBx.SetName('Bx')	 
        vtk_dataBy = vtk.vtkFloatArray()
        vtk_dataBy.SetNumberOfTuples(numpoints)
        vtk_dataBy.SetNumberOfComponents(1)
        vtk_dataBy.SetVoidArray(By, numpoints, 1)	
        vtk_dataBy.SetName('By')	 
        vtk_dataBz = vtk.vtkFloatArray()
        vtk_dataBz.SetNumberOfTuples(numpoints)
        vtk_dataBz.SetNumberOfComponents(1)
        vtk_dataBz.SetVoidArray(Bz, numpoints, 1)	
        vtk_dataBz.SetName('Bz')	 
        vtk_dataEx = vtk.vtkFloatArray()
        vtk_dataEx.SetNumberOfTuples(numpoints)
        vtk_dataEx.SetNumberOfComponents(1)
        vtk_dataEx.SetVoidArray(Ex, numpoints, 1)	
        vtk_dataEx.SetName('Ex')	 
        vtk_dataEy = vtk.vtkFloatArray()
        vtk_dataEy.SetNumberOfTuples(numpoints)
        vtk_dataEy.SetNumberOfComponents(1)
        vtk_dataEy.SetVoidArray(Ey, numpoints, 1)	
        vtk_dataEy.SetName('Ey')	 
        vtk_dataEz = vtk.vtkFloatArray()
        vtk_dataEz.SetNumberOfTuples(numpoints)
        vtk_dataEz.SetNumberOfComponents(1)
        vtk_dataEz.SetVoidArray(Ez, numpoints, 1)	
        vtk_dataEz.SetName('Ez')	 
        
        
        rtg.GetPointData().SetScalars(vtk_dataBx)
        rtg.GetPointData().AddArray(vtk_dataBy)
        rtg.GetPointData().AddArray(vtk_dataBz)
        rtg.GetPointData().AddArray(vtk_dataEx)
        rtg.GetPointData().AddArray(vtk_dataEy)
        rtg.GetPointData().AddArray(vtk_dataEz)
        
        writer =vtk.vtkXMLRectilinearGridWriter()
        writer.SetFileName(directory+"fields_cycle"+repr(num_final).rjust(8,"0")+"-"+repr(rankproc)+".vtr")
        writer.SetInput(rtg)
        writer.Write()
    ## Write the pvtr file
        if rankproc == 0:
            output_pvtr_name=directory+"fields_cycle"+repr(num_final).rjust(8,"0")+".pvtr"
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
            grid.appendChild(coordinate)
            
            p_data = pvtr.createElementNS("VTK","PPointData")
            p_data.setAttribute("Scalars","Bx")
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","Bx")
            p_data.appendChild(data_array)
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","By")
            p_data.appendChild(data_array)
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","Bz")
            p_data.appendChild(data_array)
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","Ex")
            p_data.appendChild(data_array)
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","Ey")
            p_data.appendChild(data_array)
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","Ez")
            p_data.appendChild(data_array)
    
            grid.appendChild(p_data)
            
            for i in liste_proc:
                num=repr(i) 
                filename="fields_cycle"+repr(num_final).rjust(8,"0")+"-"+repr(i)+".vtr"
                piece = pvtr.createElementNS("VTK", "Piece")
                piece.setAttribute("Extent", "0 "+ repr(Nxc)+" 0 "+repr(Nyc)+" "+repr(i*Nzc/numproc)+" "+repr((i+1)*Nzc/numproc))
                piece.setAttribute("Source",filename)
                grid.appendChild(piece)
            
            pvtr.writexml(output_pvtr_id, newl='\n')
            output_pvtr_id.close() 
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
        
        pvd.writexml(output_pvd_id, newl='\n')
        output_pvd_id.close()
    
    ###########################################################################################
    ##                       RHO                                                             ##
    ###########################################################################################
    for cycle in list_cycle:    
        ### Start reading #######
        for read_proc in list_read_proc:
            proc_filename = directory + "proc"+repr(read_proc)+".hdf" # Proc file name
            
            h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")
            list_group = h5proc_file.root.moments._f_listNodes("Group")
            
            p_coordinate=h5proc_file.root.topology.cartesian_coord.read() 
            p_coordinate[2]=p_coordinate[2]-list_read_proc[0] 
            x_start=p_coordinate[0]*(shape[0]-1) #-1 to account for the 1 node overlapping
            x_end=x_start+shape[0]
            y_start=p_coordinate[1]*(shape[1]-1)
            y_end=y_start+shape[1]
            z_start=p_coordinate[2]*(shape[2]-1)
            z_end=z_start+shape[2]
            k=0 #Counts the 1 component
            for group in list_group[:1]:
                list_array=group._f_listNodes("Array") #List of all recorded arrays in "group"
                for array in list_array:
                    num=scipy.int64(repr(array).split()[0][19:])
                    if num==cycle:
                        num_final=num
                        Fields_local[k,x_start:x_end,y_start:y_end,z_start:z_end]=array.read()
                k=k+1
            h5proc_file.close()
    #Write the .vtr file
        rtg = vtk.vtkRectilinearGrid()
        
        rtg.SetDimensions(Nxc+1,Nyc+1,Nzlocal)
        rtg.SetExtent(0,Nxc,0,Nyc,rankproc*Nzc/numproc,(rankproc+1)*Nzc/numproc)
        
        rtg.SetXCoordinates(vtkXcoordinates)
        rtg.SetYCoordinates(vtkYcoordinates)
        rtg.SetZCoordinates(vtkZcoordinates)
        
        Bx=Fields_local[0,:,:,:].swapaxes(0,2).flatten().astype("float32")
        
        vtk_dataBx = vtk.vtkFloatArray()
        vtk_dataBx.SetNumberOfTuples(numpoints)
        vtk_dataBx.SetNumberOfComponents(1)
        vtk_dataBx.SetVoidArray(Bx, numpoints, 1)	
        vtk_dataBx.SetName('Rho')	 
            
        
        rtg.GetPointData().SetScalars(vtk_dataBx)
        
        writer =vtk.vtkXMLRectilinearGridWriter()
        writer.SetFileName(directory+"rho_cycle"+repr(num_final).rjust(8,"0")+"-"+repr(rankproc)+".vtr")
        writer.SetInput(rtg)
        writer.Write()
    ## Write the pvtr file
        if rankproc == 0:
            output_pvtr_name=directory+"rho_cycle"+repr(num_final).rjust(8,"0")+".pvtr"
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
            grid.appendChild(coordinate)
            
            p_data = pvtr.createElementNS("VTK","PPointData")
            p_data.setAttribute("Scalars","Rho")
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","Rho")
            p_data.appendChild(data_array)
    
            grid.appendChild(p_data)
            
            for i in liste_proc:
                num=repr(i) 
                filename="rho_cycle"+repr(num_final).rjust(8,"0")+"-"+repr(i)+".vtr"
                piece = pvtr.createElementNS("VTK", "Piece")
                piece.setAttribute("Extent", "0 "+ repr(Nxc)+" 0 "+repr(Nyc)+" "+repr(i*Nzc/numproc)+" "+repr((i+1)*Nzc/numproc))
                piece.setAttribute("Source",filename)
                grid.appendChild(piece)
            
            pvtr.writexml(output_pvtr_id, newl='\n')
            output_pvtr_id.close() 
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
                num=repr(i).rjust(8,"0") #Must match the rjust from read_proc
                filename="rho_cycle"+num+".pvtr"
                dataSet = pvd.createElementNS("VTK", "DataSet")
                dataSet.setAttribute("timestep", str(i))
                dataSet.setAttribute("group", "")
                dataSet.setAttribute("part", "0")
                dataSet.setAttribute("file", filename)
                collection.appendChild(dataSet)
        
        pvd.writexml(output_pvd_id, newl='\n')
        output_pvd_id.close()
    
    ##########################################################################################
    ##                       CURRENTS                                                        ##
    ###########################################################################################
    
    for cycle in list_cycle:    
        Fields_local = scipy.zeros((3,Nxc+1,Nyc+1,Nzlocal),dtype="float32")
        ### Start reading #######
        for read_proc in list_read_proc:
            proc_filename = directory + "proc"+repr(read_proc)+".hdf" # Proc file name
            
            h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")
            list_group = h5proc_file.root.moments._f_listNodes("Group")
            
            p_coordinate=h5proc_file.root.topology.cartesian_coord.read() 
            p_coordinate[2]=p_coordinate[2]-list_read_proc[0] 
            x_start=p_coordinate[0]*(shape[0]-1) #-1 to account for the 1 node overlapping
            x_end=x_start+shape[0]
            y_start=p_coordinate[1]*(shape[1]-1)
            y_end=y_start+shape[1]
            z_start=p_coordinate[2]*(shape[2]-1)
            z_end=z_start+shape[2]
            for speciesgroup in list_group[1:]: #Loop on species
                list_sub_group = speciesgroup._f_listNodes("Group")
                k=0 #Counts the 3 components
                for group in list_sub_group[:3]:    #Loop on the 3 currents
                    list_array=group._f_listNodes("Array") #List of all recorded arrays in "group"
                    for array in list_array:
                        num=scipy.int64(repr(array).split()[0][28:])
                        if num==cycle:
                            num_final=num
                            Fields_local[k,x_start:x_end,y_start:y_end,z_start:z_end]+=array.read()
                    k=k+1
            h5proc_file.close()
    #Write the .vtr file
        
        rtg = vtk.vtkRectilinearGrid()
        
        rtg.SetDimensions(Nxc+1,Nyc+1,Nzlocal)
        rtg.SetExtent(0,Nxc,0,Nyc,rankproc*Nzc/numproc,(rankproc+1)*Nzc/numproc)
        
        rtg.SetXCoordinates(vtkXcoordinates)
        
        rtg.SetYCoordinates(vtkYcoordinates)
        
        rtg.SetZCoordinates(vtkZcoordinates)
        
        Bx=Fields_local[0,:,:,:].swapaxes(0,2).flatten().astype("float32")
        By=Fields_local[1,:,:,:].swapaxes(0,2).flatten().astype("float32")
        Bz=Fields_local[2,:,:,:].swapaxes(0,2).flatten().astype("float32")
        
        vtk_dataBx = vtk.vtkFloatArray()
        vtk_dataBx.SetNumberOfTuples(numpoints)
        vtk_dataBx.SetNumberOfComponents(1)
        vtk_dataBx.SetVoidArray(Bx, numpoints, 1)	
        vtk_dataBx.SetName('Jx')	 
        vtk_dataBy = vtk.vtkFloatArray()
        vtk_dataBy.SetNumberOfTuples(numpoints)
        vtk_dataBy.SetNumberOfComponents(1)
        vtk_dataBy.SetVoidArray(By, numpoints, 1)	
        vtk_dataBy.SetName('Jy')	 
        vtk_dataBz = vtk.vtkFloatArray()
        vtk_dataBz.SetNumberOfTuples(numpoints)
        vtk_dataBz.SetNumberOfComponents(1)
        vtk_dataBz.SetVoidArray(Bz, numpoints, 1)	
        vtk_dataBz.SetName('Jz')	 
            
        rtg.GetPointData().SetScalars(vtk_dataBx)
        rtg.GetPointData().AddArray(vtk_dataBy)
        rtg.GetPointData().AddArray(vtk_dataBz)
        
        writer =vtk.vtkXMLRectilinearGridWriter()
        writer.SetFileName(directory+"currents_cycle"+repr(num_final).rjust(8,"0")+"-"+repr(rankproc)+".vtr")
        writer.SetInput(rtg)
        writer.Write()
    ## Write the pvtr file
        if rankproc == 0:
            output_pvtr_name=directory+"currents_cycle"+repr(num_final).rjust(8,"0")+".pvtr"
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
            grid.appendChild(coordinate)
            
            p_data = pvtr.createElementNS("VTK","PPointData")
            p_data.setAttribute("Scalars","Jx")
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","Jx")
            p_data.appendChild(data_array)
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","Jy")
            p_data.appendChild(data_array)
    
            data_array = pvtr.createElementNS("VTK", "PDataArray")
            data_array.setAttribute("type","Float32")
            data_array.setAttribute("Name","Jz")
            p_data.appendChild(data_array)
    
            grid.appendChild(p_data)
            
            for i in liste_proc:
                num=repr(i) 
                filename="currents_cycle"+repr(num_final).rjust(8,"0")+"-"+repr(i)+".vtr"
                piece = pvtr.createElementNS("VTK", "Piece")
                piece.setAttribute("Extent", "0 "+ repr(Nxc)+" 0 "+repr(Nyc)+" "+repr(i*Nzc/numproc)+" "+repr((i+1)*Nzc/numproc))
                piece.setAttribute("Source",filename)
                grid.appendChild(piece)
            
            pvtr.writexml(output_pvtr_id, newl='\n')
            output_pvtr_id.close() 
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
        
        pvd.writexml(output_pvd_id, newl='\n')
        output_pvd_id.close()

if rankproc == 0:
    print 'Time elapsed: ', time.clock() - start_time