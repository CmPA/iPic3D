""" Convert the iPIC3D output proc* files to pvtr readable by ParaView.

    USAGE
    Set directory, data_read and cycles below in the code. Then run:  
    mpirun -n 4 python proc2vtr.py

    NOTE
    In case of 4 species assumes there are background ones and sums the corresponding rho, j, and p's.

    You may wish to change the number of processors.
    
    This is my spin-off of Arnaud's read_proc_para.py routine.
    Note, if you, e.g., work with currents, where the sum: Jx = Jx0+Jx1.
    This sum would be different when computed here comparing to what you'll get in ParaView!
    So the result could be different from what you'd get from read_proc_para.py!    

    (c) V. Olshevsky: sya@mao.kiev.ua

"""
""" Modified by MEI to read specific variables needed for generalised Ohm's law 
    (c) Maria Elena Innocenti: mariaelena.innocenti@gmail.com

    Parallelisation in Y
"""

import xml.dom.minidom
import numpy, tables, vtk, time, os
from scipy import ndimage ### for the smoothing filter, added by MEI

try:
    from mpi4py import MPI
    mpi_enabled = True
except:
    mpi_enabled = False


"""TODO: Doesn't work when the number of CPUs don't divide evenly the number of domains in Z!"""
""" You SHOULD set the following parameters: """
# Directory of the simulation data

#directory = '/home/mariaelena/iPic3D_MC/iPic3D_MC/GEM_BCOct/'
#directory = '/home/mariaelena/GIT-HUB/iPic3D-Git/T14_DoubleGEMDiego/'
#directory= '/home/mariaelena/GIT-HUB/iPic3D-Git/CollapsingTest_tanh10/'
#directory= '/home/mariaelena/GIT-HUB/iPic3D-Git/GEM_dt0_25_Lambda10_COLL/'
#directory= '/mnt/Pause_SIMS/MC/GEM_dt0_25_Lambda10_Coll_BC3/'
#directory= '/mnt/Pause_SIMS/MC/GEM_NewGasPars_2_BC2/'
#directory= '/mnt/Pause_SIMS/MC/DoubleGEM_NoColl/'
directory= '/data/leuven/300/vsc30068/MC/DoubleGEM_NoColl/'
directory= '/data/leuven/300/vsc30068/MC/MRX_NoColl/'
out_directory = directory

# Which quantities to process
# Inertia: vector components for Ohm's law
# Inertia_SC: scalar components for Ohm's law - related studies
data_read = {'Fields': 1, 'Currents': 1, 'Rho': 1, 'Pressure': 1, 'Inertia': 1, 'Inertia_SC': 0, 'MC': 0}  #NB: rhoe_MA is still missing
#data_read = {'Fields': 0, 'Currents': 0, 'Rho': 0, 'Pressure': 0, 'Inertia_1': 0, 'Inertia_2':0, 'T1':0, 'T2':0, 'T3':0, 'T4':0, 'T5':0, 'T6':1} 
# Which cycles to process
first_cycle = 0
last_cycle = 200 #8000 #1000 #800 #1900 ##138 #2000 # 200 #1200
step =  100 #100 #100


def writePVD(pvd_name, list_pvtr, list_cycle):
    """write the pvd (VTK header file) file that contains all .pvtr file names
       
       Parameters:
        pvd_name -- outpur file name
        list_pvtr -- list of pvtr file names
        list_cycle -- cycles that correspond to pvtr files. 
                      A list of (timesteps of) the same length as list_pvtr

    """    
    pvd = xml.dom.minidom.Document()                                       
    pvd_root = pvd.createElementNS("VTK", "VTKFile")                       
    pvd_root.setAttribute("type", "Collection")                            
    pvd_root.setAttribute("version", "0.1")                                
    pvd_root.setAttribute("byte_order", "LittleEndian")                    
    pvd.appendChild(pvd_root)                                              

    collection = pvd.createElementNS("VTK", "Collection")                  
    pvd_root.appendChild(collection)                                       

    for i in range(len(list_pvtr)):
        dataSet = pvd.createElementNS("VTK", "DataSet")
        dataSet.setAttribute("timestep", str(list_cycle[i]))
        dataSet.setAttribute("group", "")
        dataSet.setAttribute("part", "0")
        dataSet.setAttribute("file", list_pvtr[i])
        collection.appendChild(dataSet)

    try:
        output_pvd_id = open(pvd_name, 'w')
        pvd.writexml(output_pvd_id, newl='\n')                                 
        output_pvd_id.close()
    except:
        print 'Error writing PVD file: ', pvd_name


def writePVTR(pvtr_name, cycle, vtr_prefix, scalar_fields, vector_fields, Nxc, Nyc, Nzc, numproc):
    """Writes the Parallel VTR file that contains info about particular VTR files
       NOTE, this routine implicitly assumes 3 dimensions!
        
       Parameters:
        pvtr_name - name of the output PVTR file
        cycle - simulation cycle index
        vtr_prefix - prefix for VTR file names, e.g. 'fields'
        scalar_fields - list of names of scalar fields in the VTR. May be empty!
        vector_fields - list of names of vector fields in the VTR. May be empty!
        Nxc, Nyc, Nzc - number of grid points in each dimension.
        numproc - number of processors (corresponds to the number of VTR files)

    """
    pvtr = xml.dom.minidom.Document()                                                                                    
    pvtr_root = pvtr.createElementNS("VTK", "VTKFile")                                                                   
    pvtr_root.setAttribute("type", "PRectilinearGrid")                                                                   
    pvtr_root.setAttribute("version", "0.1")                                                                             
    pvtr_root.setAttribute("byte_order", "LittleEndian")                                                                 
    pvtr.appendChild(pvtr_root)                                                                                          

    grid = pvtr.createElementNS("VTK", "PRectilinearGrid")                                                               
    grid.setAttribute("WholeExtent", "0 " + str(Nxc) + " 0 " + str(Nyc) + " 0 " + str(Nzc))                                    
    grid.setAttribute("GhostLevel","0")                                                                                  
    pvtr_root.appendChild(grid)                                                                                          

    coordinate = pvtr.createElementNS("VTK", "PCoordinates")                                                             
    for k in range(3):                                                                                    
        data_array = pvtr.createElementNS("VTK", "PDataArray")                                                           
        data_array.setAttribute("type", "Float32")                                                                        
        data_array.setAttribute("format", "Appended")                                                                     
        coordinate.appendChild(data_array)                                                                               

    grid.appendChild(coordinate)                                                                                         

    p_data = pvtr.createElementNS("VTK", "PPointData")    
    if len(scalar_fields) > 0: p_data.setAttribute("Scalars", scalar_fields[0])    
    if len(vector_fields) > 0: p_data.setAttribute("Vectors", vector_fields[0])
    for c in scalar_fields:
        data_array = pvtr.createElementNS("VTK", "PDataArray")
        data_array.setAttribute("type", "Float32")
        data_array.setAttribute("Name", c)
        p_data.appendChild(data_array)

    for c in vector_fields:
        data_array = pvtr.createElementNS("VTK", "PDataArray")
        data_array.setAttribute("type", "Float32")
        data_array.setAttribute("NumberOfComponents", "3")
        data_array.setAttribute("Name", c)
        p_data.appendChild(data_array)

    grid.appendChild(p_data)

    for i in range(numproc):
        vtr_name = vtr_prefix + str(cycle).rjust(8,"0") + "-" + str(i) + ".vtr"
        piece = pvtr.createElementNS("VTK", "Piece")
        #TwoDims
        #piece.setAttribute("Extent", "0 " + str(Nxc) + " 0 " + str(Nyc) + " " + str(i*Nzc/numproc) + " " + str((i+1)*Nzc/numproc))
        piece.setAttribute("Extent", "0 " + str(Nxc) + " "  + str(i*Nyc/numproc) + " " + str((i+1)*Nyc/numproc) +  " 0 " + str(Nzc) )
        piece.setAttribute("Source", vtr_name)
        grid.appendChild(piece)                                                                                          

    with open(pvtr_name, 'w') as pvtr_id:
        pvtr.writexml(pvtr_id, addindent='\t', newl='\n')


def writeVTR(vtr_name, scalar_fields, vector_fields, vtkX, vtkY, vtkZ,localYrange):  #localZrange):
    """Writes a single VTR file per Python processor/variable

    Parameters:
     vtr_name - name of the VTR file
     scalar_fields - dictionary with scalar field arrays ordered [x, y, z], e.g. {'p': array[nx,ny,nz], 'rho0': array[nx,ny,nz]}
     vector_fields - dictionary with vector fields ordered [3, x, y, z], e.g. {'J': array[3,nx,ny,nz], 'B': array[3,nx,ny,nz]}
     vtkX, vtkY, vtkZ - VTR coordinates
     localZrange - local range for Z indices

    """
    Nx = vtkX.GetNumberOfTuples() - 1
    Ny = vtkY.GetNumberOfTuples() - 1
    Nz = vtkZ.GetNumberOfTuples() - 1
    numpoints = (Nx+1)*(Ny+1)*(Nz+1)
    rtg = vtk.vtkRectilinearGrid()
    #TwoDims
    #rtg.SetExtent(0, Nx, 0, Ny, localZrange[0], localZrange[1])
    rtg.SetExtent(0, Nx, localYrange[0], localYrange[1], 0, Nz)     
    rtg.SetXCoordinates(vtkX)
    rtg.SetYCoordinates(vtkY)
    rtg.SetZCoordinates(vtkZ)
    vtk_data = []
    array_list = []
    for f in scalar_fields:
        vtk_data.append(vtk.vtkFloatArray())
        vtk_data[-1].SetNumberOfTuples(numpoints)
        vtk_data[-1].SetNumberOfComponents(1)
        array_list.append(scalar_fields[f].swapaxes(0,2).flatten().astype("float32"))
        vtk_data[-1].SetVoidArray(array_list[-1], numpoints, 1)
        vtk_data[-1].SetName(f)
        if f == scalar_fields.keys()[0]:
            rtg.GetPointData().SetScalars(vtk_data[-1])
        else:
            rtg.GetPointData().AddArray(vtk_data[-1])

    for f in vector_fields:
        vtk_data.append(vtk.vtkFloatArray())
        vtk_data[-1].SetNumberOfTuples(numpoints*3)
        vtk_data[-1].SetNumberOfComponents(3)
        array_list.append(vector_fields[f].swapaxes(0,3).swapaxes(1,2).flatten().astype("float32"))
        vtk_data[-1].SetVoidArray(array_list[-1], numpoints*3, 1)
        vtk_data[-1].SetName(f)
        if f == vector_fields.keys()[0]:
            rtg.GetPointData().SetVectors(vtk_data[-1])
        else:
            rtg.GetPointData().AddArray(vtk_data[-1])

    try:
        writer = vtk.vtkXMLRectilinearGridWriter()                                                                                              
        writer.SetFileName(vtr_name)
        if (vtk.vtkVersion.GetVTKMajorVersion() < 6):
            writer.SetInput(rtg)
        else:
            writer.SetInputData(rtg)
        writer.Write()
    except:
        print 'Error writing VTR file: ', vtr_name


def openProcFile(proc_filename, read_proc0, shape):
    """Open single proc*.hdf file and return its handler
        
       Parameters
         proc_filename - name of the input proc*hdf file
         read_proc0 == list_read_proc[0], the bottom processor layer index. See below.
         shape - the shape of each individual subdomain         

       Return
         the handle of the file to read
         bounds - coordinates of the edges of the domain saved in the proc file
    """
    try :
        h5proc_file = tables.open_file(proc_filename, mode = "r", title = "Proc_file")                 
        p_coordinate = h5proc_file.root.topology.cartesian_coord.read()
    except :
        print 'Error opening proc file: ', proc_filename

    #p_coordinate[2] = p_coordinate[2] - read_proc0    
    #bounds[0] = p_coordinate[0]*(shape[0]-1) # x_start -1 to account for the 1 node overlapping
    #bounds[1] = bounds[0] + shape[0]           # x_end
    #bounds[2] = p_coordinate[1]*(shape[1]-1) # y_start
    #bounds[3] = bounds[2] + shape[1]           # y_end
    #bounds[4] = p_coordinate[2]*(shape[2]-1) # z_start
    #bounds[5] = bounds[4] + shape[2]           # z_end

    p_coordinate[1] = p_coordinate[1] - read_proc0                                                                                                       
    bounds[0] = p_coordinate[0]*(shape[0]-1) # x_start -1 to account for the 1 node overlapping                                                         
    bounds[1] = bounds[0] + shape[0]           # x_end                                                                                                  
    bounds[2] = p_coordinate[1]*(shape[1]-1) # y_start                                                                                                  
    bounds[3] = bounds[2] + shape[1]           # y_end                                                                                                  
    bounds[4] = p_coordinate[2]*(shape[2]-1) # z_start                                                                                                  
    bounds[5] = bounds[4] + shape[2]           # z_end   
    return h5proc_file, bounds
    
    
""" Main program starts here """
# Keep all MPI-related stuff in one place; can be serial
if mpi_enabled:
    try:
        comm = MPI.COMM_WORLD
        mpi_myrank = comm.Get_rank() # Proc rank
        numproc = comm.Get_size() # Total number of procs
    except:
        print "MPI library is imported, but MPI not initialized. Trying to run serial."
        mpi_myrank = 0
        numproc = 1
else:
    print "MPI is not enabled, trying to run serial."
    mpi_myrank = 0
    numproc = 1

print 'mpi_enabled? :' + str(mpi_enabled)
print 'comm.Get_size() : ' + str(comm.Get_size() )

liste_proc = range(numproc)
# Create output folder if needed
if (mpi_myrank == 0) and (not os.path.exists(out_directory)):
    os.makedirs(out_directory)
if mpi_enabled: comm.Barrier()
if mpi_myrank == 0: start_time = time.clock()

""" Constants """
str_fields = 'fields'
str_rho = 'rho'
str_currents = 'currents'
str_pressure = 'pressure'
str_moments = 'moments'
str_species = 'species_'
str_inertia = 'inertia' ## OHM specific
str_MC= 'MC'
str_inertia_SC = 'inertia_SC' ## OHM specific  
current_comp_names = ['Jx', 'Jy', 'Jz']
inertia_comp_names = ['E_MA_x', 'E_MA_y', 'E_MA_z', 'Je_MA_x', 'Je_MA_y', 'Je_MA_z', 'T1_MA_x', 'T1_MA_y', 'T1_MA_z', 'T2_MA_x', 'T2_MA_y', 'T2_MA_z', 'T3_MA_x', 'T3_MA_y', 'T3_MA_z', 'T4_MA_x', 'T4_MA_y', 'T4_MA_z', 'B_MA_x', 'B_MA_y', 'B_MA_z'] ## OHM specific
#inertia_comp_names = ['T1_MA_x', 'T1_MA_y', 'T1_MA_z', 'T2_MA_x', 'T2_MA_y', 'T2_MA_z', 'T3_MA_x', 'T3_MA_y', 'T3_MA_z', 'T4_MA_x', 'T4_MA_y', 'T4_MA_z'] ## OHM specific
inertia_SC_comp_names = ['rhoe_MA']
field_comp_names = ['Bx', 'By', 'Bz', 'Ex', 'Ey', 'Ez']
field_comp_count = len(field_comp_names)
inertia_comp_count = len(inertia_comp_names)
inertia_SC_comp_count = len(inertia_SC_comp_names)
field_vec_names = ['B', 'E']
current_vec_name = 'J'
MC_comp_names = ['C_dPx', 'C_dPy', 'C_dPz']
MC_vec_name = 'C_dP'
inertia_vec_names = ['E_MA', 'Je_MA', 'T1_MA', 'T2_MA', 'T3_MA', 'T4_MA', 'B_MA'] ## NB: i still have to extract rhoe_MA
inertia_SC_vec_name = ['rhoe_MA']
#inertia_vec_names = ['T1_MA', 'T2_MA', 'T3_MA', 'T4_MA']
press_names = ['pXX', 'pXY', 'pXZ', 'pYY', 'pYZ', 'pZZ']
press_comp_count = len(press_names)
list_cycle = numpy.arange(first_cycle, last_cycle+1, step).tolist()
bounds = numpy.zeros(6) # local domain bounds: x_start, x_end, y_start, ...

""" Gathering information and initialization """
if mpi_myrank == 0: 
    print "Reading directory " + directory
setting_filename = os.path.join(directory, "settings.hdf") # Setting file name
try:
    h5setting_file = tables.open_file(setting_filename, mode = "r", title = "Setting_file")
except:
    if mpi_myrank == 0:
        print 'Error opening the settings file, check the base directory!'
    exit()

Nprocs = h5setting_file.root.topology.Nprocs.read()[0]
XLEN = h5setting_file.root.topology.XLEN.read()[0]
YLEN = h5setting_file.root.topology.YLEN.read()[0]
ZLEN = h5setting_file.root.topology.ZLEN.read()[0]
dx = h5setting_file.root.collective.Dx.read()[0]
dy = h5setting_file.root.collective.Dy.read()[0]
dz = h5setting_file.root.collective.Dz.read()[0]
Nxc = h5setting_file.root.collective.Nxc.read()[0]
Nyc = h5setting_file.root.collective.Nyc.read()[0]
Nzc = h5setting_file.root.collective.Nzc.read()[0]
ns = h5setting_file.root.collective.Ns.read()[0] 
h5setting_file.close()
# Shape of each individual subdomain. We assume that each processor has domain of the same size!
shape = (Nxc/XLEN+1, Nyc/YLEN+1, Nzc/ZLEN+1) 
species = range(ns)

#Zbyproc = ZLEN/numproc # Minimum Nomber of Z you have to read
#reste = ZLEN - Zbyproc*numproc
#if mpi_myrank < reste:
#    Zrange = numpy.arange(mpi_myrank*(Zbyproc+1),(mpi_myrank+1)*(Zbyproc+1))
#else:
#    Zrange = numpy.arange(reste*(Zbyproc+1)+(mpi_myrank-reste)*Zbyproc, reste*(Zbyproc+1)+(mpi_myrank-reste+1)*Zbyproc)
#list_read_proc = []
#for Z in Zrange:
#    list_read_proc.extend(range(Z,Nprocs,ZLEN))
#list_read_proc = numpy.array(list_read_proc)
#zread_proc = len(list_read_proc)/(YLEN*XLEN) #Number of proc file you want to read in the z direction
#Nzlocal = zread_proc*(shape[2] - 1)
#Fields_local = numpy.zeros((field_comp_count,Nxc+1,Nyc+1,Nzlocal+1), dtype="float32")
#numproc = min(numproc, ZLEN) # To account for specific case where more processors are used here than in the Z direction of the simulation

Ybyproc = YLEN/numproc # Minimum Nomber of Z you have to read                                                                                             
reste = YLEN - Ybyproc*numproc

print "Ybyproc: " + str(Ybyproc) + ", reste: " + str(reste) + ", numproc: " + str(numproc)
if mpi_myrank < reste:
    Yrange = numpy.arange(mpi_myrank*(Ybyproc+1),(mpi_myrank+1)*(Ybyproc+1))
else:
    Yrange = numpy.arange(reste*(Ybyproc+1)+(mpi_myrank-reste)*Ybyproc, reste*(Ybyproc+1)+(mpi_myrank-reste+1)*Ybyproc)
list_read_proc = []
for Y in Yrange:
    list_read_proc.extend(range(Y,Nprocs,YLEN))
list_read_proc = numpy.array(list_read_proc)
yread_proc = len(list_read_proc)/(ZLEN*XLEN) #Number of proc file you want to read in the z direction

Nylocal = yread_proc*(shape[1] - 1)
Nzlocal= Nzc
Fields_local = numpy.zeros((field_comp_count,Nxc+1,Nylocal+1,Nzc+1), dtype="float32")
numproc = min(numproc, YLEN) # To account for specific case where more processors are used here than in the Z direction of the simulation

print 'size shape: ' + str(Nxc/XLEN+1) + " " + str( Nyc/YLEN+1) + " " + str( Nzc/ZLEN+1)
print 'size Fields local: ' + str(Nxc+1) + " " + str(Nylocal+1) + " " + str(Nzc+1)
print 'Nylocal ' + str(Nylocal)
print 'yread_proc: ' + str(yread_proc)



if len(list_read_proc) > 0:
    #numpoints = (Nxc+1)*(Nyc+1)*(Nzlocal+1)
    #RangeZlocal = mpi_myrank*Nzc/numproc, (mpi_myrank+1)*Nzc/numproc
    #Xcoordinates = numpy.arange(Nxc+1).flatten().astype('float32')*dx
    #Ycoordinates = numpy.arange(Nyc+1).flatten().astype('float32')*dy
    #Zcoordinates = (numpy.arange(Nzlocal+1).flatten().astype('float32') + list_read_proc[0]*(shape[2]-1))*dz
    #vtkXcoordinates = vtk.vtkFloatArray()
    #vtkXcoordinates.SetNumberOfComponents(1)
    #vtkXcoordinates.SetNumberOfTuples(Nxc+1)
    #vtkXcoordinates.SetVoidArray(Xcoordinates,Nxc+1,1)
    #vtkYcoordinates = vtk.vtkFloatArray()
    #vtkYcoordinates.SetNumberOfComponents(1)
    #vtkYcoordinates.SetNumberOfTuples(Nyc+1)
    #vtkYcoordinates.SetVoidArray(Ycoordinates,Nyc+1,1)
    #vtkZcoordinates = vtk.vtkFloatArray()
    #vtkZcoordinates.SetNumberOfComponents(1)
    #vtkZcoordinates.SetNumberOfTuples(Nzlocal+1)
    #vtkZcoordinates.SetVoidArray(Zcoordinates,Nzlocal+1,1)
    #
    numpoints = (Nxc+1)*(Nylocal+1)*(Nzlocal+1)
    RangeYlocal = mpi_myrank*Nyc/numproc, (mpi_myrank+1)*Nyc/numproc
    Xcoordinates = numpy.arange(Nxc+1).flatten().astype('float32')*dx
    Zcoordinates = numpy.arange(Nzc+1).flatten().astype('float32')*dz
    Ycoordinates = (numpy.arange(Nylocal+1).flatten().astype('float32') + list_read_proc[0]*(shape[1]-1))*dy
    vtkXcoordinates = vtk.vtkFloatArray()
    vtkXcoordinates.SetNumberOfComponents(1)
    vtkXcoordinates.SetNumberOfTuples(Nxc+1)
    vtkXcoordinates.SetVoidArray(Xcoordinates,Nxc+1,1)
    vtkYcoordinates = vtk.vtkFloatArray()
    vtkYcoordinates.SetNumberOfComponents(1)
    vtkYcoordinates.SetNumberOfTuples(Nylocal+1)
    vtkYcoordinates.SetVoidArray(Ycoordinates,Nylocal+1,1)
    vtkZcoordinates = vtk.vtkFloatArray()
    vtkZcoordinates.SetNumberOfComponents(1)
    vtkZcoordinates.SetNumberOfTuples(Nzc+1)
    vtkZcoordinates.SetVoidArray(Zcoordinates,Nzc+1,1)
    #
    list_proc_files = [os.path.join(directory, "proc" + repr(read_proc) + ".hdf") for read_proc in list_read_proc]

    #print 'here'
    #print str(list_proc_files)

    """ 1. Fields: magnetic and electric """
    if data_read['Fields']: 
        if mpi_myrank == 0: pvtr_fields_list = []
        for cycle in list_cycle:    
            for proc_filename in list_proc_files:
                h5proc_file, bounds = openProcFile(proc_filename, list_read_proc[0], shape)
                ## Cycle over field components
                for c in range(field_comp_count):
                    for node in h5proc_file.walk_nodes('/' + str_fields + '/'+ field_comp_names[c]):
                        if type(node) == tables.array.Array:
                            num = int(node.name[6:])
                            if num == cycle:
                                Fields_local[c, bounds[0]:bounds[1], bounds[2]:bounds[3], bounds[4]:bounds[5]] = node.read()

                h5proc_file.close()

            vector_fields = {}
            for c in range(len(field_vec_names)): 
                vector_fields[field_vec_names[c]] = Fields_local[3*c:3*(c+1), :, :, :]
            writeVTR(os.path.join(out_directory, str_fields + repr(cycle).rjust(8,"0") + "-" + repr(mpi_myrank) + ".vtr"), {}, vector_fields,
                     vtkXcoordinates, vtkYcoordinates, vtkZcoordinates, RangeYlocal)
            if mpi_myrank == 0:
                pvtr_fields_list.append(str_fields + repr(cycle).rjust(8,"0") + ".pvtr")
                writePVTR(os.path.join(out_directory, pvtr_fields_list[-1]), cycle, str_fields, [], vector_fields.keys(), Nxc, Nyc, Nzc, numproc)

        if mpi_myrank == 0:
            writePVD(os.path.join(out_directory, str_fields + '.pvd'), pvtr_fields_list, list_cycle)
            print 'Fields written'

    """ 2. Rho (Charge density) """
    if data_read['Rho']: 
        if mpi_myrank == 0: pvtr_rho_list = []
        for cycle in list_cycle:    
            Fields_local = numpy.zeros((ns,Nxc+1,Nylocal+1,Nzc+1), dtype="float32")            
            for proc_filename in list_proc_files:
                #print list_proc_files.index(proc_filename), '/', len(list_proc_files), proc_filename
                h5proc_file, bounds = openProcFile(proc_filename, list_read_proc[0], shape)
                '''
                # Read the total density first
                for node in h5proc_file.walk_nodes('/' + str_moments + '/' + str_rho):
                    if type(node) == tables.array.Array:
                        num = int(node.name[6:])
                        if num == cycle:                    
                            Fields_local[0, bounds[0]:bounds[1], bounds[2]:bounds[3], bounds[4]:bounds[5]] = node.read()
                '''
                ## Cycle over species for each specie density
                for specie in species:
                    for node in h5proc_file.walk_nodes('/' + str_moments + '/'+ str_species + str(specie) + '/' + str_rho):
                        if type(node) == tables.array.Array:
                            num = int(node.name[6:])
                            if num == cycle:
                                Fields_local[specie, bounds[0]:bounds[1], bounds[2]:bounds[3], bounds[4]:bounds[5]] = node.read()

                h5proc_file.close()

            rho_fields = {} #{str_rho: Fields_local[0,:,:,:]}

            # Consider a case of background species
            if ns == 4:
                #print 'ns' + str(ns)
                for si in range(2):
                    rho_fields[str_rho + repr(si)] = Fields_local[si,:,:,:] + Fields_local[si+2,:,:,:]
            else:
                for specie in species: 
                    rho_fields[str_rho + repr(specie)] = Fields_local[specie,:,:,:]

            # separated rho
            #for specie in species: 
            #    rho_fields[str_rho + repr(specie)] = Fields_local[specie,:,:,:]
            
            writeVTR(os.path.join(out_directory, str_rho + repr(cycle).rjust(8,"0") + "-" + repr(mpi_myrank) + ".vtr"), rho_fields, {}, vtkXcoordinates, vtkYcoordinates, vtkZcoordinates, RangeYlocal)
            if mpi_myrank == 0:            
                pvtr_rho_list.append(str_rho + repr(cycle).rjust(8,"0") + ".pvtr")
                writePVTR(os.path.join(out_directory, pvtr_rho_list[-1]), cycle, str_rho, rho_fields.keys(), [], Nxc, Nyc, Nzc, numproc)

        if mpi_myrank == 0:
            writePVD(os.path.join(out_directory, str_rho + '.pvd'), pvtr_rho_list, list_cycle)
            print 'Rho written'

    """ 3. Currents """
    if data_read['Currents']:
        if mpi_myrank == 0: pvtr_currents_list = []
        for cycle in list_cycle:    
            Fields_local = numpy.zeros((ns*3,Nxc+1,Nylocal+1,Nzc+1), dtype="float32")
            for proc_filename in list_proc_files:
                h5proc_file, bounds = openProcFile(proc_filename, list_read_proc[0], shape)
                ## Cycle over species to get current for each specie
                for specie in species:
                    for c in range(3):
                        for node in h5proc_file.walk_nodes('/' + str_moments + '/'+ str_species + repr(specie) + '/' + current_comp_names[c]):
                            if type(node) == tables.array.Array:
                                num = int(node.name[6:])
                                if num == cycle:
                                    Fields_local[specie*3 + c, bounds[0]:bounds[1], bounds[2]:bounds[3], bounds[4]:bounds[5]] = node.read()

                h5proc_file.close()

            j_fields = {}
            
            # Consider a case of background species
            if ns == 4:
                for si in range(2):
                    j_fields[current_vec_name + repr(si)] = Fields_local[3*si:3*(si+1), :, :, :] + Fields_local[3*(si+2):3*(si+3), :, :, :]
            else:
                for specie in species: 
                    j_fields[current_vec_name + repr(specie)] = Fields_local[3*specie:3*(specie+1), :, :, :]
            # separated species  
            #for specie in species: 
            #    j_fields[current_vec_name + repr(specie)] = Fields_local[3*specie:3*(specie+1), :, :, :]
            
            writeVTR(os.path.join(out_directory, str_currents + repr(cycle).rjust(8,"0") + "-" + repr(mpi_myrank) + ".vtr"), {}, j_fields,
                     vtkXcoordinates, vtkYcoordinates, vtkZcoordinates, RangeYlocal)
            if mpi_myrank == 0:
                pvtr_currents_list.append(str_currents + repr(cycle).rjust(8,"0")+".pvtr")
                writePVTR(os.path.join(out_directory, pvtr_currents_list[-1]), cycle, str_currents, [], j_fields.keys(), Nxc, Nyc, Nzc, numproc)

        if mpi_myrank == 0:
            writePVD(os.path.join(out_directory, str_currents + '.pvd'), pvtr_currents_list, list_cycle)
            print 'Currents written'

    """ 4. Pressure tensor for each specie"""
    if data_read['Pressure']: 
        if mpi_myrank == 0: pvtr_press_list = []
        for cycle in list_cycle:    
            Fields_local = numpy.zeros((ns*6, Nxc+1, Nylocal+1, Nzc+1), dtype="float32")            
            for proc_filename in list_proc_files:
                h5proc_file, bounds = openProcFile(proc_filename, list_read_proc[0], shape)
                for specie in species:
                    for c in range(press_comp_count):
                        for node in h5proc_file.walk_nodes('/' + str_moments + '/'+ str_species + str(specie) + '/' + press_names[c]):
                            if type(node) == tables.array.Array:
                                num = int(node.name[6:])
                                if num == cycle:
                                    Fields_local[press_comp_count*specie+c, bounds[0]:bounds[1], bounds[2]:bounds[3], bounds[4]:bounds[5]] = node.read()

                h5proc_file.close()

            p_fields = {}

            # A case of background species
            if ns == 4:
                for si in range(2):
                    for c in range(press_comp_count):
                        p_fields[press_names[c] + repr(si)] = Fields_local[press_comp_count*si+c,:,:,:] + Fields_local[press_comp_count*(si+2)+c,:,:,:]
            else:
                for specie in species:
                    for c in range(press_comp_count):
                        p_fields[press_names[c] + repr(specie)] = Fields_local[press_comp_count*specie+c,:,:,:]

            writeVTR(os.path.join(out_directory, str_pressure + repr(cycle).rjust(8,"0") + "-" + repr(mpi_myrank) + ".vtr"), p_fields, {}, 
                     vtkXcoordinates, vtkYcoordinates, vtkZcoordinates, RangeYlocal)
            if mpi_myrank == 0:            
                pvtr_press_list.append(str_pressure + repr(cycle).rjust(8,"0") + ".pvtr")
                writePVTR(os.path.join(out_directory, pvtr_press_list[-1]), cycle, str_pressure, p_fields.keys(), [], Nxc, Nyc, Nzc, numproc)

        if mpi_myrank == 0:
            writePVD(os.path.join(out_directory, str_pressure + '.pvd'), pvtr_press_list, list_cycle)
            print 'Pressure written'

    ### here
    """ 5. Inertia: Ohm's law component """
    if data_read['Inertia']: 
        if mpi_myrank == 0: pvtr_inertia_list = []
        for cycle in list_cycle:  
            Fields_local = numpy.zeros((inertia_comp_count,Nxc+1,Nylocal+1,Nzc+1), dtype="float32")
            for proc_filename in list_proc_files:
                h5proc_file, bounds = openProcFile(proc_filename, list_read_proc[0], shape)
                ## Cycle over field components
                for c in range(inertia_comp_count):
                    #print 'Extracting ' + str(inertia_comp_names[c])
                    for node in h5proc_file.walk_nodes('/' + str_inertia + '/'+ inertia_comp_names[c]):
                        if type(node) == tables.array.Array:
                            num = int(node.name[6:])
                            if num == cycle:
                                Fields_local[c, bounds[0]:bounds[1], bounds[2]:bounds[3], bounds[4]:bounds[5]] = node.read()

                h5proc_file.close()

            vector_inertia = {}
            for c in range(len(inertia_vec_names)): 
                vector_inertia[inertia_vec_names[c]] = Fields_local[3*c:3*(c+1), :, :, :]
            writeVTR(os.path.join(out_directory, str_inertia + repr(cycle).rjust(8,"0") + "-" + repr(mpi_myrank) + ".vtr"), {}, vector_inertia,
                     vtkXcoordinates, vtkYcoordinates, vtkZcoordinates, RangeYlocal)
            if mpi_myrank == 0:
                pvtr_inertia_list.append(str_inertia + repr(cycle).rjust(8,"0") + ".pvtr")
                writePVTR(os.path.join(out_directory, pvtr_inertia_list[-1]), cycle, str_inertia, [], vector_inertia.keys(), Nxc, Nyc, Nzc, numproc)

        if mpi_myrank == 0:
            writePVD(os.path.join(out_directory, str_inertia+ '.pvd'), pvtr_inertia_list, list_cycle)
            print 'Inertia written'
    
    ### here, scalar components for Ohm-related things
    """ 5. Inertia: Ohm's law component """
    if data_read['Inertia_SC']: 
        if mpi_myrank == 0: pvtr_inertia_SC_list = []
        for cycle in list_cycle:  
            Fields_local = numpy.zeros((inertia_SC_comp_count,Nxc+1,Nylocal+1,Nzc+1), dtype="float32")
            for proc_filename in list_proc_files:
                h5proc_file, bounds = openProcFile(proc_filename, list_read_proc[0], shape)
                ## Cycle over field components
                for c in range(inertia_SC_comp_count):
                    #print 'Extracting ' + str(inertia_SC_comp_names[c])
                    for node in h5proc_file.walk_nodes('/' + str_inertia + '/'+ inertia_SC_comp_names[c]):
                        if type(node) == tables.array.Array:
                            num = int(node.name[6:])
                            if num == cycle:
                                Fields_local[c, bounds[0]:bounds[1], bounds[2]:bounds[3], bounds[4]:bounds[5]] = node.read()
                                #print 'fucking here: ' + str(Fields_local)
                h5proc_file.close()

            vector_inertia_SC = {}
            for c in range(len(inertia_SC_comp_names)): 
                vector_inertia_SC[inertia_SC_comp_names[c]] = Fields_local[c:1*(c+1), :, :, :]

            writeVTR(os.path.join(out_directory, str_inertia_SC + repr(cycle).rjust(8,"0") + "-" + repr(mpi_myrank) + ".vtr"), vector_inertia_SC, {},
                     vtkXcoordinates, vtkYcoordinates, vtkZcoordinates, RangeYlocal)
            if mpi_myrank == 0:
                pvtr_inertia_SC_list.append(str_inertia_SC + repr(cycle).rjust(8,"0") + ".pvtr")
                #print 'vector_inertia_SC.keys(): ' + str(vector_inertia_SC.keys())
                writePVTR(os.path.join(out_directory, pvtr_inertia_SC_list[-1]), cycle, str_inertia_SC, vector_inertia_SC.keys(), [], Nxc, Nyc, Nzc, numproc)

        if mpi_myrank == 0:
            writePVD(os.path.join(out_directory, str_inertia_SC+ '.pvd'), pvtr_inertia_SC_list, list_cycle)
            print 'Inertia SC written' 

    """ MC """
    if data_read['MC']:
        if mpi_myrank == 0: pvtr_MC_list = []
        for cycle in list_cycle:    
            Fields_local = numpy.zeros((ns*3,Nxc+1,Nylocal+1,Nzc+1), dtype="float32")
            for proc_filename in list_proc_files:
                h5proc_file, bounds = openProcFile(proc_filename, list_read_proc[0], shape)
                ## Cycle over species to get MC for each specie
                for specie in species:
                    for c in range(3):
                        for node in h5proc_file.walk_nodes('/' + str_MC + '/'+ str_species + repr(specie) + '/' + MC_comp_names[c]):
                            if type(node) == tables.array.Array:
                                num = int(node.name[6:])
                                if num == cycle:
                                    Fields_local[specie*3 + c, bounds[0]:bounds[1], bounds[2]:bounds[3], bounds[4]:bounds[5]] = node.read()

                h5proc_file.close()

            MC_fields = {}
            
            # Consider a case of background species
            if ns == 4:
                for si in range(2):
                    MC_fields[MC_vec_name + repr(si)] = Fields_local[3*si:3*(si+1), :, :, :] + Fields_local[3*(si+2):3*(si+3), :, :, :]
            else:
                for specie in species: 
                    MC_fields[MC_vec_name + repr(specie)] = Fields_local[3*specie:3*(specie+1), :, :, :]
            
            #print 'MC_fields' + str(MC_fields)
            writeVTR(os.path.join(out_directory, str_MC + repr(cycle).rjust(8,"0") + "-" + repr(mpi_myrank) + ".vtr"), {}, MC_fields,
                     vtkXcoordinates, vtkYcoordinates, vtkZcoordinates, RangeYlocal)
            if mpi_myrank == 0:
                pvtr_MC_list.append(str_MC + repr(cycle).rjust(8,"0")+".pvtr")
                writePVTR(os.path.join(out_directory, pvtr_MC_list[-1]), cycle, str_MC, [], MC_fields.keys(), Nxc, Nyc, Nzc, numproc)

        if mpi_myrank == 0:
            writePVD(os.path.join(out_directory, str_MC + '.pvd'), pvtr_MC_list, list_cycle)
            print 'MC components written'


if mpi_myrank == 0: print 'Time elapsed: ', time.clock() - start_time
