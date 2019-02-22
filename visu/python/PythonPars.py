#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 11:03:15 2019

@author: innocent

Reads particle data from iPic3D sims

(c) M.E.Innocenti: mariaelena.innocenti@gmail.com

"""

import numpy as np
import h5py, tables, os
import collections   # I need it for the ordered dict
import scipy.io as sio
try:
    from mpi4py import MPI
    mpi_enabled = True
    print('mpi_enabled')
except:
    mpi_enabled = False

    
    
directory="/Users/innocent/Documents/MLMD/iPic3D/test/"
directory="/Users/innocent/Documents/ExpandingBox/iPic3D/NonEB_Lx10_resonantPars/" # UNO A CASO
directory="/Users/innocent/Documents/ExpandingBox/iPic3D/Resonance_EB_SavingOften_ByBz_EyEz_n1_Lx10_NewINIT_deltaB_002_Longer/"
directory="/Users/innocent/Documents/ExpandingBox/iPic3D/NonEB_Lx10_resonantPars_scr_R230_dt4/"
out_directory= directory + 'partFromPy/'

#NB: no error handling if file does not exist
first_cycle = 0
last_cycle = 10000 #20000 # 101 
step = 200

# which particle properties to extract
what=['x', 'y', 'z', 'q', 'u', 'v', 'w'];

cycles=range(first_cycle, last_cycle+1, step)
        
setting_filename = directory + "settings.hdf" # Setting file name
try:

    h5setting_file = tables.open_file(setting_filename, mode = "r")  
    
except:
    if mpi_myrank == 0:
        print('Error opening the settings file, check the base directory!')
    exit()
    
ns = h5setting_file.root.collective.Ns.read()[0] 
Nprocs = h5setting_file.root.topology.Nprocs.read()[0]
XLEN = h5setting_file.root.topology.XLEN.read()[0]
YLEN = h5setting_file.root.topology.YLEN.read()[0]
ZLEN = h5setting_file.root.topology.ZLEN.read()[0]



#print('ns is ', ns)

def traverse_datasets(hdf_file):

    def h5py_dataset_iterator(g, prefix=''):
        for key in g.keys():
            item = g[key]
            path = f'{prefix}/{key}'
            if isinstance(item, h5py.Dataset): # test for dataset
                yield (path, item)
            elif isinstance(item, h5py.Group): # test for group (go down)
                yield from h5py_dataset_iterator(item, path)

    with h5py.File(hdf_file, 'r') as f:
        for path, _ in h5py_dataset_iterator(f):
            yield path


# extract the (few) info that I need to read particles

         
#with h5py.File(filename, 'r') as f:
#    for dset in traverse_datasets(filename):
#        print('Path:', dset)
#        print('Shape:', f[dset].shape)
#        print('Data type:', f[dset].dtype)
        
# topology info
#with h5py.File(filename, 'r') as f:
#    for dset in traverse_datasets(filename):
#        if "particles" in dset:
#            print('Path:', dset)
#            print('Shape:', f[dset].shape)
#            print('Data type:', f[dset].dtype)
            



if mpi_enabled:
    try:
        comm = MPI.COMM_WORLD
        mpi_myrank = comm.Get_rank() # Proc rank
        numproc = comm.Get_size() # Total number of procs
    except:
        print("MPI library is imported, but MPI not initialized. Trying to run serial.")
        mpi_myrank = 0
        numproc = 1
else:
    print("MPI is not enabled, trying to run serial.")
    mpi_myrank = 0
    numproc = 1


# distribute work when parallel
my_allprocs= np.arange(mpi_myrank, XLEN*YLEN*ZLEN, numproc)
print('my mpi_myrank is ', mpi_myrank, ' my my_allprocs is ', my_allprocs)

if (not os.path.exists(out_directory) and mpi_myrank== 0 ):
    print(str(mpi_myrank) + ' is here')
    os.makedirs(out_directory)

###         

for num in my_allprocs:
    filename= directory+ 'proc' + str(num) + '.hdf'

    # create the empty dict, also create the two lists to zip for the output
    names_list=[]
    
    for wh in what:
        for sp in range(ns):
            # empty dict
            exec( wh +str(sp)+'_d=  collections.OrderedDict()')
            # list of names
            names_list.append(wh +str(sp))
            

    with h5py.File(filename, 'r') as f:
        for cyc in cycles: # ordered dict remembers the order of insertion
            for wh in what:
                for sp in range(ns):
                    gr= '/particles/species_' + str(sp) +'/' +wh +'/cycle_' + str(cyc)            
                    #print(gr)
                    test= f[gr].value  
                    exec( wh +str(sp) +'_d[' + str(cyc) +'] ='+ 'test')
                
    

    # this for the output
    out_l1=['NS','NumHDF','cycles']
    out_l2=[sp, XLEN*YLEN*ZLEN, cycles ]
    for wh in what:
        for sp in range(ns):
            out_l1.append( wh +str(sp) + '_pr' + str(num))
            str1= 'list('+ wh +str(sp) +'_d.values())'
            exec('out_l2.append (' + str1 +')')
    
    out_dict=  dict(zip(out_l1, out_l2))
    sio.savemat(out_directory+ 'part_proc'+ str(num), out_dict)
