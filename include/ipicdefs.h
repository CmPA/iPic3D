#ifndef __IPIC_DEFS_H__
#define __IPIC_DEFS_H__

// comment this out if OpenMP is not installed on your system.
#define USING_OMP

// uncomment the following line to use parallel hdf5
//#define USING_PARALLEL_HDF5

// use precprocessor to remove MPI_Barrier() calls.
#define MPI_Barrier(args...)

#endif
