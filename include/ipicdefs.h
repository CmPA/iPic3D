#ifndef __IPIC_DEFS_H__
#define __IPIC_DEFS_H__

// comment this out if OpenMP is not installed on your system.
#define USING_OMP

// uncomment the following line to use parallel hdf5
//#define USING_PARALLEL_HDF5

// use precprocessor to remove former MPI_Barrier() calls.
//#define MPI_Barrier(args...)
#define former_MPI_Barrier(args...)

#define ipicMPI_Allreduce(args...) \
  { \
    static int count=0; \
    dprint(count++); \
    MPI_Allreduce(## args); \
  }

//#define SINGLE_PRECISION_PCLS
//
// single precision does not seem to help on the MIC
#ifdef SINGLE_PRECISION_PCLS
  typedef float pfloat;
  #ifdef __MIC__
    #define VECTOR_WIDTH 16
  #else
    #define VECTOR_WIDTH 8
  #endif
#else
  #ifdef __MIC__
    #define VECTOR_WIDTH 8
  #else
    #define VECTOR_WIDTH 4
  #endif
  typedef double pfloat;
#endif

#endif
