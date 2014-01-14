#ifndef _parallel_h_
#define _parallel_h_
/*********************************
 * General header for parallelism
 * (MPI, OpenMP, and SIMD)
 *********************************/

#include "MPIdata.h"
#include "ompdefs.h"

/*! used to restrict output to a single thread of a single process */
//inline bool is_main_master_thread()
inline bool is_output_thread()
{
  return !(MPIdata::get_rank() || omp_get_thread_num());
}

#endif
