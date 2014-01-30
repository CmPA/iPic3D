#ifndef ompdefs_H
#define ompdefs_H

#include <stdio.h>
#include "asserts.h"
// the compiler sets _OPENMP if the -openmp flag is used
#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num() { return 0;}
inline int omp_get_max_threads(){ return 1;}
#define omp_set_num_threads(num_threads)
#endif

class Caller_to_SetMaxThreadsForScope{
 int max_threads;
 public:
  Caller_to_SetMaxThreadsForScope(int i)
  {
    max_threads = omp_get_max_threads();
    // omp_set_num_threads should have been
    // called omp_set_max_threads
    omp_set_num_threads(i);
  }
  ~Caller_to_SetMaxThreadsForScope()
  {
    // restore the original maximum number of threads
    omp_set_num_threads(max_threads);
  }
};

#define set_max_threads_for_scope(num_threads) \
  Caller_to_SetMaxThreadsForScope \
  instanceOfCaller_to_SetMaxThreadsForScope(num_threads);

#endif
