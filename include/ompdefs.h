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
#endif

#endif
