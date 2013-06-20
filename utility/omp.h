#ifndef omp_stubs_H
#define omp_stubs_H
#ifdef _OPENMP
  #include <omp.h>
#else // stubs for OpenMP methods that we use

inline int omp_get_thread_num() {
    return 0;
}

#endif
#endif
