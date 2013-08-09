#ifndef ompdefs_H
#define ompdefs_H

#include <stdio.h>
#include "asserts.h"
// the compiler sets _OPENMP if the -openmp flag is used
#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_thread_num() {
    return 0;
}
#endif

inline int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    #ifndef _OPENMP // USING_OMP
    assert_eq(n,1);
    #endif
    return n;
}

#endif
