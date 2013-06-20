#ifndef ompdefs_H
#define ompdefs_H

#include <stdio.h>
#include "defs.h"
#include "omp.h"
#include "asserts.h"

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
