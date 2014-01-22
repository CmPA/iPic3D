#ifndef _ipicmath_h_
#define _ipicmath_h_
#include "assert.h"

// valid if roundup power is representable.
inline int
pow2roundup (int x)
{
    assert(x>=0);
    //if (x < 0)
    //    return 0;
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return x+1;
}

// does not work if highest non-sign bit is set
inline int
pow2rounddown (int x)
{
    assert(x>=0);
    //if (x < 0)
    //    return 0;

    // set all bits below highest bit
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    // set the bit higher than the highest bit
    x++;
    // shift it down and return it
    return (x >> 1);
}

// round n up to next multiple of m
inline int roundup_to_multiple(int n, int m)
{
  //return ((n-1)/m+1)*m;
  return (n+m-1)/m*m;
}

#endif
