#ifndef Moments_H
#define Moments_H
#include "Alloc.h"

// class to accumulate node-centered species moments
// 
class Moments10
{
  private:
    arr4_double arr;
    int nx;
    int ny;
    int nz;
  public:
    void set_to_zero();

    // fetch accessors (write access)
    arr4_double fetch_arr() { return arr; }

    Moments10(int nxn, int nyn, int nzn) :
      nx(nxn),
      ny(nyn),
      nz(nzn),
      arr (nxn, nyn, nzn,10)
    {
    };
    ~Moments10(){};
};

#endif
