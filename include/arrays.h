#ifndef IPIC_ARRAYS_H
#define IPIC_ARRAYS_H
#include "Alloc.h" // variable-dimension arrays
/*
    Fixed array class developed by
      Alec Johnson

    For examples of use of this class,
    see test_arrays.cpp
*/

/*** begin FixedArray classes for use when dimensions are known at compile time. ***/
//
// These classes improve upon native fixed arrays as follows:
// - bounds-checking is performed if CHECK_BOUNDS is defined,
// - myarray(i,j) access is supported,  and
// - functions can return fixed-dimension arrays,
//   whereas since the C standard does not allow
//   one to return an array with fixed dimensions.
//
// The purpose of implementing these extensions is so that
// fixed-dimension arrays can be used in iPic3D if doing so
// yields a significant performance benefit for the choice of
// architecture and compiler.

template <class type, size_t s1>
class FixedArray1D
{
 public:
  type arr[s1];
 public:
  type& fetch(size_t n1)
  {
    check_bounds(n1,s1);
    return arr[n1];
  }
  type& operator[](size_t n1)
  {
    check_bounds(n1,s1);
    return arr[n1];
  }
};

// auxiliary class for chained operator[] dereferencing of FixedArray2D
//
template <class type, size_t s1, size_t s2>
class FixedArray2D1
{
  type (&arr)[s1][s2];
  size_t n1;
 public:
  FixedArray2D1(type (&_arr)[s1][s2], size_t _n1) :
    arr(_arr), n1(_n1) {};
  type& operator[](size_t n2)
  {
    check_bounds(n1,s1);
    check_bounds(n2,s2);
    return arr[n1][n2];
  }
};

template <class type, size_t s1, size_t s2>
class FixedArray2D
{
 public:
  type arr [s1][s2];
 public:
  type& fetch(size_t n1, size_t n2)
  {
    check_bounds(n1,s1);
    check_bounds(n2,s2);
    return arr[n1][n2];
  }
  // Chaining operator[] this way
  // does not allow bounds checking
  // and does not work beyond 2D.
  //type* operator[](size_t n1) { return arr[n1]; }
  FixedArray2D1<type,s1,s2> operator[](size_t n1)
    { return FixedArray2D1<type,s1,s2>(arr,n1); }
};

// auxiliary classes for chained operator[] dereferencing of FixedArray3D
//
template <class type, size_t s1, size_t s2, size_t s3>
class FixedArray3D1
{
  type (&arr)[s1][s2][s3];
  size_t n1, n2;
 public:
  FixedArray3D1(type (&_arr)[s1][s2][s3], size_t _n1, size_t _n2) :
    arr(_arr), n1(_n1), n2(_n2) {};
  type& operator[](size_t n3)
  {
    check_bounds(n1,s1);
    check_bounds(n2,s2);
    check_bounds(n3,s3);
    return arr[n1][n2][n3];
  }
};
//
template <class type, size_t s1, size_t s2, size_t s3>
class FixedArray3D2
{
  type (&arr)[s1][s2][s3];
  size_t n1;
 public:
  FixedArray3D2(type (&_arr)[s1][s2][s3], size_t _n1) : arr(_arr), n1(_n1) {};
  FixedArray3D1<type,s1,s2,s3> operator[](size_t n2)
    { return FixedArray3D1<type,s1,s2,s3>(arr,n1,n2); }
};

template <class type, size_t s1, size_t s2, size_t s3>
struct FixedArray3D
{
  type arr [s1][s2][s3];
 public:
  type& fetch(size_t n1, size_t n2, size_t n3)
  {
    check_bounds(n1,s1);
    check_bounds(n2,s2);
    check_bounds(n3,s3);
    return arr[n1][n2][n3];
  }
  // chained operator[] dereferencing requires
  // auxiliary array dereferencing classes,
  // since the C standard does not allow one to
  // return an array with fixed dimensions.
  FixedArray3D2<type,s1,s2,s3> operator[](size_t n1)
    { return FixedArray3D2<type,s1,s2,s3>(arr,n1); }
};

// auxiliary classes for chained operator[] dereferencing of FixedArray4D
//
template <class type, size_t s1, size_t s2, size_t s3, size_t s4>
class FixedArray4D1
{
  type (&arr)[s1][s2][s3][s4];
  size_t n1,n2,n3;
 public:
  FixedArray4D1(type(&_arr)[s1][s2][s3][s4],size_t _n1,size_t _n2,size_t _n3):
    arr(_arr), n1(_n1), n2(_n2), n3(_n3){};
  type& operator[](size_t n4) { return arr[n1][n2][n3][n4]; }
};
//
template <class type, size_t s1, size_t s2, size_t s3, size_t s4>
class FixedArray4D2
{
  type (&arr)[s1][s2][s3][s4];
  size_t n1,n2;
 public:
  FixedArray4D2(type (&_arr)[s1][s2][s3][s4], size_t _n1, size_t _n2) :
    arr(_arr), n1(_n1), n2(_n2){};
  FixedArray4D1<type,s1,s2,s3,s4> operator[](size_t n3)
    { return FixedArray4D1<type,s1,s2,s3,s4>(arr,n1,n2,n3); }
};
//
template <class type, size_t s1, size_t s2, size_t s3, size_t s4>
class FixedArray4D3
{
  type (&arr)[s1][s2][s3][s4];
  size_t n1;
 public:
  FixedArray4D3(type (&_arr)[s1][s2][s3][s4], size_t _n1) :
    arr(_arr), n1(_n1) {};
  FixedArray4D2<type,s1,s2,s3,s4> operator[](size_t n2)
    { return FixedArray4D2<type,s1,s2,s3,s4>(arr,n1,n2); }
};

template <class type, size_t s1, size_t s2, size_t s3, size_t s4>
class FixedArray4D
{
 public:
  type arr [s1][s2][s3][s4];
 public:
  type& fetch(size_t n1, size_t n2, size_t n3, size_t n4)
  {
    check_bounds(n1,s1);
    check_bounds(n2,s2);
    check_bounds(n3,s3);
    check_bounds(n4,s4);
    return arr[n1][n2][n3][n4];
  }
  FixedArray4D3<type,s1,s2,s3,s4> operator[](size_t n1)
    { return FixedArray4D3<type,s1,s2,s3,s4>(arr,n1); }
};
/*** end FixedArray classes for use when dimensions are known at compile time. ***/

#endif
