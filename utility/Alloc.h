/***************************************************************************

  -------------------

developers           : Stefano Markidis, Giovanni Lapenta

This file was rewritten by eaj (Alec Johnson)

 ***************************************************************************/
#ifndef Alloc_H
#define Alloc_H
#include "stddef.h"
#include "assert.h"
// to maintain 64-bit memory alignment
// 1 would be okay for double, but 2 is needed for int
#define ALLOC_SHIFT 2


/** subroutines for allocation and deallocation of arrays 2D, 3D, 4D */

// Methods to store array dimensions.
//
// These methods will trash memory if ALLOC_SHIFT*sizeof(type)
// is smaller than 2*sizeof(int)
//
inline void set_size(void *arr, int sz) {
  int* ptr = (int*) arr;
  ptr[-1] = sz;
}
inline void set_sizes(void *arr, int sz) {
  int* ptr = (int*) arr;
  ptr[-2] = ptr[-1] = sz;
}
inline int get_size(void * arr) {
  int* ptr = (int*) arr;
  return ptr[-1];
}
inline int get_fullsize(void * arr) {
  int* ptr = (int*) arr;
  return ptr[-2];
}

// macros to allocate arrays
//
#define newArr1(type,sz1) _new_1d_array((sz1),(type *) NULL)
#define newArr2(type,sz1,sz2) _new_2d_array((sz1),(sz2),(type *) NULL)
#define newArr3(type,sz1,sz2,sz3) _new_3d_array((sz1),(sz2),(sz3),(type *) NULL)
#define newArr4(type,sz1,sz2,sz3,sz4) \
    _new_4d_array((sz1),(sz2),(sz3),(sz4),(type *) NULL);

// methods to allocate arrays
//
template < class type > inline type * _new_1d_array(int sz1, type * dummy) {
  type *arr = new type [sz1+ALLOC_SHIFT]+ALLOC_SHIFT;
  set_sizes(arr,sz1);
  return arr;
}
template < class type > inline type ** _new_2d_array(int sz1, int sz2, type * dummy) {
  type **arr = new type *[sz1+ALLOC_SHIFT]+ALLOC_SHIFT;
  type *ptr = newArr1(type, sz1*sz2);
  for (int i = 0; i < sz1; i++)
  {
    arr[i] = ptr;
    ptr += sz2;
  }
  set_sizes(arr,sz1);
  set_size(*arr,sz2);
  return arr;
}
template < class type > inline type *** _new_3d_array(int sz1, int sz2, int sz3, type * dummy) {
  type ***arr = new type **[sz1+ALLOC_SHIFT]+ALLOC_SHIFT;
  type **ptr = newArr2(type, sz1*sz2, sz3);
  for (int i = 0; i < sz1; i++)
  {
    arr[i] = ptr;
    ptr += sz2;
  }
  set_sizes(arr,sz1);
  set_size(*arr,sz2);
  return arr;
}
template < class type > inline type **** _new_4d_array(int sz1, int sz2, int sz3, int sz4, type * dummy) {
  type ****arr = (new type ***[sz1+ALLOC_SHIFT])+ALLOC_SHIFT;
  type ***ptr = newArr3(type, sz1*sz2, sz3, sz4);
  for (int i = 0; i < sz1; i++) {
    arr[i] = ptr;
    ptr += sz2;
  }
  set_sizes(arr,sz1);
  set_size(*arr,sz2);
  return arr;
}

// methods to deallocate arrays
//
template < class type > inline void delArr1(type * arr) {
  delete[](arr-ALLOC_SHIFT);
}
template < class type > inline void delArr2(type ** arr, int sz1) {
  assert(get_size(arr)==sz1);
  delArr1(arr[0]);
  delete[](arr-ALLOC_SHIFT);
}
template < class type > inline void delArr3(type *** arr, int sz1, int sz2) {
  assert(get_size(arr)==sz1);
  delArr2(arr[0],sz2);
  delete[](arr-ALLOC_SHIFT);
}
template < class type > inline void delArr4(type **** arr, int sz1, int sz2, int sz3) {
  assert(get_size(arr)==sz1);
  delArr3(arr[0],sz2,sz3);
  delete[](arr-ALLOC_SHIFT);
}

//*********************************************
//**** classes that adopt multidimensional ****
//****  arrays and efficiently access them ****
//*********************************************

// classes to dereference arrays.
//
// The purpose of this class is to allow elements of multidimensional
// arrays to be accessed with a calculated one-dimensional index while
// using the same syntax as is used for a nested array.  This gives
// correct results, but unfortunately is too slow; evidently the
// compiler is not intelligent enough to figure out that the whole chain
// of calls to the operator[] methods and to the DerefN constructors
// reduces to computing a one-dimensional subscript used to access a
// one-dimensional array.
//
struct Deref1
{
  double* arr;
  int* sizes;
  int k;
  inline Deref1(double* arr_in, int* sizes_in, int k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline double& operator[](int idx){
    k *= sizes[0]; k += idx;
    return arr[k];
  }
};

struct Deref2
{
  double* arr;
  int* sizes;
  int k;
  inline Deref2(double* arr_in, int* sizes_in, int k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline Deref1 operator[](int idx){
    k *= sizes[0]; k += idx;
    return Deref1(arr, sizes+1, k);
  }
};

struct Deref3
{
  double* arr;
  int* sizes;
  int k;
  inline Deref3(double* arr_in, int* sizes_in, int k_in) :
    arr(arr_in),
    sizes(sizes_in),
    k(k_in) {}
  inline Deref2 operator[](int idx){
    k *= sizes[0]; k += idx;
    return Deref2(arr, sizes+1, k);
  }
};

// doubleArrN adopts an array allocated by newArrN
//
// The purpose of these classes is to provide more efficient
// and more regulated access to array elements.  The idea is to
// maintain backward compatibility while allowing us to move
// toward a proper array abstraction.
//
// proposed improvements:
// - alternative constructor that takes a list of array dimensions
// - reference counting and destructor that deallocates if appropriate
// - methods that use parallel arithmetic for omp and vectorized code

// class that can adopt array allocated by newArr1
class doubleArr1
{
  private:
    double* arr;
    int sizes[3];
  public:
    inline doubleArr1(double*** in)
    {
      arr = **in;
      sizes[0] = 2; // not used
      sizes[1] = get_size(in);
      sizes[2] = 0;
    }
    inline double operator[](int idx){
      return arr[idx]; // Deref2(arr, sizes+2, idx);
    }
};

// class that can adopt array allocated by newArr2
class doubleArr2
{
  private:
    double* arr;
    int sizes[4];
  public:
    inline doubleArr2(double*** in)
    {
      arr = **in;
      sizes[0] = 2; // not used
      sizes[1] = get_size(in);
      sizes[2] = get_size(*in);
      sizes[3] = 0;
    }
    inline Deref1 operator[](int idx){
      return Deref1(arr, sizes+2, idx);
    }
    int getidx(int n1, int n2) const
    {
        int k = n1;
        k *= sizes[2]; k += n2;
        return k;
    }
    const double& get(int n1,int n2) const
      { return arr[getidx(n1,n2)]; }
    double& fetch(int n1,int n2)
      { return arr[getidx(n1,n2)]; }
    void set(int n1,int n2, double value)
      { arr[getidx(n1,n2)] = value; }
};

// class that can adopt array allocated by newArr3
class doubleArr3
{
  private:
    double* arr;
    int sizes[5];
  public:
    inline doubleArr3(double*** in)
    {
      arr = **in;
      sizes[0] = 3; // not used
      sizes[1] = get_size(in);
      sizes[2] = get_size(*in);
      sizes[3] = get_size(**in);
      sizes[4] = 0;
    }
    inline Deref2 operator[](int idx){
      return Deref2(arr, sizes+2, idx);
    }
    int getidx(int n1, int n2, int n3) const
    {
        int k = n1;
        k *= sizes[2]; k += n2;
        k *= sizes[3]; k += n3;
        return k;
    }
    const double& get(int n1,int n2,int n3) const
      { return arr[getidx(n1,n2,n3)]; }
    double& fetch(int n1,int n2,int n3)
      { return arr[getidx(n1,n2,n3)]; }
    void set(int n1,int n2,int n3, double value)
      { arr[getidx(n1,n2,n3)] = value; }
};

// class that can adopt array allocated by newArr4
class doubleArr4
{
  private:
    double* arr;
    int sizes[6];
  public:
    doubleArr4(double**** in)
    {
      arr = ***in;
      sizes[0] = 4; // not used
      sizes[1] = get_size(in);
      sizes[2] = get_size(*in);
      sizes[3] = get_size(**in);
      sizes[4] = get_size(***in);
      sizes[5] = 0;
    }
    inline Deref3 operator[](int idx){
      return Deref3(arr, sizes+2, idx);
    }
    int getidx(int n1, int n2, int n3, int n4) const
    {
        int k = n1;
        k *= sizes[2]; k += n2;
        k *= sizes[3]; k += n3;
        k *= sizes[4]; k += n4;
        return k;
    }
    const double& get(int n1,int n2,int n3,int n4) const
      { return arr[getidx(n1,n2,n3,n4)]; }
    double& fetch(int n1,int n2,int n3,int n4)
      { return arr[getidx(n1,n2,n3,n4)]; }
    void set(int n1,int n2,int n3,int n4, double value)
      { arr[getidx(n1,n2,n3,n4)] = value; }
};

#endif
