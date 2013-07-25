#ifndef IPIC_ALLOC_H
#define IPIC_ALLOC_H
#include <cstddef> // for alignment stuff
#include "ipicdefs.h" // for CHECK_BOUNDS
#include "asserts.h" // for assert_le, assert_lt
//#include "arrays.h" // fixed-dimension arrays

/*
    Array classes developed by
      Alec Johnson,
    consolidating arrays developed by 
      Reger Ferrer, Vicenç Beltran, and Florentino Sainz
    and earlier arrays defined by
      Jorge Amaya and Stefano Markidis.

    For examples of use of this class,
    see test_arrays.cpp
*/
#define ALIGNMENT (64)
#ifdef __INTEL_COMPILER
    #define ALIGNED(X) __assume_aligned(X, ALIGNMENT)
    #define AlignedAlloc(T, NUM) \
        (T *const __restrict__)(_mm_malloc(sizeof(T)*NUM, ALIGNMENT))
    #define AlignedFree(S) (_mm_free(S))
#else
    #define ALIGNED(X)
    #define AlignedFree(S) (delete[] S)
    #define AlignedAlloc(T, NUM) (new T[NUM]) 
#endif

// Compile with -DCHECK_BOUNDS to turn on bounds checking.
//#define CHECK_BOUNDS
#ifdef CHECK_BOUNDS
  #define check_bounds(n,S) {assert_le(0, n); assert_lt(n, S);}
#else
  #define check_bounds(n,S)
#endif

/*** begin Array classes with flexible dimensions ***/

// methods to allocate arrays.
// These are a succinct equivalent of Jorge's earler methods,
// except for the use of AlignedAlloc in place of new.
//
template < class type >
inline type * newArray1(size_t sz1)
{
  type *arr = AlignedAlloc(type, sz1); // new type [sz1];
  return arr;
}
template < class type >
inline type ** newArray2(size_t sz1, size_t sz2)
{
  type **arr = AlignedAlloc(type*, sz1); // new type *[sz1];
  type *ptr = newArray1<type>(sz1*sz2);
  for (size_t i = 0; i < sz1; i++)
  {
    arr[i] = ptr;
    ptr += sz2;
  }
  return arr;
}
template < class type >
inline type *** newArray3(size_t sz1, size_t sz2, size_t sz3)
{
  type ***arr = AlignedAlloc(type**, sz1); // new type **[sz1];
  type **ptr = newArray2<type>(sz1*sz2, sz3);
  for (size_t i = 0; i < sz1; i++)
  {
    arr[i] = ptr;
    ptr += sz2;
  }
  return arr;
}
template <class type>
inline type **** newArray4(size_t sz1, size_t sz2, size_t sz3, size_t sz4)
{
  type ****arr = AlignedAlloc(type***, sz1); //(new type ***[sz1]);
  type ***ptr = newArray3<type>(sz1*sz2, sz3, sz4);
  for (size_t i = 0; i < sz1; i++) {
    arr[i] = ptr;
    ptr += sz2;
  }
  return arr;
}

// methods to deallocate arrays
//
template < class type > inline void delArray1(type * arr)
{ AlignedFree(arr); }
template < class type > inline void delArray2(type ** arr)
{ delArray1(arr[0]); AlignedFree(arr); }
template < class type > inline void delArray3(type *** arr)
{ delArray2(arr[0]); AlignedFree(arr); }
template < class type > inline void delArray4(type **** arr)
{ delArray3(arr[0]); AlignedFree(arr); }
//
// versions with dummy dimensions (for backwards compatibility)
//
template <class type> inline void delArr1(type * arr)
{ delArray1(arr); }
template <class type> inline void delArr2(type ** arr, size_t sz1)
{ delArray2(arr); }
template <class type> inline void delArr3(type *** arr, size_t sz1, size_t sz2)
{ delArray3(arr); }
template <class type> inline void delArr4(type **** arr,
  size_t sz1, size_t sz2, size_t sz3)
{ delArray3(arr); }

// classes to dereference arrays.
//
// ArrayRefN is essentially a dumbed-down version of ArrN with
// an index shift applied to the underlying array.  The purpose
// of ArrayRefN is to allow elements of multidimensional arrays
// to be accessed with a calculated one-dimensional index while
// using chained operator[] syntax (e.g. myarr[i][j]), i.e. the
// same syntax as is used for native or nested arrays.  This
// implementation is likely to be slow unless optimization is
// turned on, allowing the compiler to figure out that the whole
// chain of calls to the operator[] methods and to the ArrayRefN
// constructors reduces to computing a one-dimensional subscript
// used to access a one-dimensional array.
//
template <class type>
class ArrayRef1
{
  type* const __restrict__ arr;
  const size_t S1;
  const size_t shift;
 public:
  inline ArrayRef1(type*const arr_, size_t k, size_t s1) :
    arr(arr_), shift(k), S1(s1)
  {}
  inline type& operator[](size_t n1){
    check_bounds(n1, S1);
    ALIGNED(arr);
    return arr[shift+n1];
  }
};

template <class type>
class ArrayRef2
{
  type* const __restrict__ arr;
  const size_t shift;
  const size_t S2, S1;
 public:
  inline ArrayRef2(type*const arr_, size_t k, size_t s2, size_t s1) :
    arr(arr_), shift(k), S2(s2), S1(s1)
  {}
  inline ArrayRef1<type> operator[](size_t n2){
    check_bounds(n2,S2);
    return ArrayRef1<type>(arr, (shift+n2)*S1, S1);
  }
};

template <class type>
class ArrayRef3
{
  type* const __restrict__ arr;
  const size_t shift;
  const size_t S3, S2, S1;
 public:
  inline ArrayRef3(type*const arr_, size_t k, size_t s3, size_t s2, size_t s1) :
    arr(arr_), shift(k), S3(s3), S2(s2), S1(s1)
  {}
  inline ArrayRef2<type> operator[](size_t n3){
    check_bounds(n3, S3);
    return ArrayRef2<type>(arr, (shift+n3)*S2, S2, S1);
  }
};

// ArrN can adopt an array allocated by newArrN
//
// The purpose of these classes is to provide more efficient
// and more regulated access to array elements.  The idea is to
// maintain backward compatibility while allowing us to move
// toward a proper array abstraction.
//
// The user of ArrN is responsible for memory management.
// The ArrayN classes are the version of this class
// with automatic deallocation.
//
// Examples:
//
// Using constructor to create array:
// {
//   Arr2 arr<int>(16, 16);
//   arr[1][2] = 5;
//   arr.free();
// }
// Using ArrN to adopt an array allocated by newArrN
// {
//   int** array = newArray2<int>(16,16)
//   Arr2 arr(array,16,16); // adopt array
//   arr[1][2] = 5;
//   assert_eq(arr[1][2],array[1][2]);
//   // arr.free(); // should not do both this and next line.
//   delArray2<int>(array);
// }
//
// proposed improvements:
// - methods that use parallel arithmetic for omp and vectorized code

template <class type>
class Arr1
{
  private: // data
    type* const __restrict__ arr;
    const size_t S1;
  public:
    ~Arr1() { }
    void free() { AlignedFree(arr); }
    Arr1(size_t s1) :
      S1(s1),
      arr(AlignedAlloc(type, s1))
    { }
    Arr1(type* in,
      size_t s1) :
      S1(s1),
      arr(in)
    { }
    inline type& operator[](size_t n1){
      check_bounds(n1, S1);
      ALIGNED(arr);
      return arr[n1];
    }
    inline size_t getidx(size_t n1) const
    {
      check_bounds(n1, S1);
      return n1;
    }
    const type& get(size_t n1) const
      { ALIGNED(arr); return arr[getidx(n1)]; }
    type& fetch(size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[getidx(n1)]; }
    void set(size_t n1, type value)
      { ALIGNED(arr); arr[getidx(n1)] = value; }
};

template <class type>
class Arr2
{
  private: // data
    const size_t S2,S1;
    type* const __restrict__ arr;
  public:
    ~Arr2(){}
    void free() { AlignedFree(arr); }
    Arr2(size_t s2, size_t s1) :
      S2(s2), S1(s1),
      arr(AlignedAlloc(type, s2*s1))
    {
    }
    Arr2(type*const* in,
      size_t s2, size_t s1) :
      S2(s2), S1(s1),
      arr(*in)
    { }
    // for backwards compatibility support bracket notation
    inline ArrayRef1<type> operator[](size_t n2){
      check_bounds(n2, S2);
      return ArrayRef1<type>(arr, n2*S1, S1);
    }
    inline size_t getidx(size_t n2, size_t n1) const
    {
      check_bounds(n2, S2);
      check_bounds(n1, S1);
      return n2*S1+n1;
    }
    // I prefer "fetch" over operator() to hilight read/write access
    //type& operator()(size_t n2, size_t n1) const
    //  { ALIGNED(arr); return arr[n1+S1*n2]; }
    type& fetch(size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[getidx(n2,n1)]; }
    // better to use accessors that distinguish read from write:
    const type& get(size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[getidx(n2,n1)]; }
    void set(size_t n2,size_t n1, type value)
      { ALIGNED(arr); arr[getidx(n2,n1)] = value; }
};

template <class type>
class Arr3
{
  private: // data
    type* const __restrict__ arr;
    const size_t S3,S2,S1;
  public:
    ~Arr3(){}
    void free() { AlignedFree(arr); }
    Arr3(size_t s3, size_t s2, size_t s1) :
      S3(s3), S2(s2), S1(s1),
      arr(AlignedAlloc(type, s3*s2*s1))
    { }
    Arr3(type*const*const* in,
      size_t s3, size_t s2, size_t s1) :
      S3(s3), S2(s2), S1(s1),
      arr(**in)
    { }
    inline ArrayRef2<type> operator[](size_t n3){
      check_bounds(n3, S3);
      return ArrayRef2<type>(arr, n3*S2, S2, S1);
    }
    type* get_arr(){return arr;}
    inline size_t getidx(size_t n3, size_t n2, size_t n1) const
    {
      check_bounds(n3, S3);
      check_bounds(n2, S2);
      check_bounds(n1, S1);
      return (n3*S2+n2)*S1+n1;
    }
    //type& operator()(size_t n3, size_t n2, size_t n1) const
    //{ ALIGNED(arr); return arr[getidx(n3,n2,n1)]; }
    type& fetch(size_t n3,size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[getidx(n3,n2,n1)]; }
    const type& get(size_t n3,size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[getidx(n3,n2,n1)]; }
    void set(size_t n3,size_t n2,size_t n1, type value)
      { ALIGNED(arr); arr[getidx(n3,n2,n1)] = value; }
};

template <class type>
class Arr4
{
  private: // data
    const size_t S4,S3,S2,S1;
    type* const __restrict__ arr;
  public:
    ~Arr4(){} // nonempty destructor would kill performance
    void free() { AlignedFree(arr); }
    Arr4(size_t s4, size_t s3, size_t s2, size_t s1) :
      arr(AlignedAlloc(type, s4*s3*s2*s1)),
      S4(s4), S3(s3), S2(s2), S1(s1)
    { }
    Arr4(type*const*const*const* in,
      size_t s4, size_t s3, size_t s2, size_t s1) :
      S4(s4), S3(s3), S2(s2), S1(s1),
      arr(***in)
    { }
    inline ArrayRef3<type> operator[](size_t n4){
      check_bounds(n4, S4);
      return ArrayRef3<type>(arr, n4*S3, S3, S2, S1);
    }
    inline size_t getidx(size_t n4, size_t n3, size_t n2, size_t n1) const
    {
      check_bounds(n4, S4);
      check_bounds(n3, S3);
      check_bounds(n2, S2);
      check_bounds(n1, S1);
      return ((n4*S3+n3)*S2+n2)*S1+n1;
    }
    const type& get(size_t n4,size_t n3,size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[getidx(n4,n3,n2,n1)]; }
    type& fetch(size_t n4,size_t n3,size_t n2,size_t n1) const
      { ALIGNED(arr); return arr[getidx(n4,n3,n2,n1)]; }
    void set(size_t n4,size_t n3,size_t n2,size_t n1, type value)
      { ALIGNED(arr); arr[getidx(n4,n3,n2,n1)] = value; }
};

// Versions of array classes which automatically free memory.
//
// Note that the nonempty destructor kills performance
// unless compiling with -fno-exceptions

template <class type>
struct Array1 : public Arr1<type>
{
    ~Array1(){Arr1<type>::free();}
    Array1(size_t s1) : Arr1<type>(s1) { }
};

template <class type>
struct Array2 : public Arr2<type>
{
    ~Array2(){Arr2<type>::free();}
    Array2(size_t s2, size_t s1) : Arr2<type>(s2,s1) { }
};

template <class type>
struct Array3 : public Arr3<type>
{
    ~Array3(){Arr3<type>::free();}
    Arr3<type>& fast_accessor() { return *(Arr3<type>*)this; }
    Array3(size_t s3, size_t s2, size_t s1) : Arr3<type>(s3,s2,s1) { }
};

template <class type>
struct Array4 : public Arr4<type>
{
    ~Array4(){Arr4<type>::free();}
    Array4(size_t s4, size_t s3, size_t s2, size_t s1)
      : Arr4<type>(s4,s3,s2,s1) { }
};

// These aliases are defined for the following flexibilization purposes:
// - to avoid filling the code with template brackets
//   (i.e., to minimize explicitly template-dependent code).
// - so that they can be redefined according to the user's
//   preferred array implementation.
//
typedef Arr1<int> intArr1;
typedef Arr2<int> intArr2;
typedef Arr3<int> intArr3;
typedef Arr4<int> intArr4;
typedef Arr1<double> doubleArr1;
typedef Arr2<double> doubleArr2;
typedef Arr3<double> doubleArr3;
typedef Arr4<double> doubleArr4;
//
#define newArr4(type,sz1,sz2,sz3,sz4) newArray4<type>((sz1),(sz2),(sz3),(sz4))
#define newArr3(type,sz1,sz2,sz3) newArray3<type>((sz1),(sz2),(sz3))
#define newArr2(type,sz1,sz2) newArray2<type>((sz1),(sz2))
/*** end Array classes with flexible dimensions ***/
#endif
