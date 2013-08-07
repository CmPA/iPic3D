#ifndef IPIC_ALLOC_H
#define IPIC_ALLOC_H
#include <cstddef> // for alignment stuff
#include "asserts.h" // for assert_le, assert_lt
//#include "errors.h" // for assert_le, assert_lt
#include "arraysfwd.h"
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

    Compiler options:
    -DCHECK_BOUNDS: check bounds when performing array access
      (major performance penalty).
    -DFLAT_ARRAYS: use calculated 1d subscript to dereference
      even for arr[i][j][k] notation.
    -DCHAINED_ARRAYS: use hierarchy of pointers to dereference
      even for arr.get(i,j,k) notation.

    By default, chained pointers are used for arr[i][j][k]
    notation (unless -DCHECK_BOUNDS is turned on, in which case
    we don't care about performance anyway), and calculated 1d
    subscript is used for arr.get(i,j,k) notation.

    An alternative would have been use boost arrays.  Use of our
    own array class allows flexibility for our choice of array
    implementation, including the possibility of using boost
    for the implementation, while avoiding boost as an external
    dependency.  On some systems, it may be preferable to use
    native arrays with hard-coded dimensions; this could suit us
    well, since all arrays are approximately the same size, but
    would require a recompile when changing the maximum array size.

    Rather than using these templates directly, the typedefs
    declared in "arraysfwd.h" should be used:

    * const_array_ref3_double = const_array_ref3<double>
    * arr3_double = array_ref3<double>
    * array3_double = array3<double>

    The point is that we do not want to hard-code the fact that
    we are using templates, and we may well wish to eliminate use
    of templates in the future.  (Alternatives are to use the
    preprocessor or to have separate implementations for each
    type (double, int, possibly float) if we go to use of mixed
    precision).  Support for templates is notoriously buggy in
    compilers, particularly when it comes to inheritance, and I
    in fact had to eliminate inheriting from the base_arr class
    and use the "protected" hack below in order to get this
    code to compile on the latest intel compiler (2013) and on
    g++ 4.0 (2005); g++ 4.2 (2007) compiled (but unfortunately,
    for my g++ 4.2, iPic3D suffered from stack frame corruption.)
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

// build chained pointer hierarchy for pre-existing bottom level
//
template <class type>
inline type **** newArray4(type * in, size_t sz1, size_t sz2, size_t sz3, size_t sz4)
{
  type****arr = newArray3<type*>(sz1,sz2,sz3);
  type**arr2 = **arr;
  type *ptr = in;
  size_t szarr2 = sz1*sz2*sz3;
  for(size_t i=0;i<szarr2;i++) {
    arr2[i] = ptr;
    ptr += sz4;
  }
  return arr;
}
template <class type>
inline type *** newArray3(type * in, size_t sz1, size_t sz2, size_t sz3)
{
  type***arr = newArray2<type*>(sz1,sz2);
  type**arr2 = *arr;
  type *ptr = in;
  size_t szarr2 = sz1*sz2;
  for(size_t i=0;i<szarr2;i++) {
    arr2[i] = ptr;
    ptr += sz3;
  }
  return arr;
}
template <class type>
inline type ** newArray2(type * in, size_t sz1, size_t sz2)
{
  type**arr = newArray2<type*>(sz1);
  type**arr2 = arr;
  type *ptr = in;
  size_t szarr2 = sz1;
  for(size_t i=0;i<szarr2;i++) {
    arr2[i] = ptr;
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
  
namespace iPic3D
{
  // underlying 1-dimensional array class for arrays
  
  template <class type>
  class base_arr
  {
    private:
      size_t size;
    protected:
      type* const __restrict__ arr;
      type* get_arr()const{return arr;}
    public:
      base_arr(size_t s) : size(s), arr(AlignedAlloc(type, s)) {}
      base_arr(type* in, size_t s) : size(s), arr(in) {}
      ~base_arr(){}
      void free() { AlignedFree(arr); }
      void setall(type val){
        for(size_t i=0;i<size;i++) arr[i]=val;
      }
      //type* fetch_arr(){return arr;}
  };
  
  // classes to dereference arrays.
  //
  // ArrayGetN is essentially a dumbed-down version of ArrN with
  // an index shift applied to the underlying array.  The purpose
  // of ArrayGetN is to allow elements of multidimensional arrays
  // to be accessed with a calculated one-dimensional index while
  // using chained operator[] syntax (e.g. myarr[i][j]), i.e. the
  // same syntax as is used for native or nested arrays.  This
  // implementation is likely to be slow unless optimization is
  // turned on, allowing the compiler to figure out that the whole
  // chain of calls to the operator[] methods and to the ArrayGetN
  // constructors reduces to computing a one-dimensional subscript
  // used to access a one-dimensional array.
  //
  template <class type>
  class ArrayGet1
  {
    type* const __restrict__ arr;
    const size_t S1;
    const size_t shift;
   public:
    inline ArrayGet1(type*const arr_, size_t k, size_t s1) :
      arr(arr_), shift(k), S1(s1)
    {}
    inline type& operator[](size_t n1){
      check_bounds(n1, S1);
      ALIGNED(arr);
      return arr[shift+n1];
    }
  };
  
  template <class type>
  class ArrayGet2
  {
    type* const __restrict__ arr;
    const size_t shift;
    const size_t S2, S1;
   public:
    inline ArrayGet2(type*const arr_, size_t k, size_t s2, size_t s1) :
      arr(arr_), shift(k), S2(s2), S1(s1)
    {}
    inline ArrayGet1<type> operator[](size_t n2){
      check_bounds(n2,S2);
      return ArrayGet1<type>(arr, (shift+n2)*S1, S1);
    }
  };
  
  template <class type>
  class ArrayGet3
  {
    type* const __restrict__ arr;
    const size_t shift;
    const size_t S3, S2, S1;
   public:
    inline ArrayGet3(type*const arr_, size_t k, size_t s3, size_t s2, size_t s1) :
      arr(arr_), shift(k), S3(s3), S2(s2), S1(s1)
    {}
    inline ArrayGet2<type> operator[](size_t n3){
      check_bounds(n3, S3);
      return ArrayGet2<type>(arr, (shift+n3)*S2, S2, S1);
    }
  };
  
  // const versions
  
  template <class type>
  class const_array_get1
  {
    type* const __restrict__ arr;
    const size_t S1;
    const size_t shift;
   public:
    inline const_array_get1(type*const arr_, size_t k, size_t s1) :
      arr(arr_), shift(k), S1(s1)
    {}
    inline const type& operator[](size_t n1)const{
      check_bounds(n1, S1);
      ALIGNED(arr);
      return arr[shift+n1];
    }
  };
  
  template <class type>
  class const_array_get2
  {
    type* const __restrict__ arr;
    const size_t shift;
    const size_t S2, S1;
   public:
    inline const_array_get2(type*const arr_, size_t k, size_t s2, size_t s1) :
      arr(arr_), shift(k), S2(s2), S1(s1)
    {}
    inline const const_array_get1<type> operator[](size_t n2)const{
      check_bounds(n2,S2);
      return const_array_get1<type>(arr, (shift+n2)*S1, S1);
    }
  };
  
  template <class type>
  class const_array_get3
  {
    type* const __restrict__ arr;
    const size_t shift;
    const size_t S3, S2, S1;
   public:
    const_array_get3(type*const arr_, size_t k, size_t s3, size_t s2, size_t s1) :
      arr(arr_), shift(k), S3(s3), S2(s2), S1(s1)
    {}
    inline const const_array_get2<type> operator[](size_t n3)const{
      check_bounds(n3, S3);
      return const_array_get2<type>(arr, (shift+n3)*S2, S2, S1);
    }
  };
  
  // ArrN corresponds to multi_array_ref in the boost library.
  //
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
  //   array_ref2 arr<int>(16, 16);
  //   arr[1][2] = 5;
  //   arr.free();
  // }
  // Using ArrN to adopt an array allocated by newArrN
  // {
  //   int** array = newArray2<int>(16,16)
  //   array_ref2 arr(array,16,16); // adopt array
  //   arr[1][2] = 5;
  //   assert_eq(arr[1][2],array[1][2]);
  //   // arr.free(); // should not do both this and next line.
  //   delArray2<int>(array);
  // }
  //
  // proposed improvements:
  // - allow shifting of the base:
  //   - need "double shift" in each class
  //   - need to implement "arr3.set_bases(b1,b2,b3);"
  //     which calculates "shift".
  //   - need "const size_t b1, b2, b3;" for beginning indices
  //     to allow bounds checking.  Should not incur run-time
  //     penalty, but it so then condition on CHECK_BOUNDS.
  // - methods that use parallel arithmetic for omp and vectorized code
  
  template <class type>
  class array_ref1
  {
    private: // data
      const size_t S1;
      type* const __restrict__ arr;
    public:
      ~array_ref1() { }
      void free() { AlignedFree(arr); }
      array_ref1(size_t s1) :
        S1(s1),
        arr(AlignedAlloc(type, s1))
      { }
      array_ref1(type* in,
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
  class array_ref2
  {
    private: // data
      const size_t S2,S1;
      type* const __restrict__ arr;
    public:
      ~array_ref2(){}
      void free() { AlignedFree(arr); }
      array_ref2(size_t s2, size_t s1) :
        S2(s2), S1(s1),
        arr(AlignedAlloc(type, s2*s1))
      {
      }
      array_ref2(type*const* in,
        size_t s2, size_t s1) :
        S2(s2), S1(s1),
        arr(*in)
      { }
      // dereference via calculated index
      inline ArrayGet1<type> operator[](size_t n2){
        check_bounds(n2, S2);
        return ArrayGet1<type>(arr, n2*S1, S1);
      }
      inline size_t getidx(size_t n2, size_t n1) const
      {
        check_bounds(n2, S2);
        check_bounds(n1, S1);
        return n2*S1+n1;
      }
      type& fetch(size_t n2, size_t n1) const
        { ALIGNED(arr); return arr[n1+S1*n2]; }
      // better to use accessors that distinguish read from write:
      const type& get(size_t n2,size_t n1) const
        { ALIGNED(arr); return arr[getidx(n2,n1)]; }
      void set(size_t n2,size_t n1, type value)
        { ALIGNED(arr); arr[getidx(n2,n1)] = value; }
      //inline array_ref1<type>fetch_Arr1(){ return array_ref1<type>(arr, S1*S2); }
  };
  
  template <class type>
  class const_array_ref3 // : public base_arr<type>
  {
      //using base_arr<type>::get_arr;
      //using base_arr<type>::arr;
    protected: // data
      size_t size;
      const size_t S3,S2,S1;
      type* const __restrict__ arr;
      type*const*const*const arr3;
    public:
      ~const_array_ref3(){}
      const_array_ref3(size_t s3, size_t s2, size_t s1) :
        size(s3*s2*s1), arr(AlignedAlloc(type, size)),
        //base_arr<type>(s3*s2*s1),
        S3(s3), S2(s2), S1(s1),
        arr3(newArray3<type>(arr,s3,s2,s1))
      { }
      const_array_ref3(type*const*const* in,
        size_t s3, size_t s2, size_t s1) :
        size(s3*s2*s1), arr(**in),
        //base_arr<type>(**in, s3*s2*s1),
        S3(s3), S2(s2), S1(s1),
        arr3(in)
      { }
    #if defined(FLAT_ARRAYS) || defined(CHECK_BOUNDS)
      const const_array_get2<type> operator[](size_t n3)const{
        check_bounds(n3, S3);
        return const_array_get2<type>(arr, n3*S2, S2, S1);
      }
    #else
      // this causes operator[] to dereference via chained pointer
      operator type***(){ return (type***) arr3; }
    #endif
      void check_idx_bounds(size_t n3, size_t n2, size_t n1) const
      {
        check_bounds(n3, S3);
        check_bounds(n2, S2);
        check_bounds(n1, S1);
      }
      inline size_t getidx(size_t n3, size_t n2, size_t n1) const
        { check_idx_bounds(n3,n2,n1); return (n3*S2+n2)*S1+n1; }
    #ifdef CHAINED_ARRAYS
      const type& get(size_t n3,size_t n2,size_t n1) const
        { check_idx_bounds(n3,n2,n1); return arr3[n3][n2][n1]; }
    protected: // hack: not in const_array_ref3 due to icpc compile error
      type& fetch(size_t n3,size_t n2,size_t n1) const
        { check_idx_bounds(n3,n2,n1); return arr3[n3][n2][n1]; }
      void set(size_t n3,size_t n2,size_t n1, type value)
        { check_idx_bounds(n3,n2,n1); arr3[n3][n2][n1] = value; }
    #else
      const type& get(size_t n3,size_t n2,size_t n1) const
        { ALIGNED(arr); return arr[getidx(n3,n2,n1)]; }
    protected: // hack: not in const_array_ref3 due to icpc compile error
      type& fetch(size_t n3,size_t n2,size_t n1) const
        { ALIGNED(arr); return arr[getidx(n3,n2,n1)]; }
      void set(size_t n3,size_t n2,size_t n1, type value)
        { ALIGNED(arr); arr[getidx(n3,n2,n1)] = value; }
    #endif
  };
  
  template <class type>
  class array_ref3 : public const_array_ref3<type>
  {
      //using base_arr<type>::arr;
      //using base_arr<type>::get_arr;
      using const_array_ref3<type>::size;
      using const_array_ref3<type>::arr;
      using const_array_ref3<type>::S3;
      using const_array_ref3<type>::S2;
      using const_array_ref3<type>::S1;
      using const_array_ref3<type>::arr3;
      using const_array_ref3<type>::getidx;
    public:
      ~array_ref3(){}
      array_ref3(size_t s3, size_t s2, size_t s1) :
        const_array_ref3<type>(s3,s2,s1)
      { }
      array_ref3(type*const*const* in,
        size_t s3, size_t s2, size_t s1) :
        const_array_ref3<type>(in,s3,s2,s1)
      { }
      void free(){ delArray3<type>((type***)arr3); }
    #if defined(FLAT_ARRAYS) || defined(CHECK_BOUNDS)
      inline ArrayGet2<type> operator[](size_t n3){
        check_bounds(n3, S3);
        return ArrayGet2<type>(arr, n3*S2, S2, S1);
      }
    #else
      // this causes operator[] to dereference via chained pointer
      operator type***(){ return (type***) arr3; }
    #endif
      type& fetch(size_t n3,size_t n2,size_t n1) const
        { return const_array_ref3<type>::fetch(n3,n2,n1); }
      void set(size_t n3,size_t n2,size_t n1, type value)
        { const_array_ref3<type>::set(n3,n2,n1, value); }
      void setall(type val){
        for(size_t i=0;i<size;i++) arr[i]=val;
      }
      type*** fetch_arr3(){ return (type***) arr3; }
  };
  
  // inheriting from base_arr<type> causes problems in g++ 4.0 (2005).
  template <class type>
  class const_array_ref4 //: public base_arr<type>
  {
      //using base_arr<type>::get_arr;
    protected: // data
      size_t size;
      const size_t S4,S3,S2,S1;
      type* const __restrict__ arr;
      type*const*const*const*const arr4;
    public:
      ~const_array_ref4(){}
      const_array_ref4(size_t s4, size_t s3, size_t s2, size_t s1) :
        size(s4*s3*s2*s1), arr(AlignedAlloc(type, size)),
        //base_arr<type>(s4*s3*s2*s1),
        S4(s4), S3(s3), S2(s2), S1(s1),
        arr4(newArray4<type>(arr,s4,s3,s2,s1))
      { }
      const_array_ref4(type*const*const*const* in,
        size_t s4, size_t s3, size_t s2, size_t s1) :
        size(s4*s3*s2*s1), arr(***in),
        //base_arr<type>(***in, s4*s3*s2*s1),
        S4(s4), S3(s3), S2(s2), S1(s1),
        arr4(in)
      { }
    #if defined(FLAT_ARRAYS) || defined(CHECK_BOUNDS)
      const const_array_get3<type> operator[](size_t n4)const{
        check_bounds(n4, S4);
        return const_array_get3<type>(arr, n4*S3, S3, S2, S1);
      }
    #else
      // this causes operator[] to dereference via chained pointer
      operator type****(){ return (type****) arr4; }
    #endif
      void check_idx_bounds(size_t n4, size_t n3, size_t n2, size_t n1) const
      {
        check_bounds(n4, S4);
        check_bounds(n3, S3);
        check_bounds(n2, S2);
        check_bounds(n1, S1);
      }
      inline size_t getidx(size_t n4, size_t n3, size_t n2, size_t n1) const
        { check_idx_bounds(n4,n3,n2,n1); return ((n4*S3+n3)*S2+n2)*S1+n1; }
    #ifdef CHAINED_ARRAYS
      const type& get(size_t n4,size_t n3,size_t n2,size_t n1) const
        { ALIGNED(arr); return arr[getidx(n4,n3,n2,n1)]; }
    protected: // hack: not in const_array_ref4 due to icpc compile error
      type& fetch(size_t n4,size_t n3,size_t n2,size_t n1) const
        { ALIGNED(arr); return arr[getidx(n4,n3,n2,n1)]; }
      void set(size_t n4,size_t n3,size_t n2,size_t n1, type value)
        { ALIGNED(arr); arr[getidx(n4,n3,n2,n1)] = value; }
    #else
      const type& get(size_t n4,size_t n3,size_t n2,size_t n1) const
        { check_idx_bounds(n4,n3,n2,n1); return arr4[n4][n3][n2][n1]; }
    protected: // hack: not in const_array_ref4 due to icpc compile error
      type& fetch(size_t n4,size_t n3,size_t n2,size_t n1) const
        { check_idx_bounds(n4,n3,n2,n1); return arr4[n4][n3][n2][n1]; }
      void set(size_t n4,size_t n3,size_t n2,size_t n1, type value)
        { check_idx_bounds(n4,n3,n2,n1); arr4[n4][n3][n2][n1] = value; }
    #endif
  };
  
  template <class type>
  class array_ref4 : public const_array_ref4<type>
  {
      //using base_arr<type>::get_arr;
      using const_array_ref4<type>::arr;
      using const_array_ref4<type>::S4;
      using const_array_ref4<type>::S3;
      using const_array_ref4<type>::S2;
      using const_array_ref4<type>::S1;
      using const_array_ref4<type>::arr4;
      using const_array_ref4<type>::getidx;
    public:
      ~array_ref4(){}
      array_ref4(size_t s4, size_t s3, size_t s2, size_t s1) :
        const_array_ref4<type>(s4,s3,s2,s1)
      { }
      array_ref4(type*const*const*const* in,
        size_t s4, size_t s3, size_t s2, size_t s1) :
        const_array_ref4<type>(in,s4,s3,s2,s1)
      { }
    #if defined(FLAT_ARRAYS) || defined(CHECK_BOUNDS)
      inline ArrayGet3<type> operator[](size_t n4){
        check_bounds(n4, S4);
        return ArrayGet3<type>(arr, n4*S3, S3, S2, S1);
      }
    #else
      operator type****(){ return (type****) arr4; }
    #endif
      type& fetch(size_t n4,size_t n3,size_t n2,size_t n1) const
        { return const_array_ref4<type>::fetch(n4,n3,n2,n1); }
      void set(size_t n4,size_t n3,size_t n2,size_t n1, type value)
        { const_array_ref4<type>::set(n4,n3,n2,n1, value); }
      void free(){ delArray4<type>((type****)arr4); }
      type**** fetch_arr4(){ return (type****) arr4; }
      //bool verify_dims(size_t s4, size_t s3, size_t s2, size_t s1){
      //  if(s4==S4 && s3==S3 && s2==S2 && s1==S1) return true;
      //  Wprintf("%d==%d && %d==%d && %d==%d && %d==%d failed",
      //     s4, S4, s3, S3, s2, S2, s1, S1);
      //  return false;
      //}
  };
  
  // Versions of array classes which automatically free memory
  // (corresponding to multi_array in the boost library).
  //
  // Note that the nonempty destructor kills performance
  // unless compiling with -fno-exceptions
  
  template <class type>
  struct array1 : public array_ref1<type>
  {
      ~array1(){array_ref1<type>::free();}
      array1(size_t s1) : array_ref1<type>(s1) { }
  };
  
  template <class type>
  struct array2 : public array_ref2<type>
  {
      ~array2(){array_ref2<type>::free();}
      array2(size_t s2, size_t s1) : array_ref2<type>(s2,s1) { }
  };
  
  template <class type>
  struct array3 : public array_ref3<type>
  {
      ~array3(){array_ref3<type>::free();}
      array3(size_t s3, size_t s2, size_t s1) : array_ref3<type>(s3,s2,s1) { }
  };
  
  template <class type>
  struct array4 : public array_ref4<type>
  {
      ~array4(){array_ref4<type>::free();}
      array4(size_t s4, size_t s3, size_t s2, size_t s1)
        : array_ref4<type>(s4,s3,s2,s1) { }
  };

}

#define newArr4(type,sz1,sz2,sz3,sz4) newArray4<type>((sz1),(sz2),(sz3),(sz4))
#define newArr3(type,sz1,sz2,sz3) newArray3<type>((sz1),(sz2),(sz3))
#define newArr2(type,sz1,sz2) newArray2<type>((sz1),(sz2))
/*** end Array classes with flexible dimensions ***/
#endif
