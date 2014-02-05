/*
   Reger Ferrer
   Vicenç Beltran
   Alec Johnson
*/

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "stopwatch.h"
#include "arrays.h"
#include "Alloc.h"
#include "asserts.h"
#include "debug.h"

/**** begin Jorge Amaya's array allocation methods ****/

/*! The allocator for 4D array */
template < class type > type **** newArray4_Amaya(int sz1, int sz2, int sz3, int sz4) {

  type ****all_x;
  type ***all_y;
  type **all_z;
  type *all_r;

  all_x = new type ***[sz1];
  all_y = new type **[sz1 * sz2];
  all_z = new type *[sz1 * sz2 * sz3];
  all_r = new type[sz1 * sz2 * sz3 * sz4];

  type ****result = all_x;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
    for (int j = 0; j < sz2; j++, all_z += sz3) {
      result[i][j] = all_z;
      for (int k = 0; k < sz3; k++, all_r += sz4) {
        result[i][j][k] = all_r;
      }
    }
  }

  return result;
}

/*! Deallocator for 4D arrays */
template < class type > void delArr4_Amaya(type **** arr, int dummyx, int dummyy, int dummyz) {
  delete[]arr[0][0][0];
  delete[]arr[0][0];
  delete[]arr[0];
  delete[]arr;
}

/*! The allocator for 3D array */
template < class type > type *** newArray3_Amaya(int sz1, int sz2, int sz3) {

  type ***all_x;
  type **all_y;
  type *all_z;

  all_x = new type **[sz1];
  all_y = new type *[sz1 * sz2];
  all_z = new type[sz1 * sz2 * sz3];

  type ***result = all_x;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
    for (int j = 0; j < sz2; j++, all_z += sz3) {
      result[i][j] = all_z;
    }
  }

  return result;

}

/*! Deallocator for 3D arrays */
template < class type > void delArr3_Amaya(type *** arr, int dummyx, int dummyy) {
  delete[]arr[0][0];
  delete[]arr[0];
  delete[]arr;
}

/*! The allocator for 2D array */
template < class type > type ** newArr2_Amaya(int sz1, int sz2) {

  type **all_x;
  type *all_y;

  all_x = new type *[sz1];
  all_y = new type[sz1 * sz2];

  type **result = all_x;

  for (int i = 0; i < sz1; i++, all_y += sz2) {
    result[i] = all_y;
  }

  return result;

}

/*! Deallocator for 2D arrays */
template < class type > void delArr2_Amaya(type ** arr, int dummyx) {
  delete[]arr[0];
  delete[]arr;
}

#define newArr4_Amaya(type,sz1,sz2,sz3,sz4) newArray4_Amaya<type>((sz1),(sz2),(sz3),(sz4))
#define newArr3_Amaya(type,sz1,sz2,sz3) newArray3_Amaya<type>((sz1),(sz2),(sz3))
#define newArr2_Amaya(type,sz1,sz2) newArray2_Amaya<type>((sz1),(sz2))

/**** end Jorge Amaya's array allocation methods ****/

/****** begin (i,j) arrays from Reger Ferrer and Vicenç Beltran ******/

template <class type>
class Rank1
{
    const size_t S1;
    type  * __restrict__ const  arr;

public:

    Rank1(size_t s1) : S1(s1), arr(AlignedAlloc(type, s1)) {}

    //Rank1( const Rank1& other ) : S1( other.S1 ), arr( other.arr ) {}

    type& operator()(size_t n1) const
    {
       ALIGNED(arr);
       return arr[n1];
    }

    size_t dim1() const { return S1; }

    ~Rank1() { };
};

template <class type>
class Rank2
{
    const size_t  S1, S2;
    type * __restrict__ const arr;


public:
    Rank2(size_t s1, size_t s2) : S1(s1), S2(s2), arr(AlignedAlloc(type, s1*s2)) {}

    //Rank2( const Rank2& other ) : S1( other.S1 ), S2( other.S2 ), arr( other.arr ) {}

    type& operator()(size_t n1, size_t n2) const
    {
       check_bounds(n1,S1);
       check_bounds(n2,S2);
       ALIGNED(arr);
       return arr[n2+S2*n1];
    }
    type& fetch(size_t n1,size_t n2) const
    {
       check_bounds(n1,S1);
       check_bounds(n2,S2);
       ALIGNED(arr);
       return arr[n2+S2*n1];
    }

    size_t dim1() const { return S1; }
    size_t dim2() const { return S2; }
    
    void free() {
       AlignedFree(arr);
    }

    ~Rank2() { };
};
    
template <class type>
class Rank3
{
    const size_t S1, S2, S3;
    type *    const __restrict__ arr;


public:

    Rank3(size_t s1, size_t s2, size_t s3) : S1(s1), S2(s2), S3(s3),
    arr(AlignedAlloc(type, s1*s2*s3)) {}

    //Rank3( const Rank3& other ) : S1( other.S1 ), S2( other.S2 ), S3( other.S3 ),
    //arr( other.arr ) {}

    type& operator()(size_t n1, size_t n2, size_t n3) const
    {
       check_bounds(n1,S1);
       check_bounds(n2,S2);
       check_bounds(n3,S3);
       ALIGNED(arr);
       return arr[n3+S3*(n2+S2*n1)];
    }
    type& fetch(size_t n1, size_t n2, size_t n3) const
    {
       check_bounds(n1,S1);
       check_bounds(n2,S2);
       check_bounds(n3,S3);
       ALIGNED(arr);
       return arr[n3+S3*(n2+S2*n1)];
    }
    const type& get(size_t n1, size_t n2, size_t n3) const
    {
       check_bounds(n1,S1);
       check_bounds(n2,S2);
       check_bounds(n3,S3);
       ALIGNED(arr);
       return arr[n3+S3*(n2+S2*n1)];
    }

    ~Rank3() { }

    size_t dim1() const { return S1; }
    size_t dim2() const { return S2; }
    size_t dim3() const { return S3; }

    void free() {
       AlignedFree(arr);
    }
};

template <class type>
class Rank4
{
    const size_t S1, S2, S3, S4;
    type* __restrict__ const arr;

public:

    Rank4(size_t s1, size_t s2, size_t s3, size_t s4) : S1(s1), S2(s2), S3(s3), S4(s4),
    arr(AlignedAlloc(type, s1*s2*s3*s4)) {}

    //Rank4( const Rank4& other ) : S1( other.S1 ), S2( other.S2 ), S3( other.S3 ), S4( other.S4 ),
    //arr( other.arr ) {}

    type& operator()(size_t n1, size_t n2, size_t n3, size_t n4) const
    {
       check_bounds(n1,S1);
       check_bounds(n2,S2);
       check_bounds(n3,S3);
       check_bounds(n4,S4);
       ALIGNED(arr);
       return arr[n4+S4*(n3+S3*(n2+S2*n1))];
    }

    ~Rank4() { }

    size_t dim1() const { return S1; } 
    size_t dim2() const { return S2; }
    size_t dim3() const { return S3; }
    size_t dim4() const { return S4; }
    
    void free() { AlignedFree(arr); }

};

/******** end (i,j) arrays from Reger Ferrer and Vicenç Beltran ******/

/****** begin [i][j] arrays from Reger Ferrer and Vicenç Beltran ******/

template <class type>
class BracketRank1
{
    const size_t S1;
    type  * __restrict__ const  arr;

public:
    BracketRank1(size_t s1, void * __restrict__ const storage) : S1(s1),
         arr(reinterpret_cast<type * __restrict__ const>(storage)){}

    BracketRank1(size_t s1) : S1(s1), arr(new type[s1]){}

    type& operator[](size_t i) const
    {
        return arr[i];
    }

};

template <class type>
class BracketRank2
{
    const size_t S1, S2;
    type * __restrict__ const arr;

public:
    void free(){ delete[] arr; }
    BracketRank2(size_t s1, size_t s2, void *storage) : S1(s1), S2(s2),
         arr(reinterpret_cast<type * __restrict__ const>(storage)){}

    BracketRank2(size_t s1, size_t s2) : S1(s1), S2(s2),
         arr(new type[s1*s2]) {}

    BracketRank1<type> operator[](size_t i) const
    {
        return BracketRank1<type>(S2, arr + i * S2);
    }
    type& operator()(size_t n1, size_t n2) const
    {
       ALIGNED(arr);
       return arr[n2+S2*n1];
    }
};

/******** end [i][j] arrays from Reger Ferrer and Vicenç Beltran ******/

using namespace std;

template <class type>
void testArr2_diagonal()
{
   const int ITERS = 10000;
   const size_t dim1 = 64;
   const size_t dim2 = 64;

   BracketRank2<type> Abra(dim1, dim2);
   BracketRank2<type> Bbra(dim1, dim2);
   BracketRank2<type> Cbra(dim1, dim2);

   Rank2<type> Apar(dim1, dim2);
   Rank2<type> Bpar(dim1, dim2);
   Rank2<type> Cpar(dim1, dim2);

   FixedArray2D<type, dim1, dim2> Afix ;
   FixedArray2D<type, dim1, dim2> Bfix ;
   FixedArray2D<type, dim1, dim2> Cfix ;

   type** Aold = newArr2(type, dim1, dim2);
   type** Bold = newArr2(type, dim1, dim2);
   type** Cold = newArr2(type, dim1, dim2);

   array_ref2<type> Aarr(dim1, dim2);
   array_ref2<type> Barr(dim1, dim2);
   array_ref2<type> Carr(dim1, dim2);

   printf("Initializing data ...\n");
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Bbra[i][j] = rand();
      Cbra[i][j] = rand();
      Bpar(i,j) = Bbra[i][j];
      Cpar(i,j) = Cbra[i][j];
      Bfix.fetch(i,j) = Bbra[i][j];
      Cfix.fetch(i,j) = Cbra[i][j];
      Bold[i][j] = Bbra[i][j];
      Cold[i][j] = Cbra[i][j];
      Barr.fetch(i,j) = Bbra[i][j];
      Carr.fetch(i,j) = Cbra[i][j];
   }

   stopwatch(START);
   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Aold[i][j] = Bold[i][j] * Cold[i][j];
   }
   printf("%d ms = Total time [i][j] chained-pointer array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      //Afix.fetch(i,j) = Bfix.fetch(i,j) * Cfix.fetch(i,j);
      Afix[i][j] = Bfix[i][j] * Cfix[i][j];
   }
   printf("%d ms = Total time [i][j] fixed-dimension array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Abra[i][j] = Bbra[i][j] * Cbra[i][j];
   }
   printf("%d ms = Total time [i][j] Vicenc array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Apar.fetch(i,j) = Bpar.fetch(i,j) * Cpar.fetch(i,j);
      //Apar(i,j) = Bpar(i,j) * Cpar(i,j);
   }
   printf("%d ms = Total time (i,j) Vicenc array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      Aarr[i][j] = Barr[i][j] * Carr[i][j];
   }
   printf("%d ms = Total time [i][j] access of array_ref2\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      //Aarr(i,j) = Barr(i,j) * Carr(i,j);
      Aarr.fetch(i,j) = Barr.fetch(i,j) * Carr.fetch(i,j);
   }
   printf("%d ms = Total time (i,j) access of array_ref2\n", tv_to_ms(stopwatch(LAP)));

   for(size_t i=0; i<dim1; i++)
   for(size_t j=i; j<dim2; j++)
   {
      assert(Afix.fetch(i,j) == Abra[i][j]);
      assert(Aold[i][j] == Abra[i][j]);
      assert(Aarr.get(i,j) == Abra[i][j]);
   }

   printf("Verification done!\n");
   stopwatch(STOP);

   delArr2(Aold,dim1);
   delArr2(Bold,dim1);
   delArr2(Cold,dim1);
}

template <class type>
void testArr2()
{
   const int ITERS = 10000;
   const size_t dim1 = 64;
   const size_t dim2 = 64;

   BracketRank2<type> Abra(dim1, dim2);
   BracketRank2<type> Bbra(dim1, dim2);
   BracketRank2<type> Cbra(dim1, dim2);

   Rank2<type> Apar(dim1, dim2);
   Rank2<type> Bpar(dim1, dim2);
   Rank2<type> Cpar(dim1, dim2);

   FixedArray2D<type, dim1, dim2> Afix ;
   FixedArray2D<type, dim1, dim2> Bfix ;
   FixedArray2D<type, dim1, dim2> Cfix ;

   type** Aold = newArr2(type, dim1, dim2);
   type** Bold = newArr2(type, dim1, dim2);
   type** Cold = newArr2(type, dim1, dim2);

   array_ref2<type> Aarr(dim1, dim2);
   array_ref2<type> Barr(dim1, dim2);
   array_ref2<type> Carr(dim1, dim2);

   printf("Initializing data ...\n");
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Bbra[i][j] = rand();
      Cbra[i][j] = rand();
      Bpar(i,j) = Bbra[i][j];
      Cpar(i,j) = Cbra[i][j];
      Bfix.fetch(i,j) = Bbra[i][j];
      Cfix.fetch(i,j) = Cbra[i][j];
      Bold[i][j] = Bbra[i][j];
      Cold[i][j] = Cbra[i][j];
      Barr.fetch(i,j) = Bbra[i][j];
      Carr.fetch(i,j) = Cbra[i][j];
   }

   stopwatch(START);
   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Aold[i][j] = Bold[i][j] * Cold[i][j];
   }
   printf("%d ms = Total time [i][j] chained-pointer array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Afix.fetch(i,j) = Bfix.fetch(i,j) * Cfix.fetch(i,j);
   }
   printf("%d ms = Total time [i][j] fixed-dimension array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Abra[i][j] = Bbra[i][j] * Cbra[i][j];
   }
   printf("%d ms = Total time [i][j] Vicenc array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Apar(i,j) = Bpar(i,j) * Cpar(i,j);
   }
   printf("%d ms = Total time (i,j) Vicenc array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Aarr[i][j] = Barr[i][j] * Carr[i][j];
   }
   printf("%d ms = Total time [i][j] access of array_ref2\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      Aarr.fetch(i,j) = Barr.get(i,j) * Carr.get(i,j);
   }
   printf("%d ms = Total time (i,j) access of array_ref2\n", tv_to_ms(stopwatch(LAP)));

   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   {
      assert(Afix.fetch(i,j) == Abra[i][j]);
      assert(Aold[i][j] == Abra[i][j]);
      assert(Aarr.get(i,j) == Abra[i][j]);
   }

   printf("Verification done!\n");
   stopwatch(STOP);

   delArr2(Aold,dim1);
   delArr2(Bold,dim1);
   delArr2(Cold,dim1);
}

#define testArr3nestedFor(arg1, arg2) \
for(int t=0; t<ITERS; t++) \
for(size_t i=0; i<dim1; i++) \
for(size_t j=0; j<dim2; j++) \
for(size_t k=0; k<dim3; k++) \
{ \
   #arg1; \
} \
printf("%d ms = Total time " #arg2 "\n", tv_to_ms(stopwatch(LAP)));

template <class type>
void set_prod3(array_ref3<type> Aarr,const_arr3<type> Barr,array_ref3<type> Carr,int ITERS, size_t dim1,size_t dim2,size_t dim3)
{
   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      //Aarr[i][j][k] = Barr[i][j][k] * Carr[i][j][k];
      Aarr.fetch(i,j,k) = Barr.get(i,j,k) * Carr.get(i,j,k);
   }
   printf("%d ms = Total time [i][j][k] access of array_ref3\n", tv_to_ms(stopwatch(LAP)));
}

template <class type>
void testArr3()
{
   const int ITERS = 100;
   const size_t dim1 = 64;
   const size_t dim2 = 64;
   const size_t dim3 = 64;

   Rank3<type> Apar(dim1, dim2, dim3);
   Rank3<type> Bpar(dim1, dim2, dim3);
   Rank3<type> Cpar(dim1, dim2, dim3);

   FixedArray3D<type, dim1, dim2, dim3> Afix ;
   FixedArray3D<type, dim1, dim2, dim3> Bfix ;
   FixedArray3D<type, dim1, dim2, dim3> Cfix ;

   type*** Aold = newArr3(type, dim1, dim2, dim3);
   type*** Bold = newArr3(type, dim1, dim2, dim3);
   type*** Cold = newArr3(type, dim1, dim2, dim3);

   //array3<type> Aarr(dim1, dim2, dim3);
   //array3<type> Barr(dim1, dim2, dim3);
   //array3<type> Carr(dim1, dim2, dim3);
   array_ref3<type> Aarr(dim1, dim2, dim3);
   array_ref3<type> Barr(dim1, dim2, dim3);
   array_ref3<type> Carr(dim1, dim2, dim3);

   printf("Initializing data ...\n");
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Barr.fetch(i,j,k) = rand();
      Carr.fetch(i,j,k) = rand();
      Bpar.fetch(i,j,k) = Barr.get(i,j,k);
      Cpar.fetch(i,j,k) = Carr.get(i,j,k);
      Bfix.fetch(i,j,k) = Barr.get(i,j,k);
      Cfix.fetch(i,j,k) = Carr.get(i,j,k);
      Bold[i][j][k] = Barr.get(i,j,k);
      Cold[i][j][k] = Carr.get(i,j,k);
   }

   stopwatch(START);

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Aold[i][j][k] = Bold[i][j][k] * Cold[i][j][k];
   }
   printf("%d ms = Total time [i][j][k] chained-pointer array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Afix[i][j][k] = Bfix[i][j][k] * Cfix[i][j][k];
      //Afix.arr[i][j][k] = Bfix.arr[i][j][k] * Cfix.arr[i][j][k];
      //Afix.fetch(i,j,k) = Bfix.fetch(i,j,k) * Cfix.fetch(i,j,k);
   }
   printf("%d ms = Total time [i][j][k] fixed-dimension array\n", tv_to_ms(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Apar.fetch(i,j,k) = Bpar.fetch(i,j,k) * Cpar.fetch(i,j,k);
   }
   printf("%d ms = Total time (i,j,k) Vicenc array\n", tv_to_ms(stopwatch(LAP)));

   set_prod3(Aarr,Barr,Carr,ITERS,dim1,dim2,dim3);

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      Aarr.fetch(i,j,k) = Barr.get(i,j,k) * Carr.get(i,j,k);
   }
   printf("%d ms = Total time (i,j,k) access of array_ref3\n", tv_to_ms(stopwatch(LAP)));

   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   {
      assert_eq(Aold[i][j][k], Aarr.get(i,j,k));
      assert_eq(Apar.fetch(i,j,k), Aarr.get(i,j,k));
      assert_eq(Afix.fetch(i,j,k), Aarr.get(i,j,k));
   }

   printf("Verification done!\n");
   stopwatch(STOP);

   // automatic destructor slows array access
   // unless compiling with -fno-exceptions
   //
   Apar.free();
   Bpar.free();
   Cpar.free();
   Aarr.free();
   Barr.free();
   Carr.free();
}

template <class type>
void testArr4()
{
   // For some bizarre reason, if I comment out the code for the
   // "fbr" and "fpa" arrays below then icpc on knc2 is somehow
   // able to figure out that each iteration does the same thing
   // in the case of array_ref4, but not in the case of the chained
   // pointer or fixed-dimension arrays.  Why not?  And why
   // does this optimization occur for four-dimensional arrays
   // and not for 3- or 2-dimensional arrays?  And why is this
   // optimization no longer performed if "fbr" and "fpa" stuff
   // is included?  The times are baffling.
   const int ITERS = 1;
   const size_t dim1 = 16;
   const size_t dim2 = 16;
   const size_t dim3 = 16;
   const size_t dim4 = 16;

   FixedArray4D<type, dim1, dim2, dim3, dim4> Afix;
   FixedArray4D<type, dim1, dim2, dim3, dim4> Bfix;
   FixedArray4D<type, dim1, dim2, dim3, dim4> Cfix;

   type**** Aold = newArr4(type, dim1, dim2, dim3, dim4);
   type**** Bold = newArr4(type, dim1, dim2, dim3, dim4);
   type**** Cold = newArr4(type, dim1, dim2, dim3, dim4);

   //array4<type> Afbr(dim1, dim2, dim3, dim4);
   //array4<type> Bfbr(dim1, dim2, dim3, dim4);
   //array4<type> Cfbr(dim1, dim2, dim3, dim4);

   //array4<type> Afpa(dim1, dim2, dim3, dim4);
   //array4<type> Bfpa(dim1, dim2, dim3, dim4);
   //array4<type> Cfpa(dim1, dim2, dim3, dim4);

   array_ref4<type> Abra(dim1, dim2, dim3, dim4);
   array_ref4<type> Bbra(dim1, dim2, dim3, dim4);
   array_ref4<type> Cbra(dim1, dim2, dim3, dim4);

   array_ref4<type> Apar(dim1, dim2, dim3, dim4);
   array_ref4<type> Bpar(dim1, dim2, dim3, dim4);
   array_ref4<type> Cpar(dim1, dim2, dim3, dim4);

   printf("Initializing data ...\n");
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   for(size_t l=0; l<dim4; l++)
   {
      Bbra.fetch(i,j,k,l) = rand();
      Cbra.fetch(i,j,k,l) = rand();
      Bpar.fetch(i,j,k,l) = Bbra.get(i,j,k,l);
      Cpar.fetch(i,j,k,l) = Cbra.get(i,j,k,l);
      //Bfbr.fetch(i,j,k,l) = Bbra.get(i,j,k,l);
      //Cfbr.fetch(i,j,k,l) = Cbra.get(i,j,k,l);
      //Bfpa.fetch(i,j,k,l) = Bbra.get(i,j,k,l);
      //Cfpa.fetch(i,j,k,l) = Cbra.get(i,j,k,l);
      Bfix.fetch(i,j,k,l) = Bbra.get(i,j,k,l);
      Cfix.fetch(i,j,k,l) = Cbra.get(i,j,k,l);
      Bold[i][j][k][l] = Bbra.get(i,j,k,l);
      Cold[i][j][k][l] = Cbra.get(i,j,k,l);
   }

   stopwatch(START);

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   for(size_t l=0; l<dim4; l++)
   {
      Aold[i][j][k][l] = Bold[i][j][k][l] * Cold[i][j][k][l];
   }
   printf("%d us = Total time [i][j][k][l] chained-pointer array\n", tv_to_us(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   for(size_t l=0; l<dim4; l++)
   {
      Afix[i][j][k][l] = Bfix[i][j][k][l] * Cfix[i][j][k][l];
      //Afix.arr[i][j][k][l] = Bfix.arr[i][j][k][l] * Cfix.arr[i][j][k][l];
      //Afix.fetch(i,j,k,l) = Bfix.fetch(i,j,k,l) * Cfix.fetch(i,j,k,l);
   }
   printf("%d us = Total time [i][j][k][l] fixed-dimension array\n", tv_to_us(stopwatch(LAP)));

   //for(int t=0; t<ITERS; t++)
   //for(size_t i=0; i<dim1; i++)
   //for(size_t j=0; j<dim2; j++)
   //for(size_t k=0; k<dim3; k++)
   //for(size_t l=0; l<dim4; l++)
   //{
   //   Afbr.fetch(i,j,k,l) = Bfbr.get(i,j,k,l) * Cfbr.get(i,j,k,l);
   //}
   //printf("%d us = Total time (i,j,k,l) access of array4\n", tv_to_us(stopwatch(LAP)));

   //for(int t=0; t<ITERS; t++)
   //for(size_t i=0; i<dim1; i++)
   //for(size_t j=0; j<dim2; j++)
   //for(size_t k=0; k<dim3; k++)
   //for(size_t l=0; l<dim4; l++)
   //{
   //   Afpa.fetch(i,j,k,l) = Bfpa.get(i,j,k,l) * Cfpa.get(i,j,k,l);
   //}
   //printf("%d us = Total time (i,j,k,l) access of array4\n", tv_to_us(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   for(size_t l=0; l<dim4; l++)
   {
      Abra[i][j][k][l] = Bbra[i][j][k][l] * Cbra[i][j][k][l];
   }
   printf("%d us = Total time [i][j][k][l] access of array_ref4\n", tv_to_us(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   for(size_t l=0; l<dim4; l++)
   {
      Apar.fetch(i,j,k,l) = Bpar.get(i,j,k,l) * Cpar.get(i,j,k,l);
   }
   printf("%d us = Total time (i,j,k,l) access of array_ref4\n", tv_to_us(stopwatch(LAP)));

   for(int t=0; t<ITERS; t++)
   for(size_t i=0; i<dim1; i++)
   for(size_t j=0; j<dim2; j++)
   for(size_t k=0; k<dim3; k++)
   for(size_t l=0; l<dim4; l++)
   {
      assert_eq(Aold[i][j][k][l], Abra.get(i,j,k,l));
      assert_eq(Apar.fetch(i,j,k,l), Abra.get(i,j,k,l));
      //assert_eq(Afbr[i][j][k][l], Abra.get(i,j,k,l));
      //assert_eq(Afpa.fetch(i,j,k,l), Abra.get(i,j,k,l));
      assert_eq(Afix.fetch(i,j,k,l), Abra.get(i,j,k,l));
   }

   printf("Verification done!\n");
   stopwatch(STOP);

   Apar.free();
   Bpar.free();
   Cpar.free();
   Abra.free();
   Bbra.free();
   Cbra.free();
}

int main()
{
  //printf("=== testing array_ref2<int> (diagonal) ===\n");
  //testArr2_diagonal<int>();
  //printf("=== testing array_ref2<double> (diagonal) ===\n");
  //testArr2_diagonal<double>();
  printf("=== testing array_ref2<int> ===\n");
  testArr2<int>();
  printf("=== testing array_ref2<double> ===\n");
  testArr2<double>();
  printf("=== testing array_ref3<int> ===\n");
  testArr3<int>();
  printf("=== testing array_ref3<double> ===\n");
  testArr3<double>();
  printf("=== testing array_ref4<int> ===\n");
  testArr4<int>();
  printf("=== testing array_ref4<double> ===\n");
  testArr4<double>();
}
