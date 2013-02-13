/***************************************************************************

  -------------------

developers           : Stefano Markidis, Giovanni Lapenta

 ***************************************************************************/
#ifndef Alloc_H
#define Alloc_H
/** subroutines for allocation and deallocation of arrays 2D, 3D, 4D */
/**
  2 dimensional arrays
  */

/** The allocator for 2D array */
    template <class type>
inline type **_new_2_array(int sz1,int sz2,type *stupid)
{
    type **foo;
    //foo = new (type *)[sz1];
    foo = new type *[sz1];
    for (int i=0;i<sz1;i++) foo[i] = new type[sz2];
    return foo;
}
/** macro for allocate 2D array */
#define newArr(type,sz1,sz2) _new_2_array((sz1),(sz2),(type *) NULL)

/** deallocator for a 2D array*/
    template <class type>
inline void delArr(type **foo,int sz1)
{ for (int i=0;i<sz1;i++) delete[] foo[i]; delete[] foo; }

/**
  3 dimensional arrays
  */

/** The allocator for 3D array */
    template <class type>
inline type ***_new_3_array(int sz1,int sz2,int sz3,type *stupid)
{
    type ***foo;
    foo = new type **[sz1];
    for (int i=0;i<sz1;i++) foo[i] = newArr(type,sz2,sz3);
    return foo;
}
/** macro for allocate 3D array */
#define newArr3(type,sz1,sz2,sz3) _new_3_array((sz1),(sz2),(sz3),(type *) NULL)

/** deallocator for a 3D array*/
    template <class type>
inline void delArr3(type ***foo,int sz1,int sz2)
{ for (int i=0;i<sz1;i++) delArr(foo[i],sz2); delete[] foo; }

/**
  4 dimensional arrays
  */

/** The allocator for 4D array */
    template <class type>
inline type ****_new_4_array(int sz1,int sz2,int sz3,int sz4,type *stupid)
{
    type ****foo;
    foo = new type ***[sz1];
    for (int i=0;i<sz1;i++) foo[i] = newArr3(type,sz2,sz3,sz4);
    return foo;
}
/** macro for allocate 4D array */
#define newArr4(type,sz1,sz2,sz3,sz4) \
    _new_4_array((sz1),(sz2),(sz3),(sz4),(type *) NULL);

/** deallocator for a 4D array*/
    template <class type>
inline void delArr4(type ****foo,int sz1,int sz2,int sz3)
{ for (int i=0;i<sz1;i++) delArr3(foo[i],sz2,sz3); delete[] foo; }

#endif
