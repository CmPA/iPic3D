
#ifndef ALLOC_H
#define ALLOC_H

#include <stdlib.h>

/*! The allocator for 4D array */
template < class type > type **** _new_4_array(int sz1, int sz2, int sz3, int sz4) {

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
template < class type > void delArr4(type **** arr, int dummyx, int dummyy, int dummyz) {
  delete[]arr[0][0][0];
  delete[]arr[0][0];
  delete[]arr[0];
  delete[]arr;
}

/*! The allocator for 3D array */
template < class type > type *** _new_3_array(int sz1, int sz2, int sz3) {

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
template < class type > void delArr3(type *** arr, int dummyx, int dummyy) {
  delete[]arr[0][0];
  delete[]arr[0];
  delete[]arr;
}

/*! The allocator for 2D array */
template < class type > type ** _new_2_array(int sz1, int sz2) {

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
template < class type > void delArr2(type ** arr, int dummyx) {
  delete[]arr[0];
  delete[]arr;
}

#define newArr4(type,sz1,sz2,sz3,sz4) _new_4_array<type>((sz1),(sz2),(sz3),(sz4))
#define newArr3(type,sz1,sz2,sz3) _new_3_array<type>((sz1),(sz2),(sz3))
#define newArr2(type,sz1,sz2) _new_2_array<type>((sz1),(sz2))

#endif
