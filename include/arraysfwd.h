/* forward declaration for array classes */
#ifndef arraysfwd_h
#define arraysfwd_h

namespace iPic3D
{
  template <class T>
  class const_array_ref3;
  template <class T>
  class const_array_ref4;
  template <class T>
  class array_ref1;
  template <class T>
  class array_ref2;
  template <class T>
  class array_ref3;
  template <class T>
  class array_ref4;
  template <class T>
  class array1;
  template <class T>
  class array2;
  template <class T>
  class array3;
  template <class T>
  class array4;
}

// These aliases are defined for the following flexibilization purposes:
// - to avoid filling the code with template brackets
//   (i.e., to minimize explicitly template-dependent code).
// - so that they can be redefined according to the user's
//   preferred array implementation.
//
//typedef array_ref1<int> intArr1;
//typedef array_ref2<int> intArr2;
//typedef array_ref3<int> intArr3;
//typedef array_ref4<int> intArr4;
//typedef const_array_ref1<double> array_ref1_double;
//typedef const_array_ref2<double> array_ref2_double;
//
typedef iPic3D::const_array_ref3<double> const_arr3_double;
typedef iPic3D::const_array_ref4<double> const_arr4_double;
typedef iPic3D::array_ref1<double> array_ref1_double;
typedef iPic3D::array_ref2<double> array_ref2_double;
typedef iPic3D::array_ref3<double> array_ref3_double;
typedef iPic3D::array_ref4<double> array_ref4_double;
typedef iPic3D::array1<double> array1_double;
typedef iPic3D::array2<double> array2_double;
typedef iPic3D::array3<double> array3_double;
typedef iPic3D::array4<double> array4_double;
#endif
