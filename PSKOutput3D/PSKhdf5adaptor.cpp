
#include <mpi.h>
#include "PSKhdf5adaptor.h"

using namespace PSK;


//void HDF5OutputAdaptor::get_dataset_context(const std::string & name, std::vector < hid_t > &hid_array, std::string & dataset_name) {
//  hid_array.clear();
//
//  std::vector < std::string > name_components;
//
//  split_name(name, name_components);
//
//  // for( int i=0; i< name_components.size() ; ++i )
//  // std::cout<< i << ": <" << name_components[i] << ">\n";
//
//  int ncompx = name_components.size();
//
//  hid_array.resize(ncompx);
//
//  if (ncompx == 0) {
//    throw PSK::OutputException("HDF5OutputAdaptor::get_dataset_context()>> zero name components");
//
//  }
//  else if (ncompx == 1) {
//    hid_array[0] = _hdf5_file_id;
//    dataset_name = name_components[0];
//  }
//  else {
//    /* HDF5 error handling code from http://hdf.ncsa.uiuc.edu/HDF5/doc/Errors.html */
//    /* Save old error handler */
//    // herr_t (*old_func)(void*); // HDF 1.6
//    H5E_auto2_t old_func;        // HDF 1.8.8
//    void *old_client_data;
//    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);  // HDF 1.8
//
//    hid_array[0] = _hdf5_file_id;
//
//    for (int i = 0; i < ncompx - 1; ++i) {
//
//      // std::cout << "group open/create " << i << ": <" << name_components[i] << ">\n";
//
//      dataset_name = name_components[ncompx - 1];
//
//      /* Turn off error handling */
//      // H5Eset_auto2(NULL, NULL); // HDF 1.6
//      H5Eset_auto2(H5E_DEFAULT, 0, 0); // HDF 1.8
//      hid_array[i + 1] = H5Gopen2(hid_array[i], name_components[i].c_str(), H5P_DEFAULT);
//      /* Restore previous error handler */
//      H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);  // HDF 1.8
//
//      if (hid_array[i + 1] < 0) {
//
//        // std::cout << "group open failed \n" ;
//
//        hid_array[i + 1] = H5Gcreate2(hid_array[i], name_components[i].c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  // HDF 1.8
//
//        if (hid_array[i + 1] < 0) {
//
//          // std::cout << "group create failed \n" ;
//
//          throw PSK::OutputException("Failed to open/create group for <" + name + "> at element <" + name_components[i] + ">", "HDF5OutputAdaptor::get_dataset_context()");
//        }
//
//      }
//
//    }
//  }                             // end more than one name component
//}
//
//
///* 
// * 
// * abort if zero length
// * 
// * add leading "/" if doesn't start with "/"
// * 
// */
//std::string HDF5OutputAdaptor::purify_object_name(const std::string & objname) {
//  if (objname.length() == 0)
//    throw PSK::OutputException("Zero length tag name", "HDF5OutputAdaptor::purify_object_name()");
//
//  return objname[0] != '/' ? "/" + objname : objname;
//
//}
//
///* 
// * split name into elements
// * 
// */
//void HDF5OutputAdaptor::split_name(const std::string & name, std::vector < std::string > &elements) {
//
//  elements.clear();
//
//  int endidx = name.length() - 1;
//  int startidx;
//
//  do {
//    startidx = name.rfind("/", endidx);
//    // std::cout<< startidx << " " << endidx << "\n";
//    elements.push_back(name.substr(startidx + 1, endidx - startidx));
//    endidx = startidx - 1;
//  }
//  while (startidx > 0);
//
//  reverse(elements.begin(), elements.end());
//
//}
//
///* 
// * 
// */
//void HDF5OutputAdaptor::open(const std::string & name) {
//
//  // use H5F_ACC_TRUNC to delete any existing file
//  // or ACC_EXCL to fail if file exists
//
//  _hdf5_file_id = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//
//  if (_hdf5_file_id <= 0) {
//    PSK::OutputException e("H5FCreate fails", "HDF5OutputAdaptor::open2()");
//
//    // if using H5F_ACC_EXCL
//    // e.push( "Using H5F_ACC_EXCL: Check if file " +name + " already exists" );
//
//    throw e;
//  }
//
//  _hdf5_file_name = name;
//
//}
//
///* 
// * Opens an existing file */
//void HDF5OutputAdaptor::open_append(const std::string & name) {
//  /* HDF5 error handling code from http://hdf.ncsa.uiuc.edu/HDF5/doc/Errors.html */
//  /* Save old error handler */
//  // herr_t (*old_func)(void*); //HDF 1.6
//  H5E_auto2_t old_func;          // HDF 1.8.8
//  void *old_client_data;
//  H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);  // HDF 1.8
//
//
//  /* Turn off error handling */
//  // H5Eset_auto2(NULL, NULL); // HDF 1.6
//  H5Eset_auto2(H5E_DEFAULT, 0, 0);
//  _hdf5_file_id = H5Fcreate(name.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
//
//  /* Restore previous error handler */
//  H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);  // HDF 1.8
//
//  if (_hdf5_file_id <= 0) {
//
//    _hdf5_file_id = H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
//  }
//
//  if (_hdf5_file_id <= 0) {
//    PSK::OutputException e("H5Fopen fails", "HDF5OutputAdaptor::open_append()");
//
//    throw e;
//  }
//  _hdf5_file_name = name;
//
//}
//
//
///* 
// * 
// */
//void HDF5OutputAdaptor::close(void) {
//
//  herr_t hdf5err = H5Fclose(_hdf5_file_id);
//
//  if (hdf5err < 0) {
//    PSK::OutputException e("HDF5OutputAdaptor::close()>> H5FClose fails");
//    throw e;
//  }
//
//  _hdf5_file_name.clear();
//  _hdf5_file_id = 0;
//
//}
//
//
//
///* 
// * 
// * 
// */
//void HDF5OutputAdaptor::write(const std::string & tag, int i_value) {
//  try {
//    std::string ptag = purify_object_name(tag);
//
//    std::vector < hid_t > hid_array;
//    std::string dataset_name;
//
//    get_dataset_context(tag, hid_array, dataset_name);
//
//    hsize_t *hdf5dims = new hsize_t[1];
//    hdf5dims[0] = 1;
//
//    herr_t hdf5err = H5LTfind_dataset(hid_array[hid_array.size() - 1],
//                                      dataset_name.c_str());
//    if (hdf5err < 1) {
//      herr_t hdf5err = H5LTmake_dataset_int(hid_array[hid_array.size() - 1],
//                                            dataset_name.c_str(),
//                                            1, hdf5dims, &i_value);
//    }
//    else {
//      hid_t dataset_id = H5Dopen2(hid_array[hid_array.size() - 1], dataset_name.c_str(), H5P_DEFAULT); // HDF 1.8
//      herr_t hdf5err = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &i_value);
//      hdf5err = H5Dclose(dataset_id);
//    }
//
//    if (hdf5err < 0) {
//      PSK::OutputException e("make_dataset fails for " + tag, "HDF5OutputAdaptor::write(int)");
//      throw e;
//    }
//
//    // close groups, if any, but don't try to close the file id at [0]
//    for (int i = hid_array.size() - 1; i > 0; --i)
//      hdf5err = H5Gclose(hid_array[i]);
//
//    delete[]hdf5dims;
//
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(int)");
//    throw e;
//  }
//}
//// new long writing
//void HDF5OutputAdaptor::write(const std::string & tag, long i_value) {
//  try {
//    std::string ptag = purify_object_name(tag);
//
//    std::vector < hid_t > hid_array;
//    std::string dataset_name;
//
//    get_dataset_context(tag, hid_array, dataset_name);
//
//    hsize_t *hdf5dims = new hsize_t[1];
//    hdf5dims[0] = 1;
//
//    herr_t hdf5err = H5LTfind_dataset(hid_array[hid_array.size() - 1],
//                                      dataset_name.c_str());
//    if (hdf5err < 1) {
//      herr_t hdf5err = H5LTmake_dataset_long(hid_array[hid_array.size() - 1],
//                                             dataset_name.c_str(),
//                                             1, hdf5dims, &i_value);
//    }
//    else {
//      hid_t dataset_id = H5Dopen2(hid_array[hid_array.size() - 1], dataset_name.c_str(), H5P_DEFAULT); // HDF 1.8
//      herr_t hdf5err = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &i_value);
//      hdf5err = H5Dclose(dataset_id);
//    }
//
//    if (hdf5err < 0) {
//      PSK::OutputException e("make_dataset fails for " + tag, "HDF5OutputAdaptor::write(long)");
//      throw e;
//    }
//
//    // close groups, if any, but don't try to close the file id at [0]
//    for (int i = hid_array.size() - 1; i > 0; --i)
//      hdf5err = H5Gclose(hid_array[i]);
//
//    delete[]hdf5dims;
//
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(long)");
//    throw e;
//  }
//}
//
//// 
//void HDF5OutputAdaptor::write(const std::string & tag, const Dimens dimens, const int *i_array) {
//  try {
//    if (dimens.size() == 0) {
//      PSK::OutputException e("Zero Dimens size", "HDF5OutputAdaptor::write(int* array)");
//      throw e;
//    }
//
//    std::string ptag = purify_object_name(tag);
//
//    std::vector < hid_t > hid_array;
//    std::string dataset_name;
//
//    get_dataset_context(tag, hid_array, dataset_name);
//
//    hsize_t *hdf5dims = new hsize_t[dimens.size()];
//    for (int i = 0; i < dimens.size(); ++i)
//      hdf5dims[i] = dimens[i];
//
//    herr_t hdf5err = H5LTfind_dataset(hid_array[hid_array.size() - 1],
//                                      dataset_name.c_str());
//    if (hdf5err < 1) {
//      herr_t hdf5err = H5LTmake_dataset_int(hid_array[hid_array.size() - 1],
//                                            dataset_name.c_str(),
//                                            dimens.size(), hdf5dims, i_array);
//    }
//    else {
//      hid_t dataset_id = H5Dopen2(hid_array[hid_array.size() - 1], dataset_name.c_str(), H5P_DEFAULT); // HDF 1.8
//      herr_t hdf5err = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, i_array);
//      hdf5err = H5Dclose(dataset_id);
//    }
//
//    if (hdf5err < 0) {
//      PSK::OutputException e("make_dataset fails for " + tag, "HDF5OutputAdaptor::write(int* array)");
//      throw e;
//    }
//
//    // close groups, if any, but don't try to close the file id at [0]
//    for (int i = hid_array.size() - 1; i > 0; --i)
//      hdf5err = H5Gclose(hid_array[i]);
//
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(int* array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & tag, const Dimens dimens, const long *i_array) {
//  try {
//    if (dimens.size() == 0) {
//      PSK::OutputException e("Zero Dimens size", "HDF5OutputAdaptor::write(long* array)");
//      throw e;
//    }
//
//    std::string ptag = purify_object_name(tag);
//
//    std::vector < hid_t > hid_array;
//    std::string dataset_name;
//
//    get_dataset_context(tag, hid_array, dataset_name);
//
//    hsize_t *hdf5dims = new hsize_t[dimens.size()];
//    for (int i = 0; i < dimens.size(); ++i)
//      hdf5dims[i] = dimens[i];
//
//    herr_t hdf5err = H5LTmake_dataset_long(hid_array[hid_array.size() - 1],
//                                           dataset_name.c_str(),
//                                           dimens.size(), hdf5dims, i_array);
//
//    if (hdf5err < 0) {
//      PSK::OutputException e("make_dataset fails for " + tag, "HDF5OutputAdaptor::write(long* array)");
//      throw e;
//    }
//
//    // close groups, if any, but don't try to close the file id at [0]
//    for (int i = hid_array.size() - 1; i > 0; --i)
//      hdf5err = H5Gclose(hid_array[i]);
//
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(long* array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & tag, const Dimens dimens, const std::vector < int >&i_array) {
//  try {
//    int n = dimens.nels();
//    int *i_array_p = new int[n];
//    for (int i = 0; i < n; ++i)
//      i_array_p[i] = i_array[i];
//    write(tag, dimens, i_array_p);
//    delete[]i_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(vector<int> array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & tag, const Dimens dimens, const std::vector < long >&i_array) {
//  try {
//    int n = dimens.nels();
//    long *i_array_p = new long[n];
//    for (int i = 0; i < n; ++i)
//      i_array_p[i] = i_array[i];
//    write(tag, dimens, i_array_p);
//    delete[]i_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(vector<long> array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & objname, const Dimens dimens, const int ***i_array) {
//  if (dimens.size() != 3) {
//    PSK::OutputException e("Dimens size not 3 for object " + objname, "HDF5OutputAdaptor::write(int*** array)");
//    throw e;
//  }
//
//  try {
//    int nels = dimens.nels();
//    int *i_array_p = new int[nels];
//    const int di = dimens[0];
//    const int dj = dimens[1];
//    const int dk = dimens[2];
//    const int djk = dk * dj;
//    for (int i = 0; i < di; ++i)
//      for (int j = 0; j < dj; ++j)
//        for (int k = 0; k < dk; ++k)
//          i_array_p[i * djk + j * dk + k] = i_array[i][j][k];
//    write(objname, dimens, i_array_p);
//    delete[]i_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(int*** array)");
//    throw e;
//  }
//}
//
//
//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//// WRITE float FUNCTIONS
//// -----------------------------------------------------------------------
//
///* 
// * 
// * 
// */
//void HDF5OutputAdaptor::write(const std::string & tag, float f_value) {
//  try {
//    std::string ptag = purify_object_name(tag);
//
//    std::vector < hid_t > hid_array;
//    std::string dataset_name;
//
//    get_dataset_context(tag, hid_array, dataset_name);
//
//    hsize_t *hdf5dims = new hsize_t[1];
//    hdf5dims[0] = 1;
//
//    herr_t hdf5err = H5LTmake_dataset_float(hid_array[hid_array.size() - 1],
//                                            dataset_name.c_str(),
//                                            1, hdf5dims, &f_value);
//
//    if (hdf5err < 0) {
//      PSK::OutputException e("make_dataset fails for " + tag, "HDF5OutputAdaptor::write(float)");
//      throw e;
//    }
//
//    // close groups, if any, but don't try to close the file id at [0]
//    for (int i = hid_array.size() - 1; i > 0; --i)
//      hdf5err = H5Gclose(hid_array[i]);
//
//    delete[]hdf5dims;
//
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(float)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & tag, const Dimens dimens, const float *f_array) {
//  try {
//    if (dimens.size() == 0) {
//      PSK::OutputException e("Zero Dimens size", "HDF5OutputAdaptor::write(float* array)");
//      throw e;
//    }
//
//    std::string ptag = purify_object_name(tag);
//
//    std::vector < hid_t > hid_array;
//    std::string dataset_name;
//
//    get_dataset_context(tag, hid_array, dataset_name);
//
//    hsize_t *hdf5dims = new hsize_t[dimens.size()];
//    for (int i = 0; i < dimens.size(); ++i)
//      hdf5dims[i] = dimens[i];
//
//    herr_t hdf5err = H5LTmake_dataset_float(hid_array[hid_array.size() - 1],
//                                            dataset_name.c_str(),
//                                            dimens.size(), hdf5dims, f_array);
//
//    if (hdf5err < 0) {
//      PSK::OutputException e("make_dataset fails for " + tag, "HDF5OutputAdaptor::write(float* array)");
//      throw e;
//    }
//
//    // close groups, if any, but don't try to close the file id at [0]
//    for (int i = hid_array.size() - 1; i > 0; --i)
//      hdf5err = H5Gclose(hid_array[i]);
//
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(float* array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & tag, const Dimens dimens, const std::vector < float >&f_array) {
//  try {
//    int n = dimens.nels();
//    float *f_array_p = new float[n];
//    for (int i = 0; i < n; ++i)
//      f_array_p[i] = f_array[i];
//    write(tag, dimens, f_array_p);
//    delete[]f_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(vector<float> array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & objname, const Dimens dimens, const float ***f_array) {
//  if (dimens.size() != 3) {
//    PSK::OutputException e("Dimens size not 3 for object " + objname, "HDF5OutputAdaptor::write(float*** array)");
//    throw e;
//  }
//
//  try {
//    int nels = dimens.nels();
//    float *f_array_p = new float[nels];
//    const int di = dimens[0];
//    const int dj = dimens[1];
//    const int dk = dimens[2];
//    const int djk = dk * dj;
//    for (int i = 0; i < di; ++i)
//      for (int j = 0; j < dj; ++j)
//        for (int k = 0; k < dk; ++k)
//          f_array_p[i * djk + j * dk + k] = f_array[i][j][k];
//    write(objname, dimens, f_array_p);
//    delete[]f_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(float*** array)");
//    throw e;
//  }
//}
//
//
//// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//// WRITE double FUNCTIONS
//// -----------------------------------------------------------------------
//
///* 
// * 
// * 
// */
//void HDF5OutputAdaptor::write(const std::string & tag, double d_value) {
//  try {
//    std::string ptag = purify_object_name(tag);
//
//    std::vector < hid_t > hid_array;
//    std::string dataset_name;
//
//    get_dataset_context(tag, hid_array, dataset_name);
//
//    hsize_t *hdf5dims = new hsize_t[1];
//    hdf5dims[0] = 1;
//
//    herr_t hdf5err = H5LTfind_dataset(hid_array[hid_array.size() - 1],
//                                      dataset_name.c_str());
//    if (hdf5err < 1) {
//      herr_t hdf5err = H5LTmake_dataset_double(hid_array[hid_array.size() - 1],
//                                               dataset_name.c_str(),
//                                               1, hdf5dims, &d_value);
//    }
//    else {
//      hid_t dataset_id = H5Dopen2(hid_array[hid_array.size() - 1], dataset_name.c_str(), H5P_DEFAULT); // HDF 1.8
//      herr_t hdf5err = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &d_value);
//      hdf5err = H5Dclose(dataset_id);
//    }
//
//    if (hdf5err < 0) {
//      PSK::OutputException e("make_dataset fails for " + tag, "HDF5OutputAdaptor::write(double)");
//      throw e;
//    }
//
//    // close groups, if any, but don't try to close the file id at [0]
//    for (int i = hid_array.size() - 1; i > 0; --i)
//      hdf5err = H5Gclose(hid_array[i]);
//
//    delete[]hdf5dims;
//
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(double)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & tag, const Dimens dimens, const double *d_array) {
//  try {
//    if (dimens.size() == 0) {
//      PSK::OutputException e("Zero Dimens size", "HDF5OutputAdaptor::write(double* array)");
//      throw e;
//    }
//
//    std::string ptag = purify_object_name(tag);
//
//    std::vector < hid_t > hid_array;
//    std::string dataset_name;
//
//    get_dataset_context(tag, hid_array, dataset_name);
//
//    hsize_t *hdf5dims = new hsize_t[dimens.size()];
//    for (int i = 0; i < dimens.size(); ++i)
//      hdf5dims[i] = dimens[i];
//
//    herr_t hdf5err = H5LTfind_dataset(hid_array[hid_array.size() - 1],
//                                      dataset_name.c_str());
//    if (hdf5err < 1) {
//      herr_t hdf5err = H5LTmake_dataset_double(hid_array[hid_array.size() - 1],
//                                               dataset_name.c_str(),
//                                               dimens.size(), hdf5dims, d_array);
//    }
//    else {
//      hid_t dataset_id = H5Dopen2(hid_array[hid_array.size() - 1], dataset_name.c_str(), H5P_DEFAULT); // HDF 1.8
//      hdf5err = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, d_array);
//      hdf5err = H5Dclose(dataset_id);
//    }
//
//
//    if (hdf5err < 0) {
//      PSK::OutputException e("make_dataset fails for " + tag, "HDF5OutputAdaptor::write(double* array)");
//      throw e;
//    }
//
//    // close groups, if any, but don't try to close the file id at [0]
//    for (int i = hid_array.size() - 1; i > 0; --i)
//      hdf5err = H5Gclose(hid_array[i]);
//
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(double* array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & tag, const Dimens dimens, const std::vector < double >&d_array) {
//  try {
//    int n = dimens.nels();
//    double *d_array_p = new double[n];
//    for (int i = 0; i < n; ++i)
//      d_array_p[i] = d_array[i];
//    write(tag, dimens, d_array_p);
//    delete[]d_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(vector<double> array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & objname, const Dimens dimens, double ***d_array) {
//  if (dimens.size() != 3) {
//    PSK::OutputException e("Dimens size not 3 for object " + objname, "HDF5OutputAdaptor::write(double*** array)");
//    throw e;
//  }
//
//  try {
//    int nels = dimens.nels();
//    double *d_array_p = new double[nels];
//    const int di = dimens[0];
//    const int dj = dimens[1];
//    const int dk = dimens[2];
//    const int djk = dk * dj;
//    for (int i = 0; i < di; ++i)
//      for (int j = 0; j < dj; ++j)
//        for (int k = 0; k < dk; ++k) {
//
//          if (dk != 1)
//            d_array_p[i * djk + j * dk + k] = d_array[i + 1][j + 1][k + 1]; // I am not writing ghost cells
//          else if (dj != 1)
//            d_array_p[i * djk + j * dk] = d_array[i + 1][j + 1][0];
//          else {
//            d_array_p[i * djk + j * dk] = d_array[i + 1][0][0];
//
//          }
//        }
//    write(objname, dimens, d_array_p);
//    delete[]d_array_p;
//  }
//  catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(double*** array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & objname, const Dimens dimens, const int ns, double ****d_array) {
//  if (dimens.size() != 3) {
//    PSK::OutputException e("Dimens size not 3 for object " + objname, "HDF5OutputAdaptor::write(double**** array)");
//    throw e;
//  }
//
//  try {
//    int nels = dimens.nels();
//    double *d_array_p = new double[nels];
//    const int di = dimens[0];
//    const int dj = dimens[1];
//
//
//    const int dk = dimens[2];
//    const int djk = dk * dj;
//    for (int i = 0; i < di; ++i)
//      for (int j = 0; j < dj; ++j)
//        for (int k = 0; k < dk; ++k) {
//          if (dk != 1)
//            d_array_p[i * djk + j * dk + k] = d_array[ns][i + 1][j + 1][k + 1]; // I am not writing ghost cells
//          else if (dj != 1)
//            d_array_p[i * djk + j * dk] = d_array[ns][i + 1][j + 1][0];
//          else
//            d_array_p[i * djk + j * dk] = d_array[ns][i + 1][0][0];
//        }
//    write(objname, dimens, d_array_p);
//    delete[]d_array_p;
//  }
//  catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(double**** array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & objname, const Dimens dimens, double **d_array) {
//  if (dimens.size() != 2) {
//    PSK::OutputException e("Dimens size not 2 for object " + objname, "HDF5OutputAdaptor::write(double** array)");
//    throw e;
//  }
//
//  try {
//    int nels = dimens.nels();
//    double *d_array_p = new double[nels];
//    const int di = dimens[0];
//    const int dj = dimens[1];
//    for (int i = 0; i < di; ++i)
//      for (int j = 0; j < dj; ++j)
//        d_array_p[i * dj + j] = d_array[i + 1][j + 1];  // I am not writing ghost cells
//    write(objname, dimens, d_array_p);
//    delete[]d_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(double** array)");
//    throw e;
//  }
//}
//
//void HDF5OutputAdaptor::write(const std::string & objname, const Dimens dimens, const int ns, double ***d_array) {
//  if (dimens.size() != 2) {
//    PSK::OutputException e("Dimens size not 2 for object " + objname, "HDF5OutputAdaptor::write(double*** array)");
//    throw e;
//  }
//
//  try {
//    int nels = dimens.nels();
//    double *d_array_p = new double[nels];
//    const int di = dimens[0];
//    const int dj = dimens[1];
//    for (int i = 0; i < di; ++i)
//      for (int j = 0; j < dj; ++j)
//        d_array_p[i * dj + j] = d_array[i + 1][j + 1][ns];  // I am not writing ghost cells
//    write(objname, dimens, d_array_p);
//    delete[]d_array_p;
//  } catch(PSK::Exception & e) {
//    e.push("In HDF5OutputAdaptor::write(double*** array)");
//    throw e;
//  }
//}
