// 
// DataBuffer.h
// 
// David Burgess
// June 1999, June 2004, September 2006
// 

#ifndef KVF_DATABUFFER_H
#define KVF_DATABUFFER_H

#include <vector>
#include <string>
#include <cmath>

#include "kvfVar.h"
#include "kvfNaming.h"

namespace KVF {

  using namespace std;


  typedef enum BufferedDataMajorityE { DBUFFER_NIL_MAJOR, DBUFFER_ROW_MAJOR, DBUFFER_COL_MAJOR } BufferedDataMajority;

  /* ! \brief Description of data dimensionality in a SimpleDataBuffer.
   * 
   */

  class SimpleDBuffDescriptor:public VarTarget {
    int _num_item_elts;         // /< number of elements per item
    int _num_items;             // /< number of items
      vector < int >_dims;      // /< matrix dimensionality and size
    BufferedDataMajority _majority; // /, matrix majority
  public:
      SimpleDBuffDescriptor()
    : _num_items(0), _num_item_elts(1), _majority(DBUFFER_NIL_MAJOR) {;
    } SimpleDBuffDescriptor(const SimpleDBuffDescriptor & a)  // /< Copy ctor
    : _num_items(0), _majority(a._majority) {
      set_dims(a._dims);
    } void clear(void) {
      _num_items = 0;
      _dims.clear();
      _num_item_elts = 1;
    }

    void set_dims(const vector < int >&d);  // /< Set matrix dimensionality and size

    void set_num_items(int n)   // /< Set number of items
    {
      _num_items = n;
    }

    void set_majority(BufferedDataMajority m) // /< Set majority of matrix data.
    {
      _majority = m;
    }

    void calc_num_items(int n_elts) // /< Calculate and set number of items.
    {
      _num_items = n_elts / _num_item_elts;
    }

    void get_dims(vector < int >&d) const // /< Return matrix dimensionality.
    {
      d = _dims;
    } int num_item_elts(void) const // /< Return number of elements per item
    {
      return _num_item_elts;
    } int dims_size(void) const // /< Return number of matrix dimensions
    {
      return _dims.size();
    } int num_items(void) const // /< Return number of items
    {
      return _num_items;
    } BufferedDataMajority majority(void) // /< Return matrix majority
    {
      return _majority;
    }

    void diag_cout(void) const; // /< Diagnostic print to cout

  };

  typedef Var < SimpleDBuffDescriptor > SimpleDBuffDescriptor_var;

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ class DataBuffer
  // ----------------------------------------------------------------------------

  /* ! \brief Abstract interface class for a collection of data values
   * 
   * This is a root abstract class for the concept of a DataBuffer.
   * 
   * A DataBuffer contains or acquires data, and can be queried for that data.
   * 
   * Subclasses such as SimpleDataBuffer or TSDataBuffer have more specific behaviour.
   * 
   */
  class DataBuffer:public VarTarget {
  public:
    virtual ~ DataBuffer() {;
    } virtual string get_type_srep(void) const = 0;
    virtual void diag_print(void) const = 0;
  };

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++ class SimpleDataBuffer
  // ----------------------------------------------------------------------------

  /* ! \brief Abstract class for a DataBuffer consisting of a sequence of rectangular matrices
   * 
   * In a SimpleDataBuffer the data is organized as an ordered sequence of items. Each item is a rectangular matrix of elements, where the dimensionality, dimension sizes and majority of all matrix items is the same. The element values are accessed in sequential order ordered by item, and then element order within each item.
   * 
   * A SimpleDataBuffer object can be queried for - number of items it contains - number of elements in each item - total number of elements - the description of the data arrangement
   * 
   * A SimpleDataBuffer object can also - be rewound, ie subsequent sequential data access starts from beginning - be cleared, ie all data erased
   * 
   */

  class SimpleDataBuffer:public DataBuffer {
  public:
    virtual SimpleDBuffDescriptor_var descriptor(void) = 0; // /< return data description
    virtual void rewind(void) = 0;  // /< Rewind for sequential access
    virtual void clear(void) = 0; // /< Erase all data
    virtual int num_items(void) const = 0;  // /< Return number of items available
    virtual int num_elts(void) const = 0; // /< Return total number of elements
    virtual int num_item_elts(void) const = 0;  // /< Return number of elements per item
  };

  // ++++++++++++++++++++++++++++++++++++++++++++++++++++ class NumericDataBuffer
  // ----------------------------------------------------------------------------

  /* ! \brief Abstract class for numeric SimpleDataBuffer
   * 
   */

  class NumericDataBuffer:public SimpleDataBuffer {
  public:
    virtual bool get_data(int n, vector < double >&db, bool append = false) = 0;
    // /< Extract n elements of data
  };

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++ class StringDataBuffer
  // ----------------------------------------------------------------------------
  class StringDataBuffer:public SimpleDataBuffer {
  public:
    virtual bool get_data(int n, vector < string > &db, bool append = false) = 0;
  };


  typedef Var < DataBuffer > DataBuffer_var;
  typedef Var < NumericDataBuffer > NumericDataBuffer_var;
  typedef Var < StringDataBuffer > StringDataBuffer_var;

  // ++++++++++++++++++++++++++++++++++++++++++++++++++ class StringPtrDataBuffer
  // ----------------------------------------------------------------------------
  class DataPtrBuffer:public DataBuffer {
  public:
    DataBuffer_var data;
    StringDataBuffer_var data_ptr_name;
    string get_type_srep(void) const {
      return "DataPtrBuffer";
    } void diag_print(void) const {
      data->diag_print();
      data_ptr_name->diag_print();
    };

  };

  typedef Var < DataPtrBuffer > DataPtrBuffer_var;

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ class TSDataBuffer
  // ----------------------------------------------------------------------------


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  /* ! \brief Template class for vector based SimpleDataBuffer
   * 
   * Data elements are stored in type vector<T>.
   * 
   */
  template < class T > class BasicVectorDataBuffer {
  public:
    typedef T value_type;
  protected:
    SimpleDBuffDescriptor_var _desc;  // /< data description 
    vector < T > _v;            // /< Vector of values
    unsigned int _pos;          // /< current position used for accessing elements
  public:
    BasicVectorDataBuffer(void):_pos(0), _desc(new SimpleDBuffDescriptor()) {;
    }
    BasicVectorDataBuffer(const SimpleDBuffDescriptor & dbdesc)
    : _pos(0), _desc(new SimpleDBuffDescriptor(dbdesc)) {;
    }

    // ~BasicVectorDataBuffer(void) { cout << "~BasicVectorDataBuffer"<<endl;}

    /* ! \brief Retrieve data from data buffer
     * 
     * \param n_elts Number of elements to retrieve. \param vvals Reference to vector, in which to store values from data buffer \param append If true, then values are appended to end of vvals, otherwise vvals is cleared before values are copied.
     * 
     * \returns Boolean true if there are still values remaining in the buffer.
     * 
     * \note The position state of the buffer is used and updated. */
    bool get_data(int n_elts, vector < T > &vvals, bool append = false) {
      if (!append)
        vvals.clear();          // clear the receiving vector if append not requested
      if (_pos >= _v.size())
        return false;           // nothing to return
      int i = 0;                // number of items packed
      while (i < n_elts && _pos < _v.size()) {
        vvals.push_back(_v[_pos++]);
        i++;
      }
      return _pos < _v.size();  // return true if more to come
    }

    /* ! \brief Return data buffer description */
    SimpleDBuffDescriptor_var descriptor(void) {
      return _desc;
    }

    /* ! \brief Rewind buffer element position to start */
    void rewind(void) {
      _pos = 0;
    }

    /* ! \brief Erase all data in buffer */
    void clear(void) {
      _v.clear();
      _pos = 0;
      _desc->clear();
    }


    /* ! \brief Push supplied value as new element at end of buffer */
    void push_back(T x) {
      _v.push_back(x);
      _desc->calc_num_items(_v.size());
    }

    /* ! \brief Push supplied values as new elements at end of buffer \note The number of items is recalculated using new total number of elements. */
    void push_back(const vector < T > &x) {
      for (int i = 0; i < x.size(); i++)
        _v.push_back(x[i]);
      _desc->calc_num_items(_v.size());
    }

    /* ! \brief Push elements of supplied data buffer as new elements at end of buffer. \note The dimensionality etc. of the matrices are not used. The elements of the supplied data buffer are simply copied */
    void push_back(const BasicVectorDataBuffer < T > &x) {
      for (int i = 0; i < x._v.size(); i++)
        _v.push_back(x._v[i]);
      _desc->calc_num_items(_v.size());
    }

    /* ! \brief Reserve space of n elements for the internal vector */
    void reserve(size_t n) {
      _v.reserve(n);
    }

    /* ! \brief Set the matrix dimensionality and dimension sizes */
    void set_dims(const vector < int >&d) {
      _desc->set_dims(d);
    }

    /* ! \brief Set the matrix majority */
    void set_majority(BufferedDataMajority m) {
      _desc->set_majority(m);
    }

    /* ! \brief Return total number of elements in buffer */
    int num_elts(void) const {
      return _v.size();
    }
    /* ! \brief Return total number of elements in buffer */ int size(void) {
      return _v.size();
    }

    /* ! \brief Return number of items in buffer */
    int num_items(void) const {
      return _desc->num_items();
    }
    /* ! \brief Return number of elements per item. */ int num_item_elts(void) const {
      return _desc->num_item_elts();
    }
    /* ! \brief Return copy of i-th value in buffer. */ T operator[] (int i) const {
      return _v[i];
    }
    /* ! \brief Test item (not just element) is approximately equal. \param item_i Index of item within buffer to be tested. \param a Vector of values corresponding to the elements of the item. If there are insufficient elements in \c a, then it is used cyclically. \param precision Bound for relative precision.
     * 
     * \returns Boolean true if all elements of item \c item_i are approximately equal to corresponding elements of \c a, within given relative precision. */ bool item_approx_equal(int item_i, const vector < T > &a,
                                                                                                                                                                                     double precision) const {
      bool ok = true;
      int elt0 = item_i * _desc->num_item_elts();
      int elt1 = elt0 + _desc->num_item_elts();
      for (int i = elt0; i < elt1; i++) {
        int j = i % a.size();
        if (fabs(_v[i] - a[j]) > fabs(_v[i] * precision)) {
          ok = false;
          break;
      }} return ok;
    }

    /* ! \brief Delete i-th element in buffer. */
    void delete_elt(int i) {
      _v.erase(_v.start() + i * (_desc->num_item_elts()), _v.start() + (i + 1) * (_desc->num_item_elts()));
    }

  };


  // +++++++++++++++++++++++++++++++++++++++++++++++++++ class NumericVDataBuffer
  // ----------------------------------------------------------------------------
  class NumericVDataBuffer:public BasicVectorDataBuffer < double >, public NumericDataBuffer {

  public:
    NumericVDataBuffer(void) {;
    } NumericVDataBuffer(const SimpleDBuffDescriptor & dbdesc)
    : BasicVectorDataBuffer < double >(dbdesc) {;
    }

    SimpleDBuffDescriptor_var descriptor(void) {
      return BasicVectorDataBuffer < double >::descriptor();
    }
    void rewind(void) {
      BasicVectorDataBuffer < double >::rewind();
    }
    void clear(void) {
      BasicVectorDataBuffer < double >::clear();
    }

    bool get_data(int n, vector < double >&db, bool append);

    string get_type_srep(void) const {
      return "NumericVDataBuffer";
    } void diag_print(void) const;

    int num_items(void) const;
    int num_elts(void) const;
    int num_item_elts(void) const;
  };


  // ++++++++++++++++++++++++++++++++++++++++++++++++++++ class StringVDataBuffer
  // ----------------------------------------------------------------------------
  class StringVDataBuffer:public BasicVectorDataBuffer < string >, public StringDataBuffer {

  public:
    StringVDataBuffer(void) {;
    } StringVDataBuffer(const SimpleDBuffDescriptor & dbdesc)
    : BasicVectorDataBuffer < string > (dbdesc) {;
    }

    SimpleDBuffDescriptor_var descriptor(void);
    void rewind(void);
    void clear(void);
    bool get_data(int n, vector < string > &db, bool append);
    string get_type_srep(void) const;
    void diag_print(void) const;
    int num_items(void) const;
    int num_elts(void) const;
    int num_item_elts(void) const;
  };


  typedef Var < NumericVDataBuffer > NumericVDataBuffer_var;
  typedef Var < StringVDataBuffer > StringVDataBuffer_var;

}                               // end namespace KVF


#endif // #ifndef KVF_DATABUFFER_H
