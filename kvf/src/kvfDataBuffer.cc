//
//  DataBuffer.cc
//
//  Qdos2 - data interfacing
//
//  David Burgess
//  July 1999, June 2004, September 2006
//

#include <iostream>

#include "kvfDataBuffer.h"

namespace KVF {

/*!



*/
// update dimension info from vector<int>
// if number of dimensions is zero, then assume scalar (one element per item)
// if one or more dimensions are zero then also assume scalar!
// otherwise number of elements per item is product of all dimension sizes

void SimpleDBuffDescriptor::set_dims( const vector<int>& d )
{
  _dims=d;
  int num_dims=_dims.size();
  if(num_dims==0)
    _num_item_elts=1;
  else
  {
    _num_item_elts=1;
    for(int i=0; i<num_dims; ++i )
      _num_item_elts *= _dims[i];
    if( _num_item_elts == 0 )
      _num_item_elts = 1;
  }
}


void SimpleDBuffDescriptor::diag_cout(void) const
{
  cout << "SimpleDBuffDescriptor: ";
  cout << "(" << _num_items << ")" ;
  cout << " " << _dims.size() ;
  cout << "[";
  for(int i=0;i<_dims.size();i++)
    if(i!=_dims.size()-1)
      cout << _dims[i] << ",";
    else
      cout << _dims[i];
  cout << "]" ;
}

// +++++++++++++++++++++++++++++++++++++++++ class NumericVDataBuffer
// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// ----------------------------------------------------- XXX
void NumericVDataBuffer::diag_print(void) const
{ cout << "NumericVDataBuffer: num_elems="<< num_elts()
       <<" num_items=" << num_items() << " ";
       _desc->diag_cout(); cout << endl;
  if( num_elts()!=0)
  { for(int i=0;i<num_elts()-1;i++)
      cout << _v[i]<<", ";
    cout << _v[num_elts()-1] << endl;
  }
}

int NumericVDataBuffer::num_items(void) const
{ return BasicVectorDataBuffer<double>::num_items();}
    
int NumericVDataBuffer::num_elts(void) const
{ return BasicVectorDataBuffer<double>::num_elts(); }

int NumericVDataBuffer::num_item_elts(void) const
{ return BasicVectorDataBuffer<double>::num_item_elts(); }  

bool NumericVDataBuffer::get_data( int n, vector<double>& db, bool append )
{ return BasicVectorDataBuffer<double>::get_data( n, db, append ); }

// +++++++++++++++++++++++++++++++++++++++++ class StringVDataBuffer
// ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// ----------------------------------------------------- XXX

SimpleDBuffDescriptor_var StringVDataBuffer::descriptor(void)
{ return BasicVectorDataBuffer<string>::descriptor(); }

void StringVDataBuffer::rewind(void)
{ BasicVectorDataBuffer<string>::rewind(); }

void StringVDataBuffer::clear(void)
{ BasicVectorDataBuffer<string>::clear(); }

bool StringVDataBuffer::get_data( int n, vector<string>& db, bool append )
{ return BasicVectorDataBuffer<string>::get_data( n, db, append ); }

string StringVDataBuffer::get_type_srep(void) const
{ return "StringVDataBuffer"; }

void StringVDataBuffer::diag_print(void) const
{ cout << "StringVDataBuffer: num_elems=" << num_elts()
           <<" num_items=" << num_items()<< " ";
           _desc->diag_cout(); cout << endl;
      for(int i=0;i<num_elts();i++)
       cout <<"\"" << _v[i].c_str() << "\""<< endl;
}

int StringVDataBuffer::num_items(void) const
{ return BasicVectorDataBuffer<string>::num_items();}

int StringVDataBuffer::num_elts(void) const 
{ return BasicVectorDataBuffer<string>::num_elts();}  

int StringVDataBuffer::num_item_elts(void) const
{ return BasicVectorDataBuffer<string>::num_item_elts(); }  


} // namespace KVF
