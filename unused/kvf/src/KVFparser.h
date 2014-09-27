// 
// KVFparser.h
// 
// declarations required by KVF parser and scanner code
// 
// David Burgess
// March 2000, September 2006
// 

#include "kvfDataSource.h"

namespace KVF {

  extern int KVFyylex();        // Defined in KVFDataSource.cc
  extern int KVF_line_no;       // Defined in KVFDataSource.cc 
  extern KVFDataSource *theParser;  // Defined in KVFDataSource.cc
} using namespace KVF;
