//
//  KVFDataSource.cc
//
//  Key Value File Data Source
//
//  David Burgess
//  March 2000, June 2004, September 2006
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "kvfDataSource.h"

int KVFyyparse();                   // defined in KVFyacc.c
#define yyFlexLexer KVFyyFlexLexer
#include <FlexLexer.h>             // must be in /usr/local/include or similar

namespace KVF {

using namespace std;

int KVF_line_no;
static KVFyyFlexLexer *lexer;
KVFDataSource *theParser;

int KVFyylex()
{
  return lexer->yylex();
}

// current constructor strategy is: open file, read all data in, close file

/*! \brief Open KVF file, read all data into object, and close file
\param file_name Name of KVF file to be opened, read, and closed.

\throws KVFException, from open() call.

Parse errors are accumulated, and should be checked after object
construction.
*/
KVFDataSource::KVFDataSource(  const string& file_name )
{
  read_file( file_name );
}


KVFDataSource::~KVFDataSource( void )
{
}


void KVFDataSource::open( const string& file_name )
{
  _file_name = file_name;
  _istr.open( _file_name.c_str() );
  if( !_istr )
    throw KVFException( 
      " KVFDataSource::open: failed to open: " + _file_name );
}

/*! \brief Open, read data from, and then close a KVF file.
\param file_name Name of KVF file to be opened, read, and closed.

\throws KVFException, from open() call.

Parse errors are accumulated, and should be checked after object
construction.
*/
void KVFDataSource::read_file( const string& file_name )
{
  open( file_name );
  lexer = new KVFyyFlexLexer( &_istr );
  _parse_state = PARSE_NIL;
  KVF_line_no = 0;
  parse();        // read in the whole file now!
  close();
}

void KVFDataSource::close()
{
  delete lexer;   // no need for lexer any more
  _istr.close();   // close the file
}

/*! \brief Parse the KVF file and read/store the data.

\note Thread unsafe - KVFyyFlexLexer* lexer is static
*/

void KVFDataSource::parse(void)
{
  theParser = this;
  KVF_line_no = 0;
  KVFyyparse();
}

void KVFDataSource::set_parse_error( string s, int line_no )
{
  ostringstream ost;
  ost << "ERROR: KVFDataSource at line " << line_no << " (" << s << ")";
  _parse_errors.push_back(ost.str());
}

bool KVFDataSource::key_is_ok( string key )
{
  try
  { DataBuffer_var dbv = _data_by_key.resolve( key ); }
  catch( KVFNameNotFound& not_found )
  {
    return true;
  }
  _parse_errors.push_back(
      "ERROR: KVFDataSource: duplicated key ignored: " + key );
  return false;
}

void KVFDataSource::assign_number( const char* key, double x )
{
  if( !key_is_ok( key ) )
    return;
  NumericVDataBuffer* nvdb_p = new NumericVDataBuffer();
  nvdb_p->push_back( x );
  _data_by_key.bindx( key, nvdb_p );
  _parse_state = PARSE_NIL;
}

void KVFDataSource::assign_string( const char* key, const char* s )
{
  if( !key_is_ok( key ) )
    return;
  StringVDataBuffer* svdb_p = new StringVDataBuffer();
  svdb_p->push_back( s );
  _data_by_key.bindx( key, svdb_p );
  _parse_state = PARSE_NIL;
}

void KVFDataSource::start_nseq( double x )
{
  _nseq_buff.clear();
  _nseq_buff.push_back(x);
  _parse_state=PARSE_NSEQ;
}

void KVFDataSource::push_nseq( double x )
{
  _nseq_buff.push_back(x);
}

void KVFDataSource::assign_nseq( const char* key )
{
  if( !key_is_ok( key ) )
    return;
  NumericVDataBuffer* nvdb_p = new NumericVDataBuffer();
  reverse( _nseq_buff.begin(), _nseq_buff.end() );
  nvdb_p->push_back( _nseq_buff );
  _nseq_buff.clear();
  _data_by_key.bindx( key, nvdb_p );
  _parse_state = PARSE_NIL;
}
void KVFDataSource::start_strseq( const char *s )
{
  _strseq_buff.clear();
  _strseq_buff.push_back(s);
  _parse_state=PARSE_STRSEQ;
}

void KVFDataSource::push_strseq( const char *s )
{
  _strseq_buff.push_back(s);
}

void KVFDataSource::assign_strseq( const char* key )
{
  if( !key_is_ok( key ) )
    return;
  StringVDataBuffer* svdb_p = new StringVDataBuffer();
  reverse( _strseq_buff.begin(), _strseq_buff.end()  );
  svdb_p->push_back( _strseq_buff );
  _strseq_buff.clear();
  _data_by_key.bindx( key, svdb_p );
  _parse_state = PARSE_NIL;
}

void KVFDataSource::diag_print_data(void)
{
  NamingContext<DataBuffer_var>::recursive_iterator itr;
  BindingType bindtype;
  for( itr=_data_by_key.recursive_begin();
       itr!=_data_by_key.recursive_end(); itr++ )
  {
    cout << itr.name();
    bindtype=itr.type();
    if(bindtype==N_OBJECT)
    {
      cout << " = " << "** DataBuffer_var **" << endl;
      cout << "DataBuffer type: " << itr->second._obj->get_type_srep() << endl;
      itr->second._obj->diag_print();
    }
    else if(bindtype==N_CONTEXT)
     cout << "/" << endl;
   }
}

/*! \brief Return DataNamesTable of all data names.

Take care! The internal object DataNamesTable is only constructed
when this routine is first called. Subsequent calls do not
refresh the names table from the _data_by_key.

\todo Update needed when refresh() function added.

*/

DataNamesTable_var KVFDataSource::get_data_names(void)
{
  if( _data_names.is_nil() )
  {
    _data_names = new DataNamesTable;
    NamingContext<DataBuffer_var>::recursive_iterator itr;
    for( itr=_data_by_key.recursive_begin();
         itr!=_data_by_key.recursive_end(); itr++)
      if(itr.type()==N_OBJECT)
      {
        _data_names->bindx( itr.name(), "" );
      }
  }
  return _data_names;
}

/*! \brief Return DataBuffer var for data associated with supplied name.

\returns DataBuffer_var for data associated with given data name. Note
that the reference is to the data buffer associated with the KVFDataSource,
NOT a copy. Thus any changes to this data buffer will also affect the
KVFDatasource.

\throws KVFException if there is not such data name available.
*/
DataBuffer_var KVFDataSource::get_data( const string& data_name )
{
  DataBuffer_var dbv;
  try
  { dbv = _data_by_key.resolve( data_name ); }
  catch( KVFNameNotFound& not_found )
  {
    throw KVFException( 
     "KVFDataSource::get_data (file: " + _file_name +
     ") data name not found: " + data_name );
  }
  return dbv;
}

bool KVFDataSource::name_exists( const string& data_name )
{
  bool result = true;
  try
  { DataBuffer_var dbv = _data_by_key.resolve( data_name ); }
  catch( KVFNamingException& not_found )
  {
    result = false;
  }
  return result;
}


void KVFDataSource::get_data_names( vector<string>& names_list,
          const string& prefix, const string& separator)
{
  names_list.clear();
  NamingContext<DataBuffer_var>::recursive_iterator itr;
  string data_name;
  int prefix_length = prefix.length();
  for( itr=_data_by_key.recursive_begin();
       itr!=_data_by_key.recursive_end(); itr++)
  {
    if( itr.type() == N_OBJECT )
    {
      data_name = itr.name( separator );
      if( prefix_length == 0 )
        names_list.push_back( data_name );
      else if( data_name.substr( 0, prefix_length ) == prefix )
        names_list.push_back( data_name );
    }
  }

}

void KVFDataSource::get_data( const string& data_name, double& d_val )
{
  DataBuffer_var dbv;
  try{
    dbv = get_data( data_name );
  
    NumericDataBuffer_var n_dbv = NumericDataBuffer_var::narrow( dbv );

    if( !n_dbv.is_nil() )
    {
      n_dbv->rewind();
      if( n_dbv->num_elts() > 0 )
      {
        vector<double> d_vvals;
        n_dbv->get_data( 1, d_vvals );
        d_val = d_vvals[0];
      } else
      {
        throw KVFException(
          "NumericDataBuffer with no data elements; name <<" + data_name + ">>",
          " KVFDataSource::get_data( ... double)" );
      }
    } else
    {
      throw KVFException( 
        "Unexpected type of data (expected numeric); name <<"+data_name + ">>",
        " KVFDataSource::get_data( ... double)" );
    }
  
  } catch ( KVFException& e )
  {
    e.push_err_msg( " In function KVFDataSource::get_data( ... double)" );
    throw e;
  }
}


void KVFDataSource::get_data( const string& data_name, vector<double>& d_vvals )
{
  DataBuffer_var dbv;
  try{
    dbv = get_data( data_name );
  
    NumericDataBuffer_var n_dbv = NumericDataBuffer_var::narrow( dbv );
    
    if( !n_dbv.is_nil() )
    {
      n_dbv->rewind();
      if( n_dbv->num_elts() > 0 )
      {
        n_dbv->get_data( n_dbv->num_elts(), d_vvals );
      } else
      {
        throw KVFException(
          "NumericDataBuffer with no data elements; name <<" + data_name + ">>",
          " KVFDataSource::get_data( ... vec<double>)" );
      }
    } else
    {
      throw KVFException( 
        "Unexpected type of data (expected numeric); name <<"+data_name + ">>",
        " KVFDataSource::get_data( ... vec<double>)" );
    }
  
  } catch ( KVFException& e )
  {
    e.push_err_msg( " In function KVFDataSource::get_data( ... vec<double>)" );
    throw e;
  }
}

void KVFDataSource::get_data( const string& data_name, string& s_val )
{
  DataBuffer_var dbv;
  try{
    dbv = get_data( data_name );
  
    StringDataBuffer_var str_dbv = StringDataBuffer_var::narrow( dbv );

    if( !str_dbv.is_nil() )
    {
      str_dbv->rewind();
      if( str_dbv->num_elts() > 0 )
      {
        vector<string> str_vvals;
        str_dbv->get_data( 1, str_vvals );
        s_val = str_vvals[0];
      } else
      {
        throw KVFException(
          "StringDataBuffer_var with no data elements; name <<" + data_name + ">>",
          " KVFDataSource::get_data( ... strings)" );
      }
    } else
    {
      throw KVFException( 
        "Unexpected type of data (expected string); name <<"+data_name + ">>",
        " KVFDataSource::get_data( ... string)" );
    }
  
  } catch ( KVFException& e )
  {
    e.push_err_msg( " In function KVFDataSource::get_data( ... string)" );
    throw e;
  }
}

void KVFDataSource::get_data( const string& data_name, vector<string>& s_vvals )
{
  DataBuffer_var dbv;
  try{
    dbv = get_data( data_name );
  
    StringDataBuffer_var s_dbv = StringDataBuffer_var::narrow( dbv );
    
    if( !s_dbv.is_nil() )
    {
      s_dbv->rewind();
      if( s_dbv->num_elts() > 0 )
      {
        s_dbv->get_data( s_dbv->num_elts(), s_vvals );
      } else
      {
        throw KVFException(
          "StringDataBuffer_var with no data elements; name <<" + data_name + ">>",
          " KVFDataSource::get_data( ... vec<string>)" );
      }
    } else
    {
      throw KVFException( 
        "Unexpected type of data (expected string); name <<"+data_name + ">>",
        " KVFDataSource::get_data( ... vec<string>)" );
    }
  
  } catch ( KVFException& e )
  {
    e.push_err_msg( " In function KVFDataSource::get_data( ... vec<string>)" );
    throw e;
  }
}

/* \brief Attempt to retrieve integer data associated with supplied data name

At moment just integer truncates

\todo Correctly handle conversion double to integer

*/

void KVFDataSource::get_data( const string& data_name, int& i_val )
{
  DataBuffer_var dbv;
  try{
    dbv = get_data( data_name );
  
    NumericDataBuffer_var n_dbv = NumericDataBuffer_var::narrow( dbv );

    if( !n_dbv.is_nil() )
    {
      n_dbv->rewind();
      if( n_dbv->num_elts() > 0 )
      {
        vector<double> d_vvals;
        n_dbv->get_data( 1, d_vvals );

// more work needed ...
        i_val = static_cast<int>(d_vvals[0]);

      } else
      {
        throw KVFException(
          "NumericDataBuffer with no data elements; name <<" + data_name + ">>",
          " KVFDataSource::get_data( ... int)" );
      }
    } else
    {
      throw KVFException( 
        "Unexpected type of data (expected numeric); name <<"+data_name + ">>",
        " KVFDataSource::get_data( ... int)" );
    }
  
  } catch ( KVFException& e )
  {
    e.push_err_msg( " In function KVFDataSource::get_data( ... int)" );
    throw e;
  }
}


/* \brief Attempt to retrieve vector of 
integer data associated with supplied data name

At moment just integer truncates

\todo Correctly handle conversion double to integer

*/
void KVFDataSource::get_data( const string& data_name, vector<int>& i_vvals )
{
  DataBuffer_var dbv;
  try{
    dbv = get_data( data_name );
  
    NumericDataBuffer_var n_dbv = NumericDataBuffer_var::narrow( dbv );
    
    if( !n_dbv.is_nil() )
    {
      n_dbv->rewind();
      if( n_dbv->num_elts() > 0 )
      {
        vector<double> d_vvals;
        n_dbv->get_data( n_dbv->num_elts(), d_vvals );

// more work needed here ...
        i_vvals.resize( d_vvals.size() );
        for( int i=0; i < d_vvals.size(); ++i )
          i_vvals[i] = static_cast<int>( d_vvals[i] );
        
      } else
      {
        throw KVFException(
          "NumericDataBuffer with no data elements; name <<" + data_name + ">>",
          " KVFDataSource::get_data( ... vec<int>)" );
      }
    } else
    {
      throw KVFException( 
        "Unexpected type of data (expected numeric); name <<"+data_name + ">>",
        " KVFDataSource::get_data( ... vec<int>)" );
    }
  
  } catch ( KVFException& e )
  {
    e.push_err_msg( " In function KVFDataSource::get_data( ... vec<int>)" );
    throw e;
  }
}
/* \brief Attempt to retrieve integer data associated with supplied data name

At moment just integer truncates

\todo Correctly handle conversion double to integer

*/

void KVFDataSource::get_data( const string& data_name, float& f_val )
{
  DataBuffer_var dbv;
  try{
    dbv = get_data( data_name );
  
    NumericDataBuffer_var n_dbv = NumericDataBuffer_var::narrow( dbv );

    if( !n_dbv.is_nil() )
    {
      n_dbv->rewind();
      if( n_dbv->num_elts() > 0 )
      {
        vector<double> d_vvals;
        n_dbv->get_data( 1, d_vvals );

// more work needed ...
        f_val = static_cast<float>(d_vvals[0]);

      } else
      {
        throw KVFException(
          "NumericDataBuffer with no data elements; name <<" + data_name + ">>",
          " KVFDataSource::get_data( ... float)" );
      }
    } else
    {
      throw KVFException( 
        "Unexpected type of data (expected numeric); name <<"+data_name + ">>",
        " KVFDataSource::get_data( ... float)" );
    }
  
  } catch ( KVFException& e )
  {
    e.push_err_msg( " In function KVFDataSource::get_data( ... float)" );
    throw e;
  }
}


/* \brief Attempt to retrieve vector of 
float data associated with supplied data name

At moment just uses static cast to float

\todo Correctly handle conversion double to float

*/
void KVFDataSource::get_data( const string& data_name, vector<float>& f_vvals )
{
  DataBuffer_var dbv;
  try{
    dbv = get_data( data_name );
  
    NumericDataBuffer_var n_dbv = NumericDataBuffer_var::narrow( dbv );
    
    if( !n_dbv.is_nil() )
    {
      n_dbv->rewind();
      if( n_dbv->num_elts() > 0 )
      {
        vector<double> d_vvals;
        n_dbv->get_data( n_dbv->num_elts(), d_vvals );

// more work needed here ...
        f_vvals.resize( d_vvals.size() );
        for( int i=0; i < d_vvals.size(); ++i )
          f_vvals[i] = static_cast<float>( d_vvals[i] );
        
      } else
      {
        throw KVFException(
          "NumericDataBuffer with no data elements; name <<" + data_name + ">>",
          " KVFDataSource::get_data( ... vec<float>)" );
      }
    } else
    {
      throw KVFException( 
        "Unexpected type of data (expected numeric); name <<"+data_name + ">>",
        " KVFDataSource::get_data( ... vec<float>)" );
    }
  
  } catch ( KVFException& e )
  {
    e.push_err_msg( " In function KVFDataSource::get_data( ... vec<float>)" );
    throw e;
  }
}


} // end namespace KVF

