
#include <iostream>
#include <string>
#include "kvfDataSource.h"

using namespace std;
using namespace KVF;

main()
{
  try{

  cout << "About to open file: KVFtestdata" << endl;

  KVFDataSource kv_file( "KVFtestdata" );

  if( kv_file.num_parse_errors() != 0 )
  {
    cout << "THERE WERE ERRORS READING THE FILE!! " << endl;
    kv_file.diag_print_errors() ;
    exit(1);
  }

  cout << " DIAG PRINT DATA " << endl;
  kv_file.diag_print_data();
  cout << " DIAG PRINT DATA AGAIN!" << endl;
  kv_file.diag_print_data();

{
  cout << " ******** TEST NAMES EXIST (OR NOT) ********** " << endl;

  string data_name = "apples";
  cout << "Data name <<" << data_name << ">> ";
  if( kv_file.name_exists( data_name ) )
    cout << "EXISTS!" << endl;
  else
    cout << "DOES NOT EXIST ;-(" << endl;

  data_name = "Lists";
  cout << "Data name <<" << data_name << ">> ";
  if( kv_file.name_exists( data_name ) )
    cout << "EXISTS!" << endl;
  else
    cout << "DOES NOT EXIST ;-(" << endl;

  data_name = "Lists/Sequence";
  cout << "Data name <<" << data_name << ">> ";
  if( kv_file.name_exists( data_name ) )
    cout << "EXISTS!" << endl;
  else
    cout << "DOES NOT EXIST ;-(" << endl;

  data_name = "No_such_name_me_old_china";
  cout << "Data name <<" << data_name << ">> ";
  if( kv_file.name_exists( data_name ) )
    cout << "EXISTS!" << endl;
  else
    cout << "DOES NOT EXIST ;-(" << endl;

}

{
  cout << " ******** TEST GETTING DATA NAMES  DataNamesTable_var ********** " << endl;

  DataNamesTable_var dnames = kv_file.get_data_names();

  cout << " ===== Test of recursive iterator ======" << endl;

  DataNamesTable::recursive_iterator itr;
  for( itr=dnames->recursive_begin();itr!=dnames->recursive_end();itr++)
      cout << itr.name() << "(depth: " << itr.depth() << ")" << endl;

  cout << " ===== Test of non-recursive iterator =====" << endl;
  {
  DataNamesTable::iterator itr;
  for( itr=dnames->begin();itr!=dnames->end();itr++)
   {   cout << itr.name();
      if( itr.type()==N_CONTEXT)
      { cout << " *NameContext*"  ;
        cout <<" is empty = " << itr.naming_context()->is_empty(); }
      
      cout << endl; }
  }
  
}

{
  cout << " ******** TEST GETTING DATA NAMES vector<string> ********** " << endl;

  vector<string> names_list;

  cout << " >>>>> All data names <<<<<< " << endl;
  kv_file.get_data_names( names_list );
  for( int i=0; i< names_list.size() ; ++i )
    cout << names_list[i] << endl;

  cout << " >>>>> Select data name with prefix and change separator <<<<<< " << endl;
  kv_file.get_data_names( names_list, "a", "::" );
  for( int i=0; i< names_list.size() ; ++i )
    cout << names_list[i] << endl;

} 

{
  cout << " ******** TEST GETTING DATA BY NAME: DATA BUFFERS  ********** " << endl;

  cout << " ===== Retrieve: z_value =====" << endl;
  DataBuffer_var dbv = kv_file.get_data( "z_value" );
  dbv->diag_print();

  cout << " ===== Retrieve: apples =====" << endl;
  dbv = kv_file.get_data( "apples" );
  dbv->diag_print();

  cout << " ===== Retrieve: Lists/Sequence =====" << endl;
  dbv = kv_file.get_data( "Lists/Sequence" );
  dbv->diag_print();

  cout << " ===== Retrieve: Random/values/lottery =====" << endl;
  dbv = kv_file.get_data( "Random/values/lottery" );
  dbv->diag_print();

  cout << " ===== Retrieve non-existent object: NOTANAME ===== " << endl;
  try{ dbv = kv_file.get_data( "NOTANAME" ); }
  catch ( KVFException& e )
  {
    cout << " Caught exception! " << endl;
    e.diag_cout();
  }
}



{
  cout << " ******** TEST GETTING DATA BY NAME:  EXPLICIT NUMERIC (DOUBLE)  ********** " << endl;

  cout << " ===== get_data(  .. double ) ===== " << endl;
  
  double dval;
  kv_file.get_data( "y_value", dval );
  cout << "Found value for <<y_value>> : " << dval << endl;

  cout <<
" ===== get_data(  .. double ) If sequence then first element returned===== "
       << endl;
  kv_file.get_data( "Random/values/lottery", dval );
  cout << "Found value for <<Random/values/lottery>> : " << dval << endl;
  cout <<
   " ===== get_data(  .. double ) .. Repeat operation, same result ===== "
       << endl;
  kv_file.get_data( "Random/values/lottery", dval );
  cout << "Found value for <<Random/values/lottery>> : " << dval << endl;

  vector<double> d_vvals;
  cout << " ===== get_data(  .. vec<double> ) ===== " << endl;
  kv_file.get_data( "Random/values/lottery", d_vvals );
  cout << "Found values for <<Random/values/lottery>> : " << endl;
  for(int i=0;i<d_vvals.size();++i) cout<< i << ":  " << d_vvals[i] << endl;



  cout << " ===== Attempt get_data(... double) on string data =====" << endl;
  try
  {
    kv_file.get_data( "apples", dval );
    cout << "Found value for <<apples>> : " << dval << endl;
  } catch ( KVFException& e)
  {
    cout << "Caught exception when getting numeric data with name <<apples>>"
         << endl;
    e.diag_cout();
  }
}

{
  cout << " ******** TEST GETTING DATA BY NAME:  EXPLICIT STRING  ********** " << endl;

  cout << " ===== get_data(  .. string ) ===== " << endl;
  
  string sval;
  kv_file.get_data( "apples", sval );
  cout << "Found value for <<apples>> : " << sval << endl;

  cout <<
" ===== get_data(  .. string ) If sequence then first element returned===== "
       << endl;
  kv_file.get_data( "Lists/Sequence", sval );
  cout << "Found value for <<Lists/Sequence>> : " << sval << endl;
  cout <<
   " ===== get_data(  .. string ) .. Repeat operation, same result ===== "
       << endl;
  kv_file.get_data( "Lists/Sequence", sval );
  cout << "Found value for <<Lists/Sequence>> : " << sval << endl;

  vector<string> s_vvals;
  cout << " ===== get_data(  .. vec<string> ) ===== " << endl;
  kv_file.get_data( "Lists/Sequence", s_vvals );
  cout << "Found values for <<Lists/Sequence>> : " << endl;
  for(int i=0;i<s_vvals.size();++i) cout<< i << ":  " << s_vvals[i] << endl;


  cout << " ===== Attempt get_data(... string) on numeric data =====" << endl;
  try
  {
    kv_file.get_data( "y_value", sval );
    cout << "Found value for <<y_value>> : " << sval << endl;
  } catch ( KVFException& e)
  {
    cout << "Caught exception when getting string data with name <<y_value>>"
         << endl;
    e.diag_cout();
  }
}


{
  cout << " ******** TEST GETTING DATA:  EXPLICIT TYPE, NO SUCH NAME  ********** "
       << endl;
  cout << " ===== Attempt get_data(... double) on no such name =====" << endl;

  try
  {
    double dval;
    kv_file.get_data( "NO_SUCH_NAME", dval );
    cout << "Found value for <<y_value>> : " << dval << endl;
  } catch ( KVFException& e)
  {
    cout << "Caught exception getting numeric data with name <<NO_SUCH_NAME>>"
         << endl;
    e.diag_cout();
  }

}
  cout <<  " ************** END OF TESTS ************************" << endl;


  } catch( KVFException& e )
  {
    cout << "tKVFDataSource.cc Caught exception " << endl;
    e.diag_cout();
  }

}

