
#include <iostream>

#include "PSKOutput.h"

#include "PSKhdf5adaptor.h"



main()
{

try{

  int* iarr = new int[10];
  for(int i=0; i<10; ++i) iarr[i] = i+1;

  PSK::HDF5OutputAdaptor opa;
  
  opa.open( "test.hdf" );
  
  //opa.write( "/apples", 5 );
  
  //opa.write( "/a/b/c/d", -12345 );
  
// will fail
//  opa.write( "/a/b/c/d/e", 12345 );
  
// will fail
//  opa.write( "/a/b/c", -1 );
 
  //opa.write( "/apple_types/cox", 67 );

  //opa.write( "/apple_types/braeburn", 92 );
  //opa.write( "/apple_types/rotten/yellow", 1024 );
 
  //opa.write( "/citrus/oranges", PSK::Dimens(10), iarr );
  
  //opa.write( "/pears", PSK::Dimens(2,5), iarr );

  opa.close();
  
} catch ( PSK::Exception& e ) { e.diag_cout(); }

}
