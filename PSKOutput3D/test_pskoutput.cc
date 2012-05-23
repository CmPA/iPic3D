
#include <iostream>
#include <string>


#include "PSKOutput.h"
#include "PSKhdf5adaptor.h"

using namespace std;


template <class Toa> class myOutputAgent : public PSK::OutputAgent< Toa >
{

public:
  myOutputAgent( void ) {;}

  void output( const string& tag )
  {
  int* iarr = new int[10];
  for(int i=0; i<10; ++i) iarr[i] = i+1;
  double* darr = new double[6];
  for(int i=0; i<6; ++i) darr[i] = 100*(i+1)+i*2+0.01*i;

    cout << "myOutputAgent: hello from output()! tag: " << tag << "\n" ;
    
    cout << "Try outputting some stuff ...\n";

    if( tag == "Stage 1" )
    {
      this->output_adaptor.write( "ivalue", 10 );
      this->output_adaptor.write( "ivals_array", PSK::Dimens(10), iarr );

      float fval=123.4;
      this->output_adaptor.write( "/stage1/floats/fvalue", fval );
      double dval=9.87654321;
      this->output_adaptor.write( "/stage1/doubles/dvalue", dval );
      this->output_adaptor.write( "/stage1/doubles/dvalarray", PSK::Dimens(2,3), darr );

      
      
    } else if( tag == "Stage 2" )
    {
	    this->output_adaptor.write( "ivalue_stage2", -9999 );
      
    }
  }
};



main(){

// create the output manager which will keep track of different
//   output agents

  PSK::OutputManager< PSK::OutputAdaptor > output_mgr;

// make two instances of class myOutputAgent, one which uses
//   the HDF5 output adaptor, and another which uses a stream to std::cout

  myOutputAgent< PSK::HDF5OutputAdaptor > my_output_agent;
  
  myOutputAgent< PSK::coutOutputAdaptor > my_cout_output_agent;
  
// Add the HDF5 output agent to the Output Manager's list
  output_mgr.push_back( &my_output_agent );

// do the same again for the cout output agent
  output_mgr.push_back( &my_cout_output_agent );

// open the hdf5 file
  my_output_agent.open( "agent.hdf");

// no need to open the cout output agent

// run the actions ...

  output_mgr.output( "Stage 1" );
 
  output_mgr.output( "Stage 2" );
  
  
// close the hdf5 file  
  my_output_agent.close();

}
