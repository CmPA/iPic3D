
#include "Restart3D.h"

/** write the restart file at any RESTART_CYCLE, useful for reading intermediate results */
void writeRESTART(string SaveDirName, int myrank, int cycle, int ns, MPIdata * mpi, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles3Dcomm * part) {
  // Create an Output Manager
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr;
  // Create an Output Agent for HDF5 output
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;
  hdf5_agent.set_simulation_pointers(field, grid, vct, mpi, col);
  for (int i = 0; i < ns; ++i)
    hdf5_agent.set_simulation_pointers_part(&part[i]);

  // Add the HDF5 output agent to the Output Manager's list
  output_mgr.push_back(&hdf5_agent);
  // Print Collective informations
  stringstream ss;
  ss << myrank;

  cout << SaveDirName + "/restart" + ss.str() + ".hdf" << endl;

  hdf5_agent.open_append(SaveDirName + "/restart" + ss.str() + ".hdf");
  output_mgr.output("proc_topology ", 0);
  output_mgr.output("Eall + Ball + rhos", cycle);
  output_mgr.output("position + velocity + q + ID ", cycle, 0);
//  output_mgr.output("position + velocity + q  ", cycle, 0);
  hdf5_agent.close();

}

/** this restart function writes the last restart with the last cycle */
void writeRESTART(string SaveDirName, int myrank, int cycle, int ns, MPIdata * mpi, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles3Dcomm * part, bool fool) {
  // Create an Output Manager
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr;
  // Create an Output Agent for HDF5 output
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;
  hdf5_agent.set_simulation_pointers(field, grid, vct, mpi, col);
  for (int i = 0; i < ns; ++i)
    hdf5_agent.set_simulation_pointers_part(&part[i]);

  // Add the HDF5 output agent to the Output Manager's list
  output_mgr.push_back(&hdf5_agent);
  // Print Collective informations
  stringstream ss;
  ss << myrank;

//  cout << SaveDirName + "/restart" + ss.str() + ".hdf" << endl;
//  cout << "son qui 1";
//  MPI_Barrier(MPI_COMM_WORLD);
  hdf5_agent.open(SaveDirName + "/restart" + ss.str() + ".hdf");
//  hdf5_agent.open_append(SaveDirName + "/restart" + ss.str() + ".hdf");
//  MPI_Barrier(MPI_COMM_WORLD); 
//  cout << "son qui 2";
//  MPI_Barrier(MPI_COMM_WORLD);
  output_mgr.output("proc_topology ", 0);
  output_mgr.output("Eall + Ball + rhos", 0);
  output_mgr.output("position + velocity + q + ID", 0, 0);
//  output_mgr.output("position + velocity + q ", 0, 0);
  output_mgr.output("last_cycle", cycle);
  hdf5_agent.close();

}


/** write the restart file at any RESTART_CYCLE, useful for reading intermediate results */
void writeRESTART_ES(string SaveDirName, int myrank, int cycle, int ns, MPIdata * mpi, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles * part) {
  // Create an Output Manager
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr;
  // Create an Output Agent for HDF5 output
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;
  hdf5_agent.set_simulation_pointers(field, grid, vct, mpi, col);
  for (int i = 0; i < ns; ++i)
    hdf5_agent.set_simulation_pointers_part(&part[i]);

  // Add the HDF5 output agent to the Output Manager's list
  output_mgr.push_back(&hdf5_agent);
  // Print Collective informations
  stringstream ss;
  ss << myrank;

  cout << SaveDirName + "/restart" + ss.str() + ".hdf" << endl;

  hdf5_agent.open_append(SaveDirName + "/restart" + ss.str() + ".hdf");
  output_mgr.output("proc_topology ", 0);
  output_mgr.output("Ex + rhos", cycle);
  output_mgr.output("x + u + q ", cycle, 0);
  hdf5_agent.close();

}

/** this restart function writes the last restart with the last cycle */
void writeRESTART_ES(string SaveDirName, int myrank, int cycle, int ns, MPIdata * mpi, VCtopology3D * vct, Collective * col, Grid * grid, Field * field, Particles * part, bool fool) {
  // Create an Output Manager
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr;
  // Create an Output Agent for HDF5 output
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;
  hdf5_agent.set_simulation_pointers(field, grid, vct, mpi, col);
  for (int i = 0; i < ns; ++i)
    hdf5_agent.set_simulation_pointers_part(&part[i]);

  // Add the HDF5 output agent to the Output Manager's list
  output_mgr.push_back(&hdf5_agent);
  // Print Collective informations
  stringstream ss;
  ss << myrank;
  hdf5_agent.open(SaveDirName + "/restart" + ss.str() + ".hdf");
  output_mgr.output("proc_topology ", 0);
  output_mgr.output("Ex + rhos", 0);
  output_mgr.output("x + u + q ", 0, 0);
  output_mgr.output("last_cycle", cycle);
  hdf5_agent.close();

}
