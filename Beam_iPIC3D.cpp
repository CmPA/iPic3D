/***************************************************************************
  Parsek3D.cpp  -  Main file for 3D sim
  -------------------
begin                : Jun 2008
copyright            : (C) 2004 by Stefano Markidis, Giovanni Lapenta
 ************************************************************************** */

// MPI
#include <mpi.h>
#include "mpidata/MPIdata.h"
// Topology
#include "processtopology/VirtualTopology3D.h"
#include "processtopology/VCtopology3D.h"
// input
#include "inputoutput/CollectiveIO.h"
#include "inputoutput/Collective.h"
// grid
#include "grids/Grid.h"
#include "grids/Grid3DCU.h"
// fields
#include "fields/Field.h"
#include "fields/EMfields3D.h"
// particles
#include "particles/Particles.h"
#include "particles/Particles3Dcomm.h"
#include "particles/Particles3D.h"
// output
#include "PSKOutput3D/PSKOutput.h"
#include "PSKOutput3D/PSKhdf5adaptor.h"
#include "inputoutput/Restart3D.h"
// performance
// #include "performances/Timing.h"
// wave
// #include "perturbation/Planewave.h"
// serial ASCII output
// #include "inputoutput/SerialIO.h"
// #include "inputoutput/ParallelIO.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using std::cerr;
using std::endl;



int main(int argc, char **argv) {
  // initialize MPI environment
  // nprocs = number of processors
  // myrank = rank of tha process*/
  int nprocs, myrank, mem_avail;
  MPIdata *mpi = new MPIdata(&argc, &argv);
  nprocs = mpi->nprocs;
  myrank = mpi->rank;
  // Timing *my_clock = new Timing(myrank); // performance and timing object
  // restart_cycle = 500; // each restart_cycle it writes a RESTART FILE, that can be used to monitor the sim
  Collective *col = new Collective(argc, argv); // Every proc loads the parameters of simulation from class Collective
  bool verbose = col->getVerbose();
  int restart_cycle = col->getRestartOutputCycle();
  string SaveDirName = col->getSaveDirName();
  string RestartDirName = col->getRestartDirName();
  const int restart = col->getRestart_status();
  const int ns = col->getNs();  // get the number of particle species involved in simulation
  const int first_cycle = col->getLast_cycle() + 1; // get the last cycle from the restart
  // initialize the virtual cartesian topology 
  VCtopology3D *vct = new VCtopology3D();
  // Check if we can map the processes into a matrix ordering defined in Collective.cpp
  if (nprocs != vct->getNprocs()) {
    if (myrank == 0) {
      cerr << "Error: " << nprocs << " processes cant be mapped into " << vct->getXLEN() << "x" << vct->getYLEN() << "x" << vct->getZLEN() << " matrix: Change XLEN,YLEN, ZLEN in method VCtopology3D.init()" << endl;
      mpi->finalize_mpi();
      return (1);
    }
  }
  // We create a new communicator with a 3D virtual Cartesian topology
  vct->setup_vctopology(MPI_COMM_WORLD);
  // Print the initial settings from the INPUT file 
  if (myrank == 0) {
    mpi->Print();
    vct->Print();
    col->Print();
  }
  if (myrank == 0)
    vct->PrintMapping();
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 1)
    vct->PrintMapping();
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 4)
    vct->PrintMapping();
  MPI_Barrier(MPI_COMM_WORLD);


  // Create the local grid
  MPI_Barrier(MPI_COMM_WORLD);
  Grid3DCU *grid = new Grid3DCU(col, vct);  // Create the local grid
  EMfields3D *EMf = new EMfields3D(col, grid);  // Create Electromagnetic Fields Object
  EMf->initBEAM(vct, grid, 0.12713296, 0.12713296, 0.25426593, 0.04237765);
  // EMf->initBEAM(vct,grid,0.12713296,0.12713296,0.25426593,0.1037765);
  // Allocation of particles
  Particles3D *part = new Particles3D[ns];
  for (int i = 0; i < ns; i++)
    part[i].allocate(i, col, vct, grid);
  // Initial Condition for PARTICLES if you are not starting from RESTART
  if (restart == 0) {
    // Planewave *wave = new Planewave(col, EMf, grid, vct);
    // wave->Wave_Rotated(part); // Single Plane Wave
    for (int i = 0; i < ns; i++)
      part[i].maxwellian(grid, EMf, vct); // all the species have Maxwellian distribution in the velocity
  }
  // Initialize the output (simulation results and restart file)
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;  // Create an Output Agent for HDF5 output
  hdf5_agent.set_simulation_pointers(EMf, grid, vct, mpi, col);
  for (int i = 0; i < ns; ++i)
    hdf5_agent.set_simulation_pointers_part(&part[i]);
  output_mgr.push_back(&hdf5_agent);  // Add the HDF5 output agent to the Output Manager's list
  if (myrank == 0 & restart < 2) {
    hdf5_agent.open(SaveDirName + "/settings.hdf");
    output_mgr.output("collective + total_topology + proc_topology", 0);
    hdf5_agent.close();
    hdf5_agent.open(RestartDirName + "/settings.hdf");
    output_mgr.output("collective + total_topology + proc_topology", 0);
    hdf5_agent.close();
  }
  stringstream num_proc;
  num_proc << myrank;
  if (restart == 0) {           // new simulation from input file
    hdf5_agent.open(SaveDirName + "/proc" + num_proc.str() + ".hdf");
    output_mgr.output("proc_topology ", 0);
    hdf5_agent.close();
  }
  else {                        // restart append the results to the previous simulation 
    hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
    output_mgr.output("proc_topology ", 0);
    hdf5_agent.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  // *******************************************//
  // **** Start the Simulation! ***//
  // *******************************************//
  for (int cycle = first_cycle; cycle < (col->getNcycles() + first_cycle); cycle++) {
    if (myrank == 0 && verbose) {
      cout << "***********************" << endl;
      cout << "*   cycle = " << cycle + 1 << "        *" << endl;
      cout << "***********************" << endl;
    }
    // interpolation
    // my_clock->start_interpP2G(); // for profiling

    EMf->setZeroDensities();    // set to zero the densities
    for (int i = 0; i < ns; i++)
      part[i].interpP2G(EMf, grid, vct);  // interpolate Particles to Grid(Nodes)
    EMf->sumOverSpecies(vct);   // sum all over the species
    MPI_Barrier(MPI_COMM_WORLD);
    EMf->interpDensitiesN2C(vct, grid); // calculate densities on centers from nodes
    EMf->calculateHatFunctions(grid, vct);  // calculate the hat quantities for the implicit method
    // my_clock->stop_interpP2G(); // for profiling
    MPI_Barrier(MPI_COMM_WORLD);
    // done with output // solve Maxwell equations
    // my_clock->start_field(); // for profiling
    EMf->calculateField(grid, vct); // calculate the EM fields
    // OUTPUT to large file, called proc**
    if (cycle % (col->getFieldOutputCycle()) == 0 || cycle == first_cycle) {
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("k_energy + E_energy + B_energy", cycle);
      output_mgr.output("Eall + Ball + rhos + Jsall", cycle);
      hdf5_agent.close();
    }
    if (cycle % (col->getParticlesOutputCycle()) == 0 && col->getParticlesOutputCycle() != 1) {
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("position + velocity + q ", cycle, 1);
      hdf5_agent.close();
    }

    // my_clock->stop_field(); // for profiling
    // push the particles
    // my_clock->start_mover(); // for profiling
    for (int i = 0; i < ns; i++)  // move each species
      mem_avail = part[i].mover_PC(grid, vct, EMf); // use the Predictor Corrector scheme 
    if (mem_avail < 0) {        // not enough memory space allocated for particles: stop the simulation
      if (myrank == 0) {
        cout << "*************************************************************" << endl;
        cout << "Simulation stopped. Not enough memory allocated for particles" << endl;
        cout << "*************************************************************" << endl;
      }
      cycle = col->getNcycles() + first_cycle;  // exit from the time loop
    }
    // my_clock->stop_mover(); // for profiling
    MPI_Barrier(MPI_COMM_WORLD);
    // Output save a file for the RESTART
    if (cycle % restart_cycle == 0 && cycle != first_cycle)
      writeRESTART(RestartDirName, myrank, cycle, ns, mpi, vct, col, grid, EMf, part, 0); // without ,0 add to restart file
    // MPI_Barrier(MPI_COMM_WORLD); 


  }
  if (mem_avail == 0)           // write the restart only if the simulation finished succesfully
    writeRESTART(RestartDirName, myrank, (col->getNcycles() + first_cycle) - 1, ns, mpi, vct, col, grid, EMf, part, 0);
  // stop profiling
  // my_clock->stopTiming();
  // close MPI
  mpi->finalize_mpi();
  return (0);

}
