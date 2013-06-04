/***************************************************************************
  iPIC3D.cpp  -  Main file for 3D simulation
  -------------------
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
// diagnostics
#include "performances/Timing.h"
#include "utility/TimeTasks.h"
#include "utility/debug.h"
// wave
// #include "perturbation/Planewave.h"
// serial ASCII output
// #include "inputoutput/SerialIO.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

// implementation (should be compiled separately)
#include "utility/diagnostics.cpp"

using namespace std;
using std::cerr;
using std::endl;
using std::ofstream;

/* global variables */
MPIdata *mpi;
TimeTasks timeTasks;

int main(int argc, char **argv) {
  // initialize MPI environment
  // nprocs = number of processors
  // myrank = rank of tha process*/
  int nprocs, myrank, mem_avail;
  mpi = new MPIdata(&argc, &argv);
  nprocs = mpi->nprocs;
  myrank = mpi->rank;

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
  // initialize the central cell index
  const int nx0 = col->getNxc() / vct->getXLEN(); // get the number of cells in x for each processor
  const int ny0 = col->getNyc() / vct->getYLEN(); // get the number of cells in y for each processor
  const int nz0 = col->getNzc() / vct->getZLEN(); // get the number of cells in z for each processor
  // Print the initial settings to stdout and a file
  if (myrank == 0) {
    mpi->Print();
    vct->Print();
    col->Print();
    col->save();
  }
  // Create the local grid
  MPI_Barrier(MPI_COMM_WORLD);
  Grid3DCU *grid = new Grid3DCU(col, vct);  // Create the local grid
  EMfields3D *EMf = new EMfields3D(col, grid);  // Create Electromagnetic Fields Object
  // EMf->initGEMnoPert(vct,grid);
  // EMf->initForceFree(vct,grid);
  EMf->initGEM(vct, grid);
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
    // part[i].force_free(grid,EMf,vct); // force free
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
  // Restart
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
  double Eenergy, Benergy, TOTenergy = 0.0, TOTmomentum = 0.0;
  double *Ke = new double[ns];
  double *momentum = new double[ns];
  string cq = SaveDirName + "/ConservedQuantities.txt";
  if (myrank == 0) {
    ofstream my_file(cq.c_str());
    my_file.close();
  }
  // Distribution functions
  int nDistributionBins = 1000;
  double maxVel = 0.;
  unsigned long* VelocityDist = new unsigned long [nDistributionBins];
  string ds = SaveDirName + "/DistributionFunctions.txt";  
  if (myrank==0) {
    ofstream my_file(ds.c_str());
    my_file.close();
  }
  string cqsat = SaveDirName + "/VirtualSatelliteTraces" + num_proc.str() + ".txt";
  // if(myrank==0){
  ofstream my_file(cqsat.c_str(), fstream::binary);
  int nsat = 3;
  for (int isat = 0; isat < nsat; isat++) {
    for (int jsat = 0; jsat < nsat; jsat++) {
      for (int ksat = 0; ksat < nsat; ksat++) {
        int index1 = 1 + isat * nx0 / nsat + nx0 / nsat / 2;
        int index2 = 1 + jsat * ny0 / nsat + ny0 / nsat / 2;
        int index3 = 1 + ksat * nz0 / nsat + nz0 / nsat / 2;
        my_file << grid->getXC(index1, index2, index3) << "\t" << grid->getYC(index1, index2, index3) << "\t" << grid->getZC(index1, index2, index3) << endl;
  }}}
  my_file.close();

  // *******************************************//
  // **** Start the Simulation! ***//
  // *******************************************//
  Timing *my_clock = new Timing(myrank);
  for (int cycle = first_cycle; cycle < (col->getNcycles() + first_cycle); cycle++) {
    if (myrank == 0 && verbose) {
      cout << "***********************" << endl;
      cout << "*   cycle = " << cycle + 1 << "        *" << endl;
      cout << "***********************" << endl;
    }
    timeTasks.resetCycle();
    // interpolation
    timeTasks.start(TimeTasks::MOMENTS);
    EMf->setZeroDensities();    // set to zero the densities
    for (int i = 0; i < ns; i++)
      part[i].interpP2G(EMf, grid, vct);  // interpolate Particles to Grid(Nodes)
    EMf->sumOverSpecies(vct);   // sum all over the species
    MPI_Barrier(MPI_COMM_WORLD);
    EMf->interpDensitiesN2C(vct, grid); // calculate densities on centers from nodes
    EMf->calculateHatFunctions(grid, vct);  // calculate the hat quantities for the implicit method
    MPI_Barrier(MPI_COMM_WORLD);
    timeTasks.end(TimeTasks::MOMENTS);

    // OUTPUT to large file, called proc**
    if (cycle % (col->getFieldOutputCycle()) == 0 || cycle == first_cycle) {
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("Eall + Ball + rhos + Jsall + pressure", cycle);
      // Pressure tensor is available
      hdf5_agent.close();
    }
    if (cycle % (col->getParticlesOutputCycle()) == 0 && col->getParticlesOutputCycle() != 1) {
      hdf5_agent.open_append(SaveDirName + "/proc" + num_proc.str() + ".hdf");
      output_mgr.output("position + velocity + q ", cycle, 1);
      hdf5_agent.close();
    }
    // write the virtual satellite traces

    if (ns > 2) {
      ofstream my_file(cqsat.c_str(), fstream::app);
      for (int isat = 0; isat < nsat; isat++) {
        for (int jsat = 0; jsat < nsat; jsat++) {
          for (int ksat = 0; ksat < nsat; ksat++) {
            int index1 = 1 + isat * nx0 / nsat + nx0 / nsat / 2;
            int index2 = 1 + jsat * ny0 / nsat + ny0 / nsat / 2;
            int index3 = 1 + ksat * nz0 / nsat + nz0 / nsat / 2;
            my_file << EMf->getBx(index1, index2, index3) << "\t" << EMf->getBy(index1, index2, index3) << "\t" << EMf->getBz(index1, index2, index3) << "\t";
            my_file << EMf->getEx(index1, index2, index3) << "\t" << EMf->getEy(index1, index2, index3) << "\t" << EMf->getEz(index1, index2, index3) << "\t";
            my_file << EMf->getJxs(index1, index2, index3, 0) + EMf->getJxs(index1, index2, index3, 2) << "\t" << EMf->getJys(index1, index2, index3, 0) + EMf->getJys(index1, index2, index3, 2) << "\t" << EMf->getJzs(index1, index2, index3, 0) + EMf->getJzs(index1, index2, index3, 2) << "\t";
            my_file << EMf->getJxs(index1, index2, index3, 1) + EMf->getJxs(index1, index2, index3, 3) << "\t" << EMf->getJys(index1, index2, index3, 1) + EMf->getJys(index1, index2, index3, 3) << "\t" << EMf->getJzs(index1, index2, index3, 1) + EMf->getJzs(index1, index2, index3, 3) << "\t";
            my_file << EMf->getRHOns(index1, index2, index3, 0) + EMf->getRHOns(index1, index2, index3, 2) << "\t";
            my_file << EMf->getRHOns(index1, index2, index3, 1) + EMf->getRHOns(index1, index2, index3, 3) << "\t";
      }}}
      my_file << endl;
      my_file.close();
    }
    // done with output


    // MAXWELL'S SOLVER
    timeTasks.start(TimeTasks::FIELDS);
    EMf->calculateE(grid, vct); // calculate the E field
    timeTasks.end(TimeTasks::FIELDS);

    // PARTICLE MOVER
    timeTasks.start(TimeTasks::PARTICLES);
    for (int i = 0; i < ns; i++)  // move each species
    {
      //#pragma omp task inout(part[i]) in(grid) target_device(booster)
      mem_avail = part[i].mover_PC(grid, vct, EMf); // use the Predictor Corrector scheme 
    }
    timeTasks.end(TimeTasks::PARTICLES);
    if (mem_avail < 0) {        // not enough memory space allocated for particles: stop the simulation
      if (myrank == 0) {
        cout << "*************************************************************" << endl;
        cout << "Simulation stopped. Not enough memory allocated for particles" << endl;
        cout << "*************************************************************" << endl;
      }
      cycle = col->getNcycles() + first_cycle;  // exit from the time loop
    }
    timeTasks.start(TimeTasks::BFIELD);
    EMf->calculateB(grid, vct); // calculate the B field
    timeTasks.end(TimeTasks::BFIELD);

    // print out total time for all tasks
    timeTasks.print_cycle_times();

    // write the conserved quantities
    if (cycle % col->getDiagnosticsOutputCycle() == 0) {
      Eenergy = EMf->getEenergy();
      Benergy = EMf->getBenergy();
      TOTenergy = 0.0;
      TOTmomentum = 0.0;
      for (int is = 0; is < ns; is++) {
        Ke[is] = part[is].getKe();
        TOTenergy += Ke[is];
        momentum[is] = part[is].getP();
        TOTmomentum += momentum[is];
      }
      if (myrank == 0) {
        ofstream my_file(cq.c_str(), fstream::app);
        my_file << cycle << "\t" << "\t" << (Eenergy + Benergy + TOTenergy) << "\t" << TOTmomentum << "\t" << Eenergy << "\t" << Benergy << "\t" << TOTenergy << endl;
        my_file.close();
      }
      // Velocity distribution
      for (int is=0; is < ns; is++) {
        maxVel = part[is].getMaxVelocity();
        VelocityDist = part[is].getVelocityDistribution(nDistributionBins, maxVel);
        if (myrank == 0){
          ofstream my_file(ds.c_str(), fstream::app);
          my_file << cycle << "\t" << is << "\t" << maxVel;
          for (int i=0; i < nDistributionBins; i++)
            my_file << "\t" << VelocityDist[i];
          my_file << endl;
          my_file.close();
        }
      }           
    }
    // write the RESTART file
    if (cycle % restart_cycle == 0 && cycle != first_cycle)
      writeRESTART(RestartDirName, myrank, cycle, ns, mpi, vct, col, grid, EMf, part, 0); // without ,0 add to restart file

  }
  if (mem_avail == 0)           // write the restart only if the simulation finished succesfully
    writeRESTART(RestartDirName, myrank, (col->getNcycles() + first_cycle) - 1, ns, mpi, vct, col, grid, EMf, part, 0);

  // stop profiling
  my_clock->stopTiming();

  // deallocate
  delete[]Ke;
  delete[]momentum;
  // close MPI
  mpi->finalize_mpi();
  return (0);

}

