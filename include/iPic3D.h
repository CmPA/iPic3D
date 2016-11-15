/***************************************************************************
  iPIC3D.cpp  -  Main file for 3D simulation
  -------------------
 ************************************************************************** */

#ifndef _IPIC3D_H_
#define _IPIC3D_H_

#include "MPIdata.h"
#include "VCtopology3D.h"
#include "Collective.h"
#include "Grid3DCU.h"
#include "EMfields3D.h"
#include "Particles3D.h"
#include "Restart3D.h"
#include "Timing.h"
#include "ParallelIO.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using std::cerr;
using std::endl;
using std::ofstream;

using namespace PSK;

class MonteCarlo;


class c_Solver {
  friend class MonteCarlo;
  
  public:
  int Init(int argc, char **argv);
  void InjectBoundaryParticles();
  void GatherMoments(int cycle);
  void CalculateField();
  void CalculateBField();
  bool ParticlesMover();
  void WriteOutput(int cycle, MonteCarlo *MC);
  void WriteConserved(int cycle);
  void WriteRestart(int cycle);
  void UpdateCycleInfo(int cycle);
  void Finalize();
  void WriteCollisionDiagnostics(MonteCarlo *MCC, int cycle);
  
  inline int FirstCycle();
  inline int LastCycle();
  inline int get_myrank();
  
 private:
  MPIdata       * mpi;
  Collective    *col;
  VCtopology3D  *vct;
  Grid3DCU      *grid;
  EMfields3D    *EMf;
  Particles3D   *part;
  double        *Ke;
  double        *momentum;
  double        *Qremoved;
  unsigned long *VelocityDist;
  Timing        *my_clock;

  MonteCarlo *MC;
  
  PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
  myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;  // Create an Output Agent for HDF5 output
  
  bool verbose;
  string SaveDirName;
  string RestartDirName;
  string cqsat;
  string cq;
  string ds;
  stringstream num_proc;
  int restart_cycle;
  int restart;
  int first_cycle;
  int ns;
  int nprocs;
  int myrank;
  int mem_avail;
  int nsat;
  int nx0;
  int ny0;
  int nz0;
  int nDistributionBins;
  double Eenergy;
  double Benergy;
  double TOTenergy;
  double TOTmomentum;
};

inline int c_Solver::FirstCycle() {
  return (first_cycle);
}
inline int c_Solver::LastCycle() {
  return (col->getNcycles() + first_cycle-1);
}
inline int c_Solver::get_myrank() {
  return (myrank);
}



#endif
