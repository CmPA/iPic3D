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
#include "PetscSolver.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using std::cerr;
using std::endl;
using std::ofstream;

using namespace PSK;

namespace iPic3D {

  class c_Solver {

  public:
    int Init(int argc, char **argv);
    void InjectBoundaryParticles();
    void GatherMoments();
    void CalculateField();
    void CalculateBField();
    bool ParticlesMover();
    void WriteOutput(int cycle);
    void WriteConserved(int cycle);
    void WriteRestart(int cycle);
    void UpdateCycleInfo(int cycle);
    void Finalize();

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
    double 		  *Qtot;
    int           *totParticles;
    int           *globalTotParticles;
    double        *speciesTemp;
    double        *qom;

    PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
    myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;  // Create an Output Agent for HDF5 output

    bool verbose;
    string SaveDirName;
    string RestartDirName;
    string cqsat;
    string cq;
    string cq2;
    string ds;
    stringstream num_proc;
    int restart_cycle;
    int restart;
    int first_cycle;
    int ns;
    int nprocs;
    int myrank;
    int mem_avail;
    int nsatx;
    int nsaty;
    int nsatz;
    int nx0;
    int ny0;
    int nz0;
    int nDistributionBins;
    bool cylindrical;
    double Eenergy;
    double Benergy;
    double TOTenergy;
    double TOTmomentum;
    double x_center;
    double y_center;
    double z_center;
    double L_square;
    double L_outer;

    #ifdef __PETSC_SOLVER__
      PetscSolver *petscSolver;
    #endif

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

}

#endif
