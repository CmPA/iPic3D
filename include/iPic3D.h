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
    void GatherMoments_Init();
    void CalculateField();
    void CalculateField_ECSIM(int i);
    void interpBC2N_ECSIM(int i);
    void CalculateBField(int i);
    bool ParticlesMover();
    void WriteOutput(int cycle);
    void WriteConserved(int cycle);
    void WriteRestart(int cycle);
    void UpdateCycleInfo(int cycle);
    void Finalize();
    


    inline int FirstCycle();
    inline int LastCycle();
    inline int get_myrank();
    /*! mlmd */
    inline int get_numGrid();
    /*! end mlmd */

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

    PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
    myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;  // Create an Output Agent for HDF5 output

    bool verbose;
    string SaveDirName;
    string RestartDirName;
    string cqsat;
    string cq;
    string *ds;
    stringstream num_proc; // mlmd: my_rank, hence local to the grid
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

    /* wether to print the distribution function*/
    bool WriteDistrFun;


    /*! mlmd variables */
    /*! number of the current grid in the mlmd hierarchy */
    int numGrid;
    /*! stringstream with grid number, for output-related ops 
     analogous to stringstream num_proc*/
    stringstream num_grid_STR;

    bool MlmdVerbose;

    /* whether to perform mlmd operations; 0 or 1 -
       the structure is set in place anyway; 
       this is mostly for debugging purposes */
    int MLMD_BC;
    int MLMD_PROJECTION;
    int MLMD_ParticleREPOPULATION;
    // if false, particle repopulation is done after mover
    // if true, particle repopulation done after CG has calculated moments
    bool FluidLikeRep;
    //int MLMD_InitialInterpolation;

    // to repopulate before the mover rather than later
    bool RepopulateBeforeMover;

    /*! end mlmd variables */

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
  /*! mlmd gets*/
  inline int c_Solver::get_numGrid() {
    return (numGrid);
  }
  /*! end mlmd gets */
}

#endif
