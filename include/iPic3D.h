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

#include "EBox.h"
#include "OutputChunks.h"

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
    void CalculateField(int cycle);
    void CalculateBField();
    bool ParticlesMover();
    void WriteOutput(int cycle);
    void WriteConserved(int cycle);
    void WriteRestart(int cycle);
    void UpdateCycleInfo(int cycle);
    void Finalize();

    void WriteChunks(int cycle);
    
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
    // heat flux in the three directions, calculated summing particle by particle
    // NB: EFx, EFy, EFz saved at grid points are the heat flux DENSITIES
    double        *qx;
    double        *qy;
    double        *qz;
    double        *Qremoved;
    unsigned long *VelocityDist;

    /* velocity distribribution - vx, vy, vz */
    long double *VelocityDist_dir;
    
    Timing        *my_clock;

    EBox          *ebox;

    OutputChunks  *outputChunk;
    
    PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
    myOutputAgent < PSK::HDF5OutputAdaptor > hdf5_agent;  // Create an Output Agent for HDF5 output

    bool verbose;
    string SaveDirName;
    string RestartDirName;
    string cqsat;
    string cq;
    string ds;
    string *ds_vx;
    string *ds_vy;
    string *ds_vz;
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
    /* the cycle where to revert expansion */
    int cycle_reverseEBdir;

    /* whether to print VDF, vx, vy, vz */
    int WriteVDF_xyz; 
    
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
