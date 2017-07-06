/***************************************************************************
  Collective.h  -  Stefano Markidis, Giovanni Lapenta
  -------------------------------------------------------------------------- */


/*! Collective properties. Input physical parameters for the simulation.  Use ConfigFile to parse the input file @date Wed Jun 8 2011 @par Copyright: (C) 2011 K.U. LEUVEN @author Pierre Henri, Stefano Markidis @version 1.0 */

#ifndef Collective_H
#define Collective_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "ConfigFile.h"
#include "input_array.h"
#include "hdf5.h"

#include "Basic.h" /*! included for some basic mlmd ops */

#ifdef BATSRUS
#include "InterfaceFluid.h"
#endif

using std::cout;
using std::endl;
using std::ofstream;
using namespace std;

class Collective
#ifdef BATSRUS
: public InterfaceFluid
#endif
  {
  public:
    /*! constructor: initialize physical parameters with values */
    Collective(int argc, char **argv);
    /*! destructor */
    ~Collective();
    /*! read input file */
    void ReadInput(string inputfile);
    /*! read the restart input file from HDF5 */
    int ReadRestart(string inputfile);
    /*! Print physical parameters */
    void Print();
    /*! save setting in a file */
    void save();
    /*! get the physical space dimensions */
    int getDim();
    /*! Get object center - direction X */
    double getx_center();
    /*! Get object center - direction Y */
    double gety_center();
    /*! Get object center - direction Z */
    double getz_center();
    /*! Get object size - cubic box */
    double getL_square();
    /*! mlmd: use getNXC_mlmd insead */
    /*! Get the number of cells - direction X */
    //int getNxc();
    /*! Get the number of cells - direction Y */
    //int getNyc();
    /*! Get the number of cells - direction Z */
    //int getNzc();
    /*! end mlmd: use getNXC_mlmd instead */
    /* mlmd: use getPERIODICX_mlmd instead*/
    //bool getPERIODICX() {return (PERIODICX);};
    //bool getPERIODICY() {return (PERIODICY);};
    //bool getPERIODICZ() {return (PERIODICZ);};
    /*! mlmd: use getDx_mlmd instead */
    /*! Get the grid spacing - direction X */
    //double getDx();
    /*! Get the grid spacing - direction Y */
    //double getDy();
    /*! Get the grid spacing - direction Z */
    //double getDz();
    /*! get the light speed */
    double getC();
    /*! get the time step */
    double getDt();
    /*! get the decentering parameter */
    double getTh();
    /*! get the Smoothing value */
    double getSmooth();
    /*! get the number of time cycles */
    int getNcycles();
    /*! get the number of species */
    int getNs();
    /*! get the number of particles for different species */
    long getNp(int nspecies);
    /*! get the number of particles per cell */
    int getNpcel(int nspecies);
    /*! get the number of particles per cell - direction X */
    int getNpcelx(int nspecies);
    /*! get the number of particles per cell - direction Y */
    int getNpcely(int nspecies);
    /*! get the number of particles per cell - direction Z */
    int getNpcelz(int nspecies);
    /*! get maximum number of particles for different species */
    long getNpMax(int nspecies);
    /*! NpMax/Np is the ratio between the maximum number of particles allowed on a processor and the number of particles */
    double getNpMaxNpRatio();
    /*! get charge to mass ratio for different species */
    double getQOM(int nspecies);
    /*! get background charge for GEM challenge */
    double getRHOinit(int nspecies);
    /*! get rho injection */
    double getRHOinject(int nspecies);
    /*! get thermal velocity - X direction */
    double getUth(int nspecies);
    /*! get thermal velocity - Y direction */
    double getVth(int nspecies);
    /*! get thermal velocity - Z direction */
    double getWth(int nspecies);
    /*! get Drift velocity - Direction X */
    double getU0(int nspecies);
    /*! get Drift velocity - Direction Y */
    double getV0(int nspecies);
    /*! get Drift velocity - Direction Z */
    double getW0(int nspecies);
    /*! get the boolean value for TrackParticleID */
    bool getTrackParticleID(int nspecies);
    /*! get SaveDirName */
    string getSaveDirName();
    /*! get last_cycle */
    int getLast_cycle();
    /*! get RestartDirName */
    string getRestartDirName();

    /*! get Case type */
    string getCase();
    /*! get particle initialization type */
    string getPartInit();
    /*! get output writing method */
    string getWriteMethod();
    /*! get simulation name */
    string getSimName();
    /*! get poisson correction flag */
    string getPoissonCorrection();

    /*! get initial solution flag */
    bool getSolInit();

    /*! get Boundary Condition Particles: FaceXright */
    int getBcPfaceXright(int NG);
    /*! get Boundary Condition Particles: FaceXleft */
    int getBcPfaceXleft(int NG);
    /*! get Boundary Condition Particles: FaceYright */
    int getBcPfaceYright(int NG);
    /*! get Boundary Condition Particles: FaceYleft */
    int getBcPfaceYleft(int NG);
    /*! get Boundary Condition Particles: FaceYright */
    int getBcPfaceZright(int NG);
    /*! get Boundary Condition Particles: FaceYleft */
    int getBcPfaceZleft(int NG);

    /*! get Boundary Condition Particles: FaceXright */
    int getBcPfaceXright();
    /*! get Boundary Condition Particles: FaceXleft */
    int getBcPfaceXleft();
    /*! get Boundary Condition Particles: FaceYright */
    int getBcPfaceYright();
    /*! get Boundary Condition Particles: FaceYleft */
    int getBcPfaceYleft();
    /*! get Boundary Condition Particles: FaceYright */
    int getBcPfaceZright();
    /*! get Boundary Condition Particles: FaceYleft */
    int getBcPfaceZleft();

    /*! get Boundary Condition Electrostatic Potential: FaceXright */
    int getBcPHIfaceXright();
    /*! get Boundary Condition Electrostatic Potential:FaceXleft */
    int getBcPHIfaceXleft();
    /*! get Boundary Condition Electrostatic Potential:FaceYright */
    int getBcPHIfaceYright();
    /*! get Boundary Condition Electrostatic Potential:FaceYleft */
    int getBcPHIfaceYleft();
    /*! get Boundary Condition Electrostatic Potential:FaceYright */
    int getBcPHIfaceZright();
    /*! get Boundary Condition Electrostatic Potential:FaceYleft */
    int getBcPHIfaceZleft();

    /*! get Boundary ConditionElectric Field: FaceXright */
    int getBcEMfaceXright();
    /*! get Boundary Condition Electric Field: FaceXleft */
    int getBcEMfaceXleft();
    /*! get Boundary Condition Electric Field: FaceYright */
    int getBcEMfaceYright();
    /*! get Boundary Condition Electric Field: FaceYleft */
    int getBcEMfaceYleft();
    /*! get Boundary Condition Electric Field: FaceZright */
    int getBcEMfaceZright();
    /*! get Boundary Condition Electric Field: FaceZleft */
    int getBcEMfaceZleft();

    /*! get RESTART */
    int getRestart_status();

    /*! get the sheet thickness */
    double getDelta();
    /*! get the amplitude of the magnetic field along x */
    double getB0x();
    /*! get the amplitude of the magnetic field along y */
    double getB0y();
    /*! get the amplitude of the magnetic field along z */
    double getB0z();
    /*! get the amplitude of the magnetic field 1 along x */
    double getB1x();
    /*! get the amplitude of the magnetic field 1 along y */
    double getB1y();
    /*! get the amplitude of the magnetic field 1 along z */
    double getB1z();

    /*! get the boolean value for verbose results */
    bool getVerbose();

    /*! get the velocity of injection of the plasma from the wall */
    double getVinj();

    /*! get the converging tolerance for CG solver */
    double getCGtol();
    /*! get the converging tolerance for GMRES solver */
    double getGMREStol();
    /*! get the numbers of iteration for the PC mover */
    int getNiterMover();

    /*! output of fields */
    int getFieldOutputCycle();
    /*! output of particles */
    int getParticlesOutputCycle();
    /*! output of restart */
    int getRestartOutputCycle();
    /*! output of diagnostics */
    int getDiagnosticsOutputCycle();

    /*! Boundary condition selection for BCFace for the electric field components */
    int bcEx[6], bcEy[6], bcEz[6];
    /*! Boundary condition selection for BCFace for the magnetic field components */
    int bcBx[6], bcBy[6], bcBz[6];

    /*! get initfile */
    string getinitfile();
    /*! MLMD specific function */
    /*! MLMD gets */
    int getNgrids();
    int getgridLevel(int numgrid);
    int getparentGrid(int numgrid);
    /*! in terms of the PARENT grid */
    double getOx_P(int numgrid);
    double getOy_P(int numgrid);
    double getOz_P(int numgrid);
    /*! end: in terms of the PARENT grid */

    /*! in terms of the COARSEST grid */
    double getOx_SW(int numgrid);
    double getOy_SW(int numgrid);
    double getOz_SW(int numgrid);
    /*! end: in terms of the COARSEST grid */

    double getDx_mlmd(int numgrid);
    double getDy_mlmd(int numgrid);
    double getDz_mlmd(int numgrid);
    int getNxc_mlmd(int numgrid);
    int getNyc_mlmd(int numgrid);
    int getNzc_mlmd(int numgrid);
    double getLx_mlmd(int numgrid);
    double getLy_mlmd(int numgrid);
    double getLz_mlmd(int numgrid);

    /* return the periodicity per direction of grid with number N */
    bool getPERIODICX(int N) {return (PERIODICX_mlmd[N]);}
    bool getPERIODICY(int N) {return (PERIODICY_mlmd[N]);}
    bool getPERIODICZ(int N) {return (PERIODICZ_mlmd[N]);}

    bool getPERIODICX_P(int N) {return (PERIODICX_P_mlmd[N]);}
    bool getPERIODICY_P(int N) {return (PERIODICY_P_mlmd[N]);}
    bool getPERIODICZ_P(int N) {return (PERIODICZ_P_mlmd[N]);}

    /*! get the number of children of grid 'numgrid' */
    int getChildrenNum(int numgrid);

    /*! get the child number 'childnum' of the grid 'numgrid' 
      prints warning if childnum > ChildrenNum[ng]-1, to speed up debugging*/
    int getChildrenGrids(int numgrid, int childnum);

    /*! get the parent of grid 'numgrid'
      prints warning if numgrid=0, to speed up debugging */
    int getParentGrid(int numgrid);

    int getTopologyType();
    
    /* get the topology per grid */
    int getXLEN_mlmd(int N);
    int getYLEN_mlmd(int N);
    int getZLEN_mlmd(int N);

    /* to use only when creating topology, for comparison
       with the numGrid obtained there -
       should be the same, problems if it's not */
    int getnumGrid_clt();

    /* the rank in MPI_COMM_WORLD of the lowest-ranker core in grid n
       needed to create field parent-child communicators */
    int getLowestRankOfGrid(int n);
    /* same for highest-rank core, for particle parent-child communicator */
    int getHighestRankOfGrid(int n);

    /*! get whether to perform mlmd operations */
    bool getMLMD_BC();
    bool getMLMD_PROJECTION();
    bool getMLMD_ParticleREPOPULATION();
    bool getMLMD_InitialInterpolation();
    bool getMLMD_BCBufferArea();
    bool getMLMD_fixBCenters();

    bool getAllowPMsgResize();

    /*! returns MLMDVerbose */
    bool getMLMDVerbose();
    /*! end MLMD gets */
    /*! a first sanity check on MLMD inputs, called at the end of the constructor */
    void checkMLMDinputs();
    /*! end MLMD specific function */

  private:
    /*! inputfile */
    string inputfile;
    /*! Restart HDF5 file */
    string initfile;
    /*! light speed */
    double c;
    /*! 4 pi */
    double fourpi;
    /*! time step */
    double dt;
    /*! decentering parameter */
    double th;
    /*! Smoothing value */
    double Smooth;
    /*! number of time cycles */
    int ncycles;
    /*! physical space dimensions */
    int dim;
    /*! simulation box length - X direction */
    double Lx;
    /*! simulation box length - Y direction */
    double Ly;
    /*! simulation box length - Z direction */
    double Lz;
    /*! object center - X direction */
    double x_center;
    /*! object center - Y direction */
    double y_center;
    /*! object center - Z direction */
    double z_center;
    /*! object size - assuming a cubic box */
    double L_square;
    /*! number of cells - X direction */
    int nxc;
    /*! number of cells - Y direction */
    int nyc;
    /*! number of cells - Z direction */
    int nzc;
    /*! grid spacing - X direction */
    double dx;
    /*! grid spacing - Y direction */
    double dy;
    /*! grid spacing - Z direction */
    double dz;

    /* mlmd: periodicity values are now vectors */
    bool *PERIODICX_mlmd;     /*! Periodicity in the X direction */
    bool *PERIODICY_mlmd;     /*! Periodicity in the X direction */
    bool *PERIODICZ_mlmd;     /*! Periodicity in the X direction */

    // and for particles
    bool *PERIODICX_P_mlmd;     /*! Periodicity in the X direction */
    bool *PERIODICY_P_mlmd;     /*! Periodicity in the X direction */
    bool *PERIODICZ_P_mlmd;     /*! Periodicity in the X direction */

    /*! number of species */
    int ns;
    /*! number of particles per cell */
    int *npcel;
    /*! number of particles per cell - X direction */
    int *npcelx;
    /*! number of particles per cell - Y direction */
    int *npcely;
    /*! number of particles per cell - Z direction */
    int *npcelz;
    /*! number of particles array for different species */
    long *np;
    /*! maximum number of particles array for different species */
    long *npMax;
    /*! max number of particles */
    double NpMaxNpRatio;
    /*! charge to mass ratio array for different species */
    double *qom;
    /*! charge to mass ratio array for different species */
    double *rhoINIT;
    /*! density of injection */
    double *rhoINJECT;
    /*! thermal velocity - Direction X */
    double *uth;
    /*! thermal velocity - Direction Y */
    double *vth;
    /*! thermal velocity - Direction Z */
    double *wth;
    /*! Drift velocity - Direction X */
    double *u0;
    /*! Drift velocity - Direction Y */
    double *v0;
    /*! Drift velocity - Direction Z */
    double *w0;

    /*! Case type */
    string Case;
    /*! Fields initialization type */
    string FieldsInit;
    /*! Particle initialization type */
    string PartInit;
    /*! Output writing method */
    string wmethod;
    /*! Simulation name */
    string SimName;
    /*! Poisson correction flag */
    string PoissonCorrection;

    /*! HDF5 initial solution flag */
    bool SolInit;

    /*! TrackParticleID */
    bool *TrackParticleID;
    /*! SaveDirName */
    string SaveDirName;
    /*! RestartDirName */
    string RestartDirName;
    /*! get inputfile */
    string getinputfile();
    /*! restart_status 0 --> no restart; 1--> restart, create new; 2--> restart, append; */
    int restart_status;
    /*! last cycle */
    int last_cycle;

    /*! Boundary condition on particles per grid level -1 = mlmd 0 = exit 1 = perfect mirror 2 = riemission */
    /*! Boundary Condition Particles: FaceXright */
    int *bcPfaceXright;
    /*! Boundary Condition Particles: FaceXleft */
    int *bcPfaceXleft;
    /*! Boundary Condition Particles: FaceYright */
    int *bcPfaceYright;
    /*! Boundary Condition Particles: FaceYleft */
    int *bcPfaceYleft;
    /*! Boundary Condition Particles: FaceYright */
    int *bcPfaceZright;
    /*! Boundary Condition Particles: FaceYleft */
    int *bcPfaceZleft;


    /*! Field Boundary Condition 0 = Dirichlet Boundary Condition: specifies the valueto take pn the boundary of the domain 1 = Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain 2 = Periodic Condition */
    /*! Boundary Condition Electrostatic Potential: FaceXright */
    int bcPHIfaceXright;
    /*! Boundary Condition Electrostatic Potential:FaceXleft */
    int bcPHIfaceXleft;
    /*! Boundary Condition Electrostatic Potential:FaceYright */
    int bcPHIfaceYright;
    /*! Boundary Condition Electrostatic Potential:FaceYleft */
    int bcPHIfaceYleft;
    /*! Boundary Condition Electrostatic Potential:FaceZright */
    int bcPHIfaceZright;
    /*! Boundary Condition Electrostatic Potential:FaceZleft */
    int bcPHIfaceZleft;

    /*! Boundary Condition EM Field: FaceXright */
    int bcEMfaceXright;
    /*! Boundary Condition EM Field: FaceXleft */
    int bcEMfaceXleft;
    /*! Boundary Condition EM Field: FaceYright */
    int bcEMfaceYright;
    /*! Boundary Condition EM Field: FaceYleft */
    int bcEMfaceYleft;
    /*! Boundary Condition EM Field: FaceZright */
    int bcEMfaceZright;
    /*! Boundary Condition EM Field: FaceZleft */
    int bcEMfaceZleft;


    /*! GEM Challenge parameters */
    /*! current sheet thickness */
    double delta;
    /* Amplitude of the field */
    double B0x;
    double B0y;
    double B0z;
    double B1x;
    double B1y;
    double B1z;


    /*! boolean value for verbose results */
    bool verbose;
    /*! RESTART */
    bool RESTART1;
    /*! SOLINIT */
    bool SOLINIT1;

    /*! velocity of the injection from the wall */
    double Vinj;

    /*! CG solver stopping criterium tolerance */
    double CGtol;
    /*! GMRES solver stopping criterium tolerance */
    double GMREStol;
    /*! mover predictor correcto iteration */
    int NiterMover;

    /*! Output for field */
    int FieldOutputCycle;
    /*! Output for particles */
    int ParticlesOutputCycle;
    /*! restart cycle */
    int RestartOutputCycle;
    /*! Output for diagnostics */
    int DiagnosticsOutputCycle;

    // MLMD variables
    int Ngrids;
    int *gridLevel;
    /*! number of the mlmd parent grid */
    int *parentGrid;
    /*! number of Children each mlmd grid has*/
    int *childrenNum;
    /*! children of each mlmd grid; each grid can have more than one child;
      row: mlmd grids
      columns: number of the children grid; use ChildrenNum here */
    int **childrenGrids;
    // for when i need the numgrid BEFORE the topology is created; 
    // compared with the num created in topology
    // that number is to one to propagate
    int numGrid_clt;

    double *Ox_P; // in terms of the PARENT grid, not of the COARSEST grid
    double *Oy_P;
    double *Oz_P;
    
    double *Ox_SW; // in terms of the COARSEST grid, at the moment used only for output
    double *Oy_SW;
    double *Oz_SW;

    int *nxc_mlmd;
    int *nyc_mlmd;
    int *nzc_mlmd;
    double *dx_mlmd;
    double *dy_mlmd;
    double *dz_mlmd;
    double *Lx_mlmd;
    double *Ly_mlmd;
    double *Lz_mlmd;

    int *XLEN_mlmd;
    int *YLEN_mlmd;
    int *ZLEN_mlmd;

    // here, the lowest rank in MPI_COMM_WORLD for the different grids
    // used to create the parent/ child communicator
    int *LowestRankOfGrid;
    // analogous, for particle parent / child communicator
    int *HighestRankOfGrid;

    /*! how the grids are distributed in the cores:
      0: one piece of each grid per core
      1: each core hosts a part of a single grid only 
      --- unused at the moment */
    int TopologyType; 

    /*! picked up by other classes, true for mlmd related output */
    bool MLMDVerbose;

    /* wether to perform mlmd operations */
       
    bool MLMD_BC;
    bool MLMD_PROJECTION;
    bool MLMD_ParticleREPOPULATION;
    bool MLMD_InitialInterpolation;
    bool MLMD_BCBufferArea;
    bool MLMD_fixBCenters;

    /* wether to allow the resize of the buffers containing the particles to send CG->RG for repopulation
       false may be heavy on memory consumption
       true may be heavy on performance */
    bool AllowPMsgResize;
    // end MLMD variables
};

#endif
