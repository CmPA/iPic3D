/***************************************************************************
  Collective.h  -  Stefano Markidis, Giovanni Lapenta
  -------------------------------------------------------------------------- */


/*! Collective properties. Input physical parameters for the simulation.  Use ConfigFile to parse the input file @date Wed Jun 8 2011 @par Copyright: (C) 2011 K.U. LEUVEN @author Pierre Henri, Stefano Markidis @version 1.0 */

#ifndef Collective_H
#define Collective_H

#ifdef BATSRUS
#include "InterfaceFluid.h"
#endif


#include <math.h>
//#include <iostream>
//#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "ConfigFile.h"
#include "input_array.h"
#include "hdf5.h"
//#include "CollectiveIO.h"
using namespace std;

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

    // accessors
    //
    int getDim()const{ return (dim); }
    double getLx()const{ return (Lx); }
    double getLy()const{ return (Ly); }
    double getLz()const{ return (Lz); }
    double getx_center()const{ return (x_center); }
    double gety_center()const{ return (y_center); }
    double getz_center()const{ return (z_center); }
    double getL_square()const{ return (L_square); }
    int getNxc()const{ return (nxc); }
    int getNyc()const{ return (nyc); }
    int getNzc()const{ return (nzc); }
    int getXLEN()const{ return (XLEN); }
    int getYLEN()const{ return (YLEN); }
    int getZLEN()const{ return (ZLEN); }
    bool getPERIODICX()const{ return (PERIODICX); }
    bool getPERIODICY()const{ return (PERIODICY); }
    bool getPERIODICZ()const{ return (PERIODICZ); }
    double getDx()const{ return (dx); }
    double getDy()const{ return (dy); }
    double getDz()const{ return (dz); }
    double getC()const{ return (c); }
    double getDt()const{ return (dt); }
    double getTh()const{ return (th); }
    double getSmooth()const{ return (Smooth); }
    int getNcycles()const{ return (ncycles); }
    int getNs()const{ return (ns); }
    int getNpcel(int nspecies)const{ return (npcel[nspecies]); }
    int getNpcelx(int nspecies)const{ return (npcelx[nspecies]); }
    int getNpcely(int nspecies)const{ return (npcely[nspecies]); }
    int getNpcelz(int nspecies)const{ return (npcelz[nspecies]); }
    int getNp(int nspecies)const{ return (np[nspecies]); }
    int getNpMax(int nspecies)const{ return (npMax[nspecies]); }
    double getNpMaxNpRatio()const{ return (NpMaxNpRatio); }
    double getQOM(int nspecies)const{ return (qom[nspecies]); }
    double getRHOinit(int nspecies)const{ return (rhoINIT[nspecies]); }
    double getRHOinject(int nspecies)const { return(rhoINJECT[nspecies]); }
    double getUth(int nspecies)const{ return (uth[nspecies]); }
    double getVth(int nspecies)const{ return (vth[nspecies]); }
    double getWth(int nspecies)const{ return (wth[nspecies]); }
    double getU0(int nspecies)const{ return (u0[nspecies]); }
    double getV0(int nspecies)const{ return (v0[nspecies]); }
    double getW0(int nspecies)const{ return (w0[nspecies]); }
    int getBcPfaceXright()const{ return (bcPfaceXright); }
    int getBcPfaceXleft()const{ return (bcPfaceXleft); }
    int getBcPfaceYright()const{ return (bcPfaceYright); }
    int getBcPfaceYleft()const{ return (bcPfaceYleft); }
    int getBcPfaceZright()const{ return (bcPfaceZright); }
    int getBcPfaceZleft()const{ return (bcPfaceZleft); }
    int getBcPHIfaceXright()const{ return (bcPHIfaceXright); }
    int getBcPHIfaceXleft()const{ return (bcPHIfaceXleft); }
    int getBcPHIfaceYright()const{ return (bcPHIfaceYright); }
    int getBcPHIfaceYleft()const{ return (bcPHIfaceYleft); }
    int getBcPHIfaceZright()const{ return (bcPHIfaceZright); }
    int getBcPHIfaceZleft()const{ return (bcPHIfaceZleft); }
    int getBcEMfaceXright()const{ return (bcEMfaceXright); }
    int getBcEMfaceXleft()const{ return (bcEMfaceXleft); }
    int getBcEMfaceYright()const{ return (bcEMfaceYright); }
    int getBcEMfaceYleft()const{ return (bcEMfaceYleft); }
    int getBcEMfaceZright()const{ return (bcEMfaceZright); }
    int getBcEMfaceZleft()const{ return (bcEMfaceZleft); }
    double getDelta()const{ return (delta); }
    double getB0x()const{ return (B0x); }
    double getB0y()const{ return (B0y); }
    double getB0z()const{ return (B0z); }
    double getB1x()const{ return (B1x); }
    double getB1y()const{ return (B1y); }
    double getB1z()const{ return (B1z); }
    bool getVerbose()const{ return (verbose); }
    bool getTrackParticleID(int nspecies)const
      { return (TrackParticleID[nspecies]); }
    int getRestart_status()const{ return (restart_status); }
    string getSaveDirName()const{ return (SaveDirName); }
    string getRestartDirName()const{ return (RestartDirName); }
    string getinputfile()const{ return (inputfile); }
    string getCase()const{ return (Case); }
    string getSimName()const{ return (SimName); }
    string getWriteMethod()const{ return (wmethod); }
    string getPoissonCorrection()const{ return (PoissonCorrection); }
    int getLast_cycle()const{ return (last_cycle); }
    double getVinj()const{ return (Vinj); }
    double getCGtol()const{ return (CGtol); }
    double getGMREStol()const{ return (GMREStol); }
    int getNiterMover()const{ return (NiterMover); }
    int getFieldOutputCycle()const{ return (FieldOutputCycle); }
    int getParticlesOutputCycle()const{ return (ParticlesOutputCycle); }
    int getRestartOutputCycle()const{ return (RestartOutputCycle); }
    int getDiagnosticsOutputCycle()const{ return (DiagnosticsOutputCycle); }

    /*! Boundary condition selection for BCFace for the electric field components */
    int bcEx[6], bcEy[6], bcEz[6];
    /*! Boundary condition selection for BCFace for the magnetic field components */
    int bcBx[6], bcBy[6], bcBz[6];

  private:
    /*! inputfile */
    string inputfile;
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
    /*! number of MPI subdomains in each direction */
    int XLEN;
    int YLEN;
    int ZLEN;
    /*! periodicity in each direction */
    bool PERIODICX;
    bool PERIODICY;
    bool PERIODICZ;
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
    int *np;
    /*! maximum number of particles array for different species */
    int *npMax;
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
    /*! Output writing method */
    string wmethod;
    /*! Simulation name */
    string SimName;
    /*! Poisson correction flag */
    string PoissonCorrection;

    /*! TrackParticleID */
    bool *TrackParticleID;
    /*! SaveDirName */
    string SaveDirName;
    /*! RestartDirName */
    string RestartDirName;
    /*! restart_status 0 --> no restart; 1--> restart, create new; 2--> restart, append; */
    int restart_status;
    /*! last cycle */
    int last_cycle;

    /*! Boundary condition on particles 0 = exit 1 = perfect mirror 2 = riemission */
    /*! Boundary Condition Particles: FaceXright */
    int bcPfaceXright;
    /*! Boundary Condition Particles: FaceXleft */
    int bcPfaceXleft;
    /*! Boundary Condition Particles: FaceYright */
    int bcPfaceYright;
    /*! Boundary Condition Particles: FaceYleft */
    int bcPfaceYleft;
    /*! Boundary Condition Particles: FaceYright */
    int bcPfaceZright;
    /*! Boundary Condition Particles: FaceYleft */
    int bcPfaceZleft;


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
};
typedef Collective CollectiveIO;

#endif
