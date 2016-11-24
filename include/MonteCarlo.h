/*******************************************************************************************
  MonteCarlo.h  -  Class for MonteCarlo wrapper
  -------------------
developers: Maria Elena Innocenti
 ********************************************************************************************/

#ifndef MONTECARLO_H
#define MONTECARLO_H

//#include "Particles3Dcomm.h"
//#include "Particles3D.h"
#include "TimeTasks.h"

/* for srand */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <string.h>

class c_Solver;
class Particles3D;
class Collective;

/**
 * 
 * Class for passing/ receiving info to the MonteCarlo code
 * 
 * @date Mon July 4 2016
 * @author Maria Elena Innocenti 
 * @version 1.0
 *
 */


/* definition of the particle struct to pass to the MC part */
/*struct MC_struct{
  double X; // particle position
  double Y;
  double Z;

  double Vx; // particle velocity
  double Vy;
  double Vz;

  double qp; // the q
  long long IDp; // because there is no unsigned long

  int sp; // species of the particle 

  

  };*/
// the new struct, to start using when Wei sends the new module
struct MC_struct{
  double X; // particle position
  double Y;
  double Z;

  double Vx; // particle velocity
  double Vy;
  double Vz;

  double qp; // the q
  long long IDp; // because there is no unsigned long

  int sp_in; // when entering the MC module: species of the in particle
             // when exiting the MC module: species of the particle that entered the MC module
  int sp_out; // when entering the MC module: value does not matter
              // when exiting the MC module: species of the out particle

  int CollisionType;  // written in the MC module
                      // 0: e-ion Excitation
                      // 1: e-ion Ionisation
                      // 2: e-ion Elastic
                      // 3: ion-ion Charge Exchange
                      // 4: ion-ion elastic
  
  // when exiting MC: variation in momentum; + (-) the output species gains (loses) momentum
  // when entering MC: value does not matter
  double dP_x;
  double dP_y;
  double dP_z;

  // when exiting MC: variation in energy; + (-) the output species gains (loses) energy
  // when entering MC: value does not matter  
  double dE;


};




/* C++ / fortran */
extern "C" {
  void sumsr(double *a, double *b, double *c);

  // dt in IS
  // IndexLastSpecies: 
  // Collision
  //void MCCInitialization(double *dt, int *IndexLastSpecies, double * CollPerc_MC, int Length, char name);
  void MCCInitialization(double *dt, int *IndexLastSpecies, double * CollPerc_MC);

  // NParBefore: # par in the input vector
  // MCCBefore: input vector
  // NParAfter: # par in the output vector
  // MCCAfter: output vector
  void MCC(int *NParBefore,MC_struct* MCCBefore,int *NParAfter,MC_struct *MCCAfter);
}

class MonteCarlo {
  public:

    /** constructor */
  
  /* in: temperature (const), pressure (const),number of species
     out: collision ratio, species per species */
  MonteCarlo(c_Solver *C);
  
  /** destructor */
  ~MonteCarlo();
  /** prova */
  void Test(double prova);

  // to know wether MC operations are to be done
  string getMonteCarloPlugIn();
  
  /* wrapper for the fortran code;
     in: all the info of selected particles (also index)
     out : modified particle, created particles */
  //    void MCC();
  
  /* select the particles to send to MC, 
     based on CollPerc;
     produces the array ToMC, with size CollNop_Total,
     to send to MC */
  
  void SelectCollidingParticles(c_Solver *C);
  
  /* calls the MC code,
     deals with the results */
  void MCWrapper(c_Solver *C, int cycle);

  /* print data in the MC_struct structure */
  void printMC_struct(int index, MC_struct P);


  /* Diagnostic functions */
  /* for the 5 kinds of collisions implemented- input: v^2 in code units
     output: cross section for that type of collisions in m^2 (CS_m2) and in code units (CS_CU)*/
  /* elastic collision electron- ion*/
  void eIonsigmaElastic(double v2, double *CS_m2, double *CS_CU);
  /* electron- ion, excitation */
  void eIonigmaExc(double v2, double *CS_m2, double *CS_CU);
  /* electron- ion, ionisation */
  void eIonsigmaIz(double v2, double *CS_m2, double *CS_CU);
  /* ion-ion. charge exchange */
  void IonIonsigmaCX(double v2, double *CS_m2, double *CS_CU);
  /* ion-ion, elastic */
  void IonIonsigmaElastic(double v2, double *CS_m2, double *CS_CU);

  // this variables are public because otherwise it's just too messy

  // collision frequencies

  double **** CF_eIonElastic;
  double **** CF_eIonExc;
  double **** CF_eIonIz;

  double **** CF_IonIonCX;
  double **** CF_IonIonElastic;

  double **** CF_affected;

  // dPx, dPy, dPz
  double **** C_dPx;
  double **** C_dPy;
  double **** C_dPz;

  // dE
  double **** C_dE;

  /// for moving averages of dP, dE
  int DOI;
  int DT_counter;
  int Smooth;

  // ---
  double ***** C_dP_DO_x;
  double ***** C_dP_DO_y;
  double ***** C_dP_DO_z;

  double ***** C_dE_DO;
  // ---
  double **** C_dP_MA_x;
  double **** C_dP_MA_y;
  double **** C_dP_MA_z;

  double **** C_dE_MA;
  // ---


  
  /// end for moving averages of dP, dE


  /* Diagnostics */

  /* set collision frequencies arrays to zero */
  void setZeroMCCDiagnostics();
  void setDiagnosticsToZero();
  /* add diagnostic functions */

  void addCF_eIonElastic(double weight[][2][2], int ix, int iy, int iz, int ns);
  void addCF_eIonExc(double weight[][2][2], int ix, int iy, int iz, int ns);
  void addCF_eIonIz(double weight[][2][2], int ix, int iy, int iz, int ns);
  void addCF_IonIonCX(double weight[][2][2], int ix, int iy, int iz, int ns);
  void addCF_IonIonElastic(double weight[][2][2], int ix, int iy, int iz, int ns);

  void addCF_affected(double weight[][2][2], int ix, int iy, int iz, int ns);

  // to collect diagnostics particle by particle
  void PBP_Diagnostics(c_Solver* C, MC_struct P);
  // to finalise this diagnostic
  void PBP_Collective(c_Solver* C);


  void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology3D * vct);

  void communicateGhostP2G(VirtualTopology3D * vct);

  double ****getCF_eIonElastic();
  double ****getCF_eIonExc();
  double ****getCF_eIonIz();
  double ****getCF_IonIonCX();
  double ****getCF_IonIonElastic();
  double ****getCF_affected();

  // Oct 16

  double ****getC_dPx();
  double ****getC_dPy();
  double ****getC_dPz();

  double ****getC_dE();

  double ****getC_dP_MA_x();
  double ****getC_dP_MA_y();
  double ****getC_dP_MA_z();

  double ****getC_dE_MA();
 private:
  // read from inputfile if collisions are done at all
  string Collisions;

  // grid properties
  int ns;
  int nxn; 
  int nyn;
  int nzn;
  double invVOL;
  double xstart;
  double ystart;
  double zstart;
  double xend;
  double yend;
  double zend;
  double inv_dx;
  double inv_dy;
  double inv_dz;
  
  // end grid properties
  
  string Gas;

  /** to normalise back and forth between code units and IS**/
  double IonRAM;            // relative Atomic Mass of the Neutral- normalised to mass of proton
  double Z;                 // charge state, NOT atomic number 
  double GasN;              // gas density, in m-3
  double GasT;              // gas temperature, in K
  double PlasmaN;           // plasma density, in m-3

  // omega_pi, to calculate dt in seconds which I need to pass to the MC plug-in   
  double omega_pi_IS;
  double dt_MC; // s, needed for the MC plug-in
  double mr; // mass ratio
  // ion skin depth in m, to normalise back and forth
  double di_IS;
  /** end to normalise back and forth between code units and IS**/
  

  int prova2;
  
  
  /* percentage of the Particles which Collide; calculated in the constructore -
     [0 -> 1] */
  double *CollPerc;
  double *CollPerc_MC;
  /* number of particles which collide - recalculated every time based on current number of particles */
  long long *CollNop;

  /* total number of particles to send to MC */
  long long CollNop_Total;

  /* max number of particles that can be exchanged with the MC code */
  long long MAX_MC;

  /* particles going to the MC code;
     max size: MAX_MC
     actual size (send to MC): CollNop_Total */
  MC_struct* ToMC;

  /* particles received from the MC code;
     max size: MAX_MC
     actual size (from MC): CollRcv_Total */
  MC_struct* FromMC;


  /* total number of particles received from MC */
  // i defined it as int because it's what MCC is spitting out -
  // i would actually like a long long
  int CollRcv_Total;

  /* # of particles received from MC, for debugging purposes */
  long long *CollRcv;

  /* variables for diagnostics */
  // from v^2 in code units to energy in eV for ELECTRONS
  double factor_e;
  // from v^2 in code units to energy in eV for IONS
  double factor_i;
  /* end variables for diagnostics */

  double ReducedC; // the value of the speed of light to pass to the MC module; not necessarily the real one; real from inputfile
  double ReducedEps0;
 
};


#endif
