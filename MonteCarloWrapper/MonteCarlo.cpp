/*******************************************************************************************
 MonteCarlo.cpp  -  Wrapper to interface to the MonteCarlo code
  -------------------
developers: Maria Elena Innocenti
 ********************************************************************************************/


#include <iostream>
#include <math.h>
#include <string.h>
using namespace std;

#include "Basic.h"
#include "MPIdata.h"
#include "TimeTasks.h"

#include "Particles3D.h"

#include "MonteCarlo.h"

#include "hdf5.h"
#include <complex>

#include "iPic3D.h"
#include "Collective.h"
#include "Particles3D.h"


#include "PhyConstants.h"

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define MIN_VAL   1E-16
// particles processed together
#define P_SAME_TIME 2

/**
 * 
 * Wrapper to interface to the MonteCarlo code
 * @date Mon July 4 2016
 * @author Maria Elena Innocenti
 * @version 1.0
 *
 */

/** constructor */
MonteCarlo::MonteCarlo(c_Solver* C) { 
  //MonteCarlo::MonteCarlo(Collective *col){

  int myrank= C->vct->getCartesian_rank();
  Collisions= C->col->getCollisions();

  if (myrank == 0)
    { cout << "Collisions: " << Collisions <<endl;
    }

  if (!(Collisions == "yes")){
    return;
  }

  /* grid init */
  Smooth= C->col->getSmooth();
  nxn= C->grid->getNXN();
  nyn= C->grid->getNYN();
  nzn= C->grid->getNZN();
  invVOL = C->grid->getInvVOL();
  xstart= C->grid->getXstart();
  ystart= C->grid->getYstart();
  zstart= C->grid->getZstart();
  xend= C->grid->getXend();
  yend= C->grid->getYend();
  zend= C->grid->getZend();
  inv_dx= 1.0/C->grid->getDX();
  inv_dy= 1.0/C->grid->getDY();
  inv_dz= 1.0/C->grid->getDZ();

  /* physical init */
  ns=      C->col->getNs();
  Gas=     C->col->getGas();
  GasN=    C->col->getGasN();
  GasT=    C->col->getGasT();
  PlasmaN= GasN* C->col->getDensityRatio();
  ReducedC = C->col->getReducedC();
  ReducedEps0 = Realeps0;

  if (myrank==0 and fabs(ReducedC-Realc)/Realc > 0.01){
    cout << "We are using a reduced speed of light of c= " << ReducedC << ", reduced of a factor " << Realc/ReducedC << endl;
    cout << "I still have to decide what to do with eps0, at the moment it is the real valus eps0= " << ReducedEps0 << endl;
    
  }
  
  Z=1; // charge state can be only 1
  mr=      fabs(C->col->getQOM(0));

  if (Gas== "Hydrogen"){
    IonRAM= 1.0;
  }
  else{
    if (myrank==0){
      cout << "We only deal in Hydrogen here, aborting " << endl;
      abort();
    }
  }
  
  omega_pi_IS= sqrt(PlasmaN* Z*Z *e*e/ Realeps0 /(mP*IonRAM) ) ;
  // s
  // checked with Gianni, no 2 PI
  dt_MC= C->col->getDt()/ omega_pi_IS;
  // m
  di_IS= ReducedC/omega_pi_IS;

  // init for diagnostic
  
  // from v^2 in code units to energy in eV for ELECTRONS
  //factor_e= (mP/mr)*c*c*JtoeV/2.0;
  factor_e= (mP/mr)*ReducedC*ReducedC*JtoeV/2.0;
  // from v^2 in code units to energy in eV for IONS
  //factor_i= mP*IonRAM*c*c*JtoeV/2.0;
  factor_i= mP*IonRAM*ReducedC*ReducedC*JtoeV/2.0; 
  // end init for diagnostic

  
  if (myrank ==0){
    cout << "Gas density: " <<GasN << endl;
    cout << "Plasma density: " <<PlasmaN << endl;
    cout << "omega_pi: " << omega_pi_IS <<"s-1" <<endl;
    cout << "dt iPic: " << C->col->getDt() <<", dt_MC: " << dt_MC << endl;
    cout << "di_IS: " << di_IS << endl;
  }
  
  /// write the file MCC.txt that the fortran part reads

  int IndexLastSpecies= ns-1;

  // write in the result folder, for the record - with comments
  string SaveDirName = C->col->getSaveDirName();
  string cq = SaveDirName + "/MCC.txt";
  if (myrank == 0) {
    ofstream my_file(cq.c_str());
    my_file << mP*IonRAM << " " << GasT << " " << GasN << endl;
    my_file << "!!! Ion mass [Kg] - Gas temperature [K] - Gas density [m-3]" << endl;
    my_file << IndexLastSpecies << endl;
    my_file << "!!! Index Last Species" << endl;
    for (int is=0; is < ns; is++){
      if (C->col->getQOM(is) < 0) {//(is%2==0){ // electrons
	my_file << is%2 << "   " << mP/mr <<"   " << ReducedC <<endl; 
	my_file << "!!! Electron mass - normalisation for electron velocities (c)" << endl;
      }
      else{// ions
	my_file << is%2 << "   " << mP*IonRAM <<"   "<< ReducedC <<endl;
	my_file<< "!!! Ion mass - normalisation for ion velocities (c)" << endl;
      }
    }
    my_file << (mP*IonRAM)*ReducedC << endl;
    my_file << "!!! normalisation for moments (ion mass *c)"<<endl;
    my_file << (mP*IonRAM)*ReducedC*ReducedC << endl;
    my_file << "!!! normalisation for energy (ion mass *c *c )"<<endl;
    
    my_file.close();
  }

  // write in the current folder, fortran reads this - without comments
  string cq2 = "./MCC.txt";
  if (myrank == 0) {
    ofstream my_file(cq2.c_str());
    my_file << mP*IonRAM << " " << GasT << " " << GasN << endl;
    //my_file << "!!! Ion mass [Kg] - Gas temperature [K] - Gas density [m-3]" << endl;
    my_file << IndexLastSpecies << endl;
    //my_file << "!!! Index Last Species" << endl;
    for (int is=0; is < ns; is++){
      if (C->col->getQOM(is) < 0) {//(is%2==0){ // electrons
	my_file << is%2 << "   " << mP/mr <<"   " << ReducedC <<endl; 
	//my_file << "!!! Electron mass - normalisation for electron velocities (c)" << endl;
      }
      else{// ions
	my_file << is%2 << "   " << mP*IonRAM <<"   "<< ReducedC <<endl;
	//my_file<< "!!! Ion mass - normalisation for ion velocities (c)" << endl;
      }
    }
    my_file << (mP*IonRAM)*ReducedC << endl;
    //my_file << "!!! normalisation for moments (ion mass *c)"<<endl;
    my_file << (mP*IonRAM)*ReducedC*ReducedC << endl;
    //my_file << "!!! normalisation for energy (ion mass *c *c )"<<endl;
    
    my_file.close();
  }
  /// end the file MCC.txt that the fortran part reads


  /* end physical init */


  /* instantiation and stuff */

  // with DoubleGEM, only the upper cores will be collisional
  // so check YLEN%2 ==0, otherwise exits

  /*if (C->col->getCase()== "DoubleGEM"){
    if (! (C->col->getYLEN() % 2)){
      cout << "For some internal reasons, with collisions and DoubleGEM, i need even YLEN" <<endl;
      cout << "Aborting now... " <<endl;
      abort();
    }
    }

  if (C->col->getCase()== "DoubleGEM" and C->vct->getCoordinates(1)<  C->vct->getYLEN()/2){
    return;
    }*/


  CollPerc= new double[ns];
  CollNop= new long long[ns];
  CollRcv= new long long[ns];



  /* max number of particles that can be exchanged */
  MAX_MC= 1.0*C->col->npMax[0]; 

  /* array of particles going to MC */
  ToMC= new MC_struct[MAX_MC];

  /* array of particles coming from MC */
  FromMC= new MC_struct[MAX_MC];
  
  /* initialisation of diagnostics */
 
  // used only for electrons, but instantiated for all
  CF_eIonElastic = newArr4(double, ns, nxn, nyn, nzn);
  CF_eIonExc = newArr4(double, ns, nxn, nyn, nzn);
  CF_eIonIz = newArr4(double, ns, nxn, nyn, nzn);
  // only for ions
  CF_IonIonCX= newArr4(double, ns, nxn, nyn, nzn);
  CF_IonIonElastic= newArr4(double, ns, nxn, nyn, nzn);
  
  CF_affected = newArr4(double, ns, nxn, nyn, nzn);
  
  // Oct 16
  C_dPx = newArr4(double, ns, nxn, nyn, nzn);
  C_dPy = newArr4(double, ns, nxn, nyn, nzn);
  C_dPz = newArr4(double, ns, nxn, nyn, nzn);

  C_dE = newArr4(double, ns, nxn, nyn, nzn);

  // for moving average
  DOI= 50;

  C_dP_DO_x = newArr5(double, DOI, ns, nxn, nyn, nzn);
  C_dP_DO_y = newArr5(double, DOI, ns, nxn, nyn, nzn);
  C_dP_DO_z = newArr5(double, DOI, ns, nxn, nyn, nzn);

  C_dE_DO = newArr5(double, DOI, ns, nxn, nyn, nzn);
  //
  C_dP_MA_x = newArr4(double, ns, nxn, nyn, nzn);
  C_dP_MA_y = newArr4(double, ns, nxn, nyn, nzn);
  C_dP_MA_z = newArr4(double, ns, nxn, nyn, nzn);

  C_dE_MA = newArr4(double, ns, nxn, nyn, nzn);

  // end for moving average
  setDiagnosticsToZero(); // set to 0 only here, then accumulated
  /* end initialisation of diagnostics */


  /* end instantiation and stuff */


  /* read init */
  MCCInitialization(&dt_MC, &IndexLastSpecies, &(CollPerc[0]));

  for (int is=0; is < ns; is++){
    cout << "CollPerc, sp " << is <<": " << CollPerc[is]<< endl;
  }


  C->hdf5_agent.set_simulation_pointers_MC(this);

}
/** deallocate */
MonteCarlo::~MonteCarlo() {
  
  if (!(Collisions == "yes")){
    return;
  }


  delete[]CollPerc;
  delete[]CollNop;
  delete[]CollRcv;

  delete[]ToMC;
  delete[]FromMC;

  
}

/** prova */
void MonteCarlo::Test(double prova) {

  if (!(Collisions == "yes")){
    return;
  }

  cout << "MC: " <<prova <<endl;
  
}


void MonteCarlo::SelectCollidingParticles(c_Solver *C){
  if (!(Collisions == "yes")){
    return;
  }

  
  
  //cout << "Inside SelectCollidingParticles" <<endl;
  
  /*if (C->col->getCase()== "DoubleGEM" and C->vct->getCoordinates(1)<  C->vct->getYLEN()/2){
    return;
    }*/

  // debug
  if (1){
    for (int is=0; is < ns; is++){
      long long nop_TOTAL=0;
      long long nop_partial= C->part[is].nop;

      MPI_Allreduce(&nop_partial, &nop_TOTAL, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      if (C->vct->getCartesian_rank()==0){
	cout <<  "MONTECARLO: Bef MC, sp " << is <<", nop: " <<nop_TOTAL <<endl; }
	
      //cout << "MONTECARLO: Bef MC, proc " << C->vct->getCartesian_rank() <<" sp " << is << " has " << C->part[is].nop << " parts" << endl;
      /*long long NOP= C->part[is].nop;
      cout <<"p 0.9*nop: " << int(0.9*NOP) << "(" << C->part[is].u[int(0.9*NOP)] << " " << C->part[is].v[int(0.9*NOP)] << " " << C->part[is].w[int(0.9*NOP)] <<")" <<endl;
      cout <<"p nop-1: " <<NOP-1 << "(" << C->part[is].u[NOP-1] << " " << C->part[is].v[NOP-1] << " " << C->part[is].w[NOP-1] <<")" <<endl;*/
     
      /*      for (long long p= 0; p< C->part[is].nop; p++){
	cout << p << ": " << C->part[is].x[p] << " "  <<C->part[is].y[p] << " " <<C->part[is].z[p] << " " <<C->part[is].u[p] << " "	<<C->part[is].v[p] << " " <<C->part[is].w[p] << " " <<C->part[is].q[p] << " " << endl;
	}*/

    }
  }
  // debug
  srand(time(NULL));
  CollNop_Total=0;

  // select particles to collide
  for (int is=0; is<ns; is++){
    long long nop_old= C->part[is].nop;
    long long np_C= round(C->part[is].nop * CollPerc[is]) ; //number of particles per species to collide
    
    CollNop[is]= np_C;
    long long nplast= C->part[is].nop-1;

    /*if (CollNop[is] > round(MAX_MC/ ns)){
      cout << "You are planning on colliding " << CollNop[is] << " particles of species " << is << endl;
      cout << "The maximum number allowed per species is " << round(MAX_MC/ns) << endl;
      cout << "Aborting now... ";
      abort();
      }*/

    if (CollNop[is] > round(MAX_MC)){
      cout << "You are planning on colliding " << CollNop[is] << " particles of species " << is << endl;
      cout << "The maximum number allowed per species is " << round(MAX_MC) << endl;
      cout << "Aborting now... ";
      abort();
    }


    int TR; // index of particle to remove
    for (long long r=0; r< np_C; r++){
     TR= rand() % nplast+1; // between 0 and nop 

   
     ToMC[CollNop_Total].X = C->part[is].x[TR];
     ToMC[CollNop_Total].Y = C->part[is].y[TR];
     ToMC[CollNop_Total].Z = C->part[is].z[TR];

     ToMC[CollNop_Total].Vx = C->part[is].u[TR];
     ToMC[CollNop_Total].Vy = C->part[is].v[TR];
     ToMC[CollNop_Total].Vz = C->part[is].w[TR];

     ToMC[CollNop_Total].qp = C->part[is].q[TR];

     if (C->part[is].TrackParticleID){
       ToMC[CollNop_Total].IDp = (long long) (C->part[is].ParticleID[TR]);
     }

     ToMC[CollNop_Total].sp_in = is;

     // at this stage, I do not need to set sp_out, CollisionType, dP_x, dP_y, dP_z, dE

     //debug      
     if (0){
       if(!( CollNop_Total%100)){
	 cout << "index: "<< CollNop_Total  << " C->part[is].u[TR]: " << C->part[is].u[TR] << " ToMC.Vx "<< ToMC[CollNop_Total].Vx <<endl;
       }
     }

     CollNop_Total++;
     // end pack into the struct


     // remove particles, with update of nop, so the net cycle selects in a range going up to the new nop
     C->part[is].del_pack(TR, &nplast);
     
     if (0){
       if(!( CollNop_Total%100)){
	 cout << "TR, vx bef removing: " <<ToMC[CollNop_Total-1].Vx << " after " <<  C->part[is].u[TR] << endl;
       }
     }
 
    } // end cycle TR

    //update particle number
    C->part[is].nop= nplast+1;

    //cout << "Species  " << is << ": old particle #: " << nop_old << ", p to collide: " << CollNop[is] << ", p left: " << C->part[is].nop << endl;
    if (nop_old - (CollNop[is]+ C->part[is].nop) )
      {
      cout << "MonteCarlo::SelectCollidingParticles, sp " << is <<" : Test not passed " << endl << "Aborting now..." << endl;
      abort();
    } 
   
  } // end for on is

  int checkNop=0;
  for (int is=0; is <ns; is ++){
    checkNop+= CollNop[is];
  }

  if (!(checkNop==CollNop_Total)){
    cout << "MonteCarlo::SelectCollidingParticles failed a check" <<endl;
    cout << "Aborting now..." <<endl;
    abort();
  }

    MPI_Barrier(MPI_COMM_WORLD);
  

}


void MonteCarlo::MCWrapper(c_Solver *C, int cycle){
  if (!(Collisions == "yes")){
    return;
  }


  // this for the MA
  DT_counter= cycle;

  /*if (C->col->getCase()== "DoubleGEM" and C->vct->getCoordinates(1)<  C->vct->getYLEN()/2){
    return;
    }*/
  
  // there may be problms casting to int
  int int_MAX_MC= int(MAX_MC);
  int int_CollNop_Total= int(CollNop_Total); 


  /*cout << "Test passing info, from C++: " << endl;
  cout << "CollNop_Total: " << CollNop_Total << endl;
  cout << "First P:" <<endl;
  printMC_struct(0, ToMC[0]);
  cout << "Last P:" <<endl;
  printMC_struct(CollNop_Total-1, ToMC[CollNop_Total-1]);*/
  // there may be a problem with CollRcv_Total being a long long and getting an int
  //MCC(int *NParBefore,MC_struct* MCCBefore,int *NParAfter,MC_struct *MCCAfter);

  //MPI_Barrier(MPI_COMM_WORLD);
  //cout << "Before MCC " << endl;
  MCC(&int_CollNop_Total, &(ToMC[0]), &CollRcv_Total, &(FromMC[0]));
  
  //MPI_Barrier(MPI_COMM_WORLD);
  //cout << "After MCC "<< endl;
  
  /*cout << "Test passing info, from C++: " << endl;
  cout << "CollRcv_Total: " << CollRcv_Total << endl;
  cout << "First P:" <<endl;
  printMC_struct(0, FromMC[0]);
  cout << "Last P:" <<endl;
  printMC_struct(CollRcv_Total-1, FromMC[CollRcv_Total-1]);*/

  for (int is=0; is < ns; is++){
    CollRcv[is]=0;
  }

  setDiagnosticsToZero(); // comment if you want to accumulate
 
  //cout << "CollRcv_Total: " << CollRcv_Total <<endl;
  for (long long p=0; p< CollRcv_Total; p++){
    
    /* add the particle to the right species */
    int sp= FromMC[p].sp_out;
    if (sp <0 || sp >ns-1){
      cout << "MonteCarlo::MCWrapper: error in processing of MC results" << endl;
      cout << "Aborting now ...";
      abort();
    }

    // ok, it's a valid species, copy back
    CollRcv[sp]++;

    C->part[sp].x[C->part[sp].nop]=  FromMC[p].X;
    C->part[sp].y[C->part[sp].nop]=  FromMC[p].Y;
    C->part[sp].z[C->part[sp].nop]=  FromMC[p].Z;

    C->part[sp].u[C->part[sp].nop]= FromMC[p].Vx;
    C->part[sp].v[C->part[sp].nop]= FromMC[p].Vy;
    C->part[sp].w[C->part[sp].nop]= FromMC[p].Vz;

    C->part[sp].q[C->part[sp].nop]= FromMC[p].qp;
    if (C->part[sp].TrackParticleID){
      //C->part[sp].ParticleID[C->part[sp].nop]=(unsigned long) (FromMC[p].IDp);
      C->part[sp].ParticleID[C->part[sp].nop]=(unsigned long) (cycle); // it just means you were affected by MC at cycle 'cycle' 
    }

    // update dP and dE
    PBP_Diagnostics(C, FromMC[p]);
        
    // update particle #
    (C->part[sp].nop)++;

    int myrank= C->vct->getCartesian_rank();
    // check if nop exceeds npmax
    if (C->part[sp].nop > C->part[sp].npmax-1 ){
      cout << "MonteCarlo::MCWrapper: while copying back particles, I have exceeded npmax" <<endl;
      cout << "Species " <<sp <<endl;
      cout << "Core " << myrank << ": [" <<xstart << " - " << xend <<"]" <<": [" <<ystart << " - " << yend <<"]" <<", [" <<zstart << " - " << zend <<"]" <<endl;
      cout << "nop: " << C->part[sp].nop << ", npmax: " << C->part[sp].npmax << endl;

      cout << "Aborting now..." <<endl;
      abort();
    }

    
  } // end for (p=0; p< CollRcv_Total; p++){
  
  PBP_Collective(C);

  int checkRCV=0;
  for (int is=0; is<ns; is++){
    checkRCV+= CollRcv[is];
  }
  if (!(checkRCV== CollRcv_Total)){
    cout << "MonteCarlo::MCWrapper: some mistake in receving MC particles" << endl;
    cout << "Aborting now" << endl;
    abort();
  }
  
  // debug                                                                                                           
  if (1){
                   
  for (int is=0; is < ns; is++){
    long long nop_TOTAL=0;
    long long nop_partial= C->part[is].nop;

    MPI_Allreduce(&nop_partial, &nop_TOTAL, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if (C->vct->getCartesian_rank()==0){
      cout << "MONTECARLO: Aft MC, sp " << is <<", nop: " <<nop_TOTAL << ", npmax (per core): " << C->part[is].npmax  <<endl;}

    //cout << "MONTECARLO: Aft MC, proc " << C->vct->getCartesian_rank() <<" sp " << is << " has " << C->part[is].nop << " parts" << endl;
    /*int NOP= C->part[is].nop;
    cout <<"p 0.9*nop: " << int(0.9*NOP) << "(" << C->part[is].u[int(0.9*NOP)] << " " << C->part[is].v[int(0.9*NOP)] << " " << C->part[is].w[int(0.9*NOP)] <<")" <<endl;
    cout <<"p nop-1: " <<NOP-1 << "(" << C->part[is].u[NOP-1] << " " << C->part[is].v[NOP-1] << " " << C->part[is].w[NOP-1] <<")" <<endl;*/

    /*    for (long long p= 0; p< C->part[is].nop; p++){
      cout << p << ": " << C->part[is].x[p] << " "  <<C->part[is].y[p] << " " <<C->part[is].z[p] << " " <<C->part[is].u[p] << " "     <<C->part[is].v[p] << " " <<C->part[is].w[p] << " " <<C->part[is].q[p] << " " << endl;
      }*/
  }
  }

  MPI_Barrier(MPI_COMM_WORLD);

}

void MonteCarlo::printMC_struct(int index, MC_struct P){

  if (!(Collisions == "yes")){
    return;
  }


  cout << "index " <<index <<": " <<endl;
  cout << "P.X: " << P.X << ", P.Y: " <<P.Y <<", P.Z: " << P.Z <<endl;
  cout << "P.Vx: " << P.Vx <<", P.Vy: " <<P.Vy << " P.Vz: " <<P.Vz <<endl;
  cout << "P.qp: " <<P.qp <<", P.IDp:  " <<P.IDp << endl;
  cout << "P.sp_in: " <<P.sp_in << ", P.sp_out: " <<P.sp_out <<endl; 
  cout << "P.CollisionType: " << P.CollisionType << endl;
  cout << "P.dP_x: " << P.dP_x <<", P.dP_y: " <<P.dP_x << "P.dP_z: " <<P.dP_z << endl;
  cout << "P.dE: " <<P.dE << endl;
  return;
}
void MonteCarlo::eIonsigmaElastic(double v2, double *CS_m2, double *CS_CU){
  // from v^2 to energy in eV  double factor= (mP/mr)*c*c*JtoeV;
  
  double energy_eV= factor_e*v2;
  
  if(energy_eV < 1.0) 
    {
      if(energy_eV < 0.2) { *CS_m2= 1./pow(10.0, 19.0 +energy_eV/.11); }
      else 
	{ 
	  *CS_m2=  9.07e-19*pow(energy_eV, 1.55)*pow(energy_eV+70.0, 1.10)/pow(14.+energy_eV, 3.25); 
	}
    }
  else 
    {
      *CS_m2=  9.07e-19*pow(energy_eV, 1.55)*pow(energy_eV+70.0, 1.10)/pow(14.+energy_eV, 3.25);
    }

  *CS_CU= *CS_m2/di_IS/di_IS;
  return;
}

void MonteCarlo::eIonigmaExc(double v2, double *CS_m2, double *CS_CU){

  double energy_eV= factor_e*v2;

  if(energy_eV < 12.0) { *CS_m2= 0.0; }  
  else{
    *CS_m2= (3.85116e-19*log(energy_eV/3.4015) -4.85227e-19)/energy_eV;
  }
  
  *CS_CU= *CS_m2/di_IS/di_IS;
  return;

}

void MonteCarlo::eIonsigmaIz(double v2, double *CS_m2, double *CS_CU){

  double energy_eV= factor_e*v2;
  
  if (energy_eV < 15.76) {*CS_m2= 0.0;}
  else if (energy_eV < 79) {
    *CS_m2= 1.7155e-18*(energy_eV-15.76)/(energy_eV*energy_eV)*log(0.0634*energy_eV);}
  else {  
    *CS_m2= 2.648e-18*(energy_eV-15.76)/(energy_eV*energy_eV)*log(0.0344*energy_eV);}

  *CS_CU= *CS_m2/di_IS/di_IS;
  return;
}


void MonteCarlo::IonIonsigmaCX(double v2, double *CS_m2, double *CS_CU){
  
  double energy_eV= factor_i*v2;

  if(energy_eV > 4.0) {
    *CS_m2=(2.0e-19 +5.5e-19/sqrt(energy_eV));
  } else{
    *CS_m2=(-2.95e-19*sqrt(energy_eV) +10.65e-19);
  }

  *CS_CU= *CS_m2/di_IS/di_IS;
  return;

}

void MonteCarlo::IonIonsigmaElastic(double v2, double *CS_m2, double *CS_CU){
  
  double energy_eV= factor_i*v2;

  if(energy_eV > 4.0) {
    *CS_m2=(1.8e-19 +4.0e-19/sqrt(energy_eV));
  } else{
    *CS_m2=(-2.0e-19*sqrt(energy_eV) +7.8e-19);
  }

  *CS_CU= *CS_m2/di_IS/di_IS;
  return;
}

string MonteCarlo::getMonteCarloPlugIn(){
  return Collisions;
}



void MonteCarlo::setZeroMCCDiagnostics(){

  for (register int kk = 0; kk < ns; kk++)
    for (register int i = 0; i < nxn; i++)
      for (register int j = 0; j < nyn; j++)
        for (register int k = 0; k < nzn; k++) {
	  // electrons only
          CF_eIonElastic[kk][i][j][k] = 0.0;
	  CF_eIonExc[kk][i][j][k] = 0.0;
	  CF_eIonIz[kk][i][j][k] = 0.0;
	  // ions only
	  CF_IonIonCX[kk][i][j][k] = 0.0;
	  CF_IonIonElastic[kk][i][j][k] = 0.0;

	  CF_affected[kk][i][j][k] = 0.0;

        }

}


void MonteCarlo::communicateGhostP2G( int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology3D * vct) {
  // interpolate adding common nodes among processors                                                                                            
  communicateInterp(nxn, nyn, nzn, ns, CF_eIonElastic, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, CF_eIonExc, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, CF_eIonIz, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, CF_IonIonCX, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, CF_IonIonElastic, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, ns, CF_affected, 0, 0, 0, 0, 0, 0, vct);  


  communicateNode_P(nxn, nyn, nzn, CF_eIonElastic, ns, vct);
  communicateNode_P(nxn, nyn, nzn, CF_eIonExc, ns, vct);
  communicateNode_P(nxn, nyn, nzn, CF_eIonIz, ns, vct);
  communicateNode_P(nxn, nyn, nzn, CF_IonIonCX, ns, vct);
  communicateNode_P(nxn, nyn, nzn, CF_IonIonElastic, ns, vct);
  communicateNode_P(nxn, nyn, nzn, CF_affected, ns, vct);

}

void MonteCarlo::addCF_affected(double weight[][2][2], int X, int Y, int Z, int is){
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	CF_affected[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL;
	//CF_affected[is][X - i][Y - j][Z - k] += weight[i][j][k] ; 
}

void MonteCarlo::addCF_eIonElastic(double weight[][2][2], int X, int Y, int Z, int is) {

  double factor= ReducedC*GasN; // the c is to de-normalise the sqrt(v2), which at the moment is still in code units
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	CF_eIonElastic[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL* factor;
  //CF_eIonElastic[is][X - i][Y - j][Z - k] += weight[i][j][k] * factor;
}

void MonteCarlo::addCF_eIonExc(double weight[][2][2], int X, int Y, int Z, int is) {

  double factor= ReducedC*GasN; // the c is to de-normalise the sqrt(v2), which at the moment is still in code units
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	CF_eIonExc[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL* factor;
  //CF_eIonExc[is][X - i][Y - j][Z - k] += weight[i][j][k] * factor;
}

void MonteCarlo::addCF_eIonIz(double weight[][2][2], int X, int Y, int Z, int is) {

  double factor= ReducedC*GasN; // the c is to de-normalise the sqrt(v2), which at the moment is still in code units
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	CF_eIonIz[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL* factor;
  //CF_eIonIz[is][X - i][Y - j][Z - k] += weight[i][j][k] * factor;
}

void MonteCarlo::addCF_IonIonCX(double weight[][2][2], int X, int Y, int Z, int is) {

  double factor= ReducedC*GasN; // the c is to de-normalise the sqrt(v2), which at the moment is still in code units
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	CF_IonIonCX[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL* factor;
  //CF_IonIonCX[is][X - i][Y - j][Z - k] += weight[i][j][k] *  factor;
}

void MonteCarlo::addCF_IonIonElastic(double weight[][2][2], int X, int Y, int Z, int is) {

  double factor= ReducedC*GasN; // the c is to de-normalise the sqrt(v2), which at the moment is still in code units
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
	CF_IonIonElastic[is][X - i][Y - j][Z - k] += weight[i][j][k] * invVOL* factor;
  //CF_IonIonElastic[is][X - i][Y - j][Z - k] += weight[i][j][k] * factor;
}


double ****MonteCarlo::getCF_eIonElastic() {
  return (CF_eIonElastic);
}
double ****MonteCarlo::getCF_eIonExc() {
  return (CF_eIonExc);
}
double ****MonteCarlo::getCF_eIonIz() {
  return (CF_eIonIz);
}
double ****MonteCarlo::getCF_IonIonCX() {
  return (CF_IonIonCX);
}
double ****MonteCarlo::getCF_IonIonElastic() {
  return (CF_IonIonElastic);
}
double ****MonteCarlo::getCF_affected() {
  return (CF_affected);
}

void MonteCarlo::PBP_Diagnostics(c_Solver* C, MC_struct P){


  const int ix = 2 + int (floor((P.X - xstart) * inv_dx));
  const int iy = 2 + int (floor((P.Y - ystart) * inv_dy));
  const int iz = 2 + int (floor((P.Z - zstart) * inv_dz));
  double temp[2][2][2];
  double xi[2], eta[2], zeta[2];
  xi[0] = P.X - C->grid->getXN(ix - 1, iy, iz);
  eta[0] = P.Y - C->grid->getYN(ix, iy - 1, iz);
  zeta[0] = P.Z - C->grid->getZN(ix, iy, iz - 1);
  xi[1] = C->grid->getXN(ix, iy, iz) - P.X;
  eta[1] = C->grid->getYN(ix, iy, iz) - P.Y;
  zeta[1] = C->grid->getZN(ix, iy, iz) - P.Z;
  double weight[2][2][2];
  for (int ii = 0; ii < 2; ii++)
    for (int jj = 0; jj < 2; jj++)
      for (int kk = 0; kk < 2; kk++) {
	weight[ii][jj][kk] = P.qp * xi[ii] * eta[jj] * zeta[kk] * invVOL;
      }


  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){

        C_dPx[P.sp_out][ix - i][iy - j][iz - k] += weight[i][j][k] * invVOL* P.dP_x;
	C_dPy[P.sp_out][ix - i][iy - j][iz - k] += weight[i][j][k] * invVOL* P.dP_y;
	C_dPz[P.sp_out][ix - i][iy - j][iz - k] += weight[i][j][k] * invVOL* P.dP_z;

	C_dE[P.sp_out][ix - i][iy - j][iz - k] += weight[i][j][k] * invVOL* P.dE;
	/*C_dPx[P.sp_out][ix - i][iy - j][iz - k] += weight[i][j][k] * P.dP_x;
	C_dPy[P.sp_out][ix - i][iy - j][iz - k] += weight[i][j][k] * P.dP_y;
	C_dPz[P.sp_out][ix - i][iy - j][iz - k] += weight[i][j][k] * P.dP_z;

	C_dE[P.sp_out][ix - i][iy - j][iz - k] += weight[i][j][k] * P.dE;*/
	
      }
  return;
}

double **** MonteCarlo::getC_dPx(){
  return C_dPx;
}
double **** MonteCarlo::getC_dPy(){
  return C_dPy;
}
double **** MonteCarlo::getC_dPz(){
  return C_dPz;
}
double **** MonteCarlo::getC_dE(){
  return C_dE;
}

// for MA
double **** MonteCarlo::getC_dP_MA_x(){
  return C_dP_MA_x;
}
double **** MonteCarlo::getC_dP_MA_y(){
  return C_dP_MA_y;
}
double **** MonteCarlo::getC_dP_MA_z(){
  return C_dP_MA_z;
}
double **** MonteCarlo::getC_dE_MA(){
  return C_dE_MA;
}

void MonteCarlo::setDiagnosticsToZero(){
  eqValue(0.0, C_dPx, ns, nxn, nyn, nzn);
  eqValue(0.0, C_dPy, ns, nxn, nyn, nzn);
  eqValue(0.0, C_dPz, ns, nxn, nyn, nzn);
  eqValue(0.0, C_dE, ns, nxn, nyn, nzn);

  // for MA
  eqValue(0.0, C_dP_MA_x, ns, nxn, nyn, nzn);
  eqValue(0.0, C_dP_MA_y, ns, nxn, nyn, nzn);
  eqValue(0.0, C_dP_MA_z, ns, nxn, nyn, nzn);

  eqValue(0.0, C_dE_MA, ns, nxn, nyn, nzn);
}

void MonteCarlo::communicateGhostP2G(VirtualTopology3D * vct){
  // interpolate adding common nodes among processors                
  for (int is=0; is<ns; is++){

  communicateInterp(nxn, nyn, nzn, is, C_dPx, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, is, C_dPy, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, is, C_dPz, 0, 0, 0, 0, 0, 0, vct);
  communicateInterp(nxn, nyn, nzn, is, C_dE, 0, 0, 0, 0, 0, 0, vct);

  communicateNode_P(nxn, nyn, nzn, C_dPx, is, vct);
  communicateNode_P(nxn, nyn, nzn, C_dPy, is, vct);
  communicateNode_P(nxn, nyn, nzn, C_dPz, is, vct);
  communicateNode_P(nxn, nyn, nzn, C_dE, is, vct);

  }
}

void MonteCarlo::PBP_Collective(c_Solver* C){

  // MA of the diagnostics, species per species
  int nn= DT_counter% DOI;
  //cout << "nn: " << nn << endl;

  communicateGhostP2G(C->vct);

  for (int is=0; is<ns; is ++){
    C->EMf->smooth(Smooth, C_dPx[is], 1, C->grid, C->vct);
    C->EMf->smooth(Smooth, C_dPy[is], 1, C->grid, C->vct); 
    C->EMf->smooth(Smooth, C_dPz[is], 1, C->grid, C->vct);

    C->EMf->smooth(Smooth, C_dE[is], 1, C->grid, C->vct);

    // calculate the MA
    SMA(C_dP_MA_x[is], C_dPx[is], C_dP_DO_x[nn][is], DOI, nxn, nyn, nzn);
    SMA(C_dP_MA_y[is], C_dPy[is], C_dP_DO_y[nn][is], DOI, nxn, nyn, nzn);
    SMA(C_dP_MA_z[is], C_dPz[is], C_dP_DO_z[nn][is], DOI, nxn, nyn, nzn);

    SMA(C_dE_MA[is], C_dE[is], C_dE_DO[nn][is], DOI, nxn, nyn, nzn);

    // set the DO
    eq(C_dP_DO_x[nn][is], C_dPx[is], nxn, nyn, nzn);
    eq(C_dP_DO_y[nn][is], C_dPy[is], nxn, nyn, nzn);
    eq(C_dP_DO_z[nn][is], C_dPz[is], nxn, nyn, nzn);

    eq(C_dE_DO[nn][is], C_dE[is], nxn, nyn, nzn);
  }
    
  return;
}
