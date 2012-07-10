/***************************************************************************
                          EMfields3D.h  -  ElectroMagnetic fields definition
                             -------------------
    begin             : May 2008
    copyright         : KU Leuven
    developers        : Stefano Markidis, Giovanni Lapenta
***************************************************************************/

#ifndef EMfields3D_H
#define EMfields3D_H


#include <iostream>

#include <math.h>
#include <mpi.h>


#include "../utility/Alloc.h"
#include "../mathlib/Basic.h"
#include "../utility/TransArraySpace3D.h"
#include "../solvers/CG.h"
#include "../solvers/GMRES.h"


using std::cout;
using std::cerr;
using std::endl;



/**
*  Electromagnetic fields and sources defined for each local grid, and for an implicit maxwell's solver 
*   
* @date May 2008
* @par Copyright:
* (C) 2008 KUL
* @author Stefano Markidis, Giovanni Lapenta.
* @version 3.0
*/

class EMfields3D : public Field{
  public:
  /** constructor */
  EMfields3D(CollectiveIO *col,Grid *grid);
  /** destructor */
  ~EMfields3D();

 /** initialize the electromagnetic fields with constant values */
  void init(VirtualTopology3D *vct, Grid *grid);
  /** init beam */
  void initBEAM(VirtualTopology3D *vct, Grid *grid, double x_center, double y_center, double z_center, double radius);
  /** initialize GEM challenge */
  void initGEM(VirtualTopology3D *vct, Grid *grid);
  /** initialize GEM challenge with no Perturbation */
  void initGEMnoPert(VirtualTopology3D *vct, Grid *grid);
  /**  Init Force Free (JxB=0) */
  void initForceFree(VirtualTopology3D *vct, Grid *grid);
  /** initialized with rotated magnetic field */
  void initEM_rotate(VirtualTopology3D *vct, Grid *grid, double B, double theta);
  /** add a perturbattion to charge density */
  void AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid *grid);
  /** add a perturbattion to the EM field */
  void AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid *grid);
  
  
  /** Calculate Electric field using the implicit Maxwell solver */
  void calculateE(Grid *grid, VirtualTopology3D *vct);
  /** Image of Poisson Solver (for SOLVER)*/
  void PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology3D *vct);
  /** Image of Maxwell Solver (for Solver) */
  void MaxwellImage(double *im, double *vector, Grid *grid,VirtualTopology3D *vct);     
  /** Maxwell source term (for SOLVER) */
  void MaxwellSource(double *bkrylov, Grid *grid, VirtualTopology3D *vct);
  /** Calculate Magnetic field with the implicit solver: calculate B defined on nodes
  With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
  void calculateB(Grid *grid, VirtualTopology3D *vct);
  /** fix B on the boundary for gem challange*/
  void fixBgem(Grid *grid, VirtualTopology3D *vct);
  /** fix B on the boundary for gem challange*/
  void fixBforcefree(Grid *grid, VirtualTopology3D *vct);
  
  /** Calculate the three components of Pi(implicit pressure) cross image vector */
  void PIdot(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid *grid);
  /** Calculate the three components of mu (implicit permeattivity) cross image vector */
  void MUdot(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid *grid);
  /** Calculate rho hat, Jx hat, Jy hat, Jz hat */
  void calculateHatFunctions(Grid *grid, VirtualTopology3D *vct);
 

  /** communicate ghost for densities and interp rho from node to center */
  void interpDensitiesN2C(VirtualTopology3D *vct, Grid *grid);
  /** set to 0 all the densities fields */
  void setZeroDensities();
  /** Sum rhon over species */
  void sumOverSpecies(VirtualTopology3D *vct);
  /** Sum current over different species */
  void sumOverSpeciesJ();
  /** Smoothing after the interpolation**/
  void smooth(double value ,double ***vector,bool type, Grid *grid, VirtualTopology3D *vct);
  /** SPECIES: Smoothing after the interpolation for species fields**/
  void smooth(double value ,double ****vector,int is, bool type, Grid *grid, VirtualTopology3D *vct);
  /** smooth the electric field */
  void smoothE(double value, VirtualTopology3D *vct);
  
  /** communicate ghost for grid -> Particles interpolation */
  void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology3D *vct);
    /** add an amount of charge density to charge density field at node X,Y,Z */
  void addRho(double weight[][2][2], int X, int Y,int Z, int is);
  /** add an amount of current density - direction X to current density field at node X,Y,Z */
  void addJx(double weight[][2][2], int X, int Y, int Z,int is);
  /** add an amount of current density - direction Y to current density field at node X,Y,Z */
  void addJy(double weight[][2][2], int X, int Y, int Z, int is);
  /** add an amount of current density - direction Z to current density field at node X,Y,Z */
  void addJz(double weight[][2][2], int X, int Y, int Z, int is);
 
  /** add an amount of pressure density - direction XX to current density field at node X,Y,Z */
  void addPxx(double weight[][2][2], int X, int Y, int Z, int is);
  /** add an amount of pressure density - direction XY to current density field at node X,Y,Z */
  void addPxy(double weight[][2][2], int X, int Y, int Z,int is);
  /** add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
  void addPxz(double weight[][2][2], int X, int Y, int Z,int is);
  /** add an amount of pressure density - direction YY to current density field at node X,Y,Z */
  void addPyy(double weight[][2][2], int X, int Y, int Z, int is);
  /** add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
  void addPyz(double weight[][2][2], int X, int Y, int Z, int is);
  /** add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
  void addPzz(double weight[][2][2], int X, int Y, int Z,int is);
  
   /** adjust densities on boundaries that are not periodic */
  void adjustNonPeriodicDensities(int is, VirtualTopology3D *vct);

  
  /** Perfect conductor boundary conditions LEFT wall */
  void perfectConductorLeft(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY,double ***vectorZ,int dir,Grid *grid);
  /** Perfect conductor boundary conditions RIGHT wall */
  void perfectConductorRight(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid);
  /** Perfect conductor boundary conditions for source LEFT wall*/
  void perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);
  /** Perfect conductor boundary conditions for source RIGHT wall*/
  void perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);
  
  /** Calculate the sysceptibility tensor on the boundary */
  void sustensorRightY(double** susxy,double** susyy,double** suszy);
  void sustensorLeftY(double** susxy,double** susyy,double** suszy);
   
  /** get Potential array */
  double*** getPHI();
  /** get Electric Field component X defined on node(indexX,indexY,indexZ) */
  double &getEx(int indexX, int indexY, int indexZ) const;
  /** get Electric field X component array */
  double*** getEx() ;
   /** get Electric Field component Y defined on node(indexX,indexY,indexZ) */
  double &getEy(int indexX, int indexY, int indexZ) const;
  /** get Electric field Y component array */
  double*** getEy() ;
  /** get Electric Field component Z defined on node(indexX,indexY,indexZ) */
  double &getEz(int indexX, int indexY, int indexZ) const;
  /** get Electric field Z component array */
  double*** getEz() ;
  /** get Magnetic Field component X defined on node(indexX,indexY,indexZ) */
  double &getBx(int indexX, int indexY, int indexZ) const;
  /** get Magnetic field X component array */
  double*** getBx();
  /** get Magnetic Field component Y defined on node(indexX,indexY,indexZ) */
  double &getBy(int indexX, int indexY, int indexZ) const;
  /** get Magnetic field Y component array */
  double*** getBy();
  /** get Magnetic Field component Z defined on node(indexX,indexY,indexZ) */
  double &getBz(int indexX, int indexY, int indexZ) const;
  /** get Magnetic field Z component array */
  double*** getBz() ;
  /** get density on cell(indexX,indexY,indexZ)  */
  double &getRHOc(int indexX, int indexY, int indexZ) const;
  /** get density array on center cell */
  double*** getRHOc() ;
  /** get density on nodes(indexX,indexY,indexZ)  */
  double &getRHOn(int indexX, int indexY, int indexZ) const;
  /** get density array on nodes  */
  double*** getRHOn() ;
  /** SPECIES: get density on nodes(indexX,indexY,indexZ)*/
  double &getRHOns(int indexX, int indexY,int indexZ,int is) const;
  /** SPECIES: get density on center cell(indexX,indexY,indexZ)*/
  double &getRHOcs(int indexX, int indexY,int indexZ,int is) const;
  /** SPECIES: get density array on nodes*/
  double**** getRHOns() ;



  /** get pressure tensor XX for species */
  double**** getpXXsn() ;
  /** get pressure tensor XY for species */
  double**** getpXYsn() ;
  /** get pressure tensor XZ for species*/
  double**** getpXZsn() ;
  /** get pressure tensor YY for species */
  double**** getpYYsn() ;
  /** get pressure tensor YZ for species */
  double**** getpYZsn() ;
  /** get pressure tensor ZZ for species */
  double**** getpZZsn() ;
    
  /** get Jx(X,Y,Z)  */
  double &getJx(int indexX, int indexY, int indexZ) const;
  /** get current -Direction X */
  double*** getJx();
  /** get Jxs(X,Y,Z,is)  */
  double &getJxs(int indexX, int indexY, int indexZ, int is) const;
  /** SPECIES: get current -Direction X */
  double**** getJxs();
  /** get Jy(X,Y,Z)  */
  double &getJy(int indexX, int indexY, int indexZ) const;
  /** get current -Direction Y */
  double*** getJy();
  /** get Jys(X,Y,Z,is)  */
  double &getJys(int indexX, int indexY, int indexZ, int is) const;
  /** SPECIES: get current -Direction Y */
  double**** getJys();
  /** get Jz(X,Y,Z)  */
  double &getJz(int indexX, int indexY, int indexZ) const;
  /** get current -Direction Z */
  double*** getJz();
  /** get Jzs(X,Y,Z,is)  */
  double &getJzs(int indexX, int indexY, int indexZ, int is) const;
  /** SPECIES: get current -Direction Z */
  double**** getJzs();
  /** get the electric field energy */
  double getEenergy();
  /** get the magnetic field energy */
  double getBenergy();
  
  
  /** print electromagnetic fields info */
  void print(void) const;
  
/**********************************
// VARIABLES
//**************************** */
private:
  /** light speed */
  double c;
  /* 4*PI for normalization */
  double FourPI;
  /** time step */
  double dt;
  /** decentering parameter */
  double th;
  /** Smoothing value*/
  double Smooth;
  /** delt = c*th*dt */
  double delt;
  /** number of particles species */
  int ns;
  /** GEM challenge parameters */
  double B0x, B0y, B0z, delta;
  /** charge to mass ratio array for different species */
  double *qom;
  

  // KEEP IN MEMORY GUARD CELLS ARE INCLUDED
  /** number of cells - X direction, including + 2 (guard cells) */
  int nxc;
  /** number of nodes - X direction, including + 2 extra nodes for guard cells */
  int nxn;
  /** number of cell - Y direction, including + 2 (guard cells) */
  int nyc;
  /** number of nodes - Y direction, including + 2 extra nodes for guard cells */
  int nyn;
  /** number of cell - Z direction, including + 2 (guard cells) */
  int nzc;
  /** number of nodes - Z direction, including + 2 extra nodes for guard cells */
  int nzn;
  /** local grid boundaries coordinate  */
  double xStart, xEnd, yStart, yEnd, zStart, zEnd;
  /** grid spacing */
  double dx, dy, dz, invVOL;
  /** simulation box length - X direction   */
  double Lx;
  /** simulation box length - Y direction   */
  double Ly;
  /** simulation box length - Z direction   */
  double Lz;
  

  /** PHI: electric potential (indexX, indexY, indexZ), defined on central points between nodes*/
  double ***PHI;
  /** Ex: electric field X-component (indexX, indexY, indexZ), defined on nodes */
  double ***Ex;
  /** Exth: implicit electric field X-component (indexX, indexY, indexZ), defined on nodes */
  double ***Exth;
  /** Ey: electric field Y-component (indexX, indexY, indexZ), defined on nodes */
  double ***Ey;
  /** Eyth: implicit electric field Y-component (indexX, indexY, indexZ), defined on nodes */
  double ***Eyth;
  /** Ez: electric field Z-component (indexX, indexY, indexZ, #species), defined on nodes */
  double ***Ez;
  /** Ezth: implicit electric field Z-component (indexX, indexY, indexZ), defined on nodes */
  double ***Ezth;
  /** Bxc: magnetic field X-component (indexX, indexY, indexZ), defined on central points between nodes*/
  double ***Bxc;
  /** Byc: magnetic field Y-component (indexX, indexY, indexZ), defined on central points between nodes*/
  double ***Byc;
  /** Bzc: magnetic field Z-component (indexX, indexY, indexZ), defined on central points between nodes*/
  double ***Bzc;
  /** Bxn: magnetic field X-component (indexX, indexY, indexZ), defined on nodes*/
  double ***Bxn;
  /** Byn: magnetic field Y-component (indexX, indexY, indexZ), defined on nodes*/
  double ***Byn;
  /** Bzn: magnetic field Z-component (indexX, indexY, indexZ), defined on nodes*/
  double ***Bzn;
 
  //*************************************
  // TEMPORARY ARRAY
  //************************************
   /**some temporary arrays (for calculate hat functions)*/
  double*** tempXC;
  double*** tempYC;
  double*** tempZC;
  double*** tempXN;
  double*** tempYN;
  double*** tempZN;
  /** other temporary arrays (in MaxwellSource) */			  
  double*** tempC; 
  double*** tempX;
  double*** tempY;
  double*** tempZ;
  double*** temp2X;
  double*** temp2Y;
  double*** temp2Z;
  /** and some for MaxwellImage */
  double*** imageX; 
  double*** imageY; 
  double*** imageZ;
  double*** Dx;
  double*** Dy;
  double*** Dz;
  double*** vectX;
  double*** vectY;
  double*** vectZ;
  double*** divC;
 
 
  /*********************************************************************************
  /*************                SOURCES                                           **
  /********************************************************************************/

  /** Charge density, defined on central points of the cell */
  double*** rhoc;
  /** Charge density, defined on nodes */
  double*** rhon;
  /** Implicit charge density, defined on central points of the cell */
  double***  rhoh;
  /** SPECIES: charge density for each species, defined on nodes */
  double****  rhons;
  /** SPECIES: charge density for each species, defined on central points of the cell */
  double**** rhocs;
  /** Current density component-X, defined on nodes */
  double***  Jx;
  /** Current density component-Y, defined on nodes */
  double***  Jy;
  /** Current density component-Z, defined on nodes */
  double***  Jz;
  /** Implicit current density X-component, defined on nodes */
  double***  Jxh;
  /** Implicit current density Y-component, defined on nodes */
  double***  Jyh;
  /** Implicit current density Z-component, defined on nodes */
  double***  Jzh;
  /** SPECIES: current density component-X for species, defined on nodes */
  double**** Jxs;
  /** SPECIES: current density component-Y for species, defined on nodes */
  double****  Jys;
  /** SPECIES: current density component-Z for species, defined on nodes */
  double****  Jzs;
  
  /** SPECIES: pressure tensor component-XX, defined on nodes */
  double**** pXXsn;
  /** SPECIES: pressure tensor component-XY, defined on nodes */
  double ****pXYsn;
  /** SPECIES: pressure tensor component-XZ, defined on nodes */
  double ****pXZsn;
  /** SPECIES: pressure tensor component-XZ, defined on nodes */
  double ****pYYsn;
  /** SPECIES: pressure tensor component-YZ, defined on nodes */
  double ****pYZsn;
  /** SPECIES: pressure tensor component-ZZ, defined on nodes */
  double ****pZZsn;
  

  /** Field Boundary Condition
   0 =  Dirichlet Boundary Condition: specifies the value to take on the boundary of the domain
   1 =  Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain
   2 =  Periodic condition
  */

  /** Boundary Condition Electrostatic Potential: FaceXright */
  int bcPHIfaceXright;
  /** Boundary Condition Electrostatic Potential:FaceXleft */
  int bcPHIfaceXleft;
  /** Boundary Condition Electrostatic Potential:FaceYright */
  int bcPHIfaceYright;
  /** Boundary Condition Electrostatic Potential:FaceYleft */
  int bcPHIfaceYleft;
  /** Boundary Condition Electrostatic Potential:FaceZright */
  int bcPHIfaceZright;
  /** Boundary Condition Electrostatic Potential:FaceZleft */
  int bcPHIfaceZleft;

  /** Boundary condition for electric field
  0 = perfect conductor
  1 = magnetic mirror
  */
  /** Boundary Condition EM Field: FaceXright */
  int bcEMfaceXright;
  /** Boundary Condition EM Field: FaceXleft */
  int bcEMfaceXleft;
  /** Boundary Condition EM Field: FaceYright */
  int bcEMfaceYright;
  /** Boundary Condition EM Field: FaceYleft */
  int bcEMfaceYleft;
  /** Boundary Condition EM Field: FaceZright */
  int bcEMfaceZright;
  /** Boundary Condition EM Field: FaceZleft */
  int bcEMfaceZleft;

  
  /** GEM Challenge background ion */
  double *rhoINIT; 
  /** Drift of the species */
  bool *DriftSpecies;
  
  /** boolean for divergence cleaning */
  bool PoissonCorrection;
  /** RESTART BOOLEAN */
  int restart1 ; 
  /** String with the directory for the restart file */
  string RestartDirName;
  
   /** CG tolerance criterium for stopping iterations */
  double CGtol;
  /**  GMRES tolerance criterium for stopping iterations */
  double GMREStol;
    
};

/** constructor */
inline EMfields3D::EMfields3D(CollectiveIO *col,Grid *grid){
  nxc = grid->getNXC();
  nxn = grid->getNXN();
  nyc = grid->getNYC();
  nyn = grid->getNYN();
  nzc = grid->getNZC();
  nzn = grid->getNZN();
  dx = grid->getDX();
  dy = grid ->getDY();
  dz = grid->getDZ();
  invVOL = grid->getInvVOL();
  xStart = grid->getXstart();
  xEnd   = grid->getXend();
  yStart = grid->getYstart();
  yEnd   = grid->getYend();
  zStart = grid->getZstart();
  zEnd   = grid->getZend();
  Lx = col->getLx();
  Ly = col->getLy();
  Lz = col->getLz();
  ns  = col->getNs();
  c = col->getC();
  dt = col->getDt();
  th = col->getTh();
  delt = c*th*dt;
  PoissonCorrection = true;
  CGtol = col->getCGtol();
  GMREStol = col->getGMREStol();
  qom = new double[ns];
  for (int i=0; i < ns;i++)
    qom[i] = col->getQOM(i);
  // boundary conditions: PHI and EM fields
  bcPHIfaceXright  = col->getBcPHIfaceXright();
  bcPHIfaceXleft   = col->getBcPHIfaceXleft();
  bcPHIfaceYright  = col->getBcPHIfaceYright();
  bcPHIfaceYleft   = col->getBcPHIfaceYleft();
  bcPHIfaceZright  = col->getBcPHIfaceZright();
  bcPHIfaceZleft   = col->getBcPHIfaceZleft();

  bcEMfaceXright   = col->getBcEMfaceXright();
  bcEMfaceXleft    = col->getBcEMfaceXleft();
  bcEMfaceYright   = col->getBcEMfaceYright();
  bcEMfaceYleft    = col->getBcEMfaceYleft();
  bcEMfaceZright   = col->getBcEMfaceZright();
  bcEMfaceZleft    = col->getBcEMfaceZleft();
  // GEM challenge parameters
  B0x = col->getB0x();
  B0y = col->getB0y();
  B0z = col->getB0z();
  delta = col->getDelta();
  Smooth = col->getSmooth();
  // get the density background for the gem Challange
  rhoINIT = new double[ns];
  DriftSpecies = new bool[ns];
  for (int i=0; i < ns;i++){
	    rhoINIT[i] = col->getRHOinit(i);
	    if ((fabs(col->getW0(i))!=0) || (fabs(col->getU0(i))!=0)) // GEM and LHDI
		  DriftSpecies[i] = true;
		else
		  DriftSpecies[i] = false;
  }
  /** parameters for GEM challenge */
  FourPI =16*atan(1.0);
  /** Restart */
  restart1 = col->getRestart_status();
  RestartDirName = col->getRestartDirName();
 // arrays allocation: nodes
 Ex = newArr3(double,nxn,nyn,nzn); Ey = newArr3(double,nxn,nyn,nzn); Ez = newArr3(double,nxn,nyn,nzn);
 Exth = newArr3(double,nxn,nyn,nzn); Eyth = newArr3(double,nxn,nyn,nzn); Ezth = newArr3(double,nxn,nyn,nzn);
 Bxn  = newArr3(double,nxn,nyn,nzn);  Byn   = newArr3(double,nxn,nyn,nzn); Bzn   = newArr3(double,nxn,nyn,nzn);
 rhon  = newArr3(double,nxn,nyn,nzn);
 Jx    = newArr3(double,nxn,nyn,nzn); Jy    = newArr3(double,nxn,nyn,nzn); Jz    = newArr3(double,nxn,nyn,nzn);
 Jxh   = newArr3(double,nxn,nyn,nzn); Jyh   = newArr3(double,nxn,nyn,nzn); Jzh   = newArr3(double,nxn,nyn,nzn);
 // involving species
 rhons = newArr4(double,ns,nxn,nyn,nzn); rhocs = newArr4(double,ns,nxc,nyc,nzn);
 Jxs   = newArr4(double,ns,nxn,nyn,nzn); Jys   = newArr4(double,ns,nxn,nyn,nzn); Jzs   = newArr4(double,ns,nxn,nyn,nzn);
 pXXsn  = newArr4(double,ns,nxn,nyn,nzn); pXYsn  = newArr4(double,ns,nxn,nyn,nzn); pXZsn  = newArr4(double,ns,nxn,nyn,nzn); pYYsn  = newArr4(double,ns,nxn,nyn,nzn); pYZsn  = newArr4(double,ns,nxn,nyn,nzn); pZZsn  = newArr4(double,ns,nxn,nyn,nzn);
 // arrays allocation: central points 
 PHI    = newArr3(double,nxc,nyc,nzc);
 Bxc    = newArr3(double,nxc,nyc,nzc); Byc    = newArr3(double,nxc,nyc,nzc); Bzc    = newArr3(double,nxc,nyc,nzc);
 rhoc   = newArr3(double,nxc,nyc,nzc); rhoh   = newArr3(double,nxc,nyc,nzc);
	
 // temporary arrays
 tempXC = newArr3(double,nxc,nyc,nzc); tempYC = newArr3(double,nxc,nyc,nzc); tempZC = newArr3(double,nxc,nyc,nzc);
 
 tempXN = newArr3(double,nxn,nyn,nzn); tempYN = newArr3(double,nxn,nyn,nzn); tempZN = newArr3(double,nxn,nyn,nzn);
 tempC  = newArr3(double,nxc,nyc,nzn);
 tempX  = newArr3(double,nxn,nyn,nzn); tempY  = newArr3(double,nxn,nyn,nzn); tempZ  = newArr3(double,nxn,nyn,nzn);
 temp2X = newArr3(double,nxn,nyn,nzn); temp2Y = newArr3(double,nxn,nyn,nzn); temp2Z = newArr3(double,nxn,nyn,nzn);
 imageX = newArr3(double,nxn,nyn,nzn); imageY = newArr3(double,nxn,nyn,nzn); imageZ = newArr3(double,nxn,nyn,nzn);
 Dx = newArr3(double,nxn,nyn,nzn); Dy = newArr3(double,nxn,nyn,nzn); Dz = newArr3(double,nxn,nyn,nzn);
 vectX  = newArr3(double,nxn,nyn,nzn); vectY  = newArr3(double,nxn,nyn,nzn); vectZ  = newArr3(double,nxn,nyn,nzn);
 divC  = newArr3(double,nxc,nyc,nzc); 
}

/** Calculate Electric field with the implicit solver: the Maxwell solver method is called here */
inline void EMfields3D::calculateE(Grid *grid, VirtualTopology3D *vct){
  if (vct->getCartesian_rank() ==0)
    cout << "*** E CALCULATION ***" << endl;
  double ***divE = newArr3(double,nxc,nyc,nzc);
  double ***gradPHIX = newArr3(double,nxn,nyn,nzn);
  double ***gradPHIY = newArr3(double,nxn,nyn,nzn);
  double ***gradPHIZ = newArr3(double,nxn,nyn,nzn);

  double *xkrylov = new double[3*(nxn-2)*(nyn-2)*(nzn-2)]; // 3 E components
  double *bkrylov = new double[3*(nxn-2)*(nyn-2)*(nzn-2)]; // 3 components

  double *xkrylovPoisson = new double[(nxc-2)*(nyc-2)*(nzc-2)];
  double *bkrylovPoisson = new double[(nxc-2)*(nyc-2)*(nzc-2)];
  // set to zero all the stuff 
  eqValue (0.0, xkrylov, 3*(nxn-2)*(nyn-2)*(nzn-2)); 
  eqValue (0.0, xkrylovPoisson, (nxc-2)*(nyc-2)*(nzc-2));
  eqValue (0.0, bkrylov, 3*(nxn-2)*(nyn-2)*(nzn-2));
  eqValue (0.0, divE,nxc,nyc,nzc);
  eqValue (0.0, tempC,nxc,nyc,nzc);
  eqValue (0.0, gradPHIX,nxn,nyn,nzn);
  eqValue (0.0, gradPHIY,nxn,nyn,nzn);
  eqValue (0.0, gradPHIZ,nxn,nyn,nzn);
  // Adjust E calculating laplacian(PHI) = div(E)  -4*PI*rho DIVERGENCE CLEANING
  if (PoissonCorrection){
    if (vct->getCartesian_rank() ==0)
      cout << "*** DIVERGENCE CLEANING ***" << endl;
      grid->divN2C(divE,Ex,Ey,Ez);
      scale(tempC,rhoc,-FourPI,nxc,nyc,nzc);
      sum(divE,tempC,nxc,nyc,nzc);
      // move to krylov space  
      phys2solver(bkrylovPoisson,divE,nxc,nyc,nzc);
      // use conjugate gradient first
      if(!CG(xkrylovPoisson,(nxc-2)*(nyc-2)*(nzc-2),bkrylovPoisson,3000,CGtol,&Field::PoissonImage,grid,vct,this)){
         if (vct->getCartesian_rank() ==0) 
		  cout << "CG not Converged. Trying with GMRes. Consider to increase the number of the CG iterations" << endl;
	     eqValue(0.0,xkrylovPoisson,(nxc-2)*(nyc-2)*(nzc-2));
		 GMRES(&Field::PoissonImage,xkrylovPoisson,(nxc-2)*(nyc-2)*(nzc-2),bkrylovPoisson,20,200,GMREStol,grid,vct,this);
      }	  
	  solver2phys(PHI,xkrylovPoisson,nxc,nyc,nzc);
	  communicateCenterBC(nxc,nyc,nzc,PHI,2,2,2,2,2,2,vct);
	  // calculate the gradient
	  grid->gradC2N(gradPHIX,gradPHIY,gradPHIZ,PHI);
	  // sub
	  sub(Ex,gradPHIX,nxn,nyn,nzn);
	  sub(Ey,gradPHIY,nxn,nyn,nzn);
	  sub(Ez,gradPHIZ,nxn,nyn,nzn);
      
  }  // end of divergence cleaning 
  if (vct->getCartesian_rank() ==0)
    cout << "*** MAXWELL SOLVER ***" << endl;
  // prepare the source	
  MaxwellSource(bkrylov,grid,vct);
  phys2solver(xkrylov,Ex,Ey,Ez,nxn,nyn,nzn);
  // solver
  GMRES(&Field::MaxwellImage, xkrylov, 3*(nxn-2)*(nyn-2)*(nzn-2),bkrylov,20, 200,GMREStol, grid, vct, this);
  // move from krylov space to physical space
  solver2phys(Exth,Eyth,Ezth,xkrylov,nxn,nyn,nzn);
 
  
  
  addscale(1/th,-(1.0-th)/th,Ex,Exth,nxn,nyn,nzn);
  addscale(1/th,-(1.0-th)/th,Ey,Eyth,nxn,nyn,nzn);
  addscale(1/th,-(1.0-th)/th,Ez,Ezth,nxn,nyn,nzn);
  
  // apply to smooth to electric field 3 times
  smoothE(Smooth,vct);
  smoothE(Smooth,vct);
  smoothE(Smooth,vct);
 
  
  // communicate so the interpolation can have good values
  communicateNodeBC(nxn,nyn,nzn,Exth,1,1,1,1,1,1,vct);
  communicateNodeBC(nxn,nyn,nzn,Eyth,1,1,2,2,1,1,vct);
  communicateNodeBC(nxn,nyn,nzn,Ezth,1,1,1,1,1,1,vct);
  communicateNodeBC(nxn,nyn,nzn,Ex,1,1,1,1,1,1,vct);
  communicateNodeBC(nxn,nyn,nzn,Ey,1,1,2,2,1,1,vct);
  communicateNodeBC(nxn,nyn,nzn,Ez,1,1,1,1,1,1,vct);
 
    
  // deallocate temporary arrays
  delete[] xkrylov;
  delete[] bkrylov;
  delete[] xkrylovPoisson;
  delete[] bkrylovPoisson;
  delArr3(divE,nxc,nyc);
  delArr3(gradPHIX,nxn,nyn);
  delArr3(gradPHIY,nxn,nyn);
  delArr3(gradPHIZ,nxn,nyn);   
 
}

/** Calculate sorgent for Maxwell solver */
inline void EMfields3D::MaxwellSource(double *bkrylov, Grid *grid, VirtualTopology3D *vct){
   eqValue (0.0, tempC,nxc,nyc,nzc);
   eqValue (0.0, tempX,nxn,nyn,nzn);
   eqValue (0.0, tempY,nxn,nyn,nzn);
   eqValue (0.0, tempZ,nxn,nyn,nzn);
   eqValue (0.0, tempXN,nxn,nyn,nzn);
   eqValue (0.0, tempYN,nxn,nyn,nzn);
   eqValue (0.0, tempZN,nxn,nyn,nzn);
   eqValue (0.0, temp2X,nxn,nyn,nzn);
   eqValue (0.0, temp2Y,nxn,nyn,nzn);
   eqValue (0.0, temp2Z,nxn,nyn,nzn);
   // communicate
   communicateCenterBC(nxc,nyc,nzc,Bxc,2,2,2,2,2,2,vct);
   communicateCenterBC(nxc,nyc,nzc,Byc,1,1,1,1,1,1,vct);
   communicateCenterBC(nxc,nyc,nzc,Bzc,2,2,2,2,2,2,vct);
   // fixBforcefree(grid,vct);
   fixBgem(grid,vct);
   //prepare curl of B for known term of Maxwell solver: for the source term
   grid->curlC2N(tempXN,tempYN,tempZN,Bxc,Byc,Bzc);
   scale(temp2X,Jxh,-FourPI/c,nxn,nyn,nzn);
   scale(temp2Y,Jyh,-FourPI/c,nxn,nyn,nzn);
   scale(temp2Z,Jzh,-FourPI/c,nxn,nyn,nzn);
   
  
   sum(temp2X,tempXN,nxn,nyn,nzn);
   sum(temp2Y,tempYN,nxn,nyn,nzn);
   sum(temp2Z,tempZN,nxn,nyn,nzn);
   scale(temp2X,delt,nxn,nyn,nzn);
   scale(temp2Y,delt,nxn,nyn,nzn);
   scale(temp2Z,delt,nxn,nyn,nzn);
   
   communicateCenterBC_P(nxc,nyc,nzc,rhoh,2,2,2,2,2,2,vct);
   grid->gradC2N(tempX,tempY,tempZ,rhoh);
   
   scale(tempX,-delt*delt*FourPI,nxn,nyn,nzn);
   scale(tempY,-delt*delt*FourPI,nxn,nyn,nzn);
   scale(tempZ,-delt*delt*FourPI,nxn,nyn,nzn);
   // sum E, past values
   sum(tempX,Ex,nxn,nyn,nzn);
   sum(tempY,Ey,nxn,nyn,nzn);
   sum(tempZ,Ez,nxn,nyn,nzn);
   // sum curl(B) + jhat part
   sum(tempX,temp2X,nxn,nyn,nzn);
   sum(tempY,temp2Y,nxn,nyn,nzn);
   sum(tempZ,temp2Z,nxn,nyn,nzn);

   // Boundary condition in the known term
   // boundary condition: Xleft
   if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==0) // perfect conductor
		perfectConductorLeftS(tempX,tempY,tempZ,0);
   // boundary condition: Xright
   if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright==0) // perfect conductor
           perfectConductorRightS(tempX,tempY,tempZ,0);
	// boundary condition: Yleft
	if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==0) // perfect conductor
           perfectConductorLeftS(tempX,tempY,tempZ,1);
     // boundary condition: Yright
     if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==0) // perfect conductor
           perfectConductorRightS(tempX,tempY,tempZ,1);
     // boundary condition: Zleft
     if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==0) // perfect conductor
           perfectConductorLeftS(tempX,tempY,tempZ,2);
     // boundary condition: Zright
     if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright==0 ) // perfect conductor
           perfectConductorRightS(tempX,tempY,tempZ,2);
    	   
     // physical space -> Krylov space
     phys2solver(bkrylov,tempX,tempY,tempZ,nxn,nyn,nzn);
   
}
/** Mapping of Maxwell image to give to solver */
inline void EMfields3D::MaxwellImage(double *im, double *vector, Grid *grid, VirtualTopology3D *vct){
   eqValue (0.0, im, 3*(nxn-2)*(nyn-2)*(nzn-2));
   eqValue (0.0, imageX,nxn,nyn,nzn); eqValue (0.0, imageY,nxn,nyn,nzn); eqValue (0.0, imageZ,nxn,nyn,nzn);
   eqValue (0.0, tempX,nxn,nyn,nzn);  eqValue (0.0, tempY,nxn,nyn,nzn);  eqValue (0.0, tempZ,nxn,nyn,nzn);
   eqValue (0.0, Dx,nxn,nyn,nzn);  eqValue (0.0,Dy,nxn,nyn,nzn);  eqValue (0.0, Dz,nxn,nyn,nzn);
   // move from krylov space to physical space
   solver2phys(vectX,vectY,vectZ,vector,nxn,nyn,nzn);
   grid->lapN2N(imageX,vectX,vct);
   grid->lapN2N(imageY,vectY,vct);
   grid->lapN2N(imageZ,vectZ,vct);
   neg(imageX,nxn,nyn,nzn);
   neg(imageY,nxn,nyn,nzn);
   neg(imageZ,nxn,nyn,nzn);
   // grad(div(mu dot E(n + theta))       mu dot E(n + theta) = D
   MUdot(Dx, Dy, Dz, vectX, vectY, vectZ, grid);
   grid->divN2C(divC,Dx,Dy,Dz);
   // communicate you should put BC 
   // think about the Physics 
   //communicateCenterBC(nxc,nyc,nzc,divC,1,1,1,1,1,1,vct);
   communicateCenterBC(nxc,nyc,nzc,divC,2,2,2,2,2,2,vct); // GO with Neumann, now then go with rho
   
	
   grid->gradC2N(tempX,tempY,tempZ,divC);
   
   // -lap(E(n +theta)) -  grad(div(mu dot E(n + theta))
   sub(imageX,tempX,nxn,nyn,nzn); sub(imageY,tempY,nxn,nyn,nzn); sub(imageZ,tempZ,nxn,nyn,nzn);
   
   //scale delt*delt
   scale(imageX,delt*delt,nxn,nyn,nzn); scale(imageY,delt*delt,nxn,nyn,nzn); scale(imageZ,delt*delt,nxn,nyn,nzn);
   
   //  -lap(E(n +theta)) -  grad(div(mu dot E(n + theta)) + eps dot E(n + theta)
   sum(imageX,Dx,nxn,nyn,nzn); sum(imageY,Dy,nxn,nyn,nzn); sum(imageZ,Dz,nxn,nyn,nzn);
   sum(imageX,vectX,nxn,nyn,nzn); sum(imageY,vectY,nxn,nyn,nzn); sum(imageZ,vectZ,nxn,nyn,nzn);
   
   // boundary condition: Xleft
   if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==0) // perfect conductor
           perfectConductorLeft(imageX,imageY,imageZ,vectX,vectY,vectZ,0,grid);
   // boundary condition: Xright
   if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright==0) // perfect conductor
           perfectConductorRight(imageX,imageY,imageZ,vectX,vectY,vectZ,0,grid);
   // boundary condition: Yleft
   if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==0) // perfect conductor
           perfectConductorLeft(imageX,imageY,imageZ,vectX,vectY,vectZ,1,grid);
   // boundary condition: Yright
   if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==0) // perfect conductor
           perfectConductorRight(imageX,imageY,imageZ,vectX,vectY,vectZ,1,grid);
   // boundary condition: Zleft
   if(vct->getZleft_neighbor()==MPI_PROC_NULL && bcEMfaceZleft ==0) // perfect conductor
           perfectConductorLeft(imageX,imageY,imageZ,vectX,vectY,vectZ,2,grid);
   // boundary condition: Zright
   if(vct->getZright_neighbor()==MPI_PROC_NULL && bcEMfaceZright==0) // perfect conductor
           perfectConductorRight(imageX,imageY,imageZ,vectX,vectY,vectZ,2,grid);
   // move from physical space to krylov space
   phys2solver(im,imageX,imageY,imageZ,nxn,nyn,nzn);
     
   
}

/** Calculate PI dot (vectX, vectY, vectZ)*/
inline void EMfields3D::PIdot(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid *grid){
     double beta, edotb, omcx, omcy, omcz,denom;
     beta = .5*qom[ns]*dt/c;
     for(int i=1; i < nxn-1;i++)
         for(int j=1; j < nyn-1;j++)
          for(int k=1; k < nzn-1;k++){
                     omcx = beta*Bxn[i][j][k];
                     omcy = beta*Byn[i][j][k];
                     omcz = beta*Bzn[i][j][k];
                     edotb =  vectX[i][j][k]*omcx +  vectY[i][j][k]*omcy + vectZ[i][j][k]*omcz;
                     denom  = 1/(1.0 + omcx*omcx + omcy*omcy + omcz*omcz);
                     PIdotX[i][j][k] += (vectX[i][j][k] +(vectY[i][j][k]*omcz - vectZ[i][j][k]*omcy + edotb*omcx))*denom;
                     PIdotY[i][j][k] += (vectY[i][j][k] +(vectZ[i][j][k]*omcx - vectX[i][j][k]*omcz + edotb*omcy))*denom;
                     PIdotZ[i][j][k] += (vectZ[i][j][k] +(vectX[i][j][k]*omcy - vectY[i][j][k]*omcx + edotb*omcz))*denom;
                 }

     
}
/** Calculate MU dot (vectX, vectY, vectZ)*/
inline void EMfields3D::MUdot(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid *grid){
      double beta, edotb, omcx, omcy, omcz, denom;
      for(int i=1; i < nxn-1;i++)
        for(int j=1; j < nyn-1;j++)
          for(int k=1; k < nzn-1;k++){
               MUdotX[i][j][k] = 0.0;
               MUdotY[i][j][k] = 0.0;
               MUdotZ[i][j][k] = 0.0;
          }
      for (int is=0; is < ns; is++){
          beta = .5*qom[is]*dt/c;
          for(int i=1; i < nxn-1;i++)
             for(int j=1; j < nyn-1;j++)
               for(int k=1; k < nzn-1;k++){
                     omcx = beta*Bxn[i][j][k];
                     omcy = beta*Byn[i][j][k];
                     omcz = beta*Bzn[i][j][k];
                     edotb =  vectX[i][j][k]*omcx +  vectY[i][j][k]*omcy + vectZ[i][j][k]*omcz;
			         denom = FourPI/2*delt*dt/c*qom[is]*rhons[is][i][j][k]/(1.0 + omcx*omcx + omcy*omcy + omcz*omcz);
                     MUdotX[i][j][k] += (vectX[i][j][k] +(vectY[i][j][k]*omcz - vectZ[i][j][k]*omcy + edotb*omcx))*denom;
                     MUdotY[i][j][k] += (vectY[i][j][k] +(vectZ[i][j][k]*omcx - vectX[i][j][k]*omcz + edotb*omcy))*denom;
                     MUdotZ[i][j][k] += (vectZ[i][j][k] +(vectX[i][j][k]*omcy - vectY[i][j][k]*omcx + edotb*omcz))*denom;
                 }
                 
          }
     
}
/*Interpolation smoothing:
Smoothing (vector must already have ghost cells)
TO MAKE SMOOTH value as to be different from 1.0
type = 0 --> center based vector ; 
type = 1 --> node based vector   ;                     
*/
inline void  EMfields3D::smooth(double value,double ***vector, bool type, Grid *grid, VirtualTopology3D *vct){

 int nvolte=6;
 for (int icount=1; icount<nvolte+1; icount++)
        {

	if(value!=1.0){
		double alpha;
		int nx, ny, nz;
		switch(type){
			case (0):
				nx =grid->getNXC();
				ny =grid->getNYC();
				nz =grid->getNZC();
				communicateCenterBoxStencilBC_P(nx,ny,nz,vector,2,2,2,2,2,2,vct);
				
				break;
			case(1):
				nx = grid->getNXN();
				ny = grid->getNYN();
				nz = grid->getNZN();
				communicateNodeBoxStencilBC_P(nx,ny,nz,vector,2,2,2,2,2,2,vct);
				break;
		}
		double ***temp = newArr3(double,nx,ny,nz);
        if (icount%2==1){
        	value = 0.;
        } else {
        	value=0.5;
        }
		alpha=(1.0-value)/6;
		for (int i=1; i<nx-1; i++)
			for (int j=1; j<ny-1; j++)
				for (int k=1; k<nz-1; k++)
					temp[i][j][k] = value*vector[i][j][k]+alpha*(vector[i-1][j][k] + vector[i+1][j][k] + vector[i][j-1][k]  + vector[i][j+1][k] + vector[i][j][k-1] + vector[i][j][k+1]);
		for (int i=1; i<nx-1; i++)
			for (int j=1; j<ny-1; j++)
				for (int k=1; k<nz-1; k++)
					vector[i][j][k] = temp[i][j][k];
		delArr3(temp,nx,ny);
	}
}
}
/*Interpolation smoothing:
Smoothing (vector must already have ghost cells)
TO MAKE SMOOTH value as to be different from 1.0
type = 0 --> center based vector ; 
type = 1 --> node based vector   ;                     
*/
inline void  EMfields3D::smoothE(double value, VirtualTopology3D *vct){

 int nvolte=6;
 for (int icount=1; icount<nvolte+1; icount++)
        {
	if(value!=1.0){
		double alpha;
		communicateNodeBoxStencilBC(nxn,nyn,nzn,Ex,1,1,1,1,1,1,vct);
		communicateNodeBoxStencilBC(nxn,nyn,nzn,Ey,1,1,2,2,1,1,vct);
		communicateNodeBoxStencilBC(nxn,nyn,nzn,Ez,1,1,1,1,1,1,vct);
	
		double ***temp = newArr3(double,nxn,nyn,nzn);
        if (icount%2==1){
        	value = 0.;
        } else {
        	value=0.5;
        }
		alpha=(1.0-value)/6;
		// Exth
		for (int i=1; i<nxn-1; i++)
			for (int j=1; j<nyn-1; j++)
				for (int k=1; k<nzn-1; k++)
					temp[i][j][k] = value*Ex[i][j][k]+alpha*(Ex[i-1][j][k] + Ex[i+1][j][k] + Ex[i][j-1][k]  + Ex[i][j+1][k] + Ex[i][j][k-1] + Ex[i][j][k+1]);
		for (int i=1; i<nxn-1; i++)
			for (int j=1; j<nyn-1; j++)
				for (int k=1; k<nzn-1; k++)
					Ex[i][j][k] = temp[i][j][k];
		// Eyth
		for (int i=1; i<nxn-1; i++)
			for (int j=1; j<nyn-1; j++)
				for (int k=1; k<nzn-1; k++)
					temp[i][j][k] = value*Ey[i][j][k]+alpha*(Ey[i-1][j][k] + Ey[i+1][j][k] + Ey[i][j-1][k]  + Ey[i][j+1][k] + Ey[i][j][k-1] + Ey[i][j][k+1]);
		for (int i=1; i<nxn-1; i++)
			for (int j=1; j<nyn-1; j++)
				for (int k=1; k<nzn-1; k++)
					Ey[i][j][k] = temp[i][j][k];
		// Ezth
		for (int i=1; i<nxn-1; i++)
			for (int j=1; j<nyn-1; j++)
				for (int k=1; k<nzn-1; k++)
					temp[i][j][k] = value*Ez[i][j][k]+alpha*(Ez[i-1][j][k] + Ez[i+1][j][k] + Ez[i][j-1][k]  + Ez[i][j+1][k] + Ez[i][j][k-1] + Ez[i][j][k+1]);
		for (int i=1; i<nxn-1; i++)
			for (int j=1; j<nyn-1; j++)
				for (int k=1; k<nzn-1; k++)
					Ez[i][j][k] = temp[i][j][k];
					
					
		delArr3(temp,nxn,nyn);
	}
}
}

/* SPECIES: Interpolation smoothing
TO MAKE SMOOTH value as to be different from 1.0
type = 0 --> center based vector  
type = 1 --> node based vector                        
*/
inline void  EMfields3D::smooth(double value,double ****vector,int is, bool type, Grid *grid, VirtualTopology3D *vct){
  cout << "Smoothing for Species not implemented in 3D" << endl; 
}

/** fix the B boundary when running gem*/
inline void EMfields3D::fixBgem(Grid *grid, VirtualTopology3D *vct){
   if (vct->getYright_neighbor()==MPI_PROC_NULL){
      for (int i=0; i < nxc;i++)
	    for (int k=0; k < nzc;k++){
			Bxc[i][nyc-1][k] = B0x*tanh((grid->getYC(i,nyc-1,k) - Ly/2)/delta);
			Bxc[i][nyc-2][k] = Bxc[i][nyc-1][k];
			Bxc[i][nyc-3][k] = Bxc[i][nyc-1][k];
			Byc[i][nyc-1][k] = B0y;
			Bzc[i][nyc-1][k] = B0z;
			Bzc[i][nyc-2][k] = B0z;
			Bzc[i][nyc-3][k] = B0z;
		}
	}
	if (vct->getYleft_neighbor()==MPI_PROC_NULL){
      for (int i=0; i < nxc;i++)
	    for (int k=0; k < nzc;k++){
			Bxc[i][0][k] = B0x*tanh((grid->getYC(i,0,k) - Ly/2)/delta);
			Bxc[i][1][k] = Bxc[i][0][k];
			Bxc[i][2][k] = Bxc[i][0][k];
			Byc[i][0][k] = B0y;
			Bzc[i][0][k] = B0z;
			Bzc[i][1][k] = B0z;
			Bzc[i][2][k] = B0z;
		}
	}

}

/** fix the B boundary when running forcefree*/
inline void EMfields3D::fixBforcefree(Grid *grid, VirtualTopology3D *vct){
   if (vct->getYright_neighbor()==MPI_PROC_NULL){
      for (int i=0; i < nxc;i++)
	    for (int k=0; k < nzc;k++){
			Bxc[i][nyc-1][k] = B0x*tanh((grid->getYC(i,nyc-1,k) - Ly/2)/delta);
			Byc[i][nyc-1][k] = B0y;
			Bzc[i][nyc-1][k] = B0z/cosh((grid->getYC(i,nyc-1,k) - Ly/2)/delta);;
			Bzc[i][nyc-2][k] = B0z/cosh((grid->getYC(i,nyc-2,k) - Ly/2)/delta);;
			Bzc[i][nyc-3][k] = B0z/cosh((grid->getYC(i,nyc-3,k) - Ly/2)/delta);
		}
	}
	if (vct->getYleft_neighbor()==MPI_PROC_NULL){
      for (int i=0; i < nxc;i++)
	    for (int k=0; k < nzc;k++){
			Bxc[i][0][k] = B0x*tanh((grid->getYC(i,0,k) - Ly/2)/delta);
			Byc[i][0][k] = B0y;
			Bzc[i][0][k] = B0z/cosh((grid->getYC(i,0,k) - Ly/2)/delta);
			Bzc[i][1][k] = B0z/cosh((grid->getYC(i,1,k) - Ly/2)/delta);
			Bzc[i][2][k] = B0z/cosh((grid->getYC(i,2,k) - Ly/2)/delta);
		}
	}

}


/** adjust densities on boundaries that are not periodic */
inline void EMfields3D::adjustNonPeriodicDensities(int is, VirtualTopology3D *vct){
 if (vct->getXleft_neighbor_P()==MPI_PROC_NULL){
	   for (int i=1; i < nyn-1;i++)
	     for (int k=1; k < nzn-1;k++){
	     rhons[is][1][i][k]+=   rhons[is][1][i][k];
		 Jxs[is][1][i][k]  +=   Jxs[is][1][i][k];
		 Jys[is][1][i][k]  +=   Jys[is][1][i][k];
		 Jzs[is][1][i][k]  +=   Jzs[is][1][i][k];
		 pXXsn[is][1][i][k]  += pXXsn[is][1][i][k];
		 pXYsn[is][1][i][k]  += pXYsn[is][1][i][k];
		 pXZsn[is][1][i][k]  += pXZsn[is][1][i][k];
		 pYYsn[is][1][i][k]  += pYYsn[is][1][i][k];
		 pYZsn[is][1][i][k]  += pYZsn[is][1][i][k];
		 pZZsn[is][1][i][k]  += pZZsn[is][1][i][k];
	   }
  }
  if (vct->getYleft_neighbor_P()==MPI_PROC_NULL){
	   for (int i=1; i < nxn-1;i++)
	    for (int k=1; k < nzn-1;k++){
	     rhons[is][i][1][k]+=   rhons[is][i][1][k];
		 Jxs[is][i][1][k]  +=   Jxs[is][i][1][k];
		 Jys[is][i][1][k]  +=   Jys[is][i][1][k];
		 Jzs[is][i][1][k]  +=   Jzs[is][i][1][k];
		 pXXsn[is][i][1][k]  += pXXsn[is][i][1][k];
		 pXYsn[is][i][1][k]  += pXYsn[is][i][1][k];
		 pXZsn[is][i][1][k]  += pXZsn[is][i][1][k];
		 pYYsn[is][i][1][k]  += pYYsn[is][i][1][k];
		 pYZsn[is][i][1][k]  += pYZsn[is][i][1][k];
		 pZZsn[is][i][1][k]  += pZZsn[is][i][1][k];
	   }
  }
  if (vct->getZleft_neighbor_P()==MPI_PROC_NULL){
	   for (int i=1; i < nxn-1;i++)
	     for (int j=1; j < nyn-1;j++){
	     rhons[is][i][j][1]+=   rhons[is][i][j][1];
		 Jxs[is][i][j][1]  +=   Jxs[is][i][j][1];
		 Jys[is][i][j][1]  +=   Jys[is][i][j][1];
		 Jzs[is][i][j][1]  +=   Jzs[is][i][j][1];
		 pXXsn[is][i][j][1]  += pXXsn[is][i][j][1];
		 pXYsn[is][i][j][1]  += pXYsn[is][i][j][1];
		 pXZsn[is][i][j][1]  += pXZsn[is][i][j][1];
		 pYYsn[is][i][j][1]  += pYYsn[is][i][j][1];
		 pYZsn[is][i][j][1]  += pYZsn[is][i][j][1];
		 pZZsn[is][i][j][1]  += pZZsn[is][i][j][1];
	   }
  }
  if (vct->getXright_neighbor_P()==MPI_PROC_NULL){
	   for (int i=1; i < nyn-1;i++)
	     for (int k=1; k < nzn-1;k++){
	     rhons[is][nxn-2][i][k]+=   rhons[is][1][nxn-2][k];
		 Jxs[is][nxn-2][i][k]  +=   Jxs[is][1][nxn-2][k];
		 Jys[is][nxn-2][i][k]  +=   Jys[is][1][nxn-2][k];
		 Jzs[is][nxn-2][i][k]  +=   Jzs[is][1][nxn-2][k];
		 pXXsn[is][nxn-2][i][k]  += pXXsn[is][nxn-2][i][k];
		 pXYsn[is][nxn-2][i][k]  += pXYsn[is][nxn-2][i][k];
		 pXZsn[is][nxn-2][i][k]  += pXZsn[is][nxn-2][i][k];
		 pYYsn[is][nxn-2][i][k]  += pYYsn[is][nxn-2][i][k];
		 pYZsn[is][nxn-2][i][k]  += pYZsn[is][nxn-2][i][k];
		 pZZsn[is][nxn-2][i][k]  += pZZsn[is][nxn-2][i][k];
	   }
  }
  if (vct->getYright_neighbor_P()==MPI_PROC_NULL){
	   for (int i=1; i < nxn-1;i++)
	    for (int k=1; k < nzn-1;k++){
	     rhons[is][i][nyn-2][k]+=   rhons[is][i][nyn-2][k];
		 Jxs[is][i][nyn-2][k]  +=   Jxs[is][i][nyn-2][k];
		 Jys[is][i][nyn-2][k]  +=   Jys[is][i][nyn-2][k];
		 Jzs[is][i][nyn-2][k]  +=   Jzs[is][i][nyn-2][k];
		 pXXsn[is][i][nyn-2][k]  += pXXsn[is][i][nyn-2][k];
		 pXYsn[is][i][nyn-2][k]  += pXYsn[is][i][nyn-2][k];
		 pXZsn[is][i][nyn-2][k]  += pXZsn[is][i][nyn-2][k];
		 pYYsn[is][i][nyn-2][k]  += pYYsn[is][i][nyn-2][k];
		 pYZsn[is][i][nyn-2][k]  += pYZsn[is][i][nyn-2][k];
		 pZZsn[is][i][nyn-2][k]  += pZZsn[is][i][nyn-2][k];
	   }
  }
  if (vct->getZright_neighbor_P()==MPI_PROC_NULL){
	   for (int i=1; i < nxn-1;i++)
	     for (int j=1; j < nyn-1;j++){
	     rhons[is][i][j][nzn-2]+=   rhons[is][i][j][nzn-2];
		 Jxs[is][i][j][nzn-2]  +=   Jxs[is][i][j][nzn-2];
		 Jys[is][i][j][nzn-2]  +=   Jys[is][i][j][nzn-2];
		 Jzs[is][i][j][nzn-2]  +=   Jzs[is][i][j][nzn-2];
		 pXXsn[is][i][j][nzn-2]  += pXXsn[is][i][j][nzn-2];
		 pXYsn[is][i][j][nzn-2]  += pXYsn[is][i][j][nzn-2];
		 pXZsn[is][i][j][nzn-2]  += pXZsn[is][i][j][nzn-2];
		 pYYsn[is][i][j][nzn-2]  += pYYsn[is][i][j][nzn-2];
		 pYZsn[is][i][j][nzn-2]  += pYZsn[is][i][j][nzn-2];
		 pZZsn[is][i][j][nzn-2]  += pZZsn[is][i][j][nzn-2];
	   }
  }
}
/** Calculate Magnetic field with the implicit solver: calculate B defined on nodes
With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
inline void EMfields3D::calculateB(Grid *grid, VirtualTopology3D *vct){
   if (vct->getCartesian_rank() ==0)
    cout << "*** B CALCULATION ***" << endl;
   // calculate the curl of Eth
   grid->curlN2C(tempXC,tempYC,tempZC,Exth,Eyth,Ezth);
   // update the magnetic field
   addscale(-c*dt,1,Bxc,tempXC,nxc,nyc,nzc);
   addscale(-c*dt,1,Byc,tempYC,nxc,nyc,nzc);
   addscale(-c*dt,1,Bzc,tempZC,nxc,nyc,nzc);
    // communicate ghost 
   communicateCenterBC(nxc,nyc,nzc,Bxc,2,2,2,2,2,2,vct);
   communicateCenterBC(nxc,nyc,nzc,Byc,1,1,1,1,1,1,vct);
   communicateCenterBC(nxc,nyc,nzc,Bzc,2,2,2,2,2,2,vct);
   //fixBforcefree(grid,vct);
   fixBgem(grid,vct);
   // interpolate C2N
   grid->interpC2N(Bxn,Bxc);
   grid->interpC2N(Byn,Byc);
   grid->interpC2N(Bzn,Bzc);
    
  
   communicateNodeBC(nxn,nyn,nzn,Bxn,1,1,2,2,1,1,vct);
   communicateNodeBC(nxn,nyn,nzn,Byn,1,1,1,1,1,1,vct);
   communicateNodeBC(nxn,nyn,nzn,Bzn,1,1,2,2,1,1,vct);
  

}
/** initialize EM field with transverse electric waves 1D and rotate
anticlockwise (theta degrees)*/
inline void EMfields3D::initEM_rotate(VirtualTopology3D *vct, Grid *grid, double B, double theta){
   // initialize E and rhos on nodes
   for (int i=0; i < nxn; i++)
       for (int j=0; j <nyn; j++){
          Ex[i][j][0] =  0.0;
          Ey[i][j][0] =  0.0;
          Ez[i][j][0] =  0.0;
          Bxn[i][j][0] =  B*cos(theta*M_PI/180);
          Byn[i][j][0] =  B*sin(theta*M_PI/180);
          Bzn[i][j][0] = 0.0;
          rhons[0][i][j][0] =  0.07957747154595; // electrons: species is now first index
          rhons[1][i][j][0] =  0.07957747154595; // protons: species is now first index
       }
       // initialize B on centers
       grid->interpN2C(Bxc,Bxn);
       grid->interpN2C(Byc,Byn);
       grid->interpN2C(Bzc,Bzn);


       for (int is=0 ; is<ns; is++)
           grid->interpN2C(rhocs,is,rhons);

}
/**Add a periodic perturbation in rho exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
inline void EMfields3D::AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid *grid){

double alpha;
alpha=deltaBoB*B0/sqrt(Bx_mod*Bx_mod+By_mod*By_mod+Bz_mod*Bz_mod);

ne_mod *= alpha;
ni_mod *= alpha;
//cout<<" ne="<<ne_mod<<" ni="<<ni_mod<<" alpha="<<alpha<<endl;
for (int i=0; i < nxn; i++)
for (int j=0; j <nyn; j++){
	rhons[0][i][j][0] +=  ne_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + ne_phase);
	rhons[1][i][j][0] +=  ni_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + ni_phase);
}

for (int is=0 ; is<ns; is++)
grid->interpN2C(rhocs,is,rhons);
}


/**Add a periodic perturbation exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
inline void EMfields3D::AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid *grid){

double alpha;

alpha=deltaBoB*B0/sqrt(Bx_mod*Bx_mod+By_mod*By_mod+Bz_mod*Bz_mod);

Ex_mod *= alpha;
Ey_mod *= alpha;
Ez_mod *= alpha;
Bx_mod *= alpha;
By_mod *= alpha;
Bz_mod *= alpha;

for (int i=0; i < nxn; i++)
for (int j=0; j <nyn; j++){
	Ex[i][j][0] +=  Ex_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + Ex_phase);
	Ey[i][j][0] +=  Ey_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + Ey_phase);
	Ez[i][j][0] +=  Ez_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + Ez_phase);
	Bxn[i][j][0] +=  Bx_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + Bx_phase);
	Byn[i][j][0] +=  By_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + By_phase);
	Bzn[i][j][0] +=  Bz_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + Bz_phase);

}

// initialize B on centers
grid->interpN2C(Bxc,Bxn);
grid->interpN2C(Byc,Byn);
grid->interpN2C(Bzc,Bzn);


}


/** Calculate hat rho hat, Jx hat, Jy hat, Jz hat */
inline void EMfields3D::calculateHatFunctions(Grid *grid, VirtualTopology3D *vct){
   // smoothing
   smooth(Smooth,rhoc,0,grid,vct);
   // calculate j hat
  
  for (int is=0; is < ns; is++){
    grid->divSymmTensorN2C(tempXC,tempYC,tempZC,pXXsn,pXYsn,pXZsn,pYYsn,pYZsn,pZZsn,is);
    
	scale(tempXC,-dt/2.0,nxc,nyc,nzc);
    scale(tempYC,-dt/2.0,nxc,nyc,nzc);
    scale(tempZC,-dt/2.0,nxc,nyc,nzc);
    // communicate before interpolating
	communicateCenterBC_P(nxc,nyc,nzc,tempXC,2,2,2,2,2,2,vct);
	communicateCenterBC_P(nxc,nyc,nzc,tempYC,2,2,2,2,2,2,vct);
	communicateCenterBC_P(nxc,nyc,nzc,tempZC,2,2,2,2,2,2,vct);
	
    grid->interpC2N(tempXN,tempXC);
    grid->interpC2N(tempYN,tempYC);
    grid->interpC2N(tempZN,tempZC);
    sum(tempXN,Jxs,nxn,nyn,nzn,is);
    sum(tempYN,Jys,nxn,nyn,nzn,is);
    sum(tempZN,Jzs,nxn,nyn,nzn,is);
	//PIDOT
    PIdot(Jxh,Jyh,Jzh,tempXN,tempYN,tempZN,is,grid);
    
   }
   // smooth j
   smooth(Smooth,Jxh,1,grid,vct);
   smooth(Smooth,Jyh,1,grid,vct);
   smooth(Smooth,Jzh,1,grid,vct);
   
   // calculate rho hat = rho - (dt*theta)div(jhat)
   grid->divN2C(tempXC,Jxh,Jyh,Jzh);
   scale(tempXC,-dt*th,nxc,nyc,nzc);
   sum(tempXC,rhoc,nxc,nyc,nzc);
   eq(rhoh,tempXC,nxc,nyc,nzc);
   // communicate rhoh
   communicateCenterBC_P(nxc,nyc,nzc,rhoh,2,2,2,2,2,2,vct);
}
/** Image of Poisson Solver */
inline void EMfields3D::PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology3D *vct){
    // allocate  2 three dimensional service vectors
    double ***temp = newArr3(double,nxc,nyc,nzc);
    double ***im   = newArr3(double,nxc,nyc,nzc);
	eqValue (0.0, image,(nxc-2)*(nyc-2)*(nzc-2)); 
	eqValue (0.0, temp,nxc,nyc,nzc);
    eqValue (0.0, im,nxc,nyc,nzc);
    // move from krylov space to physical space and communicate ghost cells
	solver2phys(temp,vector,nxc,nyc,nzc);
    // calculate the laplacian
    grid->lapC2Cpoisson(im,temp,vct);
    // move from physical space to krylov space
    phys2solver(image,im,nxc,nyc,nzc);
    // deallocate temporary array and objects
    delArr3(temp,nxc,nyc);
    delArr3(im,nxc,nyc);
}
/** interpolate charge density and pressure density from node to center */
inline  void EMfields3D::interpDensitiesN2C(VirtualTopology3D *vct, Grid *grid){
   // do we need communication or not really?
   grid->interpN2C(rhoc,rhon);
}
/** communicate ghost for grid -> Particles interpolation */
inline void EMfields3D::communicateGhostP2G(int ns,int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology3D *vct){
	// interpolate adding common nodes among processors
	communicateInterp(nxn,nyn,nzn,ns,rhons,0,0,0,0,0,0,vct);
	communicateInterp(nxn,nyn,nzn,ns,Jxs,0,0,0,0,0,0,vct);
	communicateInterp(nxn,nyn,nzn,ns,Jys,0,0,0,0,0,0,vct);
	communicateInterp(nxn,nyn,nzn,ns,Jzs,0,0,0,0,0,0,vct);
	communicateInterp(nxn,nyn,nzn,ns,pXXsn,0,0,0,0,0,0,vct);
	communicateInterp(nxn,nyn,nzn,ns,pXYsn,0,0,0,0,0,0,vct);   
	communicateInterp(nxn,nyn,nzn,ns,pXZsn,0,0,0,0,0,0,vct);
	communicateInterp(nxn,nyn,nzn,ns,pYYsn,0,0,0,0,0,0,vct);
	communicateInterp(nxn,nyn,nzn,ns,pYZsn,0,0,0,0,0,0,vct);
	communicateInterp(nxn,nyn,nzn,ns,pZZsn,0,0,0,0,0,0,vct);
	// calculate the correct densities on the boundaries
	adjustNonPeriodicDensities(ns, vct);
	// put the correct values on ghost cells
	
	communicateNode_P(nxn,nyn,nzn,rhons,ns,vct);
	communicateNode_P(nxn,nyn,nzn,Jxs,ns,vct);
    communicateNode_P(nxn,nyn,nzn,Jys,ns,vct);
	communicateNode_P(nxn,nyn,nzn,Jzs,ns,vct);
	communicateNode_P(nxn,nyn,nzn,pXXsn,ns,vct);
	communicateNode_P(nxn,nyn,nzn,pXYsn,ns,vct);   
	communicateNode_P(nxn,nyn,nzn,pXZsn,ns,vct);
	communicateNode_P(nxn,nyn,nzn,pYYsn,ns,vct);
	communicateNode_P(nxn,nyn,nzn,pYZsn,ns,vct);
	communicateNode_P(nxn,nyn,nzn,pZZsn,ns,vct);
	 
}   


/** add an amount of charge density to charge density field at node X,Y */
inline void EMfields3D::addRho(double weight[][2][2], int X, int Y, int Z, int is){
    for (int i=0; i < 2; i++)
      for (int j=0; j < 2; j++)
        for(int k=0; k < 2; k++)
          rhons[is][X -i][Y -j][Z - k] += weight[i][j][k]*invVOL;
}
/** add an amount of charge density to current density - direction X to current density field on the node*/
inline void EMfields3D::addJx(double weight[][2][2], int X, int Y, int Z, int is){
	 for (int i=0; i < 2; i++)
      for (int j=0; j < 2; j++)
        for(int k=0; k < 2; k++)
            Jxs[is][X -i][Y -j][Z - k]   += weight[i][j][k]*invVOL;
}
/** add an amount of current density - direction Y to current density field on the node */
inline  void EMfields3D::addJy(double weight[][2][2], int X, int Y, int Z, int is){
	 for (int i=0; i < 2; i++)
	  for (int j=0; j < 2; j++)
        for(int k=0; k < 2; k++)
            Jys[is][X -i][Y -j][Z - k]  += weight[i][j][k]*invVOL;
}
/** add an amount of current density - direction Z to current density field on the node */
inline  void EMfields3D::addJz(double weight[][2][2], int X, int Y, int Z, int is){
	 for (int i=0; i < 2; i++)
      for (int j=0; j < 2; j++)
        for(int k=0; k < 2; k++)
			Jzs[is][X -i][Y -j][Z - k]  += weight[i][j][k]*invVOL;
}
/** add an amount of pressure density - direction XX to current density field on the node */
inline  void EMfields3D::addPxx(double weight[][2][2], int X, int Y, int Z,int is){
	for (int i=0; i < 2; i++)
      for (int j=0; j < 2; j++)
        for(int k=0; k < 2; k++)
          pXXsn[is][X -i][Y -j][Z - k]       += weight[i][j][k]*invVOL;
}
/** add an amount of pressure density - direction XY to current density field on the node*/
inline  void EMfields3D::addPxy(double weight[][2][2], int X, int Y, int Z,int is){
	for (int i=0; i < 2; i++)
        for (int j=0; j < 2; j++)
          for(int k=0; k < 2; k++)
            pXYsn[is][X -i][Y -j][Z - k]        += weight[i][j][k]*invVOL;
}
/** add an amount of pressure density - direction XZ to current density field on the node */
inline  void EMfields3D::addPxz(double weight[][2][2], int X, int Y,int Z, int is){
	 for (int i=0; i < 2; i++)
      for (int j=0; j < 2; j++)
        for(int k=0; k < 2; k++)
          pXZsn[is][X -i][Y -j][Z - k]        += weight[i][j][k]*invVOL;
}
/** add an amount of pressure density - direction YY to current density field on the node*/
inline  void EMfields3D::addPyy(double weight[][2][2], int X, int Y, int Z, int is){
	for (int i=0; i < 2; i++)
      for (int j=0; j < 2; j++)
        for(int k=0; k < 2; k++)
          pYYsn[is][X -i][Y -j][Z - k]       += weight[i][j][k]*invVOL;
}
/** add an amount of pressure density - direction YZ to current density field on the node */
inline  void EMfields3D::addPyz(double weight[][2][2], int X, int Y, int Z, int is){
	for (int i=0; i < 2; i++)
      for (int j=0; j < 2; j++)
        for(int k=0; k < 2; k++)
            pYZsn[is][X -i][Y -j][Z - k]       += weight[i][j][k]*invVOL;
}
/** add an amount of pressure density - direction ZZ to current density field on the node */
inline  void EMfields3D::addPzz(double weight[][2][2], int X, int Y, int Z, int is){
	for (int i=0; i < 2; i++)
      for (int j=0; j < 2; j++)
        for(int k=0; k < 2; k++)
            pZZsn[is][X -i][Y -j][Z - k]        += weight[i][j][k]*invVOL;
}



/** set to 0 all the densities fields */
inline  void EMfields3D::setZeroDensities(){
  for (register int i=0; i < nxn; i++)
	for (register int j=0; j < nyn; j++)
	          for (register int k=0; k < nzn; k++){
	    Jx[i][j][k]   = 0.0;
	    Jxh[i][j][k]  = 0.0;
	    Jy[i][j][k]   = 0.0;
	    Jyh[i][j][k]  = 0.0;
	    Jz[i][j][k]   = 0.0;
	    Jzh[i][j][k]  = 0.0;
		rhon[i][j][k] = 0.0;
	  }
for (register int i=0; i < nxc; i++)
	for (register int j=0; j < nyc; j++)
	 for (register int k=0; k < nzc; k++){
	   rhoc[i][j][k] = 0.0;
	   rhoh[i][j][k] = 0.0;
	}
for (register int kk=0; kk < ns; kk++)
 for (register int i=0; i < nxn; i++)
	for (register int j=0; j < nyn; j++)
	     for (register int k=0; k < nzn; k++){
		  rhons[kk][i][j][k] = 0.0;
		  Jxs[kk][i][j][k]   = 0.0;
		  Jys[kk][i][j][k]   = 0.0;
		  Jzs[kk][i][j][k]   = 0.0;
		  pXXsn[kk][i][j][k] = 0.0;
		  pXYsn[kk][i][j][k] = 0.0;
		  pXZsn[kk][i][j][k] = 0.0;
		  pYYsn[kk][i][j][k] = 0.0;
		  pYZsn[kk][i][j][k] = 0.0;
		  pZZsn[kk][i][j][k] = 0.0;
	  }  
       
}
/**SPECIES: Sum the charge density of different species on NODES*/
inline void EMfields3D::sumOverSpecies(VirtualTopology3D *vct){
	for (int is=0; is<ns; is++)
 	   for (register int i=0; i <nxn; i++)
  	      for (register int j=0; j <nyn; j++)
		    for (register int k=0; k <nzn; k++)
	           rhon[i][j][k]     +=rhons[is][i][j][k];
}

/**SPECIES: Sum current density for different species */
inline void EMfields3D::sumOverSpeciesJ(){
	for (int is=0; is < ns; is++)
		for (register int i=0; i <nxn; i++)
			for (register int j=0; j <nyn; j++)
			  for (register int k=0; j <nzn; k++) {
				Jx[i][j][k]     +=Jxs[is][i][j][k];
				Jy[i][j][k]     +=Jys[is][i][j][k];
				Jz[i][j][k]     +=Jzs[is][i][j][k];
			  }
}



/** initialize Magnetic and Electric Field with initial configuration */
inline void EMfields3D::init(VirtualTopology3D *vct, Grid *grid){
 if (restart1 ==0){
   for (int i=0; i < nxn; i++){
    for (int j=0; j < nyn; j++){
      for (int k=0; k < nzn; k++){
        for (int is=0; is < ns; is++){
		  rhons[is][i][j][k] = rhoINIT[is]/FourPI;
		}
		Ex[i][j][k] = 0.0;
        Ey[i][j][k] = 0.0;
        Ez[i][j][k] = 0.0;
        Bxn[i][j][k] = B0x;
        Byn[i][j][k] = B0y;
        Bzn[i][j][k] = B0z;
      }
    }
   }
   // initialize B on centers
   grid->interpN2C(Bxc,Bxn);
	grid->interpN2C(Byc,Byn);
	grid->interpN2C(Bzc,Bzn);
	for (int is=0 ; is<ns; is++)
		 grid->interpN2C(rhocs,is,rhons);
	} else {  // READING FROM RESTART
         	if (vct->getCartesian_rank()==0)
	    cout << "LOADING EM FIELD FROM RESTART FILE in " + RestartDirName + "/restart.hdf" << endl;
	stringstream ss;
	ss << vct->getCartesian_rank();
	string name_file = RestartDirName + "/restart" + ss.str() + ".hdf";
	 // hdf stuff 
         hid_t    file_id, dataspace;
         hid_t    datatype, dataset_id;
         herr_t   status;
	 size_t   size;
	 hsize_t     dims_out[3];           /* dataset dimensions */
	 int status_n;
	 
	 // open the hdf file
         file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
         if (file_id < 0){
           cout << "couldn't open file: " << name_file << endl;
	   cout << "RESTART NOT POSSIBLE" << endl;
	 }
	
	 dataset_id = H5Dopen(file_id,"/fields/Bx/cycle_0", H5P_DEFAULT);  // HDF 1.8.8
	 datatype  = H5Dget_type(dataset_id);  
	 size  = H5Tget_size(datatype);
	 dataspace = H5Dget_space(dataset_id);    
	 status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
	
	 
	
	 // Bxn
	double *temp_storage = new double[dims_out[0]*dims_out[1]*dims_out[2]];
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
	H5S_ALL,H5P_DEFAULT,temp_storage);
	int k=0;
	for (int i=1; i < nxn-1; i++)
	      for (int j=1; j <nyn-1; j++)
		    for (int jj=1; jj <nzn-1; jj++)
	          Bxn[i][j][jj] = temp_storage[k++];  
		 
	
	status = H5Dclose(dataset_id);
	
	// Byn
	dataset_id = H5Dopen(file_id,"/fields/By/cycle_0", H5P_DEFAULT);  // HDF 1.8.8
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
	H5S_ALL,H5P_DEFAULT,temp_storage);
	k=0;
	for (int i=1; i < nxn-1; i++)
	      for (int j=1; j <nyn-1; j++)
		    for (int jj=1; jj <nzn-1; jj++)
	          Byn[i][j][jj] = temp_storage[k++]; 
	
	status = H5Dclose(dataset_id);
	
	
	// Bzn
	dataset_id = H5Dopen(file_id,"/fields/Bz/cycle_0", H5P_DEFAULT);  // HDF 1.8.8
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
	H5S_ALL,H5P_DEFAULT,temp_storage);
	k=0;
	for (int i=1; i < nxn-1; i++)
	      for (int j=1; j <nyn-1; j++)
		    for (int jj=1; jj <nzn-1; jj++)
	          Bzn[i][j][jj] = temp_storage[k++]; 
	
	status = H5Dclose(dataset_id);
	
	
	// Ex
	dataset_id = H5Dopen(file_id,"/fields/Ex/cycle_0", H5P_DEFAULT);  // HDF 1.8.8
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
	H5S_ALL,H5P_DEFAULT,temp_storage);
	k=0;
	for (int i=1; i < nxn-1; i++)
	      for (int j=1; j <nyn-1; j++)
		    for (int jj=1; jj <nzn-1; jj++)
	          Ex[i][j][jj] = temp_storage[k++]; 
	
	status = H5Dclose(dataset_id);
	
	
	// Ey 
	dataset_id = H5Dopen(file_id,"/fields/Ey/cycle_0", H5P_DEFAULT);  // HDF 1.8.8
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
	H5S_ALL,H5P_DEFAULT,temp_storage);
	k=0;
	for (int i=1; i < nxn-1; i++)
	      for (int j=1; j <nyn-1; j++)
		   for (int jj=1; jj <nzn-1; jj++)
	        Ey[i][j][jj] = temp_storage[k++]; 
	
	status = H5Dclose(dataset_id);
	
	// Ez 
	dataset_id = H5Dopen(file_id,"/fields/Ez/cycle_0", H5P_DEFAULT);  // HDF 1.8.8
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
	H5S_ALL,H5P_DEFAULT,temp_storage);
	k=0;
	for (int i=1; i < nxn-1; i++)
	      for (int j=1; j <nyn-1; j++)
		    for (int jj=1; jj <nzn-1; jj++)
	        Ez[i][j][jj] = temp_storage[k++]; 
	
	status = H5Dclose(dataset_id);
	
	// open the charge density for species
	
	stringstream *species_name = new stringstream[ns];
	for (int is=0; is < ns;is++){ 
	   species_name[is] << is;
	   string name_dataset = "/moments/species_" + species_name[is].str() + "/rho/cycle_0";
	   dataset_id = H5Dopen(file_id,name_dataset.c_str(), H5P_DEFAULT);  // HDF 1.8.8
	   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storage);
	   k=0;
	   for (int i=1; i < nxn-1; i++)
	      for (int j=1; j <nyn-1; j++)
		    for (int jj=1; jj <nzn-1; jj++)
	          rhons[is][i][j][jj] = temp_storage[k++]; 
	   communicateNode_P(nxn, nyn,nzn, rhons,is,vct);
	   status = H5Dclose(dataset_id);
	 
	 }
	// communicate ghost
	communicateNodeBC(nxn,nyn,nzn,Bxn,1,1,2,2,1,1,vct);
    communicateNodeBC(nxn,nyn,nzn,Byn,1,1,1,1,1,1,vct);
    communicateNodeBC(nxn,nyn,nzn,Bzn,1,1,2,2,1,1,vct);
    // initialize B on centers
	grid->interpN2C(Bxc,Bxn);
	grid->interpN2C(Byc,Byn);
	grid->interpN2C(Bzc,Bzn);
	// communicate ghost
	communicateCenterBC(nxc,nyc,nzc,Bxc,2,2,2,2,2,2,vct);
    communicateCenterBC(nxc,nyc,nzc,Byc,1,1,1,1,1,1,vct);
    communicateCenterBC(nxc,nyc,nzc,Bzc,2,2,2,2,2,2,vct);
	// communicate E
	communicateNodeBC(nxn,nyn,nzn,Ex,1,1,1,1,1,1,vct);
    communicateNodeBC(nxn,nyn,nzn,Ey,1,1,2,2,1,1,vct);
    communicateNodeBC(nxn,nyn,nzn,Ez,1,1,1,1,1,1,vct);
  	for (int is=0 ; is<ns; is++)
		grid->interpN2C(rhocs,is,rhons);
    // close the hdf file
    status = H5Fclose(file_id);    
    delete[] temp_storage;
    delete[] species_name;
  }    
}

/**  initiliaze EM for GEM challange */
inline void EMfields3D::initGEM(VirtualTopology3D *vct, Grid *grid){
	// perturbation localized in X
	double pertX = 0.4;
	double xpert, ypert, exp_pert;
	if (restart1 ==0){
		// initialize
		if (vct->getCartesian_rank() ==0){
			cout << "------------------------------------------" << endl;
			cout << "Initialize GEM Challenge with Pertubation" << endl; 
			cout << "------------------------------------------" << endl;
			cout << "B0x                              = " << B0x << endl;
			cout << "B0y                              = " << B0y << endl;
			cout << "B0z                              = " << B0z << endl;
			cout << "Delta (current sheet thickness) = " << delta << endl;
			for (int i=0; i < ns; i++){
				cout << "rho species " << i <<" = " << rhoINIT[i];
				if (DriftSpecies[i])
					cout << " DRIFTING " << endl;
				else
					cout << " BACKGROUND " << endl;
			}
			cout << "-------------------------" << endl;
		}
		for (int i=0; i < nxn; i++)
			for (int j=0; j < nyn; j++)
				for (int k=0; k < nzn; k++){
					// initialize the density for species
					for (int is=0; is < ns; is++){
						if (DriftSpecies[is])
							rhons[is][i][j][k] = ((rhoINIT[is]/(cosh((grid->getYN(i,j,k)-Ly/2)/delta)*cosh((grid->getYN(i,j,k)-Ly/2)/delta))))/FourPI;
						else
							rhons[is][i][j][k] = rhoINIT[is]/FourPI;
					}
					// electric field
					Ex[i][j][k] =  0.0;
					Ey[i][j][k] =  0.0;
					Ez[i][j][k] =  0.0;
					// Magnetic field
					Bxn[i][j][k] = B0x*tanh((grid->getYN(i,j,k) - Ly/2)/delta);	
					// add the initial GEM perturbation
					// Bxn[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,k)/Lx)*sin(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly  );
					Byn[i][j][k] = B0y; // - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,k)/Lx)*cos(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly); 
					// add the initial X perturbation
					xpert = grid->getXN(i,j,k)- Lx/2;
					ypert = grid->getYN(i,j,k)- Ly/2;
					exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));
					Bxn[i][j][k] +=(B0x*pertX)*exp_pert*(-cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
														 -cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0);
					Byn[i][j][k] +=(B0x*pertX)*exp_pert*(cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
														 +sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0);
					// guide field
					Bzn[i][j][k] = B0z;
				}
					// initialize B on centers
					for (int i=0; i < nxc; i++)
						for (int j=0; j < nyc; j++)
							for (int k=0; k < nzc; k++){
								// Magnetic field
								Bxc[i][j][k] = B0x*tanh((grid->getYC(i,j,k) - Ly/2)/delta);	
								// add the initial GEM perturbation
								//Bxc[i][j][k] += (B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXC(i,j,k)/Lx)*sin(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly  );
								Byc[i][j][k] = B0y; //  - (B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXC(i,j,k)/Lx)*cos(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly); 
								// add the initial X perturbation
								xpert = grid->getXC(i,j,k)- Lx/2;
								ypert = grid->getYC(i,j,k)- Ly/2;
								exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));
								Bxc[i][j][k] +=(B0x*pertX)*exp_pert*(-cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
																	 -cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0);
								Byc[i][j][k] +=(B0x*pertX)*exp_pert*(cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
																	 +sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0);
								// guide field
								Bzc[i][j][k] = B0z;
								
							}
								for (int is=0 ; is<ns; is++)
									grid->interpN2C(rhocs,is,rhons);
	} else {
		init(vct,grid);  // use the fields from restart file
	}
}

/** initialize GEM challenge with no Perturbation */
inline void EMfields3D::initGEMnoPert(VirtualTopology3D *vct, Grid *grid){
  	if (restart1 ==0){
   
	// initialize
	if (vct->getCartesian_rank() ==0){
	 cout << "----------------------------------------------" << endl;
	 cout << "Initialize GEM Challenge without Perturbation" << endl; 
	 cout << "----------------------------------------------" << endl;
	 cout << "B0x                              = " << B0x << endl;
	 cout << "B0y                              = " << B0y << endl;
	 cout << "B0z                              = " << B0z << endl;
	 cout << "Delta (current sheet thickness) = " << delta << endl;
	 for (int i=0; i < ns; i++){
	    cout << "rho species " << i <<" = " << rhoINIT[i];
	    if (DriftSpecies[i])
              cout << " DRIFTING " << endl;
            else
              cout << " BACKGROUND " << endl;
         }
	 cout << "-------------------------" << endl;
	}
	for (int i=0; i < nxn; i++)
		for (int j=0; j < nyn; j++)
		 for (int k=0; k < nzn; k++){
		   // initialize the density for species
		   for (int is=0; is < ns; is++){
		      if (DriftSpecies[is])
			    rhons[is][i][j][k] = ((rhoINIT[is]/(cosh((grid->getYN(i,j,k)-Ly/2)/delta)*cosh((grid->getYN(i,j,k)-Ly/2)/delta))))/FourPI;
			   else
			    rhons[is][i][j][k] = rhoINIT[is]/FourPI;
		   }
		   // electric field
		   Ex[i][j][k] =  0.0;
		   Ey[i][j][k] =  0.0;
		   Ez[i][j][k] =  0.0;
		   // Magnetic field
		   Bxn[i][j][k] = B0x*tanh((grid->getYN(i,j,k) - Ly/2)/delta);	
		   Byn[i][j][k] = B0y; 
		   // guide field
		   Bzn[i][j][k] = B0z;
		}
        // initialize B on centers
	    for (int i=0; i < nxc; i++)
		  for (int j=0; j < nyc; j++)
		   for (int k=0; k < nzc; k++){
		       // Magnetic field
		       Bxc[i][j][k] = B0x*tanh((grid->getYC(i,j,k) - Ly/2)/delta);	
		       Byc[i][j][k] = B0y; 
		       // guide field
		       Bzc[i][j][k] = B0z;
		   
		   }
	    for (int is=0 ; is<ns; is++)
		 grid->interpN2C(rhocs,is,rhons);
        } else {
         init(vct,grid);  // use the fields from restart file
        }
}
/**
*
*  Init Force Free (JxB=0)
*/
inline void EMfields3D::initForceFree(VirtualTopology3D *vct, Grid *grid){
	if (restart1 ==0){
		
		// initialize
		if (vct->getCartesian_rank() ==0){
			cout << "----------------------------------------" << endl;
			cout << "Initialize Force Free with Perturbation" << endl;
			cout << "----------------------------------------" << endl;
			cout << "B0x                              = " << B0x << endl;
			cout << "B0y                              = " << B0y << endl;
			cout << "B0z                              = " << B0z << endl;
			cout << "Delta (current sheet thickness) = " << delta << endl;
			for (int i=0; i < ns; i++){
				cout << "rho species " << i <<" = " << rhoINIT[i];
			}
			cout << "Smoothing Factor = " << Smooth << endl;
			cout << "-------------------------" << endl;
		}
		for (int i=0; i < nxn; i++)
		for (int j=0; j < nyn; j++)
		 for (int k=0; k < nzn; k++){
		   // initialize the density for species
		   for (int is=0; is < ns; is++){
						rhons[is][i][j][k] = rhoINIT[is]/FourPI;
				}
				// electric field
				Ex[i][j][k] =  0.0;
				Ey[i][j][k] =  0.0;
				Ez[i][j][k] =  0.0;
				// Magnetic field
				Bxn[i][j][k] = B0x*tanh((grid->getYN(i,j,k) - Ly/2)/delta);
				// add the initial GEM perturbation
				Bxn[i][j][k] +=(B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,k)/Lx)*sin(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly  );
				Byn[i][j][k] = B0y -(B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,k)/Lx)*cos(M_PI*(grid->getYN(i,j,k)- Ly/2)/Ly);
				// guide field
				Bzn[i][j][k] = B0z/cosh((grid->getYN(i,j,k) - Ly/2)/delta);
			}
		   for (int i=0; i <nxc; i++)
			for (int j=0; j <nyc; j++)
			  for (int k=0; k <nzc; k++) {
				Bxc[i][j][k] = B0x*tanh((grid->getYC(i,j,k) - Ly/2)/delta);
				// add the perturbation
				Bxc[i][j][k] +=(B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXC(i,j,k)/Lx)*sin(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly  );
				Byc[i][j][k] = B0y -(B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXC(i,j,k)/Lx)*cos(M_PI*(grid->getYC(i,j,k)- Ly/2)/Ly);
				// guide field
				Bzc[i][j][k] = B0z/cosh((grid->getYC(i,j,k) - Ly/2)/delta);
			}
	
		for (int is=0 ; is<ns; is++)
			grid->interpN2C(rhocs,is,rhons);
	} else {
		init(vct,grid);  // use the fields from restart file
	}
}
/** Initialize the EM field with constants values or from restart*/
inline void EMfields3D::initBEAM(VirtualTopology3D *vct, Grid *grid, double x_center, double y_center, double z_center, double radius){
        double distance;
        // initialize E and rhos on nodes
        if (restart1==0){
                for (int i=0; i < nxn; i++)
                        for (int j=0; j < nyn; j++)
                          for (int k=0; k < nzn; k++){
                                Ex[i][j][k] =  0.0;
                                Ey[i][j][k] =  0.0;
                                Ez[i][j][k] =  0.0;
                                Bxn[i][j][k] =  0.0;
                                Byn[i][j][k] =  0.0;
                                Bzn[i][j][k] =  0.0;
                                distance = (grid->getXN(i,j,k) - x_center)*(grid->getXN(i,j,k)-x_center)/(radius*radius) + (grid->getYN(i,j,k)-y_center)*(grid->getYN(i,j,k)-y_center)/(radius*radius) + (grid->getZN(i,j,k)-z_center)*(grid->getZN(i,j,k)-z_center)/(4*radius*radius);
                            // plasma
                                rhons[0][i][j][k] = rhoINIT[0]/FourPI;  // initialize with constant density
                                // electrons
                                rhons[1][i][j][k] = rhoINIT[1]/FourPI;
                                // beam
                                if (distance < 1.0)
                                        rhons[2][i][j][k] = rhoINIT[2]/FourPI;
                                else
                                    rhons[2][i][j][k] = 0.0;
							}
							  // initialize B on centers
            for (int i=0; i < nxc; i++)  
                  for (int j=0; j < nyc; j++)
                   for (int k=0; k < nzc; k++){
                       // Magnetic field
                       Bxc[i][j][k] = 0.0;
                       Byc[i][j][k] = 0.0;
                       Bzc[i][j][k] = 0.0;
                       
                                        
                   }
                for (int is=0 ; is<ns; is++)
                 grid->interpN2C(rhocs,is,rhons);
        } else { // EM initialization from RESTART
                init(vct,grid);  // use the fields from restart file
        }
                           
}


/** Calculate the susceptibility on the boundary left*/
inline  void EMfields3D::sustensorLeftY(double** susxy,double** susyy,double** suszy){
          double beta, omcx, omcy, omcz, denom;
		  for(int i=0; i < nxn;i++)
		   for(int k=0; k < nzn;k++){ susxy[i][k] = 0.0; susyy[i][k] = 1.0; suszy[i][k] = 0.0;}
			for (int is=0; is < ns; is++){
				beta = .5*qom[is]*dt/c;
				for(int i=0; i < nxn;i++)
				  for(int k=0; k < nzn;k++){
					omcx = beta*Bxn[i][1][k];
					omcy = beta*Byn[i][1][k];
					omcz = beta*Bzn[i][1][k];
					denom  = FourPI/2*delt*dt/c*qom[is]*rhons[is][i][1][k]/(1.0 + omcx*omcx + omcy*omcy + omcz*omcz);
					susxy[i][k] += (omcz+omcx*omcy)*denom;
					susyy[i][k] += (1.0+omcy*omcy)*denom;
					suszy[i][k] += (omcy*omcz-omcx)*denom;
				}
			}

}
/** Calculate the susceptibility on the boundary right*/
inline  void EMfields3D::sustensorRightY(double** susxy,double** susyy,double** suszy){
		 double beta, omcx, omcy, omcz, denom;
		 for(int i=0; i < nxn;i++)
		   for(int k=0; k < nzn;k++){
				susxy[i][k] = 0.0; susyy[i][k] = 1.0; suszy[i][k] = 0.0;
			}
			for (int is=0; is < ns; is++){
				beta = .5*qom[is]*dt/c;
				for(int i=0; i < nxn;i++)
				  for(int k=0; k < nzn;k++){
					omcx = beta*Bxn[i][nyn-2][k];
					omcy = beta*Byn[i][nyn-2][k];
					omcz = beta*Bzn[i][nyn-2][k];
					denom  = FourPI/2*delt*dt/c*qom[is]*rhons[is][i][nyn-2][k]/(1.0 + omcx*omcx + omcy*omcy + omcz*omcz);
					susxy[i][k] += (omcz+omcx*omcy)*denom;
					susyy[i][k] += (1.0+omcy*omcy)*denom;
					suszy[i][k] += (omcy*omcz-omcx)*denom;
				}
			}
}

/** Perfect conductor boundary conditions: LEFT wall */
inline  void EMfields3D::perfectConductorLeft(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid){
	   double** susxy;
	   double** susyy;
	   double** suszy;
	   switch(dir){
          case 0:  // boundary condition on X-DIRECTION 
          for (int i=1; i <  nyn-1;i++)
            for (int j=1; j <  nzn-1;j++){
                imageX[1][i][j] = vectorX[1][i][j]; 
				imageY[1][i][j] = vectorY[1][i][j];
                imageZ[1][i][j] = vectorZ[1][i][j];
            }
          break;
          case 1: // boundary condition on Y-DIRECTION
		  susxy = newArr(double,nxn,nzn);
          susyy = newArr(double,nxn,nzn);
		  suszy = newArr(double,nxn,nzn);
		  sustensorLeftY(susxy, susyy, suszy);  
		  for (int i=1; i < nxn-1;i++)
              for (int j=1; j <  nzn-1;j++){ 
                 imageX[i][1][j] = vectorX[i][1][j];
				 imageY[i][1][j] = vectorY[i][1][j] - (Ey[i][1][j] - susxy[i][j]*vectorX[i][1][j] - suszy[i][j]*vectorZ[i][1][j] - Jyh[i][1][j]*dt*th*FourPI)/susyy[i][j]; 
                 imageZ[i][1][j] = vectorZ[i][1][j];
              }
          break;
          case 2: // boundary condition on Y-DIRECTION
          for (int i=1; i <  nxn-1;i++)
            for (int j=1; j <  nyn-1;j++){
                 imageX[i][j][1] = vectorX[i][j][1];
                 imageY[i][j][1] = vectorX[i][j][1];
                 imageZ[i][j][1] = vectorZ[i][j][1]; 
            }
          break;
      }
	  delArr(susxy,nxn);
	  delArr(susyy,nxn);
	  delArr(suszy,nxn);
}
/** Perfect conductor boundary conditions: RIGHT wall */
inline  void EMfields3D::perfectConductorRight(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ, int dir,Grid *grid){
      double beta, omcx, omcy, omcz, denom;
	  double** susxy;
	  double** susyy;
	  double** suszy;
	  switch(dir){
          case 0: // boundary condition on X-DIRECTION RIGHT
          
          for (int i=1; i < nyn-1;i++)
            for (int j=1; j <  nzn-1;j++){
                imageX[nxn-2][i][j] = vectorX[nxn-2][i][j];
				imageY[nxn-2][i][j] = vectorY[nxn-2][i][j];
                imageZ[nxn-2][i][j] = vectorZ[nxn-2][i][j];
            }
          break;
          case 1: // boundary condition on Y-DIRECTION RIGHT
          susxy = newArr(double,nxn,nzn);
          susyy = newArr(double,nxn,nzn);
		  suszy = newArr(double,nxn,nzn);
		  sustensorRightY(susxy, susyy, suszy);
          for (int i=1; i < nxn-1;i++)
              for (int j=1; j <  nzn-1;j++){
                 imageX[i][nyn-2][j] = vectorX[i][nyn-2][j];
                 imageY[i][nyn-2][j] = vectorY[i][nyn-2][j] - (Ey[i][nyn-2][j] - susxy[i][j]*vectorX[i][nyn-2][j] - suszy[i][j]*vectorZ[i][nyn-2][j] - Jyh[i][nyn-2][j]*dt*th*FourPI)/susyy[i][j];
				 imageZ[i][nyn-2][j] = vectorZ[i][nyn-2][j];
              }
          break;
          case 2: // boundary condition on Z-DIRECTION RIGHT
	
          for (int i=1; i < nxn-1;i++)
            for (int j=1; j < nyn-1;j++){
                 imageX[i][j][nzn-2] = vectorX[i][j][nzn-2];
                 imageY[i][j][nzn-2] = vectorY[i][j][nzn-2];
				 imageZ[i][j][nzn-2] = vectorZ[i][j][nzn-2];
            }
          break;
      }
	  delArr(susxy,nxn);
	  delArr(susyy,nxn);
	  delArr(suszy,nxn);
}
/** Perfect conductor boundary conditions for source: LEFT WALL*/
inline  void EMfields3D::perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir){
   switch(dir){
    case 0: // boundary condition on X-DIRECTION LEFT
    for (int i=1; i < nyn-1;i++)
       for (int j=1; j < nzn-1;j++){
         vectorX[1][i][j] = 0.0; 
         vectorY[1][i][j] = 0.0;
         vectorZ[1][i][j] = 0.0;
       }
    break;
    case 1: // boundary condition on Y-DIRECTION LEFT
    for (int i=1; i < nxn-1;i++)
       for (int j=1; j < nzn-1;j++){
         vectorX[i][1][j] = 0.0;
         vectorY[i][1][j] = 0.0;
         vectorZ[i][1][j] = 0.0;
       }
    break;
    case 2: // boundary condition on Z-DIRECTION LEFT
    for (int i=1; i < nxn-1;i++)
       for (int j=1; j <  nyn-1;j++){
         vectorX[i][j][1] = 0.0;
         vectorY[i][j][1] = 0.0;
         vectorZ[i][j][1] = 0.0; 
       }
    break;
   }
}

/** Perfect conductor boundary conditions for source: RIGHT WALL*/
inline  void EMfields3D::perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir){
    switch(dir){
    case 0: // boundary condition on X-DIRECTION RIGHT
    for (int i=1; i < nyn-1;i++)
       for (int j=1; j < nzn-1;j++){
         vectorX[nxn-2][i][j] = 0.0; 
         vectorY[nxn-2][i][j] = 0.0;
         vectorZ[nxn-2][i][j] = 0.0;
       }
    break;
    case 1: // boundary condition on Y-DIRECTION RIGHT
    for (int i=1; i < nxn-1;i++)
       for (int j=1; j < nzn-1;j++){
         vectorX[i][nyn-2][j] = 0.0;
         vectorY[i][nyn-2][j] = 0.0;
         vectorZ[i][nyn-2][j] = 0.0;
       }
    break;
    case 2:
    for (int i=1; i <  nxn-1;i++)
       for (int j=1; j <  nyn-1;j++){
         vectorX[i][j][nzn-2] = 0.0;
         vectorY[i][j][nzn-2] = 0.0;
         vectorZ[i][j][nzn-2] = 0.0;
       }
    break;
   }

}
/** get Potential array ***/
inline double*** EMfields3D::getPHI() {return(PHI);}
/** get Ex(X,Y,Z)  */
inline double &EMfields3D::getEx(int indexX, int indexY, int indexZ) const{
  return(Ex[indexX][indexY][indexZ]);}
/** get Electric field  component X array*/ 
inline double*** EMfields3D::getEx() {return(Ex);}
/** get Ey(X,Y,Z)  */
inline double &EMfields3D::getEy(int indexX, int indexY, int indexZ) const{
  return(Ey[indexX][indexY][indexZ]);}
/** get Electric field  component Y array*/ 
inline double*** EMfields3D::getEy() {return(Ey);}
/** get Ez(X,Y,Z)  */
inline double &EMfields3D::getEz(int indexX, int indexY, int indexZ) const{
  return(Ez[indexX][indexY][indexZ]);}
/** get Electric field  component Z array*/ 
inline double*** EMfields3D::getEz() {return(Ez);}
/** get Bx(X,Y,Z)  */
inline double &EMfields3D::getBx(int indexX, int indexY, int indexZ) const{
return(Bxn[indexX][indexY][indexZ]);}
/** get Magnetic Field component X array*/ 
inline double*** EMfields3D::getBx() {return(Bxn);}
/**  get By(X,Y,Z) */
inline double &EMfields3D::getBy(int indexX, int indexY, int indexZ) const{
  return(Byn[indexX][indexY][indexZ]);}
/** get Magnetic Field component Y array*/ 
inline double*** EMfields3D::getBy() {return(Byn);}
/**  get Bz(X,Y,Z) */
inline double &EMfields3D::getBz(int indexX, int indexY, int indexZ) const{
  return(Bzn[indexX][indexY][indexZ]);}
/** get Magnetic Field component Z array*/ 
inline double*** EMfields3D::getBz() {return(Bzn);}
/** get rhoc(X,Y,Z) */
inline double &EMfields3D::getRHOc(int indexX, int indexY, int indexZ) const{
  return(rhoc[indexX][indexY][indexZ]);}
inline double*** EMfields3D::getRHOc() {return(rhoc);}
/** get density on node(indexX,indexY,indexZ)  */
inline double &EMfields3D::getRHOn(int indexX, int indexY, int indexZ) const{ 
return(rhon[indexX][indexY][indexZ]);}
/** get density array defined on nodes*/
inline double*** EMfields3D::getRHOn() {return(rhon);}
/** get rhos(X,Y,Z) : density for species*/
inline double &EMfields3D::getRHOns(int indexX, int indexY, int indexZ, int is) const{
  return(rhons[is][indexX][indexY][indexZ]);}
/** SPECIES: get density array defined on center cells  */
inline double &EMfields3D::getRHOcs(int indexX, int indexY, int indexZ,int is) const{
  return(rhocs[is][indexX][indexY][indexZ]);}
/** get density array defined on nodes*/
inline double**** EMfields3D::getRHOns() {return(rhons);}
/** SPECIES: get pressure tensor component XX defined on nodes */
inline double**** EMfields3D::getpXXsn() {return(pXXsn);}
/** SPECIES: get pressure tensor component XY defined on nodes */
inline double**** EMfields3D::getpXYsn() {return(pXYsn);}
/** SPECIES: get pressure tensor component XZ defined on nodes */
inline double**** EMfields3D::getpXZsn() {return(pXZsn);}
/** SPECIES: get pressure tensor component YY defined on nodes */
inline double**** EMfields3D::getpYYsn() {return(pYYsn);}
/** SPECIES: get pressure tensor component YZ defined on nodes */
inline double**** EMfields3D::getpYZsn() {return(pYZsn);}
/** SPECIES: get pressure tensor component ZZ defined on nodes */
inline double**** EMfields3D::getpZZsn() {return(pZZsn);}
/** get current -Direction X */
inline double &EMfields3D::getJx(int indexX, int indexY,int indexZ) const{
  return(Jx[indexX][indexY][indexZ]);}
/** get current array X component **/ 
inline double*** EMfields3D::getJx() {return(Jx);}
/** get current -Direction Y */
inline double &EMfields3D::getJy(int indexX, int indexY, int indexZ) const{
  return(Jy[indexX][indexY][indexZ]);}
/** get current array Y component **/ 
inline double*** EMfields3D::getJy() {return(Jy);}
/** get current -Direction Z */
inline double &EMfields3D::getJz(int indexX, int indexY, int indexZ) const{
  return(Jz[indexX][indexY][indexZ]);}
/** get current array Z component **/ 
inline double*** EMfields3D::getJz() {return(Jz);}
/**SPECIES: get current array X component */
inline double**** EMfields3D::getJxs() {return(Jxs);}
/** get Jxs(X,Y,Z,is) : density for species*/
inline double &EMfields3D::getJxs(int indexX, int indexY, int indexZ, int is) const{
	return(Jxs[is][indexX][indexY][indexZ]);}
/**SPECIES: get current array Y component */
inline double**** EMfields3D::getJys() {return(Jys);}
/** get Jxs(X,Y,Z,is) : density for species*/
inline double &EMfields3D::getJys(int indexX, int indexY, int indexZ, int is) const{
	return(Jys[is][indexX][indexY][indexZ]);}
/**SPECIES: get current array Z component */
inline double**** EMfields3D::getJzs() {return(Jzs);}
/** get Jxs(X,Y,Z,is) : density for species*/
inline double &EMfields3D::getJzs(int indexX, int indexY, int indexZ, int is) const{
	return(Jzs[is][indexX][indexY][indexZ]);}


/** get the electric field energy*/
inline double EMfields3D::getEenergy(void){
   double localEenergy = 0.0;
   double totalEenergy = 0.0;
   for (int i=1; i < nxn-2; i++)
		for (int j=1; j < nyn-2; j++)
		 for (int k=1; k < nzn-2; k++)
		   localEenergy += .5*dx*dy*dz*(Ex[i][j][k]*Ex[i][j][k] + Ey[i][j][k]*Ey[i][j][k] + Ez[i][j][k]*Ez[i][j][k])/(FourPI);
		   
	MPI_Allreduce(&localEenergy,&totalEenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return(totalEenergy);

}
/** get the magnetic field energy*/
inline double EMfields3D::getBenergy(void){
   double localBenergy = 0.0;
   double totalBenergy = 0.0;
   for (int i=1; i < nxn-2; i++)
		for (int j=1; j < nyn-2; j++)
		 for (int k=1; k < nzn-2; k++)
		   localBenergy += .5*dx*dy*dz*(Bxn[i][j][k]*Bxn[i][j][k] + Byn[i][j][k]*Byn[i][j][k] + Bzn[i][j][k]*Bzn[i][j][k])/(FourPI);
   
   MPI_Allreduce(&localBenergy,&totalBenergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   return(totalBenergy);
}


/** Print info about electromagnetic field */
inline void EMfields3D::print(void) const{


}
/** destructor: deallocate arrays*/
inline EMfields3D::~EMfields3D(){
	// nodes
	delArr3(Ex,nxn,nyn);
	delArr3(Ey,nxn,nyn);
	delArr3(Ez,nxn,nyn);
	delArr3(Exth,nxn,nyn);
	delArr3(Eyth,nxn,nyn);
	delArr3(Ezth,nxn,nyn);
	delArr3(Bxn,nxn,nyn);
	delArr3(Byn,nxn,nyn);
	delArr3(Bzn,nxn,nyn);
	delArr3(rhon,nxn,nyn);
	delArr3(Jx,nxn,nyn);
	delArr3(Jy,nxn,nyn);
	delArr3(Jz,nxn,nyn);
	delArr3(Jxh,nxn,nyn);
	delArr3(Jyh,nxn,nyn);
	delArr3(Jzh,nxn,nyn);
	// nodes and species
	delArr4(rhons,ns,nxn,nyn);
	delArr4(Jxs,ns,nxn,nyn);
	delArr4(Jys,ns,nxn,nyn);
	delArr4(Jzs,ns,nxn,nyn);
	delArr4(pXXsn,ns,nxn,nyn);
	delArr4(pXYsn,ns,nxn,nyn);
	delArr4(pXZsn,ns,nxn,nyn);
	delArr4(pYYsn,ns,nxn,nyn);
	delArr4(pYZsn,ns,nxn,nyn);
	delArr4(pZZsn,ns,nxn,nyn);
	// central points
	delArr3(PHI,nxc,nyc);
	delArr3(Bxc,nxc,nyc);
	delArr3(Byc,nxc,nyc);
	delArr3(Bzc,nxc,nyc);
	delArr3(rhoc,nxc,nyc);
	delArr3(rhoh,nxc,nyc);
	// various stuff needs to be deallocated too
	delArr3(tempXC,nxc,nyc);
	delArr3(tempYC,nxc,nyc);
	delArr3(tempZC,nxc,nyc);
	delArr3(tempXN,nxn,nyn);
	delArr3(tempYN,nxn,nyn);
	delArr3(tempZN,nxn,nyn);
	delArr3(tempC,nxc,nyc);
	delArr3(tempX,nxn,nyn);
	delArr3(tempY,nxn,nyn);
	delArr3(tempZ,nxn,nyn);
	delArr3(temp2X,nxn,nyn);
	delArr3(temp2Y,nxn,nyn);
	delArr3(temp2Z,nxn,nyn);
	delArr3(imageX,nxn,nyn);
	delArr3(imageY,nxn,nyn);
	delArr3(imageZ,nxn,nyn);
	delArr3(Dx,nxn,nyn);
	delArr3(Dy,nxn,nyn);
	delArr3(Dz,nxn,nyn);
	delArr3(vectX,nxn,nyn);
	delArr3(vectY,nxn,nyn);
	delArr3(vectZ,nxn,nyn);
	delArr3(divC,nxc,nyc);
  }
#endif
