/*!************************************************************************* EMfields3D.h - ElectroMagnetic fields definition ------------------- begin : May 2008 copyright : KU Leuven developers : Stefano Markidis, Giovanni Lapenta ************************************************************************* */

#ifndef EMfields3D_H
#define EMfields3D_H


#include <iostream>
#include <sstream>

#include <math.h>
#include <mpi.h>

#include "hdf5.h"

#include "Alloc.h"
#include "Basic.h"
#include "Grid.h"
#include "TransArraySpace3D.h"
#include "CG.h"
#include "GMRES.h"
#include "Collective.h"
#include "ComNodes3D.h"
#include "ComInterpNodes3D.h"
#include "TimeTasks.h"
#include "asserts.h"
#include "BCStructure.h"

using std::cout;
using std::cerr;
using std::endl;


/*! Electromagnetic fields and sources defined for each local grid, and for an implicit maxwell's solver @date May 2008 @par Copyright: (C) 2008 KUL @author Stefano Markidis, Giovanni Lapenta. @version 3.0 */


// class to accumulate node-centered species moments
// 
class Moments {
  private:
    double invVOL;
    double ***rho;

    /** current density, defined on nodes */
    double ***Jx;
    double ***Jy;
    double ***Jz;

    /** pressure tensor components, defined on nodes */
    double ***pXX;
    double ***pXY;
    double ***pXZ;
    double ***pYY;
    double ***pYZ;
    double ***pZZ;
    int nx;
    int ny;
    int nz;
  public:
    int get_nx() const {
      return nx;
    }
    int get_ny() const {
      return ny;
    }
    int get_nz() const {
      return nz;
    }
    double get_invVOL() const {
      return invVOL;
    }
    double get_rho(int i, int j, int k) const {
      return rho[i][j][k];
    }
    double get_Jx(int i, int j, int k) const {
      return Jx[i][j][k];
    }
    double get_Jy(int i, int j, int k) const {
      return Jy[i][j][k];
    }
    double get_Jz(int i, int j, int k) const {
      return Jz[i][j][k];
    }
    double get_pXX(int i, int j, int k) const {
      return pXX[i][j][k];
    }
    double get_pXY(int i, int j, int k) const {
      return pXY[i][j][k];
    }
    double get_pXZ(int i, int j, int k) const {
      return pXZ[i][j][k];
    }
    double get_pYY(int i, int j, int k) const {
      return pYY[i][j][k];
    }
    double get_pYZ(int i, int j, int k) const {
      return pYZ[i][j][k];
    }
    double get_pZZ(int i, int j, int k) const {
      return pZZ[i][j][k];
    }
  public:
    Moments() {
    };
    Moments(int nx_, int ny_, int nz_, double invVOL_);
    ~Moments();
    void set_to_zero();
    void addRho(double weight[][2][2], int X, int Y, int Z);
    void addJx(double weight[][2][2], int X, int Y, int Z);
    void addJy(double weight[][2][2], int X, int Y, int Z);
    void addJz(double weight[][2][2], int X, int Y, int Z);

    void addPxx(double weight[][2][2], int X, int Y, int Z);
    void addPxy(double weight[][2][2], int X, int Y, int Z);
    void addPxz(double weight[][2][2], int X, int Y, int Z);
    void addPyy(double weight[][2][2], int X, int Y, int Z);
    void addPyz(double weight[][2][2], int X, int Y, int Z);
    void addPzz(double weight[][2][2], int X, int Y, int Z);
};

// construct empty instance (not zeroed)
inline Moments::Moments(int nx_, int ny_, int nz_, double invVOL_) {
  nx = nx_;
  ny = ny_;
  nz = nz_;
  invVOL = invVOL_;
  rho = newArr3(double, nx, ny, nz);
  Jx = newArr3(double, nx, ny, nz);
  Jy = newArr3(double, nx, ny, nz);
  Jz = newArr3(double, nx, ny, nz);
  pXX = newArr3(double, nx, ny, nz);
  pXY = newArr3(double, nx, ny, nz);
  pXZ = newArr3(double, nx, ny, nz);
  pYY = newArr3(double, nx, ny, nz);
  pYZ = newArr3(double, nx, ny, nz);
  pZZ = newArr3(double, nx, ny, nz);
}

inline Moments::~Moments() {
  // nodes and species
  delArr3(rho, nx, ny);
  delArr3(Jx, nx, ny);
  delArr3(Jy, nx, ny);
  delArr3(Jz, nx, ny);
  delArr3(pXX, nx, ny);
  delArr3(pXY, nx, ny);
  delArr3(pXZ, nx, ny);
  delArr3(pYY, nx, ny);
  delArr3(pYZ, nx, ny);
  delArr3(pZZ, nx, ny);
}

inline void Moments::set_to_zero() {
  // #pragma omp parallel for collapse(1)
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++) {
        rho[i][j][k] = 0.0;
        Jx[i][j][k] = 0.0;
        Jy[i][j][k] = 0.0;
        Jz[i][j][k] = 0.0;
        pXX[i][j][k] = 0.0;
        pXY[i][j][k] = 0.0;
        pXZ[i][j][k] = 0.0;
        pYY[i][j][k] = 0.0;
        pYZ[i][j][k] = 0.0;
        pZZ[i][j][k] = 0.0;
      }
}

class EMfields3D                // :public Field
{
  public:

    /*! constructor */
  /*! pre-mlmd:
    EMfields3D(Collective * col, Grid * grid); */
  /* mlmd: I need the topology in the constructore*/
  EMfields3D(Collective * col, Grid * grid, VirtualTopology3D * vct); 
    /*! destructor */
    ~EMfields3D();

    /*! initialize the electromagnetic fields with constant values */
    void init(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! init beam */
    void initBEAM(VirtualTopology3D * vct, Grid * grid, Collective *col, double x_center, double y_center, double z_center, double radius);
    /*! initialize GEM challenge */
    void initGEM(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initOriginalGEM(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initDoublePeriodicHarrisWithGaussianHumpPerturbation(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! initialize GEM challenge with dipole-like tail without perturbation */
    void initGEMDipoleLikeTailNoPert(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! initialize GEM challenge with no Perturbation */
    void initGEMnoPert(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! initialize from BATSRUS */
    void initBATSRUS(VirtualTopology3D * vct, Grid * grid, Collective * col);
    /*! Random initial field */
    void initRandomField(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! Init Force Free (JxB=0) */
    void initForceFree(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! initialized with rotated magnetic field */
    void initEM_rotate(VirtualTopology3D * vct, Grid * grid, Collective *col, double B, double theta);
    /*! initialized with Light Wave */
    void initLightWave(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initTestIntProj(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! double GEM, mlmd ready */
    void initDoubleGEM(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initDoubleGEM_CentralPerturbation(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initTestProjection(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initTestBC(VirtualTopology3D * vct, Grid * grid, Collective *col); 
    void initTestFix3B(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! inits all the fields to 0, particles only on CG, disables GMRES on RG*/
    void initMAX_Show_RG_BC(VirtualTopology3D * vct, Grid * grid, Collective *col);
    void initTestBBoundary(VirtualTopology3D * vct, Grid * grid, Collective *col);
    /*! add a perturbattion to charge density */
    void AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid * grid);
    /*! add a perturbattion to the EM field */
    void AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid * grid);
    /*! Initialise a combination of magnetic dipoles */
    void initDipole(VirtualTopology3D *vct, Grid *grid, Collective *col);
    void initDipole_2(VirtualTopology3D *vct, Grid *grid, Collective *col);
    void SetDipole_2Bext(VirtualTopology3D *vct, Grid *grid, Collective *col);
    void SetDipole_3Bext(VirtualTopology3D *vct, Grid *grid, Collective *col);

    /*! Calculate Electric field using the implicit Maxwell solver */
    void startEcalc(Grid * grid, VirtualTopology3D * vct, Collective *col);
    void endEcalc(double* xkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col);
    void endEcalc_ECSIM(double* xkrylov, VirtualTopology3D * vct, Collective *col, Grid* grid);
    /* implement divE correction, on faces (to be used in RG) */
    void divECorrection_AllFaces(double *** FX, double *** FY, double *** FZ, Grid * grid, VirtualTopology3D * vct);
    void divECorrection_OneFace(double *** FX, double *** FY, double *** FZ, Grid * grid, VirtualTopology3D * vct, int nxn, int nyn, int nzn, int DIR, int SIDE);
    void calculateE(Grid * grid, VirtualTopology3D * vct, Collective *col);
    // ECSIM
    void calculateEB_ECSIM(Grid * grid, VirtualTopology3D * vct, Collective *col, int cycle); 
    /*! Image of Poisson Solver (for SOLVER) */
    void PoissonImage(double *image, double *vector, Grid * grid, VirtualTopology3D * vct);
    /*! Image of Poisson solver (for FACE SOLVER) */
    void PoissonImage_2D(double *image, double *vector, Grid * grid, VirtualTopology3D * vct);
    void PoissonImage_2D(double *image, double *vector, Grid * grid, VirtualTopology3D * vct, int nxc, int nyc, int nzc);
    /*! Image of Maxwell Solver (for Solver) */
    void MaxwellImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct);
    /*! Maxwell source term (for SOLVER) */
    void MaxwellSource(double *bkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col);
    /* Maxwell source term for ECSIM solver */
    void MaxwellSource_ECSIM(double *bkrylov, Grid * grid, VirtualTopology3D * vct, Collective *col);
    /* Image of Maxwell Solver for ECSIM */
    void MaxwellImage_ECSIM(double *im, double *vector, Grid * grid, VirtualTopology3D * vct, Collective *col);
    /* interpolate B field centers to nodes - for ECSIM */
    void centers2nodesB(VirtualTopology3D * vct, Grid *grid, Collective * col);
    /*! Impose a constant charge inside a spherical zone of the domain */
    void ConstantChargePlanet(Grid * grid, VirtualTopology3D * vct, double R, double x_center, double y_center, double z_center);
    /*! Impose a constant charge in the OpenBC boundaries */
    void ConstantChargeOpenBC(Grid * grid, VirtualTopology3D * vct);
    /*! Impose a constant charge in the OpenBC boundaries */
    void ConstantChargeOpenBCv2(Grid * grid, VirtualTopology3D * vct);
    /*! Calculate Magnetic field with the implicit solver: calculate B defined on nodes With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
    void calculateB(Grid * grid, VirtualTopology3D * vct, Collective *col, int cycle);
    /*! fix B on the boundary for gem challange */
    void fixBgem(Grid * grid, VirtualTopology3D * vct);
    /*! fix B on the boundary for gem challange */
    void fixBforcefree(Grid * grid, VirtualTopology3D * vct);
    /* fix B ghost cell after receive B node BC */
    void fixBghostCells_Left(int Dir, int NCells);
    void fixBghostCells_Right(int Dir, int NCells);
    /*! Calculate the three components of Pi(implicit pressure) cross image vector */
    void PIdot(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid * grid);
    /*! Calculate the three components of mu (implicit permeattivity) cross image vector */
    void MUdot(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid * grid, VirtualTopology3D * vct);
    /*! Calculate rho hat, Jx hat, Jy hat, Jz hat */
    void calculateHatFunctions(Grid * grid, VirtualTopology3D * vct);

    void UpdateRHOcs(Grid * grid);
    void SetLambda  (Grid * grid, VirtualTopology3D * vct);
    double ***GetLambda();

    /*! communicate ghost for densities and interp rho from node to center */
    void interpDensitiesN2C(VirtualTopology3D * vct, Grid * grid);
    /*! set to 0 all the densities fields */
    void setZeroDensities();
    /*! Sum rhon over species */
    void sumOverSpecies(VirtualTopology3D * vct);
    /*! Sum current over different species */
    void sumOverSpeciesJ();
    /*! Smoothing after the interpolation* */
    void smooth(double value, double ***vector, int type, Grid * grid, VirtualTopology3D * vct);
    /*! SPECIES: Smoothing after the interpolation for species fields* */
    void smooth(double value, double ****vector, int is, int type, Grid * grid, VirtualTopology3D * vct);
    /*! smooth the electric field */
    void smoothE(double value, VirtualTopology3D * vct, Collective *col);
    /*! without communication, to use in applyProjection */
    void smoothE_NoComm(double value, VirtualTopology3D * vct, Collective *col);
    /* smoothing applied specifically to external faces */
    void smoothFaces(double value, double *** F, Grid * grid, VirtualTopology3D * vct, int type);
    void smoothFaces(double value, double **** F, Grid * grid, VirtualTopology3D * vct, int is, int type);
    /* external faces + interior */
    void smoothAll(double value, double **** F, Grid * grid, VirtualTopology3D * vct, int is, int type);
    void smoothAll(double value, double *** F, Grid * grid, VirtualTopology3D * vct, int type);
    /*! communicate ghost for grid -> Particles interpolation */
    void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct, Grid * grid);
    /*! add accumulated moments to the moments for a given species */
    void addToSpeciesMoments(const Moments & in, int is);
    /*! add an amount of charge density to charge density field at node X,Y,Z */
    void addRho(double weight[][2][2], int X, int Y, int Z, int is);
    /*! for debugging purposes */
    void addRho(double weight[][2][2], int X, int Y, int Z, int is, VirtualTopology3D * vct, double xp, double yp, double zp);
    /*! add an amount of current density - direction X to current density field at node X,Y,Z */
    void addJx(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of current density - direction Y to current density field at node X,Y,Z */
    void addJy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of current density - direction Z to current density field at node X,Y,Z */
    void addJz(double weight[][2][2], int X, int Y, int Z, int is);

    /*! add an amount of pressure density - direction XX to current density field at node X,Y,Z */
    void addPxx(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction XY to current density field at node X,Y,Z */
    void addPxy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
    void addPxz(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction YY to current density field at node X,Y,Z */
    void addPyy(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
    void addPyz(double weight[][2][2], int X, int Y, int Z, int is);
    /*! add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
    void addPzz(double weight[][2][2], int X, int Y, int Z, int is);

    /*! adjust densities on boundaries that are not periodic */
    void adjustNonPeriodicDensities(int is, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct);


    /*! Perfect conductor boundary conditions LEFT wall */
    void perfectConductorLeft(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid);
    /*! Perfect conductor boundary conditions RIGHT wall */
    void perfectConductorRight(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid);
    /*! Perfect conductor boundary conditions for source LEFT wall */
    void perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir);
    /*! Perfect conductor boundary conditions for source RIGHT wall */
    void perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir);

    /*! Calculate the sysceptibility tensor on the boundary */
    void sustensorX(double **susxx, double **susxy, double **susxz, int N);
    void sustensorY (double **susyx, double **susyy, double **susyz, int N);
    void sustensorZ(double **suszx, double **suszy, double **suszz, int N);

    /*! get Potential array */
    double ***getPHI();
    /*! get Electric Field component X defined on node(indexX,indexY,indexZ) */
    double &getEx(int indexX, int indexY, int indexZ) const;
    /*! get Electric field X component array */
    double ***getEx();
    double ***getExth();
    /*! get Electric field X component cell array without the ghost cells */
    double ***getExc();
    /*! get Electric Field component Y defined on node(indexX,indexY,indexZ) */
    double &getEy(int indexX, int indexY, int indexZ) const;
    /*! get Electric field Y component array */
    double ***getEy();
    double ***getEyth();
    /*! get Electric field Y component cell array without the ghost cells */
    double ***getEyc();
    /*! get Electric Field component Z defined on node(indexX,indexY,indexZ) */
    double &getEz(int indexX, int indexY, int indexZ) const;
    /*! get Electric field Z component array */
    double ***getEz();
    double ***getEzth();
    /*! get Electric field Z component cell array without the ghost cells */
    double ***getEzc();
    /*! get Magnetic Field component X defined on node(indexX,indexY,indexZ) */
    double &getBx(int indexX, int indexY, int indexZ) const;
    /*! get Magnetic field X component array */
    double ***getBx();
    /*! get Magnetic field X component cell array without the ghost cells */
    double getBxc(int i, int j, int k);
    /*! get Magnetic field X component cell array without the ghost cells */
    double ***getBxc();
    /*! get Magnetic Field component Y defined on node(indexX,indexY,indexZ) */
    double &getBy(int indexX, int indexY, int indexZ) const;
    /*! get Magnetic field Y component array */
    double ***getBy();
    /*! get Magnetic field Y component cell array without the ghost cells */
    double getByc(int i, int j, int k);
    /*! get Magnetic field Y component cell array without the ghost cells */
    double ***getByc();
    /*! get Magnetic Field component Z defined on node(indexX,indexY,indexZ) */
    double &getBz(int indexX, int indexY, int indexZ) const;
    /*! get Magnetic field Z component array */
    double ***getBz();
    /*! get Magnetic field Z component cell array without the ghost cells */
    double getBzc(int i, int j, int k);
    /*! get Magnetic field Z component cell array without the ghost cells */
    double ***getBzc();
    /*! get density on cell(indexX,indexY,indexZ) */
    double &getRHOc(int indexX, int indexY, int indexZ) const;
    /*! get density array on center cell */
    double ***getRHOc();
    /*! get density on nodes(indexX,indexY,indexZ) */
    double &getRHOn(int indexX, int indexY, int indexZ) const;
    /*! get density array on nodes */
    double ***getRHOn();
    /*! SPECIES: get density on nodes(indexX,indexY,indexZ) */
    double &getRHOns(int indexX, int indexY, int indexZ, int is) const;
    /*! SPECIES: get density on center cell(indexX,indexY,indexZ) */
    double &getRHOcs(int indexX, int indexY, int indexZ, int is) const;
    /*! SPECIES: get density on center cell(indexX,indexY,indexZ) */
    double ****& getRHOcs();
    /*! SPECIES: get density array on nodes */
    double ****getRHOns();
    /*! SPECIES: get density array on nodes */
    double ***& getRHOns(int is);
    /*! SPECIES: get density array on cells without the ghost cells */
    double ***getRHOcs(int is);
    double ***& getRHOcs(int is, int dummy);

    /** get Magnetic Field component X defined on node(indexX,indexY,indexZ) */
    double &getBx_ext(int indexX, int indexY, int indexZ) const;
    /** get Magnetic Field component Y defined on node(indexX,indexY,indexZ) */
    double &getBy_ext(int indexX, int indexY, int indexZ) const;
    /** get Magnetic Field component Z defined on node(indexX,indexY,indexZ) */
    double &getBz_ext(int indexX, int indexY, int indexZ) const;

    /** get Magnetic Field component X */
    double ***getBx_ext();
    /** get Magnetic Field component Y */
    double ***getBy_ext();
    /** get Magnetic Field component Z */
    double ***getBz_ext();

    double getRHOINIT(int is, int i, int j, int k);
    
    double ***&getBxTot();
    double ***&getByTot();
    double ***&getBzTot();

    void UpdateFext(int cycle);
    void UpdateCycle(int cycle);
    double getFext();

    /*! get pressure tensor XX for species */
    double ****getpXXsn();
    /*! get pressure tensor XY for species */
    double ****getpXYsn();
    /*! get pressure tensor XZ for species */
    double ****getpXZsn();
    /*! get pressure tensor YY for species */
    double ****getpYYsn();
    /*! get pressure tensor YZ for species */
    double ****getpYZsn();
    /*! get pressure tensor ZZ for species */
    double ****getpZZsn();

    /*! get Jx(X,Y,Z) */
    double &getJx(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction X */
    double ***getJx();
    /*! get Jxs(X,Y,Z,is) */
    double &getJxs(int indexX, int indexY, int indexZ, int is) const;
    /*! SPECIES: get current -Direction X */
    double ****getJxs();
    /*! SPECIES: get current X component for species is in all cells except ghost */
    double ***getJxsc(int is);
    /*! get Jy(X,Y,Z) */
    double &getJy(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction Y */
    double ***getJy();
    /*! get Jys(X,Y,Z,is) */
    double &getJys(int indexX, int indexY, int indexZ, int is) const;
    /*! SPECIES: get current -Direction Y */
    double ****getJys();
    /*! SPECIES: get current Y component for species is in all cells except ghost */
    double ***getJysc(int is);
    /*! get Jz(X,Y,Z) */
    double &getJz(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction Z */
    double ***getJz();
    /*! get Jzs(X,Y,Z,is) */
    double &getJzs(int indexX, int indexY, int indexZ, int is) const;
    /*! SPECIES: get current -Direction Z */
    double ****getJzs();
    /*! SPECIES: get current Z component for species is in all cells except ghost */
    double ***getJzsc(int is);

    double ***&getJxs(int is);
    double ***&getJys(int is);
    double ***&getJzs(int is);

    /*! get the electric field energy */
    /*! mlmd: i need the communicator also
      double getEenergy(); */
    double getEenergy(MPI_Comm Comm);
    /*! get the magnetic field energy */
    /*! mlmd: i need the communicator also
      double getBenergy(); */
    double getBenergy(MPI_Comm Comm); 

    /*! print electromagnetic fields info */
    void print(void) const;

    // OpenBC
    void updateInfoFields(Grid *grid,VirtualTopology3D *vct,Collective *col);
    bool getDriftSpecies(int is){
      return DriftSpecies[is];
    }
    double getDelta(){
      return delta;
    }

    /*! mlmd specific functions */
    /* children: calculate which points they need to communicate to which parent core
       parent grid: recceives the info
       here, fieldBC communication is set up */
    void initWeightBC(Grid *grid, VirtualTopology3D *vct);
    void initWeightBC_InitialInterpolation(Grid *grid, VirtualTopology3D *vct);
    void initWeightBC_InitialInterpolation_Phase1(Grid *grid, VirtualTopology3D *vct);
    void initWeightProj(Grid *grid, VirtualTopology3D *vct);
    void initWeightProj_Phase1(Grid *grid, VirtualTopology3D *vct);
    // used in initWeightProj_Phase1, initWeightBCBuffer_Phase1, copied form particles
    void Explore3DAndCommit(Grid *grid, int i_s, int i_e, int j_s, int j_e, int k_s, int k_e, RGBC_struct *RGBC_Info, int *numMsg, int *MaxSizeMsg, VirtualTopology3D * vct, char  FACE );
    void sendProjection(Grid *grid, VirtualTopology3D *vct);
    void receiveProjection(Grid *grid, VirtualTopology3D *vct);
    void applyProjection(Grid *grid, VirtualTopology3D *vct, Collective *col);
    void TestProjection(Grid *grid, VirtualTopology3D *vct);

    /* to create the MPI_Datatype associate to RGBC_struct */
    void MPI_RGBC_struct_commit();

    /* different phases of initWeightBC */
    /* performs the real BC calculations; to be reused for ghost (which =-1, *_Ghost)
       and active nodes (which= 0; *Active) */
    void initWeightBC_Phase1(Grid *grid, VirtualTopology3D *vct, RGBC_struct * RGBC_Info, int* RG_numBCMessages, int which);
    // for the inital interpolation
    void initWeightBC_InitialInterpolation_Phase1(Grid *grid, VirtualTopology3D *vct, RGBC_struct * RGBC_Info, int* RG_numBCMessages);
    void initWeightBCBuffer_Phase1(Grid *grid, VirtualTopology3D *vct, RGBC_struct *RGBC_Info, int *numMsg, int *MaxSizeMsg);
    void initWeightBCfix3B_Phase1(Grid *grid, VirtualTopology3D *vct, RGBC_struct *RGBC_Info, int *numMsg, int *MaxSizeMsg);

    /* return the factor for the BC multiplication */
    double BufferFactor(char Dir, int X, int Y, int Z, int BufLenX, int BufLenY, int BufLenZ);

    /* phase 2a of initWeightBC:                                                                          
   core 0 of each child grid receives all the messages to send to the corresponding coarse grid           
   a level-wide message structure is built */
    void initWeightBC_Phase2a(Grid *grid, VirtualTopology3D *vct, RGBC_struct *RGBC_Info_LevelWide, int * RG_numBCMessages_LevelWide, RGBC_struct *RGBC_Info, int RG_numBCMessages, int which);

    /* phase 2b of initWeightBC:
       core 0 of the child grid assembles & sends messages for all CG cores*/
    void initWeightBC_Phase2b(Grid *grid, VirtualTopology3D *vct, RGBC_struct *RGBC_Info_LevelWide,int RG_numBCMessages_LevelWide, int which);

    /* phase 2c, only for the parent     
       each CG core will receive a message per child grid  */
    void initWeightBC_Phase2c(Grid *grid, VirtualTopology3D *vct, RGBC_struct ** CG_Info, int* CG_numBCMessages, int which);
    /* end different phases of initWeightBC */

    /* Trim functions: RGBC_Info_** and CG_Info_** are initially allocated very large for
       internal reasons;
       once initWeightBC has finished, we want to free to unused memory */
    void Trim_RGBC_Vectors(VirtualTopology3D *vct);
    /* same for interpolation vectors */
    void Trim_RGBC_Vectors_II(VirtualTopology3D *vct);
    /* same for projection */
    void Trim_Proj_Vectors(VirtualTopology3D *vct);
    /* mlmd: BC related functions */
    /* sendBC: coarse grids sends BCs to the refined grids */
    void sendBC(Grid *grid, VirtualTopology3D *vct);
    void sendInitialInterpolation(Grid *grid, VirtualTopology3D *vct);
    /* receiveBC: refined grids receive BCs from the coarse grids */
    void receiveBC(Grid *grid, VirtualTopology3D *vct);
    void receiveInitialInterpolation(Grid *grid, VirtualTopology3D *vct);
    void ApplyInitialInterpolation(VirtualTopology3D *vct, Grid * grid);
    void DeallocateII();
    /* Initial interpolation Ã*/
    void sendBC_InitialInterpolation(Grid *grid, VirtualTopology3D *vct);
    void receiveBC_InitialInterpolation(Grid *grid, VirtualTopology3D *vct);
    /* end Initial interpolation */
    /* the RG sets the received BC 
       Fx, Fy, Fz: the field where BC are applied
       Fx_BC, Fy_BC, Fz_BC: the BCs
       RGBC_Info: the RGBC_Info struct (Active or Ghost)
       RG_numBCMessages: number of messages (Active or Ghost) */
    void setBC_Nodes(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages);
    void setBC_Nodes_TwoLess(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages);
    /* for when I need to put the BC on the image */
    void setBC_NodesImage(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***vectX, double ***vectY, double ***vectZ, double ***Fz, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages);
    /* this function is exactly the same as the previous setBC_NodesImage 
 only differences: unnecessary inputs deleted & vector size included
 (this just for the safety check, to be deleted in later versions)*/
    void setBC_NodesImage(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double ***vectX, double ***vectY, double ***vectZ, RGBC_struct * RGBC_Info, int RG_numBCMessages, int nx, int ny, int nz);

    void setBC_Nodes_RENORM(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages);
    /* for when I need to put the BC on the image */
    void setBC_NodesImage_RENORM(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***vectX, double ***vectY, double ***vectZ, double ***Fz, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages);
    /* average native RG and interpolated CG solution */
    void AverageBC_Nodes(VirtualTopology3D * vct, double ***Fx, double ***Fy, double ***Fz, double **Fx_BC, double **Fy_BC, double **Fz_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages);

    /* set average between CG BC and local intermediate solution as BC */
    void averageBC_BufferSource(VirtualTopology3D * vct, double *** SX, double *** SY, double *** SZ, double *** vectX, double *** vectY, double *** vectZ, double ** Exth_Buffer_BC, double ** Eyth_Buffer_BC, double ** Ezth_Buffer_BC, RGBC_struct * RGBC_Info_Buffer, int RG_numBCMessages_Buffer);
    void setBC_Nodes(VirtualTopology3D * vct, double **Fx_BC, double **Fy_BC, double **Fz_BC, double **Fa_BC, RGBC_struct * RGBC_Info, int RG_numBCMessages);
    
    /* set MaxwellSource BC in mlmd; NB: you need to set the actives */
    void MLMDSourceLeft(double ***vectorX, double ***vectorY, double ***vectorZ, int dir);
    void MLMDSourceRight(double ***vectorX, double ***vectorY, double ***vectorZ, int dir);
        
    /* CG core builds the BC msg for a RG core */
    /* ch: number of the child
       RGInfo: the message that you got from the refined grid, which holds all the info you need
       Msg: stuff to send to the RG; this is already the info POINT BY POINT
       NP: # of points to send
       Msg: the msg where to pack this
     */
    void buildBCMsg(VirtualTopology3D *vct, Grid * grid, int ch, RGBC_struct RGInfo,  int NP , double *Msg);
    void buildIIMsg(VirtualTopology3D *vct, Grid * grid, int ch, RGBC_struct RGInfo, int Size );
    void buildFix3BMsg(VirtualTopology3D *vct, Grid * grid, int ch, RGBC_struct RGInfo, int Size, double *Msg );

    /* CG_Info: the entry with the info for that msg
       ch: number of child
       which: active or ghost */
    void sendOneBC(VirtualTopology3D * vct, Grid * grid,  RGBC_struct CG_Info, int ch, int which);
    /* end mlmd: BC related functions */
    /* a barrier on both parent and child side of the field communicator, to prevent messages from crossing */
    void MPI_Barrier_ParentChild(VirtualTopology3D* vct);    
    /* copy moment vector, species ns */
    void copyMoments(Grid * grid, VirtualTopology3D * vct, double ***P_rho, double ***P_Jx, double ***P_Jy, double ***P_Jz, double ***P_pxx, double ***P_pxy, double ***P_pxz, double ***P_pyy, double ***P_pyz, double ***P_pzz, int is);

    /* some ops done after all init **/
    void PostInit();
    /* this for testing purposes */

    void TESTGhost( double **** vec);
    void TESTGhost (double **** dest, double **** source); 
    /* to enforce divB= 0 on B centers */
    void correctDivB(Grid * grid, VirtualTopology3D * vct);
    /* Image for the divB=0 correction */
    void BPoissonImage(double *image, double *vector, Grid * grid, VirtualTopology3D * vct);
    /*! end mlmd specific functions */
    /* ********************************* // VARIABLES ********************************* */
  private:
    /*! light speed */
    double c;
    /* 4*PI for normalization */
    double FourPI;
    /*! time step */
    double dt;
    /*! decentering parameter */
    double th;
    /*! Smoothing value */
    double Smooth;
    /* wether to smooth */
    bool SmoothGrid;
    /*! delt = c*th*dt */
    double delt;
    /*! number of particles species */
    int ns;
    /*! GEM challenge parameters */
    double B0x, B0y, B0z, delta;
    /** Earth Model parameters */
    double B1x, B1y, B1z;
    /*! charge to mass ratio array for different species */
    double *qom;
    /*! Boundary electron speed */
    double ue0, ve0, we0;

    // KEEP IN MEMORY GUARD CELLS ARE INCLUDED
    /*! mlmd: these values are for the local grid */

    /*! number of cells - X direction, including + 2 (guard cells) */
    int nxc;
    /*! number of nodes - X direction, including + 2 extra nodes for guard cells */
    int nxn;
    /*! number of cell - Y direction, including + 2 (guard cells) */
    int nyc;
    /*! number of nodes - Y direction, including + 2 extra nodes for guard cells */
    int nyn;
    /*! number of cell - Z direction, including + 2 (guard cells) */
    int nzc;
    /*! number of nodes - Z direction, including + 2 extra nodes for guard cells */
    int nzn;
    /*! local grid boundaries coordinate */
    double xStart, xEnd, yStart, yEnd, zStart, zEnd;
    /*! grid spacing */
    double dx, dy, dz, invVOL;
    /*! simulation box length - X direction */
    double Lx;
    /*! simulation box length - Y direction */
    double Ly;
    /*! simulation box length - Z direction */
    double Lz;

    /* mlmd: i need dx, dy, dz of the children */
    double *dx_Ch;
    double *dy_Ch;
    double *dz_Ch;
    
    /** parent grid spacing **/
    double DxP, DyP, DzP;
    /* RF, used in split */
    double RFx, RFy, RFz;

    /*! end mlmd: these values are for the local grid */

    /** source center - X direction   */
    double x_center;
    /** source center - Y direction   */
    double y_center;
    /** source center - Z direction   */
    double z_center;
    /** Characteristic length */
    double L_square;

    /*! PHI: electric potential (indexX, indexY, indexZ), defined on central points between nodes */
    double ***PHI;
    /*! Ex: electric field X-component (indexX, indexY, indexZ), defined on nodes */
    double ***Ex;
    /*! Exth: implicit electric field X-component (indexX, indexY, indexZ), defined on nodes */
    double ***Exth;
    /*! Ey: electric field Y-component (indexX, indexY, indexZ), defined on nodes */
    double ***Ey;
    /*! Eyth: implicit electric field Y-component (indexX, indexY, indexZ), defined on nodes */
    double ***Eyth;
    /*! Ez: electric field Z-component (indexX, indexY, indexZ, #species), defined on nodes */
    double ***Ez;
    /*! Ezth: implicit electric field Z-component (indexX, indexY, indexZ), defined on nodes */
    double ***Ezth;
    /*! Bxc: magnetic field X-component (indexX, indexY, indexZ), defined on central points between nodes */
    double ***Bxc;
    /*! Byc: magnetic field Y-component (indexX, indexY, indexZ), defined on central points between nodes */
    double ***Byc;
    /*! Bzc: magnetic field Z-component (indexX, indexY, indexZ), defined on central points between nodes */
    double ***Bzc;
    /*! Bxn: magnetic field X-component (indexX, indexY, indexZ), defined on nodes */
    double ***Bxn;
    /*! Byn: magnetic field Y-component (indexX, indexY, indexZ), defined on nodes */
    double ***Byn;
    /*! Bzn: magnetic field Z-component (indexX, indexY, indexZ), defined on nodes */
    double ***Bzn;

    // *************************************
    // TEMPORARY ARRAY
    // ************************************
    /*!some temporary arrays (for calculate hat functions) */
    double ***tempXC;
    double ***tempYC;
    double ***tempZC;
    double ***tempXN;
    double ***tempYN;
    double ***tempZN;
    /*! other temporary arrays (in MaxwellSource) */
    double ***tempC;
    double ***tempX;
    double ***tempY;
    double ***tempZ;
    double ***temp2X;
    double ***temp2Y;
    double ***temp2Z;
    /*! and some for MaxwellImage */
    double ***imageX;
    double ***imageY;
    double ***imageZ;
    double ***imageEX;
    double ***imageEY;
    double ***imageEZ;
    double ***imageBX;
    double ***imageBY;
    double ***imageBZ;
    double ***Dx;
    double ***Dy;
    double ***Dz;
    double ***vectX;
    double ***vectY;
    double ***vectZ;
    double ***divC;
    double ***arr;


    // *******************************************************************************
    // *********** SOURCES **
    // *******************************************************************************

    /*! Charge density, defined on central points of the cell */
    double ***rhoc;
    /*! Charge density, defined on nodes */
    double ***rhon;
    /*! Implicit charge density, defined on central points of the cell */
    double ***rhoh;
    /*! SPECIES: charge density for each species, defined on nodes */
    double ****rhons;
    /*! SPECIES: charge density for each species, defined on central points of the cell */
    double ****rhocs;
    /*! Current density component-X, defined on nodes */
    double ***Jx;
    /*! Current density component-Y, defined on nodes */
    double ***Jy;
    /*! Current density component-Z, defined on nodes */
    double ***Jz;
    /*! Implicit current density X-component, defined on nodes */
    double ***Jxh;
    /*! Implicit current density Y-component, defined on nodes */
    double ***Jyh;
    /*! Implicit current density Z-component, defined on nodes */
    double ***Jzh;
    /*! SPECIES: current density component-X for species, defined on nodes */
    double ****Jxs;
    /*! SPECIES: current density component-Y for species, defined on nodes */
    double ****Jys;
    /*! SPECIES: current density component-Z for species, defined on nodes */
    double ****Jzs;
    /*! External magnetic field component-X, defined on nodes */
    double***  Bx_ext;
    /*! External magnetic field component-Y, defined on nodes */
    double***  By_ext;
    /*! External magnetic field component-Z, defined on nodes */
    double***  Bz_ext;
    /*! External current field component-X, defined on nodes */
    double***  Jx_ext;
    /*! External current field component-Y, defined on nodes */
    double***  Jy_ext;
    /*! External current field component-Z, defined on nodes */
    double***  Jz_ext;

    double Fext;

    int currentCycle;

    /*! SPECIES: pressure tensor component-XX, defined on nodes */
    double ****pXXsn;
    /*! SPECIES: pressure tensor component-XY, defined on nodes */
    double ****pXYsn;
    /*! SPECIES: pressure tensor component-XZ, defined on nodes */
    double ****pXZsn;
    /*! SPECIES: pressure tensor component-XZ, defined on nodes */
    double ****pYYsn;
    /*! SPECIES: pressure tensor component-YZ, defined on nodes */
    double ****pYZsn;
    /*! SPECIES: pressure tensor component-ZZ, defined on nodes */
    double ****pZZsn;


    /*! Field Boundary Condition 0 = Dirichlet Boundary Condition: specifies the value to take on the boundary of the domain 1 = Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain 2 = Periodic condition */

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

    /*! Boundary condition for electric field 0 = perfect conductor 1 = magnetic mirror */
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


    /*! GEM Challenge background ion */
    double *rhoINIT;
    double *rhoINJECT;
    /*! Drift of the species */
    bool *DriftSpecies;

    /*! boolean for divergence cleaning */
    bool PoissonCorrection;
    bool PoissonCorrection_RGFace;
    /*! RESTART BOOLEAN */
    int restart1;
    /*! String with the directory for the restart file */
    string RestartDirName;
    /*! Case */
    string Case;

    /*! CG tolerance criterium for stopping iterations */
    double CGtol;
    /*! GMRES tolerance criterium for stopping iterations */
    double GMREStol;

    /*! Temporal damping parameter */
    double*** Lambda;

    // OpenBC implementation

    injInfoFields *injFieldsLeft, *injFieldsRight, *injFieldsTop, *injFieldsBottom, *injFieldsFront, *injFieldsRear;

    injInfoFields* get_InfoFieldsLeft();
    injInfoFields* get_InfoFieldsTop();
    injInfoFields* get_InfoFieldsBottom();
    injInfoFields* get_InfoFieldsFront();
    injInfoFields* get_InfoFieldsRear();
    injInfoFields* get_InfoFieldsRight();

    void BoundaryConditionsB(double ***vectorX, double ***vectorY, double ***vectorZ,int nx, int ny, int nz,Grid *grid, VirtualTopology3D *vct);
    void BoundaryConditionsE(double ***vectorX, double ***vectorY, double ***vectorZ,int nx, int ny, int nz,Grid *grid, VirtualTopology3D *vct);
    void BoundaryConditionsEImage(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ,int nx, int ny, int nz, VirtualTopology3D *vct,Grid *grid);    
    void BoundaryConditionsECSIMImage(double*** imageNX, double*** imageNY, double*** imageNZ, double*** imageCX, double*** imageCY, double*** imageCZ, double*** vectNX, double*** vectNY, double*** vectNZ, double***vectCX, double *** vectCY, double***vectCZ, int nxn, int nyn, int nzn, VirtualTopology3D * vct, Grid * grid);
    void BoundaryConditionsECSIMSource(double *** vectNX, double *** vectNY, double *** vectNZ, double *** vectCX, double *** vectCY, double *** vectCZ, VirtualTopology3D *vct, Grid *grid);


    /*! mlmd specific variables */
    /*! MLMDVerbose: when true, MLMD-related output */
    bool MLMDVerbose;
    /*! grid number in the mlmd grid hierarchy */
    int numGrid;

    /* wether to perform these operations */
    bool MLMD_BC;
    bool MLMD_PROJECTION;
    bool ParticleREPOPULATION;
    bool MLMD_BCBufferArea;
    // to test the new CG-RG mixing at RG level
    bool MLMD_fixBCenters;
    /* interpolates B cells BEFORE solving for E n+th -
       predates on MLMD_fixBCenters because I'm lazy */
    bool MLMD_InterpolateOldBCell;

    /* number of fields I am sending as BC: Ex, Ey, Ez, Exth, Eyth, Ezth, Bxn, Byn, Bzn */
    int NumF;
    int Numfix3B;

    /* coordinates of the origin on the PARENT grid */
    //    double Ox, Oy, Oz;

    /*! number of children in the mlmd hierarchy */ 
    int numChildren;

    /* Number of BC messages for the RG core -
       size of the local RGBC_struct; valid only for RGs;
       set in initWeightBC 
       vectors repeated for ghost nodes and for the one after that, 
       with different names
    */

    // RG_MaxMsgSize is used as a max size for receivinc BC msg
    // value varies core per core and decided at initWeightBC
    int RG_MaxMsgSize;
    RGBC_struct * RGBC_Info_Ghost;
    int RG_numBCMessages_Ghost;

    RGBC_struct * RGBC_Info_Active;
    int RG_numBCMessages_Active;

    /* in case MLMD_BCBufferArea= true */
    int RG_MaxMsgBufferSize;

    /* the length over which to do the buffering */
    int BufX;
    int BufY;
    int BufZ;

    // direction of the BC Buffer - (L)eft / (R)ight - (T)op / (B)ottom - (F)ront / (b)ack
    char * DirBuffer;

    RGBC_struct * RGBC_Info_Buffer;
    int RG_numBCMessages_Buffer;
    /* end in case MLMD_BCBufferArea= true */

    /* fixBghost */
    int RG_numBCMessages_fix3B;
    int RG_Maxfix3BMsgSize;
    RGBC_struct * RGBC_Info_fix3B;

    int CG_MaxSizeFix3BMsg;

    RGBC_struct ** CG_Info_Fix3B;
    int *CG_numBCMessages_Fix3B;
    double * CGMsgFix3B;

    /* end fixBghost */

    // rows: [0 - RG_numBCMessages_Active]
    // columns: [0 - RG_MaxMsgSize]

    double **Ex_Active_BC;
    double **Ey_Active_BC;
    double **Ez_Active_BC;

    double **Exth_Active_BC;
    double **Eyth_Active_BC;
    double **Ezth_Active_BC;

    double **Bxn_Active_BC;
    double **Byn_Active_BC;
    double **Bzn_Active_BC;

    // rows: [0 - RG_numBCMessages_Ghost]
    // columns: [0 - RG_MaxMsgSize]
    double **Ex_Ghost_BC;
    double **Ey_Ghost_BC;
    double **Ez_Ghost_BC;

    double **Exth_Ghost_BC;
    double **Eyth_Ghost_BC;
    double **Ezth_Ghost_BC;
    
    double **Bxn_Ghost_BC;
    double **Byn_Ghost_BC;
    double **Bzn_Ghost_BC;

    /* in case MLMD_BCBufferArea= true */
    // [RG_numBCMessages_Buffer][RG_MaxMsgBufferSize]
    double **Ex_Buffer_BC;
    double **Ey_Buffer_BC;
    double **Ez_Buffer_BC;

    double **Exth_Buffer_BC;
    double **Eyth_Buffer_BC;
    double **Ezth_Buffer_BC;
    
    double **Bxn_Buffer_BC;
    double **Byn_Buffer_BC;
    double **Bzn_Buffer_BC;
    /* end in case MLMD_BCBufferArea= true */
    
    /* in case MLMD_fixBCenters= 1 **/
    //[RG_numBCMessages_fix3B][RG_Maxfix3BMsgSize]
    double ** Bxc_fix3B_BC;
    double ** Byc_fix3B_BC;
    double ** Bzc_fix3B_BC;
    /* end in case MLMD_fixBCenters= 1 **/

    // this for ghost and active [RG_MaxMsgSize *NumF]
    double * RGMsg;
    // this for buffer [RG_MaxMsgBufferSize *NumF]
    double * RGMsgBuffer;
    // this for fix3B [RG_Maxfix3BMsgSize *Num3B]
    double * RGMsgfix3B;

    /* Number of BC messages for the CG core - 
       the rows are the children, in the same order as the communicators in the
       CommToChildren communicator in vct;
       in columns, the number of messages to the cores of that grid
       each elements is a RGBC_struct -
       int *CG_numBCMessages_*** gives the number of message for each child grid -
       
       structure repeated for the ghost nodes and for the active ones just following 
       CG_MaxSizeMsg gives you a max size for the msg to send; Size changes core by core and decided during initWeight
    */

    int CG_MaxSizeMsg;
    //[CG_MaxSizeMsg * NumF]
    double *CGMsg; // used to build msgs

    RGBC_struct ** CG_Info_Ghost;
    int *CG_numBCMessages_Ghost;

    //[numChildren][MAX_RG_numBCMessages], then resized to [numChildren][max(CG_numBCMessages_Ghost[ch)+1]
    RGBC_struct ** CG_Info_Active;
    int *CG_numBCMessages_Active;

    /* in case MLMD_BCBufferArea= true */
    //[CG_MaxSizeBufferMsg * NumF] 
    double *CGMsgBuffer; // used to build msgs

    int CG_MaxSizeBufferMsg;
    RGBC_struct ** CG_Info_Buffer;
    int *CG_numBCMessages_Buffer;
    /* end in case MLMD_BCBufferArea= true */

    /* tags for send/receive of BCs*/
    int TAG_BC_GHOST;
    int TAG_BC_ACTIVE;
    int TAG_II;
    int TAG_BC_BUFFER;
    int TAG_BC_FIX3B;

    /* max vector dimensions */
    int MAX_RG_numBCMessages;
    int MAX_size_LevelWide;

    // projection

    /* tags for send/receive of proj*/
    int TAG_PROJ;

    // RG side
    int RG_numProjMessages;
    // [MAX_RG_numBCMessages], later trimmed
    RGBC_struct * RGProj_Info;
    
    // largest numbers of point that RG has to send as proj, set in initWeightProj_Phase1
    int size_RG_ProjMsg;
    // each projection msg RG sends; [(size_RG_ProjMsg+1)*3]
    double * RG_ProjMsg;
    
    // CG side
    // [numChildren]
    int *CG_numProjMessages;
    
    // [numChildren][MAX_RG_numBCMessages], later trimmed
    RGBC_struct ** CGProj_Info;

    // largest numbers of point that CG receives as proj, set in initWeightProj
    int size_CG_ProjMsg;
    // each projection msg CG receives; [(size_CG_ProjMsg+1)*3] 
    double * CG_ProjMsg;
    // per core, this is max(CG_numProjMessages); used to instantiate the weights
    int Max_CG_numProjMessages;
    
    // [numChildren][Max_CG_numProjMessages][size_CG_ProjMsg]
    double ***ProjWeight000;
    double ***ProjWeight001;
    double ***ProjWeight010;
    double ***ProjWeight011;
    double ***ProjWeight100;
    double ***ProjWeight101;
    double ***ProjWeight110;
    double ***ProjWeight111;
    double ***ProjIX;
    double ***ProjIY;
    double ***ProjIZ;
    // end projection

    // the vectors where CG accumulates info on projection before applying it to fields
    // instantiated only if at least one CG_numProjMessages >0
    
    // [numChildren][nxn][nzn][nzn]
    double ****ExthProjSt;
    double ****EythProjSt;
    double ****EzthProjSt;

    double ****DenProjSt;

    /* if you, as a CG core, need to apply projection 
       (at any level; just checked to see wether ApplyProjection can be skipped) */
    bool ApplyProjection; 

    // instantiated only if ApplyProjection= true, used to regenerated En+1 after receiving projection
    double ***Ex_n;
    double ***Ey_n;
    double ***Ez_n;

    // instantiated only if CG_numBCMessages_Ghost[any] >0 or CG_numBCMessages_Active[any]>0 --> BSNeeded= true
    bool BSNeeded;
    double ***Ex_BS;
    double ***Ey_BS;
    double ***Ez_BS;

    // stuff for initial interpolation
    // CG side
    int CG_MaxSizeMsg_II;

    //[numChildren][MAX_RG_numBCMessages], then resized to [numChildren][max(CG_numBCMessages_Ghost[ch)+1] 
    RGBC_struct ** CG_Info_II;
    int *CG_numBCMessages_II;
    double *CGMsg_II;
    // RG side
    RGBC_struct * RGBC_Info_II;
    int RG_numBCMessages_II;
    int RG_MaxMsgSize_II;


    double **Ex_II;
    double **Ey_II;
    double **Ez_II;

    double **Bxn_II;
    double **Byn_II;
    double **Bzn_II;

    double **rho0_II;
    double **rho1_II;
    double **rho2_II;
    double **rho3_II;

    double *RGMsg_II;

    int NumF_II;

    /* MPI Datatype associated to RGBC_struct; init in MPI_RGBC_struct_commit  */
    MPI_Datatype MPI_RGBC_struct;
    /*! end mlmd specidic variables */

    // instantiated only if SmoothRGBC= true and 
    double *** Ex_TH_SRGBC;
    double *** Ey_TH_SRGBC;
    double *** Ez_TH_SRGBC;

    double *** Ex_SRGBC;
    double *** Ey_SRGBC;
    double *** Ez_SRGBC;

    // at the moment, this is used only to test the fluid repopulation method
    double **** RHOINIT;

    bool NewSmoothing;
    bool RepopulateBeforeMover;

    /* used in function related to divB cleaning for the RG */
    int iStart_BP; int iEnd_BP; int lenX_BP;
    int jStart_BP; int jEnd_BP; int lenY_BP;
    int kStart_BP; int kEnd_BP; int lenZ_BP;

    double sumEzBC;
    
    bool SmoothFaces;
};


typedef EMfields3D Field;

#endif
