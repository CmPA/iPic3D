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
#include "CollectiveIO.h"
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
    EMfields3D(CollectiveIO * col, Grid * grid);
    /*! destructor */
    ~EMfields3D();

    /*! initialize the electromagnetic fields with constant values */
    void init(VirtualTopology3D * vct, Grid * grid);
    /*! init beam */
    void initBEAM(VirtualTopology3D * vct, Grid * grid, double x_center, double y_center, double z_center, double radius);
    /*! initialize GEM challenge */
    void initGEM(VirtualTopology3D * vct, Grid * grid);
    void initOriginalGEM(VirtualTopology3D * vct, Grid * grid);
    void initDoublePeriodicHarrisWithGaussianHumpPerturbation(VirtualTopology3D * vct, Grid * grid);
    /*! initialize GEM challenge with dipole-like tail without perturbation */
    void initGEMDipoleLikeTailNoPert(VirtualTopology3D * vct, Grid * grid);
    /*! initialize GEM challenge with no Perturbation */
    void initGEMnoPert(VirtualTopology3D * vct, Grid * grid);
    /** Random initial field */
    void initRandomField(VirtualTopology3D * vct, Grid * grid);
    /*! Init Force Free (JxB=0) */
    void initForceFree(VirtualTopology3D * vct, Grid * grid);
    /*! initialized with rotated magnetic field */
    void initEM_rotate(VirtualTopology3D * vct, Grid * grid, double B, double theta);
    /*! add a perturbattion to charge density */
    void AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid * grid);
    /*! add a perturbattion to the EM field */
    void AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid * grid);


    /*! Calculate Electric field using the implicit Maxwell solver */
    void calculateE(Grid * grid, VirtualTopology3D * vct);
    /*! Image of Poisson Solver (for SOLVER) */
    void PoissonImage(double *image, double *vector, Grid * grid, VirtualTopology3D * vct);
    /*! Image of Maxwell Solver (for Solver) */
    void MaxwellImage(double *im, double *vector, Grid * grid, VirtualTopology3D * vct);
    /*! Maxwell source term (for SOLVER) */
    void MaxwellSource(double *bkrylov, Grid * grid, VirtualTopology3D * vct);
    /*! Impose a constant charge inside a spherical zone of the domain */
    void ConstantChargePlanet(Grid * grid, VirtualTopology3D * vct, double R, double x_center, double y_center, double z_center);
    /*! Impose a constant charge in the OpenBC boundaries */
    void ConstantChargeOpenBC(Grid * grid, VirtualTopology3D * vct);
    /*! Calculate Magnetic field with the implicit solver: calculate B defined on nodes With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
    void calculateB(Grid * grid, VirtualTopology3D * vct);
    /*! fix B on the boundary for gem challange */
    void fixBgem(Grid * grid, VirtualTopology3D * vct);
    /*! fix B on the boundary for gem challange */
    void fixBforcefree(Grid * grid, VirtualTopology3D * vct);

    /*! Calculate the three components of Pi(implicit pressure) cross image vector */
    void PIdot(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid * grid);
    /*! Calculate the three components of mu (implicit permeattivity) cross image vector */
    void MUdot(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid * grid);
    /*! Calculate rho hat, Jx hat, Jy hat, Jz hat */
    void calculateHatFunctions(Grid * grid, VirtualTopology3D * vct);


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
    void smoothE(double value, VirtualTopology3D * vct);

    /*! communicate ghost for grid -> Particles interpolation */
    void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology3D * vct);
    /*! add accumulated moments to the moments for a given species */
    void addToSpeciesMoments(const Moments & in, int is);
    /*! add an amount of charge density to charge density field at node X,Y,Z */
    void addRho(double weight[][2][2], int X, int Y, int Z, int is);
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
    void adjustNonPeriodicDensities(int is, VirtualTopology3D * vct);


    /*! Perfect conductor boundary conditions LEFT wall */
    void perfectConductorLeft(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid);
    /*! Perfect conductor boundary conditions RIGHT wall */
    void perfectConductorRight(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir, Grid * grid);
    /*! Perfect conductor boundary conditions for source LEFT wall */
    void perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir);
    /*! Perfect conductor boundary conditions for source RIGHT wall */
    void perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ, int dir);

    /*! Calculate the sysceptibility tensor on the boundary */
    void sustensorRightX(double **susxy, double **susyy, double **suszy);
    void sustensorLeftX (double **susxy, double **susyy, double **suszy);
    void sustensorRightY(double **susxy, double **susyy, double **suszy);
    void sustensorLeftY (double **susxy, double **susyy, double **suszy);
    void sustensorRightZ(double **susxy, double **susyy, double **suszy);
    void sustensorLeftZ (double **susxy, double **susyy, double **suszy);

    /*! get Potential array */
    double ***getPHI();
    /*! get Electric Field component X defined on node(indexX,indexY,indexZ) */
    double &getEx(int indexX, int indexY, int indexZ) const;
    /*! get Electric field X component array */
    double ***getEx();
    /*! get Electric field X component cell array without the ghost cells */
    double ***getExc(Grid3DCU *grid);
    /*! get Electric Field component Y defined on node(indexX,indexY,indexZ) */
    double &getEy(int indexX, int indexY, int indexZ) const;
    /*! get Electric field Y component array */
    double ***getEy();
    /*! get Electric field Y component cell array without the ghost cells */
    double ***getEyc(Grid3DCU *grid);
    /*! get Electric Field component Z defined on node(indexX,indexY,indexZ) */
    double &getEz(int indexX, int indexY, int indexZ) const;
    /*! get Electric field Z component array */
    double ***getEz();
    /*! get Electric field Z component cell array without the ghost cells */
    double ***getEzc(Grid3DCU *grid);
    /*! get Magnetic Field component X defined on node(indexX,indexY,indexZ) */
    double &getBx(int indexX, int indexY, int indexZ) const;
    /*! get Magnetic field X component array */
    double ***getBx();
    /*! get Magnetic field X component cell array without the ghost cells */
    double ***getBxc();
    /*! get Magnetic Field component Y defined on node(indexX,indexY,indexZ) */
    double &getBy(int indexX, int indexY, int indexZ) const;
    /*! get Magnetic field Y component array */
    double ***getBy();
    /*! get Magnetic field Y component cell array without the ghost cells */
    double ***getByc();
    /*! get Magnetic Field component Z defined on node(indexX,indexY,indexZ) */
    double &getBz(int indexX, int indexY, int indexZ) const;
    /*! get Magnetic field Z component array */
    double ***getBz();
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
    /*! SPECIES: get density array on nodes */
    double ****getRHOns();
    /*! SPECIES: get density array on cells without the ghost cells */
    double ***getRHOcs(int is);



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
    double ***getJxsc(Grid3DCU *grid, int is);
    /*! get Jy(X,Y,Z) */
    double &getJy(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction Y */
    double ***getJy();
    /*! get Jys(X,Y,Z,is) */
    double &getJys(int indexX, int indexY, int indexZ, int is) const;
    /*! SPECIES: get current -Direction Y */
    double ****getJys();
    /*! SPECIES: get current Y component for species is in all cells except ghost */
    double ***getJysc(Grid3DCU *grid, int is);
    /*! get Jz(X,Y,Z) */
    double &getJz(int indexX, int indexY, int indexZ) const;
    /*! get current -Direction Z */
    double ***getJz();
    /*! get Jzs(X,Y,Z,is) */
    double &getJzs(int indexX, int indexY, int indexZ, int is) const;
    /*! SPECIES: get current -Direction Z */
    double ****getJzs();
    /*! SPECIES: get current Z component for species is in all cells except ghost */
    double ***getJzsc(Grid3DCU *grid, int is);
    /*! get the electric field energy */
    double getEenergy();
    /*! get the magnetic field energy */
    double getBenergy();


    /*! print electromagnetic fields info */
    void print(void) const;

    // OpenBC
    void updateInfoFields(Grid *grid,VirtualTopology3D *vct,CollectiveIO *col);

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
    /*! delt = c*th*dt */
    double delt;
    /*! number of particles species */
    int ns;
    /*! GEM challenge parameters */
    double B0x, B0y, B0z, delta;
    /*! charge to mass ratio array for different species */
    double *qom;
    /*! Boundary electron speed */
    double ue0, ve0, we0;


    // KEEP IN MEMORY GUARD CELLS ARE INCLUDED
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
    double ***Dx;
    double ***Dy;
    double ***Dz;
    double ***vectX;
    double ***vectY;
    double ***vectZ;
    double ***divC;


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
    /*! Drift of the species */
    bool *DriftSpecies;

    /*! boolean for divergence cleaning */
    bool PoissonCorrection;
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

};

typedef EMfields3D Field;

#endif
