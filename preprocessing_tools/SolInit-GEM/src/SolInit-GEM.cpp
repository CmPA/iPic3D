#include "H5hut-io.h"
#include "Fields.h"

#include <sstream>
#include <string>
#include <cmath>
#include <mpi.h>

#define FourPI 12.56637

double getNodeR(int r, double dr, int coorr, int Lpr) {return (coorr*Lpr + (r-1)*dr);}
double getCellR(int r, double dr, int coorr, int Lpr) {return (coorr*Lpr + (r-1)*dr + 0.5*dr);}

int main (int argc, char *argv[])
{

  /* ---------------------- */
  /* Simulation parametters */
  /* ---------------------- */

  int       nspec    = 2;      /*!< Number of species */
  int       ntcx     = 64;     /*!< Number of cells in X */
  int       ntcy     = 64;     /*!< Number of cells in Y */
  int       ntcz     = 64;      /*!< Number of cells in Z */
  int       nppc     = 125;    /*!< Number of particles per cell */
  double    Lx       = 10.0;   /*!< Domain size in X */
  double    Ly       = 10.0;   /*!< Domain size in Y */
  double    Lz       = 10.0;   /*!< Domain size in Z */
  double    B0x      = 0.0097; /*!< Background B field in X */
  double    B0y      = 0.0;    /*!< Background B field in Y */
  double    B0z      = 0.0;    /*!< Background B field in Z */

  /* --------- */
  /* Variables */
  /* --------- */

  double    dx;
  double    dy;
  double    dz;
  int       rank, size;
  int       *perio;
  int       *dimns;
  int       *coord;
  int       MPI_ndim = 2;
  int       cycle    = 0;
  int       reorder  = 1;
  int       ghostc_x = 2;
  int       ghostc_y = 2;
  int       ghostc_z = 2;
  int       ntxn;
  int       ntyn;
  int       ntzn;
  double    *L;
  Fields    EMf;
  H5hutpart *part;
  H5output  file;
  MPI_Comm  CART_COMM;

  /* -------------------------- */
  /* Initialize MPI environment */
  /* -------------------------- */

  perio = new int[MPI_ndim];
  dimns = new int[MPI_ndim];
  coord = new int[MPI_ndim];

  std::cout << " Call MPI..." << std::endl;

  MPI_Init (&argc, &argv);
  MPI_Comm_size  (MPI_COMM_WORLD, &size);
  if (MPI_ndim==3) MPI_Dims_create(size, 3, dimns);
  else       { MPI_Dims_create(size, 2, dimns); dimns[2]=1;}
  MPI_Cart_create(MPI_COMM_WORLD, MPI_ndim, dimns, perio, reorder, &CART_COMM);
  MPI_Comm_rank  (CART_COMM,&rank);
  MPI_Cart_coords(CART_COMM, rank, MPI_ndim, coord);

  /* ---------------------------------------------- */
  /* Verify that the domain has a good partitioning */
  /* ---------------------------------------------- */

  if      (ntcx%dimns[0]!=0) { std::cout << " ERROR: Wrong partitioning" << std::endl; abort(); }
  else if (ntcy%dimns[1]!=0) { std::cout << " ERROR: Wrong partitioning" << std::endl; abort(); }
  else if (MPI_ndim==3) {
       if (ntcz%dimns[2]!=0) { std::cout << " ERROR: Wrong partitioning" << std::endl; abort(); }
  }

  /* ------------------------- */
  /* Set simulation conditions */
  /* ------------------------- */

  ntxn = ntcx + 1;
  ntyn = ntcy + 1;
  ntzn = ntcz + 1;

  int nxc  = (ntcx / dimns[0]) + ghostc_x;
  int nyc  = (ntcy / dimns[1]) + ghostc_y;
  int nzc  = (ntcz / dimns[2]) + ghostc_z;
  int nxn  = nxc + 1;
  int nyn  = nyc + 1;
  int nzn  = nzc + 1;

  dx = Lx / ntcx;
  dy = Ly / ntcy;
  dz = Lz / ntcz;

  /* -------------------------------------- */
  /* Allocate fields in the cells and nodes */
  /* -------------------------------------- */

  std::cout << " Init fields..." << std::endl;

  EMf.Init(nspec, nxc, nyc, nzc);

  double rhoINIT = 1.0;
  double rho     = 0.0;
  double delta   = 0.5;

  /* NODE DATA */
  for(int i=0; i<nxn; i++) {
    for(int j=0; j<nyn; j++) {
      for(int k=0; k<nzn; k++) {

        // Electric field
        EMf.SetEn(i, j, k, 0.0, 0.0, 0.0);

        // Magnetic field
        double Bx;
        double By;
        double Bz;

        double rx = getNodeR(i, dx, coord[0], Lx/dimns[0]);
        double ry = getNodeR(j, dy, coord[1], Ly/dimns[1]);

        Bx = B0x * tanh((ry - Ly/2) / delta);
        By = B0y;
        Bz = B0z;

        // Perturbation:
        double pertX    = 0.4;
        double xpert    = rx - Lx/2;
        double ypert    = ry - Ly/2;
        double exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));

        Bx += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * 
                                           cos(M_PI * ypert / 10.0 / delta) * 
                                           2.0 * ypert / delta 
                                         - cos(M_PI * xpert / 10.0 / delta) * 
                                           sin(M_PI * ypert / 10.0 / delta) * 
                                           M_PI / 10.0);
        By += (B0x * pertX) * exp_pert * ( cos(M_PI * xpert / 10.0 / delta) * 
                                           cos(M_PI * ypert / 10.0 / delta) *
                                           2.0 * xpert / delta
                                         + sin(M_PI * xpert / 10.0 / delta) *
                                           cos(M_PI * ypert / 10.0 / delta) *
                                           M_PI / 10.0);

        EMf.SetBn(i, j, k, Bx, By, Bz);

      }
    }
  }

  for (int s=0; s<nspec; s++) {
    for(int i=0; i<nxn; i++) {
      for(int j=0; j<nyn; j++) {
        for(int k=0; k<nzn; k++) {
          double rx = getNodeR(i, dx, coord[0], Lx/dimns[0]);
          double ry = getNodeR(j, dy, coord[1], Ly/dimns[1]);
          rho = ((rhoINIT / ( cosh((ry - Ly/2) / delta) * 
                              cosh((ry - Ly/2) / delta)))) / FourPI;
          EMf.SetRhon(s, i, j, k, rho);
        }
      }
    }
  }

  /* CELL DATA */
  for(int i=0; i<nxc; i++) {
    for(int j=0; j<nyc; j++) {
      for(int k=0; k<nzc; k++) {

        // Electric field
        EMf.SetEc(i, j, k, 0.0, 0.0, 0.0);

        // Magnetic field
        double Bx;
        double By;
        double Bz;

        double rx = getCellR(i, dx, coord[0], Lx/dimns[0]);
        double ry = getCellR(j, dy, coord[1], Ly/dimns[1]);

        Bx = B0x * tanh((ry - Ly/2) / delta);
        By = B0y;
        Bz = B0z;

        // Perturbation:
        double pertX    = 0.4;
        double xpert    = rx - Lx/2;
        double ypert    = ry - Ly/2;
        double exp_pert = exp(-(xpert / delta) * (xpert / delta) - (ypert / delta) * (ypert / delta));

        Bx += (B0x * pertX) * exp_pert * (-cos(M_PI * xpert / 10.0 / delta) * 
                                           cos(M_PI * ypert / 10.0 / delta) * 
                                           2.0 * ypert / delta 
                                         - cos(M_PI * xpert / 10.0 / delta) * 
                                           sin(M_PI * ypert / 10.0 / delta) * 
                                           M_PI / 10.0);
        By += (B0x * pertX) * exp_pert * ( cos(M_PI * xpert / 10.0 / delta) * 
                                           cos(M_PI * ypert / 10.0 / delta) *
                                           2.0 * xpert / delta
                                         + sin(M_PI * xpert / 10.0 / delta) *
                                           cos(M_PI * ypert / 10.0 / delta) *
                                           M_PI / 10.0);

        EMf.SetBc(i, j, k, Bx, By, Bz);

      }
    }
  }

  for (int s=0; s<nspec; s++) {
    for(int i=0; i<nxc; i++) {
      for(int j=0; j<nyc; j++) {
        for(int k=0; k<nzc; k++) {
          double rx = getCellR(i, dx, coord[0], Lx/dimns[0]);
          double ry = getCellR(j, dy, coord[1], Ly/dimns[1]);
          rho = ((rhoINIT / ( cosh((ry - Ly/2) / delta) * 
                              cosh((ry - Ly/2) / delta)))) / FourPI;
          EMf.SetRhoc(s, i, j, k, rho);
        }
      }
    }
  }

  std::cout << " fields initialized!" << std::endl;

  /* ------------------ */
  /* Allocate particles */
  /* ------------------ */

  std::cout << " Set particles..." << std::endl;

  part = new H5hutpart[nspec];

  srand(rank*10);
  for (int i=0; i < nspec; i++){
    double chg = -pow((-1),i);
    part[i].init(nxc*nyc*nzc*nppc, chg, Lx, Ly, Lz);
  }

  std::cout << " Particles initialized!" << std::endl;

  /* ------------------------------------------- */
  /* Open/Write/Close fields file using H5hut-io */
  /* ------------------------------------------- */

  std::cout << " Setting file name and cycle..." << std::endl;
  file.SetNameCycle("GEMtest", cycle);

  std::cout << " Open file for fields..." << std::endl;
  std::cout << " Total Nodes: " << ntxn << " " << ntyn << " " << ntzn << std::endl;
  file.OpenFieldsFile("Node", nspec, ntxn, ntyn, ntzn, coord, dimns, CART_COMM);

  file.WriteFields(EMf.getEx(), "Ex", nxn, nyn, nzn);
  file.WriteFields(EMf.getEy(), "Ey", nxn, nyn, nzn);
  file.WriteFields(EMf.getEz(), "Ez", nxn, nyn, nzn);
  file.WriteFields(EMf.getBx(), "Bx", nxn, nyn, nzn);
  file.WriteFields(EMf.getBy(), "By", nxn, nyn, nzn);
  file.WriteFields(EMf.getBz(), "Bz", nxn, nyn, nzn);

  for (int is=0; is<nspec; is++) {
    std::stringstream  ss;
    ss << is;
    std::string s_is = ss.str();
    file.WriteFields(EMf.getRHOns(is), "rho_"+ s_is, nxn, nyn, nzn);
  }
  file.CloseFieldsFile();

//  file.OpenFieldsFile("Cell", nspec, ntcx, ntcy, ntcz, coord, dimns, CART_COMM);
//
//  file.WriteFields(EMf.getExc(), "Exc", nxc, nyc, nzc);
//  file.WriteFields(EMf.getEyc(), "Eyc", nxc, nyc, nzc);
//  file.WriteFields(EMf.getEzc(), "Ezc", nxc, nyc, nzc);
//  file.WriteFields(EMf.getBxc(), "Bxc", nxc, nyc, nzc);
//  file.WriteFields(EMf.getByc(), "Byc", nxc, nyc, nzc);
//  file.WriteFields(EMf.getBzc(), "Bzc", nxc, nyc, nzc);
//
//  for (int is=0; is<nspec; is++) {
//    std::stringstream  ss;
//    ss << is;
//    std::string s_is = ss.str();
//    file.WriteFields(EMf.getRHOcs(is), "rhoc_"+ s_is, nxc, nyc, nzc);
//  }
//  file.CloseFieldsFile();

  std::cout << " file closed!" << std::endl;

  /* ---------------------------------------------- */
  /* Open/Write/Close particles file using H5hut-io */
  /* ---------------------------------------------- */

  std::cout << " Open file for particles..." << std::endl;
  file.OpenPartclFile(nspec, CART_COMM);

  std::cout << " Write particles for each species..." << std::endl;
  for (int i=0; i<nspec; i++) {
    file.WriteParticles(i, part[i].getNP(),
                           part[i].getQ(),
                           part[i].getX(),
                           part[i].getY(),
                           part[i].getZ(),
                           part[i].getU(),
                           part[i].getV(),
                           part[i].getW(),
                           CART_COMM);
  }

  std::cout << " Close particles file..." << std::endl;
  file.ClosePartclFile();

  std::cout << " file closed!" << std::endl;

  /* ------------------ */
  /* Deallocate vectors */
  /* ------------------ */

  delete [] perio;
  delete [] dimns;
  delete [] coord;
  delete [] part;

  MPI_Finalize();

  return 0;
}
