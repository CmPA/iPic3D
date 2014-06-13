#include "H5hut-io.h"
#include "Fields.h"

#include <sstream>
#include <string>
#include <mpi.h>

int main (int argc, char *argv[])
{
  int       rank, size;
  int       *perio;
  int       *dimns;
  int       *coord;
  int       ndim     = 2;
  int       cycle    = 0;
  int       reorder  = 1;
  int       nspec    = 2;
  int       ntcx     = 64;
  int       ntcy     = 64;
  int       ntcz     = 1;
  int       ntxn;
  int       ntyn;
  int       ntzn;
  int       nppc     = 27;
  double    Lx       = 10.0;
  double    Ly       = 10.0;
  double    Lz       = 1.0;
  double    *L;
  Fields    EMf;
  H5hutpart *part;
  H5output  file;
  MPI_Comm  CART_COMM;

  /* -------------------------- */
  /* Initialize MPI environment */
  /* -------------------------- */

  perio = new int[ndim];
  dimns = new int[ndim];
  coord = new int[ndim];

  std::cout << " Call MPI..." << std::endl;

  MPI_Init (&argc, &argv);
  MPI_Comm_size  (MPI_COMM_WORLD, &size);
  if (ndim==3) MPI_Dims_create(size, 3, dimns);
  else       { MPI_Dims_create(size, 2, dimns); dimns[2]=1;}
  MPI_Cart_create(MPI_COMM_WORLD, ndim, dimns, perio, reorder, &CART_COMM);
  MPI_Comm_rank  (CART_COMM,&rank);
  MPI_Cart_coords(CART_COMM, rank, ndim, coord);

  /* ------------------------- */
  /* Set simulation conditions */
  /* ------------------------- */

  ntxn = ntcx + 1;
  ntyn = ntcy + 1;
  ntzn = ntcz + 1;

  int nxc  = (ntcx / dimns[0]) + 2;
  int nyc  = (ntcy / dimns[1]) + 2;
  int nzc  = (ntcz / dimns[2]) + 2;
  int nxn  = nxc + 1;
  int nyn  = nyc + 1;
  int nzn  = nzc + 1;

  /* ---------------------------------------------- */
  /* Verify that the domain has a good partitioning */
  /* ---------------------------------------------- */

  if      (nxc%dimns[0]!=0) { std::cout << " ERROR: Wrong partitioning" << std::endl; abort(); }
  else if (nyc%dimns[1]!=0) { std::cout << " ERROR: Wrong partitioning" << std::endl; abort(); }
  else if (ndim==3) {
       if (nzc%dimns[2]!=0) { std::cout << " ERROR: Wrong partitioning" << std::endl; abort(); }
  }

  /* -------------------------------------- */
  /* Allocate fields in the cells and nodes */
  /* -------------------------------------- */

  std::cout << " Init fields..." << std::endl;

  EMf.Init(nspec, nxc, nyc, nzc);

  for(int i=0; i<nxn; i++) {
    for(int j=0; j<nyn; j++) {
      for(int k=0; k<nzn; k++) {
        EMf.SetBn(i, j, k, 0.0, 0.0, 0.1);
        EMf.SetEn(i, j, k, 0.0, 0.0, 0.0);
      }
    }
  }

  for (int s=0; s<nspec; s++) {
    for(int i=0; i<nxn; i++) {
      for(int j=0; j<nyn; j++) {
        for(int k=0; k<nzn; k++) {
          double chg = -1;
          if (s!=0) chg = 1;
          EMf.SetRhon(s, i, j, k, chg);
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

  /* ----------------------------------- */
  /* -- Performing HDF5 fields output -- */

  std::cout << " Setting file name and cycle..." << std::endl;
  file.SetNameCycle("Test", cycle);

//--- SAVE IN THE CELLS:
//  file.OpenFieldsFile("Cell", nspec, ntcx, ntcy, ntcz, coord, dimns, CART_COMM);
//  file.WriteFields(EMf.getExc(), "Ex", nxc, nyc, nzc);
//  file.WriteFields(EMf.getEyc(), "Ey", nxc, nyc, nzc);
//  file.WriteFields(EMf.getEzc(), "Ez", nxc, nyc, nzc);
//  file.WriteFields(EMf.getBxc(), "Bx", nxc, nyc, nzc);
//  file.WriteFields(EMf.getByc(), "By", nxc, nyc, nzc);
//  file.WriteFields(EMf.getBzc(), "Bz", nxc, nyc, nzc);
//  for (int is=0; is<nspec; is++) {
//    std::stringstream  ss;
//    ss << is;
//    std::string s_is = ss.str();
//    file.WriteFields(EMf.getRHOcs(is), "rho_"+ s_is, nxc, nyc, nzc);
//  }
//  file.CloseFieldsFile();
//--- SAVE IN THE NODES:
  std::cout << " Open file for fields..." << std::endl;
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
  std::cout << " file closed!" << std::endl;
//--- END SAVE FIELDS

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
  /* ----------------------------------- */

  delete [] perio;
  delete [] dimns;
  delete [] coord;
  delete [] part;
  MPI_Finalize();
  return 0;
}
