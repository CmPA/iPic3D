/* ---------------------------------------------------------- */
//  Author     : Jorge AMAYA
//  Date       : January 2013
//  Affiliation: Center for mathematical Plasma-Astrophysics
//               Mathematics Department
//               KU Leuven
/* ---------------------------------------------------------- */

#include <mpi.h>

#include "H5hut-io.h"
#include "Fields.h" 

int main (int argc, char *argv[])
{
  int       rank, size;
  int       *perio;
  int       *dimns;
  int       *coord;
  int       cycle0   = 0;
  int       tsave    = 100;
  int       ncycle   = 200;
  int       ndim     = 3;
  int       reorder  = 1;
  int       nspec    = 2;
  int       ntcx     = 64;
  int       ntcy     = 64;
  int       ntcz     = 64;
  int       ntxn;
  int       ntyn;
  int       ntzn;
  double    Lx       = 10.0;
  double    Ly       = 10.0;
  double    Lz       = 10.0;
  double    L[3];
  MPI_Comm  CART_COMM;
  Fields    EMf;
  H5input   file;
  H5output  outfile;
  H5hutpart *myPart;

  /* -------------------------- */
  /* Initialize MPI environment */
  /* -------------------------- */

  perio = new int[ndim];
  dimns = new int[ndim];
  coord = new int[ndim];

  MPI_Init (&argc, &argv);
  MPI_Comm_size  (MPI_COMM_WORLD, &size);
  if (ndim==3) MPI_Dims_create(size, 3, dimns);
  else         MPI_Dims_create(size, 2, dimns);
  MPI_Cart_create(MPI_COMM_WORLD, ndim, dimns, perio, reorder, &CART_COMM);
  MPI_Comm_rank  (CART_COMM,&rank);
  MPI_Cart_coords(CART_COMM, rank, ndim, coord);

  /* ------------------------- */
  /* Set simulation conditions */
  /* ------------------------- */

  ntxn = ntcx + 1;
  ntyn = ntcy + 1;
  ntzn = ntcz + 1;

  int nxc = (ntcx / dimns[0]) + 2;
  int nyc = (ntcy / dimns[1]) + 2;
  int nzc = (ntcz / dimns[2]) + 2;
  int nxn = nxc + 1;
  int nyn = nyc + 1;
  int nzn = nzc + 1;

  L[0] = Lx;
  L[1] = Ly;
  L[2] = Lz;

  /* ------------------ */
  /* Read the particles */
  /* ------------------ */

  file.SetNameCycle("GEMtest", cycle0);

  myPart = new H5hutpart[nspec];

  file.OpenPartclFile(nspec);
  file.ReadParticles(rank, size, dimns, L, CART_COMM);

  for (int s = 0; s < nspec; s++){
    myPart[s].memalloc(file.GetNp(s));

    file.DumpPartclX(myPart[s].getXref(), s);
    file.DumpPartclY(myPart[s].getYref(), s);
    file.DumpPartclZ(myPart[s].getZref(), s);
    file.DumpPartclU(myPart[s].getUref(), s);
    file.DumpPartclV(myPart[s].getVref(), s);
    file.DumpPartclW(myPart[s].getWref(), s);
    file.DumpPartclQ(myPart[s].getQref(), s);
  }
  file.ClosePartclFile();

  /* --------------- */
  /* Read the fields */
  /* --------------- */

  EMf.Init(nspec, nxc, nyc, nzc);

  file.OpenFieldsFile("Node", nspec, ntxn, ntyn, ntzn, coord, dimns, CART_COMM);
  file.ReadFields(EMf.getEx(), "Ex", nxn, nyn, nzn);
  file.ReadFields(EMf.getEy(), "Ey", nxn, nyn, nzn);
  file.ReadFields(EMf.getEz(), "Ez", nxn, nyn, nzn);
  file.ReadFields(EMf.getBx(), "Bx", nxn, nyn, nzn);
  file.ReadFields(EMf.getBy(), "By", nxn, nyn, nzn);
  file.ReadFields(EMf.getBz(), "Bz", nxn, nyn, nzn);
  file.CloseFieldsFile();

  /* -------------------------------------- */
  /* Empty simulation with multiple outputs */
  /* -------------------------------------- */

  for (int t = cycle0; t <= cycle0+ncycle; t++) {

    std::cout << " Cycle: " << t << std::endl;
    /* here goes the computation */

    if (t%tsave==0) {

      outfile.SetNameCycle("Output", t);

      /* ----------- */
      /* Save fields */
      /* ----------- */

      outfile.OpenFieldsFile("Node", nspec, ntxn, ntyn, ntzn, coord, dimns, CART_COMM);
      outfile.WriteFields(EMf.getEx(), "Ex", nxn, nyn, nzn);
      outfile.WriteFields(EMf.getEy(), "Ey", nxn, nyn, nzn);
      outfile.WriteFields(EMf.getEz(), "Ez", nxn, nyn, nzn);
      outfile.WriteFields(EMf.getBx(), "Bx", nxn, nyn, nzn);
      outfile.WriteFields(EMf.getBy(), "By", nxn, nyn, nzn);
      outfile.WriteFields(EMf.getBz(), "Bz", nxn, nyn, nzn);
      outfile.CloseFieldsFile();

      /* -------------- */
      /* Save particles */
      /* -------------- */

      outfile.OpenPartclFile(nspec, CART_COMM);
      for (int i=0; i<nspec; i++) {
        outfile.WriteParticles(i, myPart[i].getNP(),
                                  myPart[i].getQ(),
                                  myPart[i].getX(),
                                  myPart[i].getY(),
                                  myPart[i].getZ(),
                                  myPart[i].getU(),
                                  myPart[i].getV(),
                                  myPart[i].getW(),
                                  CART_COMM);
      }
      outfile.ClosePartclFile();
    }
  }

  std::cout << rank << " : Freeing memory..." << std::endl;

  delete [] myPart;
  delete [] perio;
  delete [] dimns;
  delete [] coord;

  std::cout << rank << " : Finalizing MPI" << std::endl;

  MPI_Finalize();

  return 0;
}
