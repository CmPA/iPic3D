
#include <mpi.h>
#include "WriteOutputParallel.h"

void WriteOutputParallel(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle){

  string       grpname;
  string       dtaname;

  stringstream filenmbr;
  string       filename;

  /* ------------------- */
  /* Setup the file name */
  /* ------------------- */

  filenmbr << setfill('0') << setw(5) << cycle;
  filename = col->getSaveDirName() + "/" + col->getSimName() + "_" + filenmbr.str() + ".h5";

  /* ---------------------------------------------------------------------------- */
  /* Define the number of cells in the globa and local mesh and set the mesh size */
  /* ---------------------------------------------------------------------------- */

  int nxc = grid->getNXC();
  int nyc = grid->getNYC();
  int nzc = grid->getNZC();

  int    dglob[3] = { col ->getNxc()  , col ->getNyc()  , col ->getNzc()   };
  int    dlocl[3] = { nxc-2,            nyc-2,            nzc-2 };
  double L    [3] = { col ->getLx ()  , col ->getLy ()  , col ->getLz ()   };

  /* --------------------------------------- */
  /* Declare and open the parallel HDF5 file */
  /* --------------------------------------- */

  PHDF5fileClass outputfile(filename, 3, vct->getCoordinates(), vct->getComm());

  outputfile.CreatePHDF5file(L, dglob, dlocl, false);

  /* ------------------------ */
  /* Write the Electric field */
  /* ------------------------ */

  grpname = "Fields";
  dtaname = "Ex";
  outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getExc(grid), nxc-2, nyc-2, nzc-2);

  grpname = "Fields";
  dtaname = "Ey";
  outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getEyc(grid), nxc-2, nyc-2, nzc-2);

  grpname = "Fields";
  dtaname = "Ez";
  outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getEzc(grid), nxc-2, nyc-2, nzc-2);

  /* ------------------------ */
  /* Write the Magnetic field */
  /* ------------------------ */

  grpname = "Fields";
  dtaname = "Bx";
  outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getBxc(), nxc-2, nyc-2, nzc-2);

  grpname = "Fields";
  dtaname = "By";
  outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getByc(), nxc-2, nyc-2, nzc-2);

  grpname = "Fields";
  dtaname = "Bz";
  outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getBzc(), nxc-2, nyc-2, nzc-2);

  /* ----------------------------------------------- */
  /* Write the Charge Density field for each species */
  /* ----------------------------------------------- */

  for (int is = 0; is < col->getNs(); is++)
  {
    stringstream snmbr;
    snmbr << is;

    grpname = "Fields";
    dtaname = "Rho_" + snmbr.str();
    outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getRHOcs(grid, is), nxc-2, nyc-2, nzc-2);
  }

  /* ---------------------------------------- */
  /* Write the Current field for each species */
  /* ---------------------------------------- */

  for (int is = 0; is < col->getNs(); is++)
  {
    stringstream snmbr;
    snmbr << is;

    grpname = "Fields";
    dtaname = "Jx_" + snmbr.str();
    outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getJxsc(grid, is), nxc-2, nyc-2, nzc-2);

    grpname = "Fields";
    dtaname = "Jy_" + snmbr.str();
    outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getJysc(grid, is), nxc-2, nyc-2, nzc-2);

    grpname = "Fields";
    dtaname = "Jz_" + snmbr.str();
    outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getJzsc(grid, is), nxc-2, nyc-2, nzc-2);
  }

  outputfile.ClosePHDF5file();

}
