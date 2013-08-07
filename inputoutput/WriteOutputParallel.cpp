
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

  array3_double arr3(nxc-2,nyc-2,nzc-2);

  grpname = "Fields";
  dtaname = "Ex";
  EMf->getExc(arr3,grid);
  outputfile.WritePHDF5dataset(grpname, dtaname, arr3, nxc-2, nyc-2, nzc-2);

  grpname = "Fields";
  dtaname = "Ey";
  EMf->getEyc(arr3,grid);
  outputfile.WritePHDF5dataset(grpname, dtaname, arr3, nxc-2, nyc-2, nzc-2);

  grpname = "Fields";
  dtaname = "Ez";
  EMf->getEzc(arr3,grid);
  outputfile.WritePHDF5dataset(grpname, dtaname, arr3, nxc-2, nyc-2, nzc-2);

  /* ------------------------ */
  /* Write the Magnetic field */
  /* ------------------------ */

  grpname = "Fields";
  dtaname = "Bx";
  EMf->getBxc(arr3);
  outputfile.WritePHDF5dataset(grpname, dtaname, arr3, nxc-2, nyc-2, nzc-2);

  grpname = "Fields";
  dtaname = "By";
  EMf->getByc(arr3);
  outputfile.WritePHDF5dataset(grpname, dtaname, arr3, nxc-2, nyc-2, nzc-2);

  grpname = "Fields";
  dtaname = "Bz";
  EMf->getBzc(arr3);
  outputfile.WritePHDF5dataset(grpname, dtaname, arr3, nxc-2, nyc-2, nzc-2);

  /* ----------------------------------------------- */
  /* Write the Charge Density field for each species */
  /* ----------------------------------------------- */

  for (int is = 0; is < col->getNs(); is++)
  {
    stringstream snmbr;
    snmbr << is;

    grpname = "Fields";
    dtaname = "Rho_" + snmbr.str();
    EMf->getRHOcs(arr3,grid, is);
    EMf->getRHOcs(arr3, grid, is);
    outputfile.WritePHDF5dataset(grpname, dtaname, arr3, nxc-2, nyc-2, nzc-2);
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
    EMf->getJxsc(arr3, grid, is);
    outputfile.WritePHDF5dataset(grpname, dtaname, arr3, nxc-2, nyc-2, nzc-2);

    grpname = "Fields";
    dtaname = "Jy_" + snmbr.str();
    EMf->getJysc(arr3, grid, is);
    outputfile.WritePHDF5dataset(grpname, dtaname, arr3, nxc-2, nyc-2, nzc-2);

    grpname = "Fields";
    dtaname = "Jz_" + snmbr.str();
    EMf->getJzsc(arr3, grid, is);
    outputfile.WritePHDF5dataset(grpname, dtaname, arr3, nxc-2, nyc-2, nzc-2);
  }

  outputfile.ClosePHDF5file();

}
