
#include "WriteOutputParallel.h"

void WriteOutputParallel(Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle){

  string       grpname;
  string       dtaname;

  stringstream filenmbr;
  string       filename;

  filenmbr << setfill('0') << setw(5) << cycle;
  filename = col->getSaveDirName() + "/" + col->getSimName() + "_" + filenmbr.str() + ".h5";

  int    dglob[3] = { col ->getNxc()  , col ->getNyc()  , col ->getNzc()   };
  int    dlocl[3] = { grid->getNXC()-2, grid->getNYC()-2, grid->getNZC()-2 };
  double L    [3] = { col ->getLx ()  , col ->getLy ()  , col ->getLz ()   };

  PHDF5fileClass outputfile(filename, 3, vct->getCoordinates(), vct->getComm());

  outputfile.CreatePHDF5file(L, dglob, dlocl, false);

  grpname = "Fields";
  dtaname = "Bx";
  outputfile.WritePHDF5dataset(grpname, dtaname, EMf->getBxc(), grid->getNXC()-2, grid->getNYC()-2, grid->getNZC()-2);

  outputfile.ClosePHDF5file();

}
