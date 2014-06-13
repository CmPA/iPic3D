#include <mpi.h>
#include <fstream>

#include "ParallelIO.h"
#include "MPIdata.h"

/*! Function used to write the EM fields using the parallel HDF5 library */
void WriteOutputParallel(Grid3DCU *grid, EMfields3D *EMf, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle){

#ifdef PHDF5
  string       grpname;
  string       dtaname;

  stringstream filenmbr;
  string       filename;

  bool         bp;

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

  bp = false;
  if (col->getParticlesOutputCycle() > 0) bp = true;

  outputfile.CreatePHDF5file(L, dglob, dlocl, bp);

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

#else  
  cout << " ERROR: The input file request the use of the Parallel HDF5 functions, but the code has been compiled using the sequential HDF5 library. " << endl;
  cout << "        Recompile the code using the parallel HDF5 options or change the input file options. " << endl;
  abort();
#endif

}

/*! Function to write the EM fields using the H5hut library. */
void WriteFieldsH5hut(int nspec, Grid3DCU *grid, EMfields3D *EMf, CollectiveIO *col, VCtopology3D *vct, int cycle){

#ifdef USEH5HUT

  H5output file;


  /* ---------------- */
  /* Write the fields */
  /* ---------------- */

  string filename = col->getSaveDirName() + "/" + col->getSimName();

  file.SetNameCycle(filename, cycle);

  file.OpenFieldsFile("Node", nspec, col->getNxc()+1, col->getNyc()+1, col->getNzc()+1, vct->getCoordinates(), vct->getDivisions(), vct->getComm());

  file.WriteFields(EMf->getEx(), "Ex", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getEy(), "Ey", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getEz(), "Ez", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getBx(), "Bx", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getBy(), "By", grid->getNXN(), grid->getNYN(), grid->getNZN());
  file.WriteFields(EMf->getBz(), "Bz", grid->getNXN(), grid->getNYN(), grid->getNZN());

  for (int is=0; is<nspec; is++) {
    stringstream  ss;
    ss << is;
    string s_is = ss.str();
    file.WriteFields(EMf->getRHOns(is), "rho_"+ s_is, grid->getNXN(), grid->getNYN(), grid->getNZN());
  }

  file.CloseFieldsFile();

//--- SAVE FIELDS IN THE CELLS:
//
//  file.OpenFieldsFile("Cell", nspec, col->getNxc(), col->getNyc(), col->getNzc(), vct->getCoordinates(), vct->getDivisions(), vct->getComm());
//
//  file.WriteFields(EMf->getExc(), "Exc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getEyc(), "Eyc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getEzc(), "Ezc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getBxc(), "Bxc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getByc(), "Byc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  file.WriteFields(EMf->getBzc(), "Bzc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//
//  for (int is=0; is<nspec; is++) {
//    stringstream  ss;
//    ss << is;
//    string s_is = ss.str();
//    file.WriteFields(EMf->getRHOcs(is), "rhoc_"+ s_is, grid->getNXC(), grid->getNYC(), grid->getNZC());
//  }
//
//  file.CloseFieldsFile();
//
//--- END SAVE FIELDS IN THE CELLS.

#else  
  cout << " ERROR: The input file request the use of the H5hut library, but the code has been compiled using other method. " << endl;
  cout << "        Recompile the code using the H5hut options or change the input file. " << endl;
  abort();
#endif

}

/*! Function to write the particles using the H5hut library. */
void WritePartclH5hut(int nspec, Grid3DCU *grid, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle){

#ifdef USEH5HUT

  H5output file;

  string filename = col->getSaveDirName() + "/" + col->getSimName();

  file.SetNameCycle(filename, cycle);

  /* ------------------- */
  /* Write the particles */
  /* ------------------- */

  file.OpenPartclFile(nspec, vct->getComm());
  for (int i=0; i<nspec; i++){
    file.WriteParticles(i, part[i].getNOP(),
                           part[i].getQall(),
                           part[i].getXall(),
                           part[i].getYall(),
                           part[i].getZall(),
                           part[i].getUall(),
                           part[i].getVall(),
                           part[i].getWall(),
                           vct->getComm());
  }
  file.ClosePartclFile();

#else  
  cout << " ERROR: The input file request the use of the H5hut library, but the code has been compiled using other method. " << endl;
  cout << "        Recompile the code using the H5hut options or change the input file. " << endl;
  abort();
#endif

}

void ReadPartclH5hut(int nspec, Particles3Dcomm *part, Collective *col, VCtopology3D *vct, Grid3DCU *grid){
#ifdef USEH5HUT

  H5input infile;
  double L[3] = {col->getLx(), col->getLy(), col->getLz()};

  infile.SetNameCycle(col->getinitfile(), col->getLast_cycle());
  infile.OpenPartclFile(nspec);

  infile.ReadParticles(vct->getCartesian_rank(), vct->getNproc(), vct->getDivisions(), L, vct->getComm());

  for (int s = 0; s < nspec; s++){
    part[s].allocate(s, infile.GetNp(s), col, vct, grid);

    infile.DumpPartclX(part[s].getXref(), s);
    infile.DumpPartclY(part[s].getYref(), s);
    infile.DumpPartclZ(part[s].getZref(), s);
    infile.DumpPartclU(part[s].getUref(), s);
    infile.DumpPartclV(part[s].getVref(), s);
    infile.DumpPartclW(part[s].getWref(), s);
    infile.DumpPartclQ(part[s].getQref(), s);
  }
  infile.ClosePartclFile();

//--- TEST PARTICLE LECTURE:
//  for (int s = 0; s < nspec; s++){
//    for (int n = 0; n < part[s].getNOP(); n++){
//      double ix = part[s].getX(n);
//      double iy = part[s].getY(n);
//      double iz = part[s].getZ(n);
//      if (ix<=0 || iy<=0 || iz <=0) {
//        cout << " ERROR: This particle has negative position. " << endl;
//        cout << "        n = " << n << "/" << part[s].getNOP();
//        cout << "       ix = " << ix;
//        cout << "       iy = " << iy;
//        cout << "       iz = " << iz;
//      }
//    }
//  }
//--- END TEST

#endif
}

void ReadFieldsH5hut(int nspec, EMfields3D *EMf, Collective *col, VCtopology3D *vct, Grid3DCU *grid){
#ifdef USEH5HUT

  H5input infile;

  infile.SetNameCycle(col->getinitfile(), col->getLast_cycle());

  infile.OpenFieldsFile("Node", nspec, col->getNxc()+1,
                                       col->getNyc()+1,
                                       col->getNzc()+1,
                                       vct->getCoordinates(),
                                       vct->getDivisions(),
                                       vct->getComm());

  infile.ReadFields(EMf->getEx(), "Ex", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getEy(), "Ey", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getEz(), "Ez", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getBx(), "Bx", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getBy(), "By", grid->getNXN(), grid->getNYN(), grid->getNZN());
  infile.ReadFields(EMf->getBz(), "Bz", grid->getNXN(), grid->getNYN(), grid->getNZN());

  for (int is = 0; is < nspec; is++){
    std::stringstream  ss;
    ss << is;
    std::string s_is = ss.str();
    infile.ReadFields(EMf->getRHOns(is), "rho_"+s_is, grid->getNXN(), grid->getNYN(), grid->getNZN());
  }

  infile.CloseFieldsFile();

  // initialize B on centers
  MPI_Barrier(MPI_COMM_WORLD);

  // Comm ghost nodes for B-field
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx(), col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy(), col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz(), col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  grid->interpN2C(EMf->getBxc(), EMf->getBx());
  grid->interpN2C(EMf->getByc(), EMf->getBy());
  grid->interpN2C(EMf->getBzc(), EMf->getBz());

  // Comm ghost cells for B-field
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx(), col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy(), col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz(), col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  // communicate E
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEx(), col->bcBx[0],col->bcBx[1],col->bcBx[2],col->bcBx[3],col->bcBx[4],col->bcBx[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEy(), col->bcBy[0],col->bcBy[1],col->bcBy[2],col->bcBy[3],col->bcBy[4],col->bcBy[5], vct);
  communicateNodeBC(grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEz(), col->bcBz[0],col->bcBz[1],col->bcBz[2],col->bcBz[3],col->bcBz[4],col->bcBz[5], vct);

  for (int is = 0; is < nspec; is++)
    grid->interpN2C(EMf->getRHOcs(), is, EMf->getRHOns());

//---READ FROM THE CELLS:
//
//  infile.OpenFieldsFile("Cell", nspec, col->getNxc(),
//                                       col->getNyc(),
//                                       col->getNzc(),
//                                       vct->getCoordinates(),
//                                       vct->getDivisions(),
//                                       vct->getComm());
//
//  infile.ReadFields(EMf->getExc(), "Exc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getEyc(), "Eyc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getEzc(), "Ezc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getBxc(), "Bxc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getByc(), "Byc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//  infile.ReadFields(EMf->getBzc(), "Bzc", grid->getNXC(), grid->getNYC(), grid->getNZC());
//
//  for (int is = 0; is < nspec; is++){
//    std::stringstream  ss;
//    ss << is;
//    std::string s_is = ss.str();
//    infile.ReadFields(EMf->getRHOcs(is, 0), "rhoc_"+s_is, grid->getNXC(), grid->getNYC(), grid->getNZC());
//  }
//
//  infile.CloseFieldsFile();
//
//  // initialize B on nodes
//  grid->interpC2N(EMf->getBx(), EMf->getBxc());
//  grid->interpC2N(EMf->getBy(), EMf->getByc());
//  grid->interpC2N(EMf->getBz(), EMf->getBzc());
//
//  for (int is = 0; is < nspec; is++)
//    grid->interpC2N(EMf->getRHOns(), is, EMf->getRHOcs());
//
//---END READ FROM THE CELLS

#endif
}
