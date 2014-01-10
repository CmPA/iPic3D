
#include "WriteOutputH5hut.h"

void WriteOutputH5hut(int nspec, Grid3DCU *grid, EMfields3D *EMf, Particles3Dcomm *part, CollectiveIO *col, VCtopology3D *vct, int cycle){

#ifdef USEH5HUT

  H5output file;


  /* ---------------- */
  /* Write the fields */
  /* ---------------- */

  string filename = col->getSaveDirName() + "/" + col->getSimName();

  file.SetBaseName(filename);

  file.OpenFieldsCellFile(nspec, col->getNxc(), col->getNyc(), col->getNzc(), cycle, vct->getCoordinates(), vct->getDivisions(), vct->getComm());

  file.WriteFieldsCell(EMf->getExc(), "Ex", grid->getNXC(), grid->getNYC(), grid->getNZC());
  file.WriteFieldsCell(EMf->getEyc(), "Ey", grid->getNXC(), grid->getNYC(), grid->getNZC());
  file.WriteFieldsCell(EMf->getEzc(), "Ez", grid->getNXC(), grid->getNYC(), grid->getNZC());
  file.WriteFieldsCell(EMf->getBxc(), "Bx", grid->getNXC(), grid->getNYC(), grid->getNZC());
  file.WriteFieldsCell(EMf->getByc(), "By", grid->getNXC(), grid->getNYC(), grid->getNZC());
  file.WriteFieldsCell(EMf->getBzc(), "Bz", grid->getNXC(), grid->getNYC(), grid->getNZC());

  for (int is=0; is<nspec; is++) {
    stringstream  ss;
    ss << is;
    string s_is = ss.str();
    file.WriteFieldsCell(EMf->getRHOcs(is), "rho_"+ s_is, grid->getNXC(), grid->getNYC(), grid->getNZC());
    file.WriteFieldsCell(EMf->getJxsc(is) , "Jx_" + s_is, grid->getNXC(), grid->getNYC(), grid->getNZC());
    file.WriteFieldsCell(EMf->getJysc(is) , "Jy_" + s_is, grid->getNXC(), grid->getNYC(), grid->getNZC());
    file.WriteFieldsCell(EMf->getJzsc(is) , "Jz_" + s_is, grid->getNXC(), grid->getNYC(), grid->getNZC());
  }

  file.CloseFieldsCellFile();

  /* ------------------- */
  /* Write the particles */
  /* ------------------- */

  file.OpenPartclFile(nspec, cycle, vct->getComm());
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
