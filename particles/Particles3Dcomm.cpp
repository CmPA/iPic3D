/*******************************************************************************************
  Particles3Dcomm.cpp  -  Class for particles of the same species, in a 2D space and 3component velocity
  -------------------
developers: Stefano Markidis, Giovanni Lapenta.
 ********************************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include "VirtualTopology3D.h"
#include "VCtopology3D.h"
#include "Collective.h"
#include "ComParticles3D.h"
#include "Alloc.h"
#include "Basic.h"
#include "BcParticles.h"
#include "Grid.h"
#include "Grid3DCU.h"
#include "Field.h"
#include "MPIdata.h"

#include "Particles3Dcomm.h"

#include "hdf5.h"
#include <vector>
#include <complex>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define MIN_VAL   1E-32

#define EXIT_VAL 1E-14
/**
 * 
 * Class for particles of the same species, in a 2D space and 3component velocity
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */

/** constructor */
Particles3Dcomm::Particles3Dcomm() {
  // see allocate(int species, Collective* col, VirtualTopology3D* vct, Grid* grid)

}
/** deallocate particles */
Particles3Dcomm::~Particles3Dcomm() {
  delete[]x;
  delete[]y;
  delete[]z;
  delete[]u;
  delete[]v;
  delete[]w;
  delete[]q;
  // deallocate buffers
  delete[]b_X_RIGHT;
  delete[]b_X_LEFT;
  delete[]b_Y_RIGHT;
  delete[]b_Y_LEFT;
  delete[]b_Z_RIGHT;
  delete[]b_Z_LEFT;

  // deallocate mlmd-related stuff
  delete[]CommToChild_P;
  delete[]Rank_CommToChildren_P;

  delete[]RGPBC_Info;
  delArr3(PCGMsg, numChildren, MaxNumMsg);
  delArr2(nopPCGMsg, numChildren);
  delArr2(PRGMsg, RG_numPBCMessages);
  delete[]nopPRGMsg;
  delete[]PRGMsgArrived;
  delete[]PRGMsg_General;

  delete CRP_ToCoreH;
  delArr2(H_CRP_Msg, XLEN*YLEN*ZLEN);
  delete[]H_num_CRP_nop;
  delete[]H_CRP_General;
  delete[]H_CRP_cores;
}
/** constructors fo a single species*/
void Particles3Dcomm::allocate(int species, long long initnpmax, Collective * col, VirtualTopology3D * vct, Grid * grid) {

  /*! use this for system-wide mlmd output */
  MLMDVerbose= col->getMLMDVerbose();
  MPI_Comm_rank(MPI_COMM_WORLD, &SpokePerson);
  /*! end use this for system-wide mlmd output */
  
  /*! mlmd */
  numGrid = grid->getNumGrid();
  /*! end mlmd */
  
  // info from collectiveIO
  ns = species;
  npcel = col->getNpcel(species);
  npcelx = col->getNpcelx(species);
  npcely = col->getNpcely(species);
  npcelz = col->getNpcelz(species);

  XLEN= vct->getXLEN();
  YLEN= vct->getYLEN();
  ZLEN= vct->getZLEN();

  // This if is necessary to restart with H5hut-io
  if (initnpmax==0){
    /*! pre-mlmd
    long ncproc = int(col->getNxc()/col->getXLEN()) *
                  int(col->getNyc()/col->getYLEN()) *
                  int(col->getNzc()/col->getZLEN()); */
    
    long ncproc = int(col->getNxc_mlmd(numGrid)/XLEN) *
      int(col->getNyc_mlmd(numGrid)/YLEN) *
      int(col->getNzc_mlmd(numGrid)/ZLEN);

    nop   = ncproc * npcel;
    npmax = nop * col->getNpMaxNpRatio();
    
  }
  else {
    npmax = initnpmax*col->getNpMaxNpRatio();
    nop   = initnpmax;
  }

  /*if (vct->getCartesian_rank()== XLEN*YLEN*ZLEN-1){
    npmax= npmax *4;
    }*/

  rhoINIT   = col->getRHOinit(species);
  rhoINJECT = col->getRHOinject(species);

  qom = col->getQOM(species);
  uth = col->getUth(species);
  vth = col->getVth(species);
  wth = col->getWth(species);
  u0 = col->getU0(species);
  v0 = col->getV0(species);
  w0 = col->getW0(species);
  dt = col->getDt();
  Lx = col->getLx_mlmd(numGrid);
  Ly = col->getLy_mlmd(numGrid);
  Lz = col->getLz_mlmd(numGrid);
  dx = grid->getDX();  /*! already local to mlmd*/
  dy = grid->getDY();
  dz = grid->getDZ();
  int P= vct->getParentGridNum();
  DxP= grid->getDx_mlmd(P);
  DyP= grid->getDy_mlmd(P);
  DzP= grid->getDz_mlmd(P);
  RFx= DxP/dx;
  RFy= DyP/dy;
  RFz= DzP/dz;
  delta = col->getDelta();
  TrackParticleID = col->getTrackParticleID(species);
  c = col->getC();
  // info for mover
  NiterMover = col->getNiterMover();
  // velocity of the injection from the wall
  Vinj = col->getVinj();
  Ninj = col->getRHOinject(species);
  // info from Grid
  xstart = grid->getXstart();
  xend = grid->getXend();
  ystart = grid->getYstart();
  yend = grid->getYend();
  zstart = grid->getZstart();
  zend = grid->getZend();

  xStart_GC= grid->getxStart_GC();
  yStart_GC= grid->getyStart_GC();
  zStart_GC= grid->getzStart_GC();
  xEnd_GC= grid->getxEnd_GC();
  yEnd_GC= grid->getyEnd_GC();
  zEnd_GC= grid->getzEnd_GC();

  nxn = grid->getNXN();
  nyn = grid->getNYN();
  nzn = grid->getNZN();
  invVOL = grid->getInvVOL();
  // info from VirtualTopology3D
  cVERBOSE = vct->getcVERBOSE();

  // boundary condition for particles
  bcPfaceXright = col->getBcPfaceXright();
  bcPfaceXleft = col->getBcPfaceXleft();
  bcPfaceYright = col->getBcPfaceYright();
  bcPfaceYleft = col->getBcPfaceYleft();
  bcPfaceZright = col->getBcPfaceZright();
  bcPfaceZleft = col->getBcPfaceZleft();

  // //////////////////////////////////////////////////////////////
  // ////////////// ALLOCATE ARRAYS /////////////////////////
  // //////////////////////////////////////////////////////////////
  // positions
  x = new double[npmax];
  y = new double[npmax];
  z = new double[npmax];
  // velocities
  u = new double[npmax];
  v = new double[npmax];
  w = new double[npmax];
  // charge
  q = new double[npmax];
  // ID
  if (TrackParticleID) {
    ParticleID = new unsigned long[npmax];
    BirthRank[0] = vct->getCartesian_rank();
    if (vct->getNprocs() > 1)
      BirthRank[1] = (int) ceil(log10((double) (vct->getNprocs())));  // Number of digits needed for # of process in ID
    else
      BirthRank[1] = 1;
    if (BirthRank[1] + (int) ceil(log10((double) (npmax))) > 10 && BirthRank[0] == 0) {
      cerr << "Error: can't Track particles in Particles3Dcomm::allocate" << endl;
      cerr << "Unsigned long 'ParticleID' cannot store all the particles" << endl;
      return;
    }
  }
  // BUFFERS
  // the buffer size should be decided depending on number of particles
  // the buffer size should be decided depending on number of particles
  if (TrackParticleID)
    nVar = 8;
  else
    nVar = 7;
  buffer_size = (int) (.05 * nop * nVar + 1); // max: 5% of the particles in the processors is going out
  buffer_size_small = (int) (.01 * nop * nVar + 1); // max 1% not resizable 

  b_X_RIGHT = new double[buffer_size];
  b_X_RIGHT_ptr = b_X_RIGHT;    // alias to make the resize
  b_X_LEFT = new double[buffer_size];
  b_X_LEFT_ptr = b_X_LEFT;      // alias to make the resize
  b_Y_RIGHT = new double[buffer_size];
  b_Y_RIGHT_ptr = b_Y_RIGHT;    // alias to make the resize
  b_Y_LEFT = new double[buffer_size];
  b_Y_LEFT_ptr = b_Y_LEFT;      // alias to make the resize
  b_Z_RIGHT = new double[buffer_size];
  b_Z_RIGHT_ptr = b_Z_RIGHT;    // alias to make the resize
  b_Z_LEFT = new double[buffer_size];
  b_Z_LEFT_ptr = b_Z_LEFT;      // alias to make the resize

  // if RESTART is true initialize the particle in allocate method
  restart = col->getRestart_status();
  if (restart != 0) {
    if (vct->getCartesian_rank() == 0 && ns == 0)
      cout << "LOADING PARTICLES FROM RESTART FILE in " + col->getRestartDirName() + "/restart.hdf" << endl;
    stringstream ss;
    ss << vct->getCartesian_rank();
    string name_file = col->getRestartDirName() + "/restart" + ss.str() + ".hdf";
    // hdf stuff 
    hid_t file_id, dataspace;
    hid_t datatype, dataset_id;
    herr_t status;
    size_t size;
    hsize_t dims_out[1];        /* dataset dimensions */
    int status_n;

    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0) {
      cout << "couldn't open file: " << name_file << endl;
      cout << "RESTART NOT POSSIBLE" << endl;
    }

    stringstream species_name;
    species_name << ns;
    // the cycle of the last restart is set to 0
    string name_dataset = "/particles/species_" + species_name.str() + "/x/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    datatype = H5Dget_type(dataset_id);
    size = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id); /* dataspace handle */
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);

    // get how many particles there are on this processor for this species
    status_n = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    nop = dims_out[0];          // this the number of particles on the processor!
    // get x
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
    // close the data set
    status = H5Dclose(dataset_id);

    // get y
    name_dataset = "/particles/species_" + species_name.str() + "/y/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y);
    status = H5Dclose(dataset_id);

    // get z
    name_dataset = "/particles/species_" + species_name.str() + "/z/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, z);
    status = H5Dclose(dataset_id);

    // get u
    name_dataset = "/particles/species_" + species_name.str() + "/u/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, u);
    status = H5Dclose(dataset_id);
    // get v
    name_dataset = "/particles/species_" + species_name.str() + "/v/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
    status = H5Dclose(dataset_id);
    // get w
    name_dataset = "/particles/species_" + species_name.str() + "/w/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, w);
    status = H5Dclose(dataset_id);
    // get q
    name_dataset = "/particles/species_" + species_name.str() + "/q/cycle_0";
    dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, q);
    status = H5Dclose(dataset_id);
    // ID 
    if (TrackParticleID) {
      // herr_t (*old_func)(void*); // HDF 1.6
      H5E_auto2_t old_func;      // HDF 1.8.8
      void *old_client_data;
      H5Eget_auto2(H5E_DEFAULT, &old_func, &old_client_data);  // HDF 1.8.8
      /* Turn off error handling */
      // H5Eset_auto(NULL, NULL); // HDF 1.6
      H5Eset_auto2(H5E_DEFAULT, 0, 0); // HDF 1.8
      name_dataset = "/particles/species_" + species_name.str() + "/ID/cycle_0";
      dataset_id = H5Dopen2(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8

      // H5Eset_auto(old_func, old_client_data); // HDF 1.6
      H5Eset_auto2(H5E_DEFAULT, old_func, old_client_data);
      if (dataset_id > 0)
        status = H5Dread(dataset_id, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, ParticleID);
      else {
        for (register long long counter = 0; counter < nop; counter++)
          ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];
      }
    }
    // close the hdf file
    status = H5Fclose(file_id);
  }

  // here, mlmd stuff which does not necessarily need to be at the beginning
  numChildren= vct->getNumChildren();

  // communicators to parent/ children for this species
  CommToParent_P= vct->getCommToParent_P(ns);
  Rank_CommToParent_P= vct->getRank_CommToParent_P(ns);
  CommToChild_P= new MPI_Comm[numChildren];
  Rank_CommToChildren_P= new int[numChildren];
  for (int ch=0; ch<numChildren; ch++) {
    CommToChild_P[ch]= vct->getCommToChild_P(ch, ns);
    Rank_CommToChildren_P[ch]= vct->getRank_CommToChildren_P(ch, ns);
  }


  int MaxGridCoreN= vct->getMaxGridCoreN();
  MAX_RG_numPBCMessages= (int) (MaxGridCoreN*6+1);
  MAX_RG_numPBCMessages_LevelWide= MAX_RG_numPBCMessages*4;

  // sizes of the PCGMsg values set here (but allocated only if needed) 
  // to have the send/ receive vectors with the same size, I cook up a number based 
  //   on the coarsest grid 
  int MaxNopMLMD= npcelx*npcely*npcelz *( col->getNxc_mlmd(0)/ col->getXLEN_mlmd(0) )*( col->getNyc_mlmd(0)/ col->getYLEN_mlmd(0) )*(col->getNzc_mlmd(0)/ col->getZLEN_mlmd(0)) *4*4 ; 

  sizeCG_PBCMsg = MaxNopMLMD; // CG side
  sizeRG_PBCMsg = MaxNopMLMD; // RG side
  MAXsizePBCMsg = 20*MaxNopMLMD;

  // wether to allow resize of repopulated particle buffers
  AllowPMsgResize= col->getAllowPMsgResize();
  
  // end sizes of the PCGMsg values set here 

  // here, to be able to used RGPBC_struct as an MPI_Datatype
  MPI_RGPBC_struct_commit();
  // here, to be able to use the RepP_struct as an MPI_Datatype
  MPI_RepP_struct_commit();
  /* commit the structure for repopulated particle exchange within the RG */
  MPI_CRP_struct_commit();

  // PRA related initialisation
  // the # of PRA cells is a tmp variable
  //   the 'visible' variables are the index at which PRA starts/ ends 
  
  int PRACells = int(RFx); 
  if (RFy > RFx) PRACells= int (RFy);

  // added
  if (bcPfaceXleft== -1 )
    PRACells= 4;

  // Index at which the PRA starts/ ends
  // 
  PRA_XLeft_Start= 0;
  PRA_XLeft_End= 0+PRACells;
  PRA_XRight_End= nxn-1;
  PRA_XRight_Start= nxn-1 -PRACells ;

  PRA_YLeft_Start= 0;
  PRA_YLeft_End= 0+PRACells;
  PRA_YRight_End= nyn-1;
  PRA_YRight_Start= nyn-1 -PRACells ;

  PRA_ZLeft_Start= 0;
  PRA_ZLeft_End= 0+PRACells;
  PRA_ZRight_End= nzn-1;
  PRA_ZRight_Start= nzn-1 -PRACells ;

  /*// for extreme testing
  PRA_XLeft_Start= 0;
  PRA_XLeft_End= ceil(nxn/2);
  PRA_XRight_End= nxn-1;
  PRA_XRight_Start= ceil(nxn/2);

  PRA_YLeft_Start= 0;
  PRA_YLeft_End= ceil(nyn/2);
  PRA_YRight_End= nyn-1;
  PRA_YRight_Start= ceil(nyn/2);

  PRA_ZLeft_Start= 0;
  PRA_ZLeft_End= ceil(nzn/2);
  PRA_ZRight_End= nzn-1;
  PRA_ZRight_Start= ceil(nzn/2);*/

  /*if (numGrid ==1){
    cout << "nxn " << nxn << " nyn " << nyn << " nzn " << nzn << endl;
    cout << "PRA_XLeft_Start " << PRA_XLeft_Start << " PRA_XLeft_End " << PRA_XLeft_End << " PRA_XRight_End " << PRA_XRight_End << " PRA_XRight_Start " << PRA_XRight_Start <<endl; 
    cout << "PRA_YLeft_Start " << PRA_YLeft_Start << " PRA_YLeft_End " << PRA_YLeft_End << " PRA_YRight_End " << PRA_YRight_End << " PRA_YRight_Start " << PRA_YRight_Start <<endl; 
    cout << "PRA_ZLeft_Start " << PRA_ZLeft_Start << " PRA_ZLeft_End " << PRA_ZLeft_End << " PRA_ZRight_End " << PRA_ZRight_End << " PRA_ZRight_Start " << PRA_ZRight_Start <<endl;
    } */

  // corresponding coordinates, relative to the grid (NOT to the core)
  // when checking on these coords, you do not need to check if the core is a boundary core
  Coord_XLeft_Start  = -dx;
  Coord_XLeft_End    = (PRACells-1)*dx;
  Coord_XRight_End   = Lx+dx;
  Coord_XRight_Start = Lx+dx -dx*(PRACells);

  Coord_YLeft_Start  = -dy;
  Coord_YLeft_End    = (PRACells-1)*dy;
  Coord_YRight_End   = Ly+dy;
  Coord_YRight_Start = Ly+dy -dy*(PRACells);

  Coord_ZLeft_Start  = -dz;
  Coord_ZLeft_End    = (PRACells-1)*dz;
  Coord_ZRight_End   = Lz+dz;
  Coord_ZRight_Start = Lz+dz -dz*(PRACells);

  // for extreme testing
  /*Coord_XLeft_Start  = -dx;
  Coord_XLeft_End    = Lx/2;
  Coord_XRight_End   = Lx+dx;
  Coord_XRight_Start = Lx/2;

  Coord_YLeft_Start  = -dy;
  Coord_YLeft_End    = Ly/2;
  Coord_YRight_End   = Ly+dy;
  Coord_YRight_Start = Ly/2;

  Coord_ZLeft_Start  = -dz;
  Coord_ZLeft_End    = Lz/2;
  Coord_ZRight_End   = Lz+dz;
  Coord_ZRight_Start = Lz/2;*/

  /*if (vct->getCartesian_rank() ==0 and numGrid ==1){
    cout << "Coord_XLeft_Start: " << Coord_XLeft_Start <<" Coord_XLeft_End: " << Coord_XLeft_End << " Coord_XRight_End " << Coord_XRight_End << " Coord_XRight_Start " << Coord_XRight_Start <<" Lx " << Lx<< endl;
    cout << "Coord_YLeft_Start: " << Coord_YLeft_Start <<" Coord_YLeft_End: " << Coord_YLeft_End << " Coord_YRight_End " << Coord_YRight_End << " Coord_YRight_Start " << Coord_YRight_Start <<" Ly " << Ly<< endl;
    cout << "Coord_ZLeft_Start: " << Coord_ZLeft_Start <<" Coord_ZLeft_End: " << Coord_ZLeft_End << " Coord_ZRight_End " << Coord_ZRight_End << " Coord_ZRight_Start " << Coord_ZRight_Start <<" Lz " << Lz<< endl;
    }*/


  PRA_PAdded=0;
  MAX_np_CRP= 10*npmax; // exagerated number
  size_CRP=  nop;
  // //FOR TEST:
  // nvDistLoc = 3;
  // vDist     = new c_vDist[nvDistLoc];

  // double vR    = 2 * sqrt(dx*dx + dy*dy + dz*dz);
  // double vFact = 2.0;

  // vDist[0].init(species, 3.0 , 11.00, 0.02625, 256, 256, 256, vR, vFact, col, grid);
  // vDist[1].init(species, 6.36, 11.00, 0.02625, 256, 256, 256, vR, vFact, col, grid);
  // vDist[2].init(species, 10.0, 15.0 , 0.02625, 256, 256, 256, vR, vFact, col, grid);
  // //END FOR TEST

}

/** Initialie arrays for velocity distributions in 3D */
void c_vDist::init(int ispec, double vX, double vY, double vZ, int bi, int bj, int bk, double vR, double vFact, Collective * col, Grid * grid){
  vDistRad   = vR;

  dovDist3D = false;
  if (vX > grid->getXstart() && vX < grid->getXend())
  if (vY > grid->getYstart() && vY < grid->getYend())
  if (vZ > grid->getZstart() && vZ < grid->getZend())
    dovDist3D = true;

  if (dovDist3D) {
    vDistLoc_x = vX;
    vDistLoc_y = vY;
    vDistLoc_z = vZ;
    cout << " :: x,y,z= " << vDistLoc_x << " " << vDistLoc_y << " " << vDistLoc_z << endl;
    nBins_i = bi;
    nBins_j = bj;
    nBins_k = bk;
    vDist3D = newArr3(unsigned long, nBins_i, nBins_j, nBins_k);
    for (int i = 0; i < nBins_i; i++)
    for (int j = 0; j < nBins_j; j++)
    for (int k = 0; k < nBins_k; k++) {
      vDist3D[i][j][k] = 0;
    }
    vBinBeg_i = col->getU0(ispec) - vFact * col->getUth(ispec);
    vBinEnd_i = col->getU0(ispec) + vFact * col->getUth(ispec);
    vBinBeg_j = col->getV0(ispec) - vFact * col->getVth(ispec);
    vBinEnd_j = col->getV0(ispec) + vFact * col->getVth(ispec);
    vBinBeg_k = col->getW0(ispec) - vFact * col->getWth(ispec);
    vBinEnd_k = col->getW0(ispec) + vFact * col->getWth(ispec);
    dv_i      = (vBinEnd_i - vBinBeg_i) / double(nBins_i);
    dv_j      = (vBinEnd_j - vBinBeg_j) / double(nBins_j);
    dv_k      = (vBinEnd_k - vBinBeg_k) / double(nBins_k);
  }
}

void c_vDist::add(double x, double y, double z, double u, double v, double w){

  double r2 = (x-vDistLoc_x)*(x-vDistLoc_x) + (y-vDistLoc_y)*(y-vDistLoc_y) + (z-vDistLoc_z)*(z-vDistLoc_z);
  cout << " --+-- r2, vDistRad = " << r2 << " " << vDistRad << endl;

  if (r2 < vDistRad*vDistRad) {

    int i = int((u-vBinBeg_i)/dv_i);
    int j = int((v-vBinBeg_j)/dv_j);
    int k = int((w-vBinBeg_k)/dv_k);

    if (i < 0) i = 0; else if (i >= nBins_i) i = nBins_i;
    if (j < 0) j = 0; else if (j >= nBins_j) j = nBins_j;
    if (k < 0) k = 0; else if (k >= nBins_k) k = nBins_k;

    vDist3D[i][j][k] += 1;
  }
  
}

/** Add velocity distributions in 3D */
void Particles3Dcomm::Add_vDist3D(){

  for (int i=0; i<nvDistLoc; i++) {

    if (vDist[i].get_doVdist()) {
      for (long long p=0; p<nop ;p++) {
        cout << p << "/" << nop;
        vDist[i].add(x[p], y[p], z[p], u[p], v[p], w[p]);
      }
    }

  }
}

/** Print velocity distributions in 3D */
void Particles3Dcomm::Write_vDist3D(string SaveDirName){

  for (int n=0; n<nvDistLoc; n++) {

    if (vDist[n].get_doVdist()) {

      ofstream myfile;
      stringstream ss;
      ss << SaveDirName << "/vDist_" << n << ".vtk";
      string filename = ss.str();
      myfile.open(filename.c_str(), ios::trunc);

      myfile << "# vtk DataFile Version 1.0" << endl;
      myfile << "Electric Field from Parsek" << endl;
      myfile << "ASCII" << endl;
      myfile << "DATASET STRUCTURED_POINTS" << endl;
      myfile << "DIMENSIONS " << vDist[n].get_dim_i() << " " << vDist[n].get_dim_j() << " " << vDist[n].get_dim_k() << endl;
      myfile << "ORIGIN " << vDist[n].get_vBinBeg_i() << " " << vDist[n].get_vBinBeg_j() << " " << vDist[n].get_vBinBeg_k() << endl;
      myfile << "SPACING " << vDist[n].get_dvi() << " " << vDist[n].get_dvj() << " " << vDist[n].get_dvk() << endl;
      myfile << "POINT_DATA " << vDist[n].get_ntotBins() << endl;
      myfile << "VECTORS E float" << endl;
      myfile << "LOOKUP_TABLE default" << endl;

      for (int i=0; i<vDist[n].get_nBinsi(); i++)
        for (int j=0; j<vDist[n].get_nBinsj(); j++)
          for (int k=0; k<vDist[n].get_nBinsk(); k++)
            myfile << vDist[n].get(i, j, k);

      myfile.close();
      
    }

  }
}

/** calculate the weights given the position of particles 0,0,0 is the left,left, left node */
void Particles3Dcomm::calculateWeights(double weight[][2][2], double xp, double yp, double zp, int ix, int iy, int iz, Grid * grid) {
  double xi[2], eta[2], zeta[2];
  xi[0] = xp - grid->getXN(ix - 1, iy, iz);
  eta[0] = yp - grid->getYN(ix, iy - 1, iz);
  zeta[0] = zp - grid->getZN(ix, iy, iz - 1);
  xi[1] = grid->getXN(ix, iy, iz) - xp;
  eta[1] = grid->getYN(ix, iy, iz) - yp;
  zeta[1] = grid->getZN(ix, iy, iz) - zp;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        weight[i][j][k] = xi[i] * eta[j] * zeta[k] * invVOL;
}


/** Interpolation Particle --> Grid */
void Particles3Dcomm::interpP2G(Field * EMf, Grid * grid, VirtualTopology3D * vct) {
  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
  const double nxn = grid->getNXN();
  const double nyn = grid->getNYN();
  const double nzn = grid->getNZN();

  for (int i=0; i< nop; i++){
    if (z[i]==1.00208){
      cout << "z== 1.00208 appeared in interpP2G" << endl;
    }
  }

  //#pragma omp parallel
  {
    //Moments speciesMoments(nxn,nyn,nzn,invVOL);
    //speciesMoments.set_to_zero();
    //#pragma omp for
    for (register long long i = 0; i < nop; i++)
    {
      const int ix = 2 + int (floor((x[i] - xstart) * inv_dx));
      const int iy = 2 + int (floor((y[i] - ystart) * inv_dy));
      const int iz = 2 + int (floor((z[i] - zstart) * inv_dz));
      double temp[2][2][2];
      double xi[2], eta[2], zeta[2];
      xi[0] = x[i] - grid->getXN(ix - 1, iy, iz);
      eta[0] = y[i] - grid->getYN(ix, iy - 1, iz);
      zeta[0] = z[i] - grid->getZN(ix, iy, iz - 1);
      xi[1] = grid->getXN(ix, iy, iz) - x[i];
      eta[1] = grid->getYN(ix, iy, iz) - y[i];
      zeta[1] = grid->getZN(ix, iy, iz) - z[i];
      double weight[2][2][2];
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++) {
            weight[ii][jj][kk] = q[i] * xi[ii] * eta[jj] * zeta[kk] * invVOL;
          }
      //weight[0][0][0] = q[i] * xi[0] * eta[0] * zeta[0] * invVOL;
      //weight[0][0][1] = q[i] * xi[0] * eta[0] * zeta[1] * invVOL;
      //weight[0][1][0] = q[i] * xi[0] * eta[1] * zeta[0] * invVOL;
      //weight[0][1][1] = q[i] * xi[0] * eta[1] * zeta[1] * invVOL;
      //weight[1][0][0] = q[i] * xi[1] * eta[0] * zeta[0] * invVOL;
      //weight[1][0][1] = q[i] * xi[1] * eta[0] * zeta[1] * invVOL;
      //weight[1][1][0] = q[i] * xi[1] * eta[1] * zeta[0] * invVOL;
      //weight[1][1][1] = q[i] * xi[1] * eta[1] * zeta[1] * invVOL;
      // add charge density
      //EMf->addRho(weight, ix, iy, iz, ns);
      EMf->addRho(weight, ix, iy, iz, ns, vct, x[i], y[i], z[i]);
      // add current density - X
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = u[i] * weight[ii][jj][kk];
      EMf->addJx(temp, ix, iy, iz, ns);
      // add current density - Y
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = v[i] * weight[ii][jj][kk];
      EMf->addJy(temp, ix, iy, iz, ns);
      // add current density - Z
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = w[i] * weight[ii][jj][kk];
      EMf->addJz(temp, ix, iy, iz, ns);
      // Pxx - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = u[i] * u[i] * weight[ii][jj][kk];
      EMf->addPxx(temp, ix, iy, iz, ns);
      // Pxy - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = u[i] * v[i] * weight[ii][jj][kk];
      EMf->addPxy(temp, ix, iy, iz, ns);
      // Pxz - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = u[i] * w[i] * weight[ii][jj][kk];
      EMf->addPxz(temp, ix, iy, iz, ns);
      // Pyy - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = v[i] * v[i] * weight[ii][jj][kk];
      EMf->addPyy(temp, ix, iy, iz, ns);
      // Pyz - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = v[i] * w[i] * weight[ii][jj][kk];
      EMf->addPyz(temp, ix, iy, iz, ns);
      // Pzz - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = w[i] * w[i] * weight[ii][jj][kk];
      EMf->addPzz(temp, ix, iy, iz, ns);
    }
    // change this to allow more parallelization after implementing array class
    //#pragma omp critical
    //EMf->addToSpeciesMoments(speciesMoments,ns);
  }
  // communicate contribution from ghost cells 

  EMf->communicateGhostP2G(ns, bcPfaceXright, bcPfaceXleft, bcPfaceYright, bcPfaceYleft, bcPfaceZright, bcPfaceZleft, vct);
}

/** communicate buffers */
int Particles3Dcomm::communicate(VirtualTopology3D * ptVCT) {
  // allocate buffers
  MPI_Status status;
  int new_buffer_size;
  int npExitingMax;
  // variable for memory availability of space for new particles
  int avail, availALL, avail1, avail2, avail3, avail4, avail5, avail6;
  for (int i = 0; i < buffer_size; i++) {
    b_X_RIGHT[i] = MIN_VAL;
    b_X_LEFT[i] = MIN_VAL;
    b_Y_RIGHT[i] = MIN_VAL;
    b_Y_LEFT[i] = MIN_VAL;
    b_Z_RIGHT[i] = MIN_VAL;
    b_Z_LEFT[i] = MIN_VAL;
  }
  npExitXright = 0, npExitXleft = 0, npExitYright = 0, npExitYleft = 0, npExitZright = 0, npExitZleft = 0, npExit = 0, rightDomain = 0;
  long long np_current = 0, nplast = nop - 1;

  double xMin, yMin, zMin;
  double xMax, yMax, zMax;

  while (np_current < nplast+1){
     
    xMin=0; yMin=0; zMin=0;
    xMax=Lx; yMax=Ly; zMax=Ly;

    // BC on particles
    if (x[np_current] < xMin && ptVCT->getXleft_neighbor_P() == MPI_PROC_NULL)
      BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft);
    else if (x[np_current] > xMax && ptVCT->getXright_neighbor_P() == MPI_PROC_NULL)
      BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft); 
    if (y[np_current] < yMin && ptVCT->getYleft_neighbor_P() == MPI_PROC_NULL)  // check it here
      BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft);
    else if (y[np_current] > yMax && ptVCT->getYright_neighbor_P() == MPI_PROC_NULL) //check it here
      BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft); 
    if (z[np_current] < zMin && ptVCT->getZleft_neighbor_P() == MPI_PROC_NULL)  // check it here
      BCpart(&z[np_current],&w[np_current],&u[np_current],&v[np_current],Lz,wth,uth,vth,bcPfaceZright,bcPfaceZleft);
    else if (z[np_current] > zMax && ptVCT->getZright_neighbor_P() == MPI_PROC_NULL) //check it here
      BCpart(&z[np_current],&w[np_current],&u[np_current],&v[np_current],Lz,wth,uth,vth,bcPfaceZright,bcPfaceZleft);

    // mlmd: delete only particles that get out of the grid completely
    if (bcPfaceXleft <0 and CommToParent_P!= MPI_COMM_NULL){
      xMin= Coord_XLeft_Start; xMax= Coord_XRight_End;
      yMin= Coord_YLeft_Start; yMax= Coord_YRight_End;
      zMin= Coord_ZLeft_Start; zMax= Coord_ZRight_End;
    }else{
      // so, if periodic, particles are not deleted
      xMin=-Lx; xMax=2*Lx;
      yMin=-Ly; yMax=2*Ly;
      zMin=-Lz; zMax=2*Lz;
    }
    

    if (x[np_current] < xMin or x[np_current]> xMax or y[np_current] < yMin or y[np_current]> yMax or z[np_current] < zMin or z[np_current]> zMax){
      // particle to delete
      del_pack(np_current,&nplast);
    }else if (x[np_current] < xstart && ptVCT->getXleft_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitXleft+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferXleft(b_X_LEFT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      npExitXleft++;
      } 
    
    else if (x[np_current] > xend && ptVCT->getXright_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitXright+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferXright(b_X_RIGHT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      npExitXright++;
    }
    
    else  if (y[np_current] < ystart && ptVCT->getYleft_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitYleft+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferYleft(b_Y_LEFT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      npExitYleft++;
    }
    
    else if (y[np_current] > yend && ptVCT->getYright_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitYright+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferYright(b_Y_RIGHT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      npExitYright++;
    }
    else if (z[np_current] < zstart && ptVCT->getZleft_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitZleft+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferZleft(b_Z_LEFT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      
      npExitZleft++;
    } 
    
    else if (z[np_current] > zend && ptVCT->getZright_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitZright+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferZright(b_Z_RIGHT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      
      npExitZright++;
    
    } // end else you have to move particle
    else {
      // particle ok
      // particle is still in the domain, procede with the next particle
      np_current++;
    }
    
  }
  

  nop = nplast + 1;
  npExitingMax = 0;
  // calculate the maximum number of particles exiting from this domain
  // use this value to check if communication is needed
  // and to resize the buffer
  npExitingMax = maxNpExiting();
  // broadcast the maximum number of particles exiting for sizing the buffer and to check if communication is really needed
  /*! mlmd: i need the communicator also */
  //npExitingMax = reduceMaxNpExiting(npExitingMax);
  npExitingMax = reduceMaxNpExiting(npExitingMax, ptVCT->getCommGrid()); 

  /*****************************************************/
  /* SEND AND RECEIVE MESSAGES */
  /*****************************************************/

  new_buffer_size = npExitingMax * nVar + 1;

  if (new_buffer_size > buffer_size) {
    cout << "resizing the receiving buffer" << endl;
    resize_buffers(new_buffer_size);
  }

  if (npExitingMax > 0) {
  
    communicateParticles(new_buffer_size, b_X_LEFT, b_X_RIGHT, b_Y_LEFT, b_Y_RIGHT, b_Z_LEFT, b_Z_RIGHT, ptVCT);

    // UNBUFFERING
    /*! mlmd: need the communicator also */
    avail1 = unbuffer(b_X_RIGHT, ptVCT->getCommGrid());
    avail2 = unbuffer(b_X_LEFT, ptVCT->getCommGrid());
    avail3 = unbuffer(b_Y_RIGHT, ptVCT->getCommGrid());
    avail4 = unbuffer(b_Y_LEFT, ptVCT->getCommGrid());
    avail5 = unbuffer(b_Z_RIGHT, ptVCT->getCommGrid());
    avail6 = unbuffer(b_Z_LEFT, ptVCT->getCommGrid());

    // if one of these numbers is negative than there is not enough space for particles
    avail = avail1 + avail2 + avail3 + avail4 + avail5 + avail6;
    /*! mlmd: i need the communicator also */
    //availALL = reduceNumberParticles(avail);
    availALL = reduceNumberParticles(avail, ptVCT->getCommGrid());
    if (availALL < 0)
      return (-1);              // too many particles coming, save data nad stop simulation

  }
  
  nop_EndCommunicate= nop;

  return(0);

}

int Particles3Dcomm::communicateAfterMover(VirtualTopology3D * ptVCT) {
  // allocate buffers

  if (! (CommToParent_P!= MPI_COMM_NULL and bcPfaceXleft <0)) return 0;

  MPI_Status status;
  int new_buffer_size;
  int npExitingMax;
  // variable for memory availability of space for new particles
  int avail, availALL, avail1, avail2, avail3, avail4, avail5, avail6;
  for (int i = 0; i < buffer_size; i++) {
    b_X_RIGHT[i] = MIN_VAL;
    b_X_LEFT[i] = MIN_VAL;
    b_Y_RIGHT[i] = MIN_VAL;
    b_Y_LEFT[i] = MIN_VAL;
    b_Z_RIGHT[i] = MIN_VAL;
    b_Z_LEFT[i] = MIN_VAL;
  }
  npExitXright = 0, npExitXleft = 0, npExitYright = 0, npExitYleft = 0, npExitZright = 0, npExitZleft = 0, npExit = 0, rightDomain = 0;
  long long np_current = 0, nplast = nop - 1;

  double xMin, yMin, zMin;
  double xMax, yMax, zMax;

  while (np_current < nplast+1){
     
    xMin=0; yMin=0; zMin=0;
    xMax=Lx; yMax=Ly; zMax=Ly;

    // BC on particles
    if (x[np_current] < xMin && ptVCT->getXleft_neighbor_P() == MPI_PROC_NULL)
      BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft);
    else if (x[np_current] > xMax && ptVCT->getXright_neighbor_P() == MPI_PROC_NULL)
      BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft); 
    if (y[np_current] < yMin && ptVCT->getYleft_neighbor_P() == MPI_PROC_NULL)  // check it here
      BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft);
    else if (y[np_current] > yMax && ptVCT->getYright_neighbor_P() == MPI_PROC_NULL) //check it here
      BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft); 
    if (z[np_current] < zMin && ptVCT->getZleft_neighbor_P() == MPI_PROC_NULL)  // check it here
      BCpart(&z[np_current],&w[np_current],&u[np_current],&v[np_current],Lz,wth,uth,vth,bcPfaceZright,bcPfaceZleft);
    else if (z[np_current] > zMax && ptVCT->getZright_neighbor_P() == MPI_PROC_NULL) //check it here
      BCpart(&z[np_current],&w[np_current],&u[np_current],&v[np_current],Lz,wth,uth,vth,bcPfaceZright,bcPfaceZleft);

    if (bcPfaceXleft==-1 and CommToParent_P!= MPI_COMM_NULL){
      xMin= Coord_XLeft_End; xMax= Coord_XRight_Start;
      yMin= Coord_YLeft_End; yMax= Coord_YRight_Start;
      zMin= Coord_ZLeft_End; zMax= Coord_ZRight_Start;
    }
    else if (bcPfaceXleft==-2 and CommToParent_P!= MPI_COMM_NULL){
      xMin= Coord_XLeft_Start; xMax= Coord_XRight_End;
      yMin= Coord_YLeft_Start; yMax= Coord_YRight_End;
      zMin= Coord_ZLeft_Start; zMax= Coord_ZRight_End;
    }
    else{
      // so, if periodic, particles are not deleted
      xMin=-Lx; xMax=2*Lx;
      yMin=-Ly; yMax=2*Ly;
      zMin=-Lz; zMax=2*Lz;
    }
    

    if (x[np_current] < xMin or x[np_current]> xMax or y[np_current] < yMin or y[np_current]> yMax or z[np_current] < zMin or z[np_current]> zMax){
      // particle to delete
      del_pack(np_current,&nplast);

    }else if (x[np_current] < xstart && ptVCT->getXleft_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitXleft+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferXleft(b_X_LEFT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      npExitXleft++;
      } 
    
    else if (x[np_current] > xend && ptVCT->getXright_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitXright+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferXright(b_X_RIGHT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      npExitXright++;
    }
    
    else  if (y[np_current] < ystart && ptVCT->getYleft_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitYleft+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferYleft(b_Y_LEFT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      npExitYleft++;
    }
    
    else if (y[np_current] > yend && ptVCT->getYright_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitYright+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferYright(b_Y_RIGHT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      npExitYright++;
    }
    else if (z[np_current] < zstart && ptVCT->getZleft_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitZleft+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferZleft(b_Z_LEFT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      
      npExitZleft++;
    } 
    
    else if (z[np_current] > zend && ptVCT->getZright_neighbor_P() != MPI_PROC_NULL){
      // check if there is enough space in the buffer before putting in the particle
      if(((npExitZright+1)*nVar)>=buffer_size){
	resize_buffers((int) (buffer_size*2)); 
      }
      // put it in the communication buffer
      bufferZright(b_Z_RIGHT,np_current,ptVCT);
      // delete the particle and pack the particle array, the value of nplast changes
      del_pack(np_current,&nplast);
      
      npExitZright++;
    
    } // end else you have to move particle
    else {
      // particle ok
      // particle is still in the domain, procede with the next particle
      np_current++;
    }
    
  }
  

  nop = nplast + 1;
  npExitingMax = 0;
  // calculate the maximum number of particles exiting from this domain
  // use this value to check if communication is needed
  // and to resize the buffer
  npExitingMax = maxNpExiting();
  // broadcast the maximum number of particles exiting for sizing the buffer and to check if communication is really needed
  /*! mlmd: i need the communicator also */
  //npExitingMax = reduceMaxNpExiting(npExitingMax);
  npExitingMax = reduceMaxNpExiting(npExitingMax, ptVCT->getCommGrid()); 

  /*****************************************************/
  /* SEND AND RECEIVE MESSAGES */
  /*****************************************************/

  new_buffer_size = npExitingMax * nVar + 1;

  if (new_buffer_size > buffer_size) {
    cout << "resizing the receiving buffer" << endl;
    resize_buffers(new_buffer_size);
  }

  if (npExitingMax > 0) {
  
    communicateParticles(new_buffer_size, b_X_LEFT, b_X_RIGHT, b_Y_LEFT, b_Y_RIGHT, b_Z_LEFT, b_Z_RIGHT, ptVCT);

    // UNBUFFERING
    /*! mlmd: need the communicator also */
    avail1 = unbuffer(b_X_RIGHT, ptVCT->getCommGrid());
    avail2 = unbuffer(b_X_LEFT, ptVCT->getCommGrid());
    avail3 = unbuffer(b_Y_RIGHT, ptVCT->getCommGrid());
    avail4 = unbuffer(b_Y_LEFT, ptVCT->getCommGrid());
    avail5 = unbuffer(b_Z_RIGHT, ptVCT->getCommGrid());
    avail6 = unbuffer(b_Z_LEFT, ptVCT->getCommGrid());

    // if one of these numbers is negative than there is not enough space for particles
    avail = avail1 + avail2 + avail3 + avail4 + avail5 + avail6;
    /*! mlmd: i need the communicator also */
    //availALL = reduceNumberParticles(avail);
    availALL = reduceNumberParticles(avail, ptVCT->getCommGrid());
    if (availALL < 0)
      return (-1);              // too many particles coming, save data nad stop simulation

  }
  
  nop_EndCommunicate= nop;

  return(0);

}



/** resize the buffers */
void Particles3Dcomm::resize_buffers(int new_buffer_size) {
  cout << "RESIZING FROM " << buffer_size << " TO " << new_buffer_size << endl;
  // resize b_X_LEFT
  double *temp = new double[buffer_size];

  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_X_LEFT_ptr[i];
  // delete[] b_X_LEFT_ptr;
  delete[]b_X_LEFT;
  b_X_LEFT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_X_LEFT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_X_LEFT[i] = MIN_VAL;

  // resize b_X_RIGHT 
  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_X_RIGHT_ptr[i];
  // delete[] b_X_RIGHT_ptr;
  delete[]b_X_RIGHT;
  b_X_RIGHT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_X_RIGHT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_X_RIGHT[i] = MIN_VAL;

  // resize b_Y_RIGHT
  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_Y_RIGHT_ptr[i];
  // delete[] b_Y_RIGHT_ptr;
  delete[]b_Y_RIGHT;
  b_Y_RIGHT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_Y_RIGHT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_Y_RIGHT[i] = MIN_VAL;

  // resize b_Y_LEFT
  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_Y_LEFT_ptr[i];
  // delete[] b_Y_LEFT_ptr;
  delete[]b_Y_LEFT;
  b_Y_LEFT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_Y_LEFT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_Y_LEFT[i] = MIN_VAL;

  // resize b_Z_RIGHT
  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_Z_RIGHT_ptr[i];
  // delete[] b_Z_RIGHT_ptr;
  delete[]b_Z_RIGHT;
  b_Z_RIGHT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_Z_RIGHT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_Z_RIGHT[i] = MIN_VAL;

  // resize b_Z_LEFT
  for (int i = 0; i < buffer_size; i++)
    temp[i] = b_Z_LEFT_ptr[i];
  // delete[] b_Z_LEFT_ptr;
  delete[]b_Z_LEFT;
  b_Z_LEFT = new double[new_buffer_size];
  for (int i = 0; i < buffer_size; i++)
    b_Z_LEFT[i] = temp[i];
  for (int i = buffer_size; i < new_buffer_size; i++)
    b_Z_LEFT[i] = MIN_VAL;

  delete[]temp;

  b_X_RIGHT_ptr = b_X_RIGHT;
  b_Y_RIGHT_ptr = b_Y_RIGHT;
  b_Z_RIGHT_ptr = b_Z_RIGHT;
  b_Y_LEFT_ptr = b_Y_LEFT;
  b_X_LEFT_ptr = b_X_LEFT;
  b_Z_LEFT_ptr = b_Z_LEFT;

  buffer_size = new_buffer_size;
}
/** Put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
void Particles3Dcomm::bufferXleft(double *b_, long long np_current, VirtualTopology3D * vct) {
  if (x[np_current] < 0 and vct->getPERIODICX_P())
    b_[npExitXleft * nVar] = x[np_current] + Lx;  // this applies to the the leftmost processor
  else
    b_[npExitXleft * nVar] = x[np_current];
  b_[npExitXleft * nVar + 1] = y[np_current];
  b_[npExitXleft * nVar + 2] = z[np_current];
  b_[npExitXleft * nVar + 3] = u[np_current];
  b_[npExitXleft * nVar + 4] = v[np_current];
  b_[npExitXleft * nVar + 5] = w[np_current];
  b_[npExitXleft * nVar + 6] = q[np_current];
  if (TrackParticleID)
    b_[npExitXleft * nVar + 7] = ParticleID[np_current];
}
/** put a particle exiting to X-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
void Particles3Dcomm::bufferXright(double *b_, long long np_current, VirtualTopology3D * vct) {
  if (x[np_current] > Lx and vct->getPERIODICX_P())
    b_[npExitXright * nVar] = x[np_current] - Lx; // this applies to the right most processor
  else
    b_[npExitXright * nVar] = x[np_current];
  b_[npExitXright * nVar + 1] = y[np_current];
  b_[npExitXright * nVar + 2] = z[np_current];
  b_[npExitXright * nVar + 3] = u[np_current];
  b_[npExitXright * nVar + 4] = v[np_current];
  b_[npExitXright * nVar + 5] = w[np_current];
  b_[npExitXright * nVar + 6] = q[np_current];
  if (TrackParticleID)
    b_[npExitXright * nVar + 7] = ParticleID[np_current];
}
/** put a particle exiting to Y-LEFT in the bufferYLEFT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferYleft(double *b_, long long np_current, VirtualTopology3D * vct) {
  b_[npExitYleft * nVar] = x[np_current];
  if (y[np_current] < 0 and vct->getPERIODICY_P())
    b_[npExitYleft * nVar + 1] = y[np_current] + Ly;
  else
    b_[npExitYleft * nVar + 1] = y[np_current];
  b_[npExitYleft * nVar + 2] = z[np_current];
  b_[npExitYleft * nVar + 3] = u[np_current];
  b_[npExitYleft * nVar + 4] = v[np_current];
  b_[npExitYleft * nVar + 5] = w[np_current];
  b_[npExitYleft * nVar + 6] = q[np_current];
  if (TrackParticleID)
    b_[npExitYleft * nVar + 7] = ParticleID[np_current];
}
/** put a particle exiting to Y-RIGHT in the bufferYRIGHT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferYright(double *b_, long long np_current, VirtualTopology3D * vct) {
  b_[npExitYright * nVar] = x[np_current];
  if (y[np_current] > Ly and vct->getPERIODICY_P())
    b_[npExitYright * nVar + 1] = y[np_current] - Ly;
  else
    b_[npExitYright * nVar + 1] = y[np_current];
  b_[npExitYright * nVar + 2] = z[np_current];
  b_[npExitYright * nVar + 3] = u[np_current];
  b_[npExitYright * nVar + 4] = v[np_current];
  b_[npExitYright * nVar + 5] = w[np_current];
  b_[npExitYright * nVar + 6] = q[np_current];
  if (TrackParticleID)
    b_[npExitYright * nVar + 7] = ParticleID[np_current];
}
/** put a particle exiting to Z-LEFT in the bufferZLEFT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferZleft(double *b_, long long np_current, VirtualTopology3D * vct) {
  b_[npExitZleft * nVar] = x[np_current];
  b_[npExitZleft * nVar + 1] = y[np_current];
  if (z[np_current] < 0 and vct->getPERIODICZ_P())
    b_[npExitZleft * nVar + 2] = z[np_current] + Lz;
  else
    b_[npExitZleft * nVar + 2] = z[np_current];
  b_[npExitZleft * nVar + 3] = u[np_current];
  b_[npExitZleft * nVar + 4] = v[np_current];
  b_[npExitZleft * nVar + 5] = w[np_current];
  b_[npExitZleft * nVar + 6] = q[np_current];
  if (TrackParticleID)
    b_[npExitZleft * nVar + 7] = ParticleID[np_current];
}
/** put a particle exiting to Z-RIGHT in the bufferZRIGHT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferZright(double *b_, long long np_current, VirtualTopology3D * vct) {
  b_[npExitZright * nVar] = x[np_current];
  b_[npExitZright * nVar + 1] = y[np_current];
  if (z[np_current] > Lz and vct->getPERIODICZ_P())
    b_[npExitZright * nVar + 2] = z[np_current] - Lz;
  else
    b_[npExitZright * nVar + 2] = z[np_current];
  b_[npExitZright * nVar + 3] = u[np_current];
  b_[npExitZright * nVar + 4] = v[np_current];
  b_[npExitZright * nVar + 5] = w[np_current];
  b_[npExitZright * nVar + 6] = q[np_current];
  if (TrackParticleID)
    b_[npExitZright * nVar + 7] = ParticleID[np_current];
}
/** This unbuffer the last communication */
/*! mlmd: i need the communicator also */
/*!int Particles3Dcomm::unbuffer(double *b_) { */
int Particles3Dcomm::unbuffer(double *b_, MPI_Comm Comm) {
  long long np_current = 0;
  // put the new particles at the end of the array, and update the number of particles
  while (b_[np_current * nVar] != MIN_VAL) {
    x[nop] = b_[nVar * np_current];
    y[nop] = b_[nVar * np_current + 1];
    z[nop] = b_[nVar * np_current + 2];
    u[nop] = b_[nVar * np_current + 3];
    v[nop] = b_[nVar * np_current + 4];
    w[nop] = b_[nVar * np_current + 5];
    q[nop] = b_[nVar * np_current + 6];
    if (TrackParticleID)
      ParticleID[nop] = (unsigned long) b_[nVar * np_current + 7];
    np_current++;
    // these particles need further communication
    if (x[nop] < xStart_GC || x[nop] > xEnd_GC || y[nop] < yStart_GC || y[nop] > yEnd_GC || z[nop] < zStart_GC || z[nop] > zEnd_GC)
    //if (x[nop] < xstart || x[nop] > xend || y[nop] < ystart || y[nop] > yend || z[nop] < zstart || z[nop] > zend)
      rightDomain++;            // the particle is not in the domain
    nop++;
    if (nop > npmax) {
      cout << "Number of particles in the domain " << nop << " and maxpart = " << npmax << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
      return (-1);              // end the simulation because you dont have enough space on the array
    }
  }
  return (0);                   // everything was fine
}
/** Delete the a particle from the array and pack the the array, update the number of 
 * particles that are exiting
 * For deleting the particle from the array take the last particle and put it
 * in the position of the particle you want to delete
 * @param np = the index of the particle that must be deleted
 * @param nplast = the index of the last particle in the array
 */
void Particles3Dcomm::del_pack(long long np_current, long long *nplast) {
  x[np_current] = x[*nplast];
  y[np_current] = y[*nplast];
  z[np_current] = z[*nplast];
  u[np_current] = u[*nplast];
  v[np_current] = v[*nplast];
  w[np_current] = w[*nplast];
  q[np_current] = q[*nplast];
  if (TrackParticleID)
    ParticleID[np_current] = ParticleID[*nplast];
  npExit++;
  (*nplast)--;
}
/** method to calculate how many particles are out of right domain */
int Particles3Dcomm::isMessagingDone(VirtualTopology3D * ptVCT) {
  int result = 0;
  /*! mlmd: i need the communicator also */
  //result = reduceNumberParticles(rightDomain);
  result = reduceNumberParticles(rightDomain, ptVCT->getCommGrid());
  if (result > 0 && cVERBOSE && ptVCT->getCartesian_rank() == 0)
    cout << "Further Comunication: " << result << " particles not in the right domain" << endl;
  return (result);

}
/** calculate the maximum number exiting from this domain */
int Particles3Dcomm::maxNpExiting() {
  int maxNp = 0;
  if (npExitXright > maxNp)
    maxNp = npExitXright;
  if (npExitXleft > maxNp)
    maxNp = npExitXleft;
  if (npExitYright > maxNp)
    maxNp = npExitYright;
  if (npExitYleft > maxNp)
    maxNp = npExitYleft;
  if (npExitZright > maxNp)
    maxNp = npExitZright;
  if (npExitZleft > maxNp)
    maxNp = npExitZleft;
  return (maxNp);
}
/** return X-coordinate of particle array */
double *Particles3Dcomm::getXall()  const {
  return (x);
}
/** return Y-coordinate  of particle array */
double *Particles3Dcomm::getYall()  const {
  return (y);
}
/** return Z-coordinate  of particle array*/
double *Particles3Dcomm::getZall()  const {
  return (z);
}
/** get X-velocity of particle with label indexPart */
double *Particles3Dcomm::getUall()  const {
  return (u);
}
/** get Y-velocity of particle with label indexPart */
double *Particles3Dcomm::getVall()  const {
  return (v);
}
/**get Z-velocity of particle with label indexPart */
double *Particles3Dcomm::getWall()  const {
  return (w);
}
/**get ID of particle with label indexPart */
unsigned long *Particles3Dcomm::getParticleIDall()  const {
  return (ParticleID);
}
/**get charge of particle with label indexPart */
double *Particles3Dcomm::getQall()  const {
  return (q);
}
/** return X-coordinate of particle array as reference */
double *& Particles3Dcomm::getXref() {
  return (x);
}
/** return Y-coordinate  of particle array as reference */
double *& Particles3Dcomm::getYref() {
  return (y);
}
/** return Z-coordinate  of particle array as reference */
double *& Particles3Dcomm::getZref() {
  return (z);
}
/** get X-velocity of particle with label indexPart as reference */
double *& Particles3Dcomm::getUref() {
  return (u);
}
/** get Y-velocity of particle with label indexPart as reference */
double *& Particles3Dcomm::getVref() {
  return (v);
}
/**get Z-velocity of particle with label indexPart as reference */
double *& Particles3Dcomm::getWref() {
  return (w);
}
/**get charge of particle with label indexPart as reference */
double *& Particles3Dcomm::getQref() {
  return (q);
}
/** return X-coordinate of particle with index indexPart */
double Particles3Dcomm::getX(long long indexPart)  const {
  return (x[indexPart]);
}
/** return Y-coordinate  of particle with index indexPart */
double Particles3Dcomm::getY(long long indexPart)  const {
  return (y[indexPart]);
}
/** return Y-coordinate  of particle with index indexPart */
double Particles3Dcomm::getZ(long long indexPart)  const {
  return (z[indexPart]);
}
/** get u (X-velocity) of particle with label indexPart */
double Particles3Dcomm::getU(long long indexPart)  const {
  return (u[indexPart]);
}
/** get v (Y-velocity) of particle with label indexPart */
double Particles3Dcomm::getV(long long indexPart)  const {
  return (v[indexPart]);
}
/**get w (Z-velocity) of particle with label indexPart */
double Particles3Dcomm::getW(long long indexPart)  const {
  return (w[indexPart]);
}
/**get ID of particle with label indexPart */
unsigned long Particles3Dcomm::getParticleID(long long indexPart)  const {
  return (ParticleID[indexPart]);
}
/**get charge of particle with label indexPart */
double Particles3Dcomm::getQ(long long indexPart)  const {
  return (q[indexPart]);
}
/** return the number of particles */
long long Particles3Dcomm::getNOP()  const {
  return (nop);
}
/** return the Kinetic energy */
/*! mlmd: i need communicator also */
//double Particles3Dcomm::getKe() {
double Particles3Dcomm::getKe(MPI_Comm Comm) {
  double localKe = 0.0;
  double totalKe = 0.0;
  for (register long long i = 0; i < nop; i++)
    localKe += .5 * (q[i] / qom) * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
  MPI_Allreduce(&localKe, &totalKe, 1, MPI_DOUBLE, MPI_SUM, Comm);
  return (totalKe);
}
/** return the total momentum */
/*! mlmd: i need communicator also */
//double Particles3Dcomm::getP() {
double Particles3Dcomm::getP(MPI_Comm Comm) {
  double localP = 0.0;
  double totalP = 0.0;
  for (register long long i = 0; i < nop; i++)
    localP += (q[i] / qom) * sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
  MPI_Allreduce(&localP, &totalP, 1, MPI_DOUBLE, MPI_SUM, Comm);
  return (totalP);
}

/** return the highest kinetic energy */
/*! mlmd: i need communicator also */
//double Particles3Dcomm::getMaxVelocity() {
double Particles3Dcomm::getMaxVelocity(MPI_Comm Comm) { 
  double localVel = 0.0;
  double maxVel = 0.0;
  for (long long i = 0; i < nop; i++)
    localVel = max(localVel, sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]));
  MPI_Allreduce(&localVel, &maxVel, 1, MPI_DOUBLE, MPI_MAX, Comm);
  return (maxVel);
}


/** get energy spectrum */
/*! mlmd: i need the communicator also */
//unsigned long *Particles3Dcomm::getVelocityDistribution(int nBins, double maxVel) {
unsigned long *Particles3Dcomm::getVelocityDistribution(int nBins, double maxVel, MPI_Comm Comm) {
  unsigned long *f = new unsigned long[nBins];
  for (int i = 0; i < nBins; i++)
    f[i] = 0;
  double Vel = 0.0;
  double dv = maxVel / nBins;
  int bin = 0;
  for (long long i = 0; i < nop; i++) {
    Vel = sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
    bin = int (floor(Vel / dv));
    if (bin >= nBins)
      f[nBins - 1] += 1;
    else
      f[bin] += 1;
  }

  MPI_Allreduce(MPI_IN_PLACE, f, nBins, MPI_LONG_LONG, MPI_SUM, Comm);

  return f;
}


/** print particles info */
void Particles3Dcomm::Print(VirtualTopology3D * ptVCT) const {
  cout << endl;
  cout << "Number of Particles: " << nop << endl;
  cout << "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << "," << ptVCT->getCoordinates(2) << ")" << endl;
  cout << "Xin = " << xstart << "; Xfin = " << xend << endl;
  cout << "Yin = " << ystart << "; Yfin = " << yend << endl;
  cout << "Zin = " << zstart << "; Zfin = " << zend << endl;
  cout << "Number of species = " << ns << endl;
  for (long long i = 0; i < nop; i++)
    cout << "Particles #" << i << " x=" << x[i] << " y=" << y[i] << " z=" << z[i] << " u=" << u[i] << " v=" << v[i] << " w=" << w[i] << endl;
  cout << endl;
}
/** print just the number of particles */
void Particles3Dcomm::PrintNp(VirtualTopology3D * ptVCT)  const {
  cout << endl;
  cout << "Number of Particles of species " << ns << ": " << nop << endl;
  cout << "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << "," << ptVCT->getCoordinates(2) << ")" << endl;
  cout << endl;
}

void Particles3Dcomm::initWeightPBC(Grid * grid, VirtualTopology3D * vct){

  
  // Phase 1: as a child
  // Phase 1a: RG cores which need BCs flag themselves and identify the number of RG cells they will need BC for and the CG which will send
  // Phase 1b: they all send the msg to the highest ranking core
  // Phase 1c: highest ranking core assembles the msg and sends and handshake msg to each CG core

  // Phase 2: CG side
  // Phase 2a: each GC core receives the msg and stores appropriate variables

  RG_numPBCMessages= 0;
  int rank_local= vct->getCartesian_rank();  // on the grid communicator
  int HighestRank= XLEN*YLEN*ZLEN-1;
  MPI_Status status;
  int TAG_CG_RG= ns; // i need to put it here to be visible from both CG and RG

  /* phase 1: as a child */
  if (CommToParent_P != MPI_COMM_NULL) { // meaning: you are a child AND you want to receive PBC
    
    RGPBC_Info = new RGPBC_struct[MAX_RG_numPBCMessages];
    
    initWeightPBC_Phase1(grid, vct, RGPBC_Info, &RG_numPBCMessages);
    
    
    /**** check begins ****/
    // checks, aborting if checks fails
    int PG= vct->getParentGridNum();
    int localRank= vct->getCartesian_rank();
    double xmin, xmax, ymin, ymax, zmin, zmax;
    

    for (int i=0; i< RG_numPBCMessages; i++){
      int CG= RGPBC_Info[i].CG_core;
      double MsgLimsXMin= RGPBC_Info[i].CG_x_first;
      double MsgLimsXMax= RGPBC_Info[i].CG_x_first+ dx*(RGPBC_Info[i].np_x-1);

      double MsgLimsYMin= RGPBC_Info[i].CG_y_first;
      double MsgLimsYMax= RGPBC_Info[i].CG_y_first+ dy*(RGPBC_Info[i].np_y-1);

      double MsgLimsZMin= RGPBC_Info[i].CG_z_first;
      double MsgLimsZMax= RGPBC_Info[i].CG_z_first+ dz*(RGPBC_Info[i].np_z-1);

      grid->getParentLimits(vct, CG, &xmin, &xmax, &ymin, &ymax, &zmin, &zmax);
	
      /*cout <<"G" <<numGrid <<"R" << localRank << " Checking msg "<< i <<" after phase 1" << endl;
      cout <<"G" <<numGrid <<"R" << localRank  << " Msg from RG core " <<RGPBC_Info[i].RG_core <<" to CG " << RGPBC_Info[i].CG_core << "(PC comm)" << endl;
      cout <<"G" <<numGrid <<"R" << localRank << " Msg limits: X=" << MsgLimsXMin <<"-" << MsgLimsXMax <<" Y: " << MsgLimsYMin << "-" << MsgLimsYMax << " Z=" << MsgLimsZMin <<"-" <<MsgLimsZMax  << " RGPBC_Info[i].np_z " << RGPBC_Info[i].np_z << endl;
      cout <<"G" <<numGrid <<"R" << localRank << " Core limits: X=" << xmin << "-" << xmax << " Y= " <<ymin <<"-" << ymax <<" Z= " <<zmin <<"-" <<zmax << endl;*/
      if (MsgLimsXMin < xmin or MsgLimsXMax > xmax or MsgLimsYMin < ymin or MsgLimsYMax > ymax or MsgLimsZMin < zmin or MsgLimsZMax > zmax){
	cout <<"G" <<numGrid <<"R" << localRank <<" Msg " << i << " we have a problem in initWeightPBC, aborting ... " << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
      
    }
    /**** check ends ****/

    // now i can instantiate the vector for receiving PBC 
    // (but only if needed)
    if (RG_numPBCMessages>0){
      PRGMsg= newArr2(RepP_struct, RG_numPBCMessages, sizeRG_PBCMsg);
      nopPRGMsg= new int[RG_numPBCMessages];
      PRGMsgArrived= new bool [RG_numPBCMessages];
      PRGMsg_General= new RepP_struct[sizeRG_PBCMsg];

      // these are vectors needed to exchange repopulated particles inside the RG 
      if (rank_local==HighestRank){ 
	H_CRP_General= new CRP_struct[size_CRP];
	H_CRP_General_ptr= H_CRP_General; // for resize
	H_CRP_Msg= newArr2(CRP_struct, XLEN*YLEN*ZLEN, size_CRP);
	H_CRP_Msg_ptr= H_CRP_Msg;
	H_num_CRP_nop= new int[XLEN*YLEN*ZLEN];
	H_CRP_cores= new int[XLEN*YLEN*ZLEN];
	for (int ii=0; ii< XLEN*YLEN*ZLEN; ii++){
	  H_CRP_cores[ii]=0;
	  H_num_CRP_nop[ii]=0;
	}
	num_H_CRP_cores=0;
      }
      else{
	CRP_ToCoreH= new CRP_struct[size_CRP];
	CRP_ToCoreH_ptr= CRP_ToCoreH;
      }
    }
  } // end if (CommToParent_P != MPI_COMM_NULL) { 


  // you are a child grid which wants to receive BC;
  // all cores have to do this, even if RG_numPBCMessages== 0
  if (CommToParent_P != MPI_COMM_NULL) {
    int TAG_1b1c= ns;
    // phase 1b: all RG cores (but itself) send their msgs to highest ranking core 
    if (rank_local < HighestRank){
      /* send one message more; the last message has -1 in the RG_core   
	 to signal end of 'valid' messages */
      int dest= XLEN*YLEN*ZLEN-1;
      MPI_Send(RGPBC_Info, RG_numPBCMessages+1, MPI_RGPBC_struct, HighestRank, TAG_1b1c, vct->getComm());
    } // end phase 1b
   
    //Phase 1c: highest ranking core assembles the msg and sends and handshake msg to each CG core
    if (rank_local==HighestRank) {
      RGPBC_Info_LevelWide = new RGPBC_struct[MAX_RG_numPBCMessages_LevelWide];
      RG_numPBCMessages_LevelWide= 0;

      // first, highest ranking core copies its own msg into the level-wide structure
      for (int m=0; m< RG_numPBCMessages; m++){
	RGPBC_Info_LevelWide[m]= RGPBC_Info[m];
      }
      RG_numPBCMessages_LevelWide= RG_numPBCMessages;

      // then, it receives msgs from all the other cores in the grid
      RGPBC_struct * buffer_rcv;
      buffer_rcv= new RGPBC_struct[MAX_RG_numPBCMessages];
      
      int src;
      // i am receiveing XLEN*YLEN*ZLEN-1 --> HighestRank msg
      for (int c=0; c< HighestRank; c++){
	int recv=0;
	MPI_Recv(buffer_rcv, MAX_RG_numPBCMessages, MPI_RGPBC_struct, MPI_ANY_SOURCE, TAG_1b1c, vct->getComm(), &status);
	src= status.MPI_SOURCE;
	int recv_ThisMsg=0;
	while(buffer_rcv[recv].RG_core != -1){
	  RGPBC_Info_LevelWide[RG_numPBCMessages_LevelWide]= buffer_rcv[recv];
	  RG_numPBCMessages_LevelWide++;

	  recv++;

	  if (RG_numPBCMessages_LevelWide== MAX_RG_numPBCMessages_LevelWide){
	    cout << "initWeightPBC: MAX_RG_numPBCMessages_LevelWide is " << MAX_RG_numPBCMessages_LevelWide <<", but you need more..."<< endl << "Aborting...";
	    MPI_Abort(MPI_COMM_WORLD, -1);
	  }
	  
	  recv_ThisMsg++;
	} // end while
	if (recv_ThisMsg>0){
	  H_CRP_cores[src]= 1;
	}
	
      } // end for (int c=0; c< HighestRank; c++){

      num_H_CRP_cores=0;
      for (int i=0; i<XLEN*YLEN*ZLEN; i++){
	num_H_CRP_cores+= H_CRP_cores[i];
      }

      delete[]buffer_rcv;

      // now send the msg to the coarse grid
      int parentGrid= vct->getParentGridNum();
      int ParentCoreNum= vct->getXLEN(parentGrid) * vct->getYLEN(parentGrid) *vct->getZLEN(parentGrid);
      RGPBC_struct ** RGPBC_Info_ToCGCore= newArr2(RGPBC_struct, ParentCoreNum, MAX_RG_numPBCMessages);
      int * RG_numPBCMessages_ToCGCore = new int[ParentCoreNum];

      for (int c=0; c< ParentCoreNum; c++){
	RG_numPBCMessages_ToCGCore[c]=0;
      }
      
      // find the parent core where msg has to be sent
      for (int m=0; m< RG_numPBCMessages_LevelWide; m++){
	int where= RGPBC_Info_LevelWide[m].CG_core;
	RGPBC_Info_ToCGCore[where][RG_numPBCMessages_ToCGCore[where]]=RGPBC_Info_LevelWide[m];
	RG_numPBCMessages_ToCGCore[where]++;

	if (RG_numPBCMessages_ToCGCore[where] == MAX_RG_numPBCMessages){
	  cout << "Too many msgs to be sent for PBC from RG to CG, aborting..."<< endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);
	}
      } // end for (int m=0; m< RG_numPBCMessages_LevelWide; m++)
      
      // this is to stop the reading when CG core receives
      for (int cg=0; cg< ParentCoreNum; cg++){
	RGPBC_Info_ToCGCore[cg][RG_numPBCMessages_ToCGCore[cg]].RG_core= -1;
      } 
      // this is the send; send +1 in the # of msgs
      // rank as a child in the PC communicator
      int rankAsChild;
      MPI_Comm_rank(CommToParent_P, &rankAsChild);

      //cout << "ParentCoreNum= "<< ParentCoreNum<< endl;
      for (int cg=0; cg< ParentCoreNum; cg++){
	MPI_Send(&(RGPBC_Info_ToCGCore[cg][0]), RG_numPBCMessages_ToCGCore[cg]+1, MPI_RGPBC_struct, cg, TAG_CG_RG, CommToParent_P);
	//cout << "R " << rankAsChild << " on the PC communicator has just sent a msg to core " << cg << endl; 
      }

      delete[]RGPBC_Info_LevelWide;
      
      delete[]RGPBC_Info_ToCGCore;
      delete[]RG_numPBCMessages_ToCGCore;
      
    } // end if if (rank_local==HighestRank) 
  } //  if (CommToParent_P != MPI_COMM_NULL)

  // Phase 2: CG receives the handshake and reacts
  // NB: for all children, I have to check if they want PBC
  // if no, do not allocate stuff

  if (numChildren>0){
    CG_Info= newArr2(RGPBC_struct, numChildren, MAX_RG_numPBCMessages);
    CG_numPBCMessages= new int[numChildren];
  }

  for (int ch=0; ch <  numChildren; ch++){
    CG_numPBCMessages[ch]=0;

    RGPBC_struct * CG_buffer=  new RGPBC_struct[MAX_RG_numPBCMessages];
    if (CommToChild_P[ch] != MPI_COMM_NULL){ // it may happen that a particular child does not want PC
      
      int ChildHighestRank;
      MPI_Comm_size(CommToChild_P[ch], &ChildHighestRank); ChildHighestRank--;
      int PCrank= Rank_CommToChildren_P[ch];
      
      //cout << "Grid " << numGrid << ", local rank on PC comm " << PCrank << " is trying to receive from ch " << ch << endl;

      // receive one msg from HighestRank on the PC comm

      // to put ChildHighestRank as source is just a precaution, MPI_ANY_Source should work
      MPI_Recv(CG_buffer, MAX_RG_numPBCMessages, MPI_RGPBC_struct, ChildHighestRank, TAG_CG_RG, CommToChild_P[ch], &status);

      //cout << "Grid " << numGrid << ", local rank on PC comm " << PCrank << " has recevied from ch " << ch << " tag " << TAG_CG_RG << endl;

      
      // process update the msg structure
      while (CG_buffer[CG_numPBCMessages[ch]].RG_core!= -1){
	//cout << "INSIDE WHILE: ch is " << ch << " CG_numPBCMessages[ch] is " << CG_numPBCMessages[ch] << " CG_buffer[CG_numPBCMessages[ch]].RG_core is " <<CG_buffer[CG_numPBCMessages[ch]].RG_core << endl;
	CG_Info[ch][CG_numPBCMessages[ch]]= CG_buffer[CG_numPBCMessages[ch]];
	CG_numPBCMessages[ch]++;
      }

    } // end if (vct->getCommToChild_P()!= MPI_COMM_NULL)

    delete[]CG_buffer;
    
  } // end for (int ch=0; ch< numChildren; ch++)

  // coarse grids allocate vectors
  if (numChildren>0){
    MaxNumMsg=0;
    for (int i=0; i< numChildren; i++){
      if (CG_numPBCMessages[i]> MaxNumMsg) MaxNumMsg= CG_numPBCMessages[i];
    }
    
    if (MaxNumMsg >0){ // otherwise, i am not sending PBC and i don't need to allocate
      PCGMsg= newArr3(RepP_struct, numChildren, MaxNumMsg, sizeCG_PBCMsg);
      PCGMsg_ptr= PCGMsg; // for the resize
      nopPCGMsg = newArr2(int, numChildren, MaxNumMsg);
      
      //cout << "Grid " << numGrid << " core " << vct->getCartesian_rank() << " nopPCGMsg has sizes " << numChildren << " x " <<MaxNumMsg << endl;  
    }
  } // end if (numChildren>0){


  // some prints for checks
  
  /*if (numGrid >0){
    cout << "Grid " <<numGrid << " R " << rank_local << ": RG_numPBCMessages: " << RG_numPBCMessages <<endl;
    if (rank_local== HighestRank){
      cout << "H_CRP_cores is" << endl;
      for (int i=0; i< XLEN*YLEN*ZLEN; i++){
	cout <<"i: " << i << ", H_CRP_cores[i]: " << H_CRP_cores[i] << endl;
      }
      cout << "Grid " << numGrid << ", num_H_CRP_cores: " << num_H_CRP_cores <<endl; 
    }
    }*/

  // RG prints its messages
  /*if (numGrid >0){
    cout << "Grid " << numGrid <<" R " << vct->getCartesian_rank() << " total msg: " <<RG_numPBCMessages<< endl;
    for (int m=0; m< RG_numPBCMessages; m++){
      cout << "msg " << m << " of " << RG_numPBCMessages << endl;
      //cout << "RG points: " << RGPBC_Info[m].np_x << "/ " nxn <<" - " << RGPBC_Info[m].np_y <<"/ " nyn << " - " << RGPBC_Info[m].np_z << "/ "<<nzn << endl;
      cout << "CG_x_first-Ox --> CG_x_first-Ox + (RGPBC_Info[m].np_x-1)*dx " << RGPBC_Info[m].CG_x_first- grid->getOx() <<"-->" << RGPBC_Info[m].CG_x_first- grid->getOx() + (RGPBC_Info[m].np_x-1)*dx <<" -2dx "  << -2*dx <<" Lx+2dx " << Lx+2*dx << " Ox " << grid->getOx()  << endl;
      cout << "CG_y_first-Oy --> CG_y_first-Oy + (RGPBC_Info[m].np_y-1)*dy " << RGPBC_Info[m].CG_y_first- grid->getOy() <<"-->" << RGPBC_Info[m].CG_y_first- grid->getOy() + (RGPBC_Info[m].np_y-1)*dy << " -2dy " << -2*dy << " Ly+2dy " << Ly+2*dy << " Oy " << grid->getOy() << endl;
      cout << "CG_z_first-Oz --> CG_z_first-Oz + (RGPBC_Info[m].np_z-1)*dz " << RGPBC_Info[m].CG_z_first- grid->getOz() <<"-->" << RGPBC_Info[m].CG_z_first- grid->getOz() + (RGPBC_Info[m].np_z-1)*dz << " -2dz " << -2*dz << " Lz+2dz " << Lz+2*dz << " Oz " << grid->getOz() << endl;
      
    }
    }*/
  
  // create communicators involving only the cores of the grid involved in particle BC communication
  if (true){
  MPI_Barrier(vct->getComm());

  // as a child, creation of a communicator of cores exchanging PBC
  int color; 
  int key=0; // so rank is assigned based on the 'old' rank
  if (CommToParent_P != MPI_COMM_NULL and RG_numPBCMessages>0) color=1; else color= MPI_UNDEFINED;
  MPI_Comm_split(vct->getComm(), color, key, &COMM_RG_PBCSubset_P);

  MPI_Barrier(vct->getComm()); // not sure it's needed
  if (numChildren >0){
    COMM_CG_PBCSubset_P= new MPI_Comm[numChildren];
    CGSide_CGLeader_PBCSubset= new int[numChildren];
  }
  // as a parent, creation of a communicator of cores sending PBC
  for (int ch=0; ch< numChildren; ch++){
    if (vct->getCommToChild_P(ch, ns) != MPI_COMM_NULL and CG_numPBCMessages[ch]> 0) color=1; else color= MPI_UNDEFINED;
    MPI_Comm_split(vct->getComm(), color, key, COMM_CG_PBCSubset_P+ch);
    MPI_Barrier(vct->getComm()); // not sure it's needed
  }

  /* the CG core EXCHANGING PBC which will be leader in the communication to the RG throught the 'usual' Parent-Child communicator
     (which most probably will pass through another CG core) */
  for (int ch=0; ch< numChildren; ch++){
    int rankPC= vct->getRank_CommToChildren_P(ch, ns);
    if (COMM_CG_PBCSubset_P[ch] != MPI_COMM_NULL)
      MPI_Allreduce(&rankPC, CGSide_CGLeader_PBCSubset+ ch, 1, MPI_INT, MPI_MIN, COMM_CG_PBCSubset_P[ch]);
  }

  // build the list of RG cores to notify of the resize
  if (AllowPMsgResize==true){
    // build the list of RG cores to notify of the resize; each core has to receive only one msg
    // hence this

    bool Needed= false; // am i the spokeperson to any child?
    int Max=0;
    for (int ch=0; ch< numChildren; ch++){
      if (COMM_CG_PBCSubset_P[ch] == MPI_COMM_NULL) continue;
      int S; 
      MPI_Comm_size(vct->getCommToChild_P(ch, ns), &S);
      if (S> Max) Max=S;
      if (CGSide_CGLeader_PBCSubset[ch]== vct->getRank_CommToChildren_P(ch, ns))
	Needed = Needed or true;
    }
    if (Needed){
      numRcv= new int[numChildren];
      RcvList= newArr2(int, numChildren, Max);
    }

    for (int ch=0; ch< numChildren; ch++){
      if (COMM_CG_PBCSubset_P[ch] == MPI_COMM_NULL) continue;
      
      int size;
      MPI_Comm_size(vct->getCommToChild_P(ch, ns), &size);
      bool *tmp= new bool[size];
      bool *tmp_2= new bool[size];
      for (int ss=0; ss< size; ss++){tmp[ss]=false; tmp_2[ss]=false; }

      // list of RG cores this CG cores communicates with
      for (int m=0; m< CG_numPBCMessages[ch]; m++){
	tmp[CG_Info[ch][m].RG_core]= true;
      }
      // now i have to share it
      MPI_Allreduce(tmp, tmp_2, size, MPI::BOOL, MPI_LOR, COMM_CG_PBCSubset_P[ch]);
	
      // i am the spokeperson
      if (vct->getRank_CommToChildren_P(ch, ns) == CGSide_CGLeader_PBCSubset[ch]){
	cout << "Grid " << numGrid << " R "<< vct->getCartesian_rank() << "is spokeperson " << endl;
	numRcv[ch]=0;
	for (int i=0; i< size; i++){
	  if (tmp_2[i] == true){ RcvList[ch][numRcv[ch]]=i;  numRcv[ch]++;}
	} // end for (int i=0; i< size; i++){

      }// end if (vct->getRank_CommToChildren_P(nc, ns) == CGSide_CGLeader_PBCSubset[ch]){
    } // end for (int ch=0; ch< numChildren; ch++){
  } // end if (AllowPMsgResize==true)
  
}
  
}

/** commit the RGPBC structure for initial handshake between coarse and refined grids **/
void Particles3Dcomm::MPI_RGPBC_struct_commit(){

  /* struct RGPBC_struct {  // when changing this, change MPI_RGPBC_struct_commit also   
    int np_x;
    int np_y;
    int np_z;
    double CG_x_first;
    double CG_y_first;
    double CG_z_first;
    int CG_core;
    int RG_core;
    int MsgID;
    };*/


  RGPBC_struct *a;

  MPI_Datatype type[9]={MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT};
  int blocklen[9]={1,1,1,1,1,1,1,1,1};

  // displacement in bytes      
  MPI_Aint disp[9];

  // np_*      
  disp[0]= (MPI_Aint) &(a->np_x) - (MPI_Aint)a ;
  disp[1]= (MPI_Aint) &(a->np_y) - (MPI_Aint)a ;
  disp[2]= (MPI_Aint) &(a->np_z) - (MPI_Aint)a ;

  // CG_*_first                                            
  disp[3]= (MPI_Aint) &(a->CG_x_first) - (MPI_Aint)a ;
  disp[4]= (MPI_Aint) &(a->CG_y_first) - (MPI_Aint)a ;
  disp[5]= (MPI_Aint) &(a->CG_z_first) - (MPI_Aint)a ;

  // the cores                                                                                          
  disp[6]= (MPI_Aint) &(a->CG_core) - (MPI_Aint)a ;
  disp[7]= (MPI_Aint) &(a->RG_core) - (MPI_Aint)a ;

  // the msg id           
  disp[8]= (MPI_Aint) &(a->MsgID) - (MPI_Aint)a ;

  MPI_Type_create_struct(9, blocklen, disp, type, &MPI_RGPBC_struct);
  MPI_Type_commit(&MPI_RGPBC_struct); 

}

/* commit the structure for the particle CG/RG exchange as MPI_Datatype */
void Particles3Dcomm::MPI_RepP_struct_commit(){

  /*struct RepP_struct{
    // position    
    double x;
    double y;
    double z;
    // velocity    
    double u;
    double v;
    double q;
    // charge    
    double q;
    // ID               
    unsigned long ID;
    }; */

  RepP_struct *a;

  MPI_Datatype type[8]={MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UNSIGNED_LONG};

  int blocklen[8]={1,1,1,1,1,1,1,1};

  // displacement in bytes      
  MPI_Aint disp[8];

  // position
  disp[0]= (MPI_Aint) &(a->x) - (MPI_Aint)a ;
  disp[1]= (MPI_Aint) &(a->y) - (MPI_Aint)a ;
  disp[2]= (MPI_Aint) &(a->z) - (MPI_Aint)a ;

  // velocity
  disp[3]= (MPI_Aint) &(a->u) - (MPI_Aint)a ;
  disp[4]= (MPI_Aint) &(a->v) - (MPI_Aint)a ;
  disp[5]= (MPI_Aint) &(a->w) - (MPI_Aint)a ;

  // charge
  disp[6]= (MPI_Aint) &(a->q) - (MPI_Aint)a ;

  // ID
  disp[7]= (MPI_Aint) &(a->ID) - (MPI_Aint)a ;


  MPI_Type_create_struct(8, blocklen, disp, type, &MPI_RepP_struct);
  MPI_Type_commit(&MPI_RepP_struct); 

}

/* commit the structure for repopulated particle exchange within the RG */
void Particles3Dcomm::MPI_CRP_struct_commit(){

  /* struct CRP_struct{
  // position         
  double x;
  double y;
  double z;
  // velocity                    
  double u;
  double v;
  double w;
  // charge               
  double q;
  // ID                
  unsigned long ID;
  // rank of destination- on the local RG communicator              
  int Destination;
}; */

  CRP_struct *a;

  MPI_Datatype type[9]={MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UNSIGNED_LONG, MPI_INT};

  int blocklen[9]={1,1,1,1,1,1,1,1,1};

  // displacement in bytes      
  MPI_Aint disp[9];

  // position
  disp[0]= (MPI_Aint) &(a->x) - (MPI_Aint)a ;
  disp[1]= (MPI_Aint) &(a->y) - (MPI_Aint)a ;
  disp[2]= (MPI_Aint) &(a->z) - (MPI_Aint)a ;

  // velocity
  disp[3]= (MPI_Aint) &(a->u) - (MPI_Aint)a ;
  disp[4]= (MPI_Aint) &(a->v) - (MPI_Aint)a ;
  disp[5]= (MPI_Aint) &(a->w) - (MPI_Aint)a ;

  // charge
  disp[6]= (MPI_Aint) &(a->q) - (MPI_Aint)a ;

  // ID
  disp[7]= (MPI_Aint) &(a->ID) - (MPI_Aint)a ;

  // Destination
  disp[8]= (MPI_Aint) &(a->Destination) - (MPI_Aint)a ;

  MPI_Type_create_struct(9, blocklen, disp, type, &MPI_CRP_struct);
  MPI_Type_commit(&MPI_CRP_struct);
  
}


/* Phase 1: RG cores build their side of the map for PBC */
void Particles3Dcomm::initWeightPBC_Phase1(Grid *grid, VirtualTopology3D * vct, RGPBC_struct *RGPBC_Info, int *RG_numPBCMessages){

  // NB: when build the PBC msg on the CG side, I have to pay attention not to replicate the info; that was not a problem with fields, but it is now
  
  int SW_rank=vct->getSystemWide_rank();

  string FACE;

  bool DIR_0=true;
  bool DIR_1=true;
  bool DIR_2=true;

  int MS= nxn; if (nyn>MS) MS= nyn; if (nzn>MS) MS= nzn;

  int i_s, i_e;
  int j_s, j_e;
  int k_s, k_e;
  // careful when copying from the initWeight in fields: getXXX_neighbor_P !!!
  // (so i can have different periodicities in fields and particles)

  // this is the bottom face
  if (vct->getCoordinates(2)==0 && vct->getZleft_neighbor_P()== MPI_PROC_NULL && DIR_2){

    FACE= "BOTTOM";
    i_s= 0-1; i_e= nxn-1+1;
    j_s= 0-1; j_e= nyn-1+1;
    k_s= PRA_ZLeft_Start-1; k_e= PRA_ZLeft_End+1;

    /*cout << "nxn " << nxn << " nyn " << nyn << " nzn " << nzn << endl;
  cout << "PRA_XLeft_Start " << PRA_XLeft_Start << " PRA_XLeft_End " << PRA_XLeft_End << " PRA_XRight_End " << PRA_XRight_End << " PRA_XRight_Start " << PRA_XRight_Start <<endl; 
  cout << "PRA_YLeft_Start " << PRA_YLeft_Start << " PRA_YLeft_End " << PRA_YLeft_End << " PRA_YRight_End " << PRA_YRight_End << " PRA_YRight_Start " << PRA_YRight_Start <<endl; 
  cout << "PRA_ZLeft_Start " << PRA_ZLeft_Start << " PRA_ZLeft_End " << PRA_ZLeft_End << " PRA_ZRight_End " << PRA_ZRight_End << " PRA_ZRight_Start " << PRA_ZRight_Start <<endl; 
  cout <<"FACE " << FACE << " i_s " << i_s << " i_e " << i_e << " j_s " << j_s << " j_e " << j_e << " k_s " << k_s << " k_e " <<k_e <<endl;*/
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    //cout <<"FACE " << FACE << endl;
  } // end bottom face
  
  // this is the top face
  if (vct->getCoordinates(2) ==ZLEN-1 && vct->getZright_neighbor_P() == MPI_PROC_NULL && DIR_2){
    
    FACE= "TOP";
    i_s=0-1; i_e= nxn-1+1;
    j_s=0-1; j_e= nyn-1+1;
    k_s= PRA_ZRight_Start-1; k_e= PRA_ZRight_End +1;

    /*cout << "nxn " << nxn << " nyn " << nyn << " nzn " << nzn << endl;
  cout << "PRA_XLeft_Start " << PRA_XLeft_Start << " PRA_XLeft_End " << PRA_XLeft_End << " PRA_XRight_End " << PRA_XRight_End << " PRA_XRight_Start " << PRA_XRight_Start <<endl; 
  cout << "PRA_YLeft_Start " << PRA_YLeft_Start << " PRA_YLeft_End " << PRA_YLeft_End << " PRA_YRight_End " << PRA_YRight_End << " PRA_YRight_Start " << PRA_YRight_Start <<endl; 
  cout << "PRA_ZLeft_Start " << PRA_ZLeft_Start << " PRA_ZLeft_End " << PRA_ZLeft_End << " PRA_ZRight_End " << PRA_ZRight_End << " PRA_ZRight_Start " << PRA_ZRight_Start <<endl; 
  cout <<"FACE " << FACE << " i_s " << i_s <<" i_e "<< i_e << " j_s " << j_s << " j_e " << j_e << " k_s " << k_s <<" k_e "<<k_e <<endl;*/
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    //cout<<"FACE " << FACE << endl;
  } // end top face

  // this is the left face
  if (vct->getCoordinates(0) ==0  && vct->getXleft_neighbor_P() == MPI_PROC_NULL && DIR_0){

    FACE= "LEFT";
    i_s= PRA_XLeft_Start-1; i_e= PRA_XLeft_End +1;
    j_s= 0-1; j_e= nyn-1+1;
    k_s= 0-1; k_e= nzn-1+1;

    /*cout << "nxn " << nxn << " nyn " << nyn << " nzn " << nzn << endl;
  cout << "PRA_XLeft_Start " << PRA_XLeft_Start << " PRA_XLeft_End " << PRA_XLeft_End << " PRA_XRight_End " << PRA_XRight_End << " PRA_XRight_Start " << PRA_XRight_Start <<endl; 
  cout << "PRA_YLeft_Start " << PRA_YLeft_Start << " PRA_YLeft_End " << PRA_YLeft_End << " PRA_YRight_End " << PRA_YRight_End << " PRA_YRight_Start " << PRA_YRight_Start <<endl; 
  cout << "PRA_ZLeft_Start " << PRA_ZLeft_Start << " PRA_ZLeft_End " << PRA_ZLeft_End << " PRA_ZRight_End " << PRA_ZRight_End << " PRA_ZRight_Start " << PRA_ZRight_Start <<endl; 
  cout <<"FACE " << FACE << " i_s " << i_s <<" i_e "<< i_e << " j_s " << j_s << " j_e " << j_e << " k_s " << k_s <<" k_e "<<k_e <<endl;*/
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    //cout<<"FACE " << FACE << endl;
  } // end left face

  // this is the right face
  if (vct->getCoordinates(0) ==XLEN-1 && vct->getXright_neighbor_P() == MPI_PROC_NULL && DIR_0){

    FACE= "RIGHT";
    i_s= PRA_XRight_Start-1; i_e= PRA_XRight_End+1;
    j_s= 0-1; j_e= nyn-1+1;
    k_s= 0-1; k_e= nzn-1+1;

    /*cout << "nxn " << nxn << " nyn " << nyn << " nzn " << nzn << endl;
  cout << "PRA_XLeft_Start " << PRA_XLeft_Start << " PRA_XLeft_End " << PRA_XLeft_End << " PRA_XRight_End " << PRA_XRight_End << " PRA_XRight_Start " << PRA_XRight_Start <<endl; 
  cout << "PRA_YLeft_Start " << PRA_YLeft_Start << " PRA_YLeft_End " << PRA_YLeft_End << " PRA_YRight_End " << PRA_YRight_End << " PRA_YRight_Start " << PRA_YRight_Start <<endl; 
  cout << "PRA_ZLeft_Start " << PRA_ZLeft_Start << " PRA_ZLeft_End " << PRA_ZLeft_End << " PRA_ZRight_End " << PRA_ZRight_End << " PRA_ZRight_Start " << PRA_ZRight_Start <<endl; 
  cout <<"FACE " << FACE << " i_s " << i_s <<" i_e "<< i_e << " j_s " << j_s << " j_e " << j_e << " k_s " << k_s <<" k_e "<<k_e <<endl;*/
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    //cout<<"FACE " << FACE << endl;
  } // end right face
  
  // this is the front face
  if (vct->getCoordinates(1) ==0 && vct->getYleft_neighbor_P() == MPI_PROC_NULL && DIR_1){ 

    FACE= "FRONT";
    i_s= 0-1; i_e= nxn-1+1;
    j_s= PRA_YLeft_Start-1; j_e= PRA_YLeft_End+1;
    k_s= 0-1; k_e= nzn-1+1;

    /*cout << "nxn " << nxn << " nyn " << nyn << " nzn " << nzn << endl;
  cout << "PRA_XLeft_Start " << PRA_XLeft_Start << " PRA_XLeft_End " << PRA_XLeft_End << " PRA_XRight_End " << PRA_XRight_End << " PRA_XRight_Start " << PRA_XRight_Start <<endl; 
  cout << "PRA_YLeft_Start " << PRA_YLeft_Start << " PRA_YLeft_End " << PRA_YLeft_End << " PRA_YRight_End " << PRA_YRight_End << " PRA_YRight_Start " << PRA_YRight_Start <<endl; 
  cout << "PRA_ZLeft_Start " << PRA_ZLeft_Start << " PRA_ZLeft_End " << PRA_ZLeft_End << " PRA_ZRight_End " << PRA_ZRight_End << " PRA_ZRight_Start " << PRA_ZRight_Start <<endl; 
  cout <<"FACE " << FACE << " i_s " << i_s <<" i_e "<< i_e << " j_s " << j_s << " j_e " << j_e << " k_s " << k_s <<" k_e "<<k_e <<endl;*/
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    //cout<<"FACE " << FACE << endl;
  } // end front face

  // this is the back face
  if (vct->getCoordinates(1) == YLEN-1 && vct->getYright_neighbor_P() == MPI_PROC_NULL && DIR_1 ){

    FACE= "BACK";
    i_s= 0-1; i_e= nxn-1+1;
    j_s= PRA_YRight_Start-1; j_e= PRA_YRight_End+1;
    k_s= 0-1; k_e= nzn-1+1;
    
    /*cout << "nxn " << nxn << " nyn " << nyn << " nzn " << nzn << endl;
  cout << "PRA_XLeft_Start " << PRA_XLeft_Start << " PRA_XLeft_End " << PRA_XLeft_End << " PRA_XRight_End " << PRA_XRight_End << " PRA_XRight_Start " << PRA_XRight_Start <<endl; 
  cout << "PRA_YLeft_Start " << PRA_YLeft_Start << " PRA_YLeft_End " << PRA_YLeft_End << " PRA_YRight_End " << PRA_YRight_End << " PRA_YRight_Start " << PRA_YRight_Start <<endl; 
  cout << "PRA_ZLeft_Start " << PRA_ZLeft_Start << " PRA_ZLeft_End " << PRA_ZLeft_End << " PRA_ZRight_End " << PRA_ZRight_End << " PRA_ZRight_Start " << PRA_ZRight_Start <<endl; 
  cout <<"FACE " << FACE << " i_s " << i_s <<" i_e "<< i_e << " j_s " << j_s << " j_e " << j_e << " k_s " << k_s <<" k_e "<<k_e <<endl;*/
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    //cout<<"FACE " << FACE << endl;
  } // end back face
  
  /* for further use, i need to set the RG_core field of the first unused slot to -1  
     but DO NOT MODIFY THE NUMBER OF MSGs;                             
     I will just send a +1 */

  //cout << "R" <<SW_rank <<" RG_numPBCMessages= " <<*RG_numPBCMessages <<endl;
  RGPBC_Info[*RG_numPBCMessages].RG_core= -1;
  RGPBC_Info[*RG_numPBCMessages].CG_core= -1;
 

}

void Particles3Dcomm::Explore3DAndCommit(Grid *grid, int i_s, int i_e, int j_s, int j_e, int k_s, int k_e, int *numMsg, int *MaxSizeMsg, VirtualTopology3D * vct ){

  // policy:
  // explore Z dir
  // for on the number of cores found there: explore Y dir
  // for on the number of cores found there: explore X dir
  // finally, commit  NB: all faces should have the same c

  int MS= nxn; if (nyn>MS) MS= nyn; if (nzn>MS) MS= nzn;
  int rank_CTP= vct->getRank_CommToParent();
  /*******************************************************************/
  // DIR1: starting point, in CG coordinates, per core                       
  double *Dir1_SPXperC= new double[MS];
  double *Dir1_SPYperC= new double[MS];
  double *Dir1_SPZperC= new double[MS];
  // DIR1: number of Refined grid point in this direction, per core     
  int *Dir1_NPperC= new int[MS];  // this does not need to be this big, but anyway   
  // DIR1: core ranks in the CommToParent communicator              
  int *Dir1_rank= new int [MS]; // this does not need to be this big, but anyway      
  int *Dir1_IndexFirstPointperC= new int [MS];
  int Dir1_Ncores=0;

  // DIR2: starting point, in CG coordinates, per core                                  
  double *Dir2_SPXperC= new double[MS];
  double *Dir2_SPYperC= new double[MS];
  double *Dir2_SPZperC= new double[MS];
  // DIR2: number of Refined grid point in this direction, per core     
  int *Dir2_NPperC= new int[MS];  // this does not need to be this big, but anyway   
  // DIR2: core ranks in the CommToParent communicator              
  int *Dir2_rank= new int [MS]; // this does not need to be this big, but anyway      
  int *Dir2_IndexFirstPointperC= new int [MS];
  int Dir2_Ncores=0;

  // DIR3: starting point, in CG coordinates, per core                     
  double *Dir3_SPXperC= new double[MS];
  double *Dir3_SPYperC= new double[MS];
  double *Dir3_SPZperC= new double[MS];
  // DIR3: number of Refined grid point in this direction, per core     
  int *Dir3_NPperC= new int[MS];  // this does not need to be this big, but anyway   
  // DIR3: core ranks in the CommToParent communicator              
  int *Dir3_rank= new int [MS]; // this does not need to be this big, but anyway      
  int *Dir3_IndexFirstPointperC= new int [MS];
  int Dir3_Ncores=0;
  /*******************************************************************/


  // variables which i do not need to update for particle BCs
  string FACE= "nn";

  // for debug 
  int PG= vct->getParentGridNum();
  // end for debug

  double CG_C_i_s= grid->getXN_P(i_s, j_s, k_s);
  double CG_C_j_s= grid->getYN_P(i_s, j_s, k_s);
  double CG_C_k_s= grid->getZN_P(i_s, j_s, k_s);
  double CG_C_i_e= grid->getXN_P(i_e, j_e, k_e);
  double CG_C_j_e= grid->getYN_P(i_e, j_e, k_e);
  double CG_C_k_e= grid->getZN_P(i_e, j_e, k_e);

  /*cout << "Inside explore and commit: " << endl;
  cout << "i_s= " <<i_s <<" -> " << CG_C_i_s << endl;
  cout << "i_e= " <<i_e <<" -> " << CG_C_i_e <<endl;

  cout << "j_s= " <<j_s <<" -> " << CG_C_j_s << endl;
  cout << "j_e= " <<j_e <<" -> " << CG_C_j_e <<endl;

  cout << "k_s= " <<k_s <<" -> " << CG_C_k_s << endl;
  cout << "k_e= " <<k_e <<" -> " << CG_C_k_e <<endl;*/

  // Z dir / Dir 1
  grid->RGBCExploreDirection(vct, FACE, 2, k_s, k_e, i_s, j_s, Dir1_SPXperC, Dir1_SPYperC, Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);

  
  /*// debug
  cout << "G"<< numGrid <<"R" << rank_CTP <<": Z dir, " << Dir1_Ncores << " cores " << endl;
  for (int ii=0; ii< Dir1_Ncores; ii++){
    cout << "G"<< numGrid <<"R" << rank_CTP << " " << ii <<": " << Dir1_rank[ii] << " first coord: " <<Dir1_SPZperC[ii] << ", n points: " << Dir1_NPperC[ii] << " which means up to " << Dir1_SPZperC[ii]+(Dir1_NPperC[ii] -1)*dz << endl;
    cout << "In terms of RG coordinates: " << Dir1_SPZperC[ii] - grid->getOz() << " - " << Dir1_SPZperC[ii]+(Dir1_NPperC[ii] -1)*dz  - grid->getOz() << endl;
  }
  // end debug*/
  
  for (int n=0; n<Dir1_Ncores; n++){ // it will find again the core in Dir 1, but it will also explore Dir 2 
    // Y dir / Dir 2
    grid->RGBCExploreDirection(vct, FACE, 1, j_s, j_e, i_s, Dir1_IndexFirstPointperC[n],  Dir2_SPXperC, Dir2_SPYperC, Dir2_SPZperC, Dir2_NPperC, Dir2_rank, &Dir2_Ncores, Dir2_IndexFirstPointperC); 
    /* // debug
    cout << "G"<< numGrid <<"R" << rank_CTP <<": Y dir, " << Dir2_Ncores << " cores " << endl;
    for (int ii=0; ii< Dir2_Ncores; ii++){
      cout << "G"<< numGrid <<"R" << rank_CTP << " " << ii <<": " << Dir2_rank[ii] << " first coord: " <<Dir2_SPYperC[ii] << ", n points: " << Dir2_NPperC[ii] << " which means up to " << Dir2_SPYperC[ii]+(Dir2_NPperC[ii] -1)*dy << endl;
    }
    // end debug */
    
    for (int m=0; m< Dir2_Ncores; m++){ //it will find again the core in Dir 1, but it will also explore Dir 2 
      // X dir / Dir 3
      grid->RGBCExploreDirection(vct, FACE, 0, i_s, i_e, Dir2_IndexFirstPointperC[m],  Dir1_IndexFirstPointperC[n], Dir3_SPXperC, Dir3_SPYperC, Dir3_SPZperC, Dir3_NPperC, Dir3_rank, &Dir3_Ncores, Dir3_IndexFirstPointperC);

      /* // debug
      cout << "G"<< numGrid <<"R" << rank_CTP <<": X dir, " << Dir3_Ncores << " cores " << endl;
      for (int ii=0; ii< Dir3_Ncores; ii++){
	cout << "G"<< numGrid <<"R" << rank_CTP << " " << ii <<": " << Dir3_rank[ii] << " first coord: " <<Dir3_SPXperC[ii] << ", n points: " << Dir3_NPperC[ii] << " which means up to " << Dir3_SPXperC[ii]+(Dir3_NPperC[ii] -1)*dx << endl;
      }
      // end debug*/

	// now commit msg
	for (int NN=0; NN< Dir3_Ncores; NN++){
	  Assign_RGBC_struct_Values(RGPBC_Info + (*numMsg), Dir3_NPperC[NN], Dir2_NPperC[m], Dir1_NPperC[n], Dir3_SPXperC[NN], Dir3_SPYperC[NN], Dir3_SPZperC[NN], Dir3_rank[NN], rank_CTP, *numMsg);
	  (*numMsg)++;

	  int tmp= Dir3_NPperC[NN]*Dir2_NPperC[m]*Dir1_NPperC[n];
	  if (tmp > *MaxSizeMsg) (*MaxSizeMsg)= tmp;
	  
	} // end for (int NN=0; NN< Dir3_Ncores; NN++){
    } // end  for (int m=0; m< Dir2_Ncores; m++){
  } // end for (int n=0; n<Dir1_Ncores; n++){

  delete[]Dir1_SPXperC;
  delete[]Dir1_SPYperC;
  delete[]Dir1_SPZperC;
  delete[]Dir1_NPperC;
  delete[]Dir1_rank;
  delete[]Dir1_IndexFirstPointperC;

  delete[]Dir2_SPXperC;
  delete[]Dir2_SPYperC;
  delete[]Dir2_SPZperC;
  delete[]Dir2_NPperC;
  delete[]Dir2_rank;
  delete[]Dir2_IndexFirstPointperC;

  delete[]Dir3_SPXperC;
  delete[]Dir3_SPYperC;
  delete[]Dir3_SPZperC;
  delete[]Dir3_NPperC;
  delete[]Dir3_rank;
  delete[]Dir3_IndexFirstPointperC;
}
/* add one handshake msg to the list */
void Particles3Dcomm::Assign_RGBC_struct_Values(RGPBC_struct *s, int np_x_tmp, int np_y_tmp, int np_z_tmp, double CG_x_first_tmp, double CG_y_first_tmp, double CG_z_first_tmp, int CG_core_tmp, int RG_core_tmp, int MsgID_tmp ) {

  s->np_x = np_x_tmp;
  s->np_y = np_y_tmp;
  s->np_z = np_z_tmp;

  s->CG_x_first = CG_x_first_tmp;
  s->CG_y_first = CG_y_first_tmp;
  s->CG_z_first = CG_z_first_tmp;

  s->CG_core = CG_core_tmp;
  
  s->RG_core = RG_core_tmp;

  s->MsgID = MsgID_tmp;

  return;
}

/* build and send particle BC msg -- CG to RG */
void Particles3Dcomm::SendPBC(Grid* grid, VirtualTopology3D * vct){
  if (numChildren==0) return;

  MPI_Request request;
  MPI_Status status;
  
  for (int ch=0; ch < numChildren; ch++){
    if (CommToChild_P[ch] != MPI_COMM_NULL and CG_numPBCMessages[ch]>0){ // this child wants PBC and this core participates

      // build all the PBC msgs that this core has to send 
      buildPBCMsg(grid, vct, ch);

      if (AllowPMsgResize){
	// the CG cores involved in PBC agree on the biggest buffer size
	int MaxGrid_sizeCG_PBCMsg;
	MPI_Allreduce(&sizeCG_PBCMsg, &MaxGrid_sizeCG_PBCMsg, 1, MPI_INT, MPI_MAX, COMM_CG_PBCSubset_P[ch]);
	
	// CGSide_CGLeader_PBCSubset[ch] sends down to all RG grid cores involved in the communication
	if (vct->getRank_CommToChildren_P(ch, ns) == CGSide_CGLeader_PBCSubset[ch]){
	  for (int i= 0; i< numRcv[ch]; i++){
	    MPI_Isend(&MaxGrid_sizeCG_PBCMsg, 1, MPI_INT, RcvList[ch][i], 200, CommToChild_P[ch], &request);
	    MPI_Wait(&request, &status);
	  }
	}
      } // end if (AllowPMsgResize){

      // now send; send an extra msg to signal the end of the meaningful part
      // cycle on all the msgs this core has to send

      /*if (false){
	for (int m=0; m< CG_numPBCMessages[ch]; m++){
	  for (int ii=0; ii< nopPCGMsg[ch][m]+1; ii++){
	    cout <<"CHECK CG, ch: " << ch << "m: "<< m << " ii: " << ii << " PCGMsg[ch][m][jj].q: "<< PCGMsg[ch][m][ii].q <<endl;
	  }
	}
	}*/

      for (int m=0; m< CG_numPBCMessages[ch]; m++){
	int dest= CG_Info[ch][m].RG_core;
	int tag= CG_Info[ch][m].MsgID;

	MPI_Isend(PCGMsg[ch][m], nopPCGMsg[ch][m]+1, MPI_RepP_struct, dest, tag, CommToChild_P[ch], &request );
	MPI_Wait(&request, &status);
      } // end for (int m=0; m< CG_numPBCMessages[ch]; m++)

    } //  if (CommToChild_P[ch] != MPI_COMM_NULL){ 
    
  } // end for (int ch=0; ch < numChildren; ch++){ 

    /*MPI_Barrier(vct->getComm());
  int RR= vct->getCartesian_rank();
  if (RR==0)
  cout << "Grid " << numGrid << " has finished sending PBC ns  " <<ns << endl;*/
  
  return;
}

/* RG receives PBC msg and acts accordingly */
void Particles3Dcomm::ReceivePBC(Grid* grid, VirtualTopology3D * vct){
  
  PRA_PAdded=0;


  /* this method is at the moment particularly inefficient to catch eventual bugs_
     improve after some testing */

  int RR= vct->getCartesian_rank();

  /*if (numGrid >0){
    //cout << "Grid " << numGrid << " core " << RR << " of " << vct->getXLEN()*vct->getYLEN()*vct->getZLEN() << " has started receiving PBC for sp " << ns << endl;

    cout << "GRID " << numGrid << " CORE " << RR << " IS INSIDE ReceivePBC" << endl;
    }*/
 
  // first, check if if i am a child and if i want to receive PBC
  // if not, exits
  if (CommToParent_P!= MPI_COMM_NULL and RG_numPBCMessages >0){

    MPI_Status status;

    if (AllowPMsgResize){ // to do before anybody has started receiving, so i don't have to copy info
      int NEW_sizePBCMsg;
      // as it is set now, each RG core receives a msg, so just do a rcv 
      MPI_Recv(&NEW_sizePBCMsg, 1, MPI_INT, MPI_ANY_SOURCE, 200, CommToParent_P, &status);

      // in case, resize
      if (NEW_sizePBCMsg > sizeRG_PBCMsg){
	resize_RG_MLMD_buffers(NEW_sizePBCMsg);
      }
    } // end if (AllowPMsgResize){ // to do before anybody has started receiving, so i don't have to copy info
    
    for (int i=0; i<RG_numPBCMessages; i++ ){
      PRGMsgArrived[i]= false;
      nopPRGMsg[i]=0;
    }
    
    int count, src, MsgID;;
    for (int i=0; i< RG_numPBCMessages; i++){
      MPI_Recv(PRGMsg_General, sizeRG_PBCMsg, MPI_RepP_struct, MPI_ANY_SOURCE, MPI_ANY_TAG, CommToParent_P, &status);
      MPI_Get_count( &status, MPI_RepP_struct, &count );
      
      /*cout << "Received msg " << i << " of " << RG_numPBCMessages << endl;
	cout<< "G" << numGrid  << "R"<< RR << " ns " << ns << " has received one PBC msg from core " << status.MPI_SOURCE<< " count is " << count << endl;*/


      /*for (int jj=0; jj< count ; jj++){
	cout << "msg " << i << " jj "<< jj << " PRGMsg_General[jj].q: " << PRGMsg_General[jj].q << endl;
	}*/

      // match the message to the RG_info
      src= status.MPI_SOURCE;
      MsgID= status.MPI_TAG;
      //cout << "src: " << src <<" MsgID: " << MsgID << endl;
      bool found= false;
      for (int mm=0; mm< RG_numPBCMessages; mm++){
	if (MsgID == RGPBC_Info[mm].MsgID and src == RGPBC_Info[mm].CG_core and PRGMsgArrived[mm]== false){ // message matched
	  // update the PRGMsgArrived 
	  PRGMsgArrived[mm]=true;
	  found= true;
	  // update the structure
	  //memcpy ( void * destination, const void * source, size_t num );
	  memcpy(PRGMsg[mm], PRGMsg_General, count*sizeof(RepP_struct));
	  break;
	} // end msg matched

      } // end for (int mm=0; mm< RG_numPBCMessages; mm++){
      if (found== false){
	cout << "In ReceivePBC, G" << numGrid <<"C" <<vct->getCartesian_rank() << " could not match a PBC msg from C" << src <<" (PC comm) with tag " << MsgID << ", aborting ..." << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      } 
      
    } // end for (int i=0; i< RG_numPBCMessages; i++){

    // here, PBC have been successfully received
    // now I have to apply them (split the particles )
    ApplyPBC(vct, grid);

  } // end if (CommToParent_P!= MPI_COMM_NULL and RG_numPBCMessages >0){
  
  nop_AfterReceivePBC= nop;

  //cout << "Grid " << numGrid <<" R " << vct->getCartesian_rank() << " ns " << ns << " added " << PRA_PAdded << " particles, deleted " << PRA_deleted << endl;

  // here, check how many PBC particles have been added at grid level
  if (false){
    int ToTPRA_PAdded;
    int TotPRA_deleted;
    MPI_Allreduce(&PRA_PAdded, &ToTPRA_PAdded, 1, MPI_INT, MPI_SUM, vct->getComm());
    MPI_Allreduce(&PRA_deleted, &TotPRA_deleted, 1, MPI_INT, MPI_SUM, vct->getComm());
    if (RR== XLEN*YLEN*ZLEN-1){
      cout << "Grid " << numGrid << " ns " << ns << " added " << ToTPRA_PAdded << " particles, deleted " << TotPRA_deleted << " (this before communicateRepopulatedParticles)" << endl;
    }
  } // end debug stuff

  if (CommToParent_P!= MPI_COMM_NULL) {
    communicateRepopulatedParticles(grid, vct);
  }

  /*if (numGrid>0)
    cout <<"Grid " << numGrid << " R " << vct->getCartesian_rank() << " ns " << ns << " nop " << nop << " at the end of ReceivePBC "<< endl;*/
}

void Particles3Dcomm::buildPBCMsg(Grid* grid, VirtualTopology3D * vct, int ch){
  
  /* returns the number - not the level or the order in the children vector - of the child grid n */
  int childNum= vct->getChildGridNum(ch);
  // these are the dx, dy, dz of the child
  double dx= grid->getDx_mlmd(childNum);
  double dy= grid->getDy_mlmd(childNum);
  double dz= grid->getDz_mlmd(childNum);

  double x_min, x_max, y_min, y_max, z_min, z_max;
  
  if (CG_numPBCMessages[ch]>0){ // this particular core has to send BC

    // pre
    for (int m=0; m< CG_numPBCMessages[ch]; m++){
      nopPCGMsg[ch][m]= 0;
    } 

    // core
    for (int p=0; p<nop; p++){
      // here, i have to check all the msgs in else if, to make sure that I am sending 
      // this particular particle only to one core
      
      for (int n=0; n< CG_numPBCMessages[ch]; n++){ 
	// send an extra dx; this is needed in case one refined grid core superimpose to more coarse grid cores-- do not touch this!!!
	x_min= CG_Info[ch][n].CG_x_first -dx ;
	x_max= CG_Info[ch][n].CG_x_first+ dx*(CG_Info[ch][n].np_x);
	y_min= CG_Info[ch][n].CG_y_first -dy;
	y_max= CG_Info[ch][n].CG_y_first+ dy*(CG_Info[ch][n].np_y);
	z_min= CG_Info[ch][n].CG_z_first -dz;
	z_max= CG_Info[ch][n].CG_z_first+ dz*(CG_Info[ch][n].np_z);
	
	if (x_min <= x[p] && x_max >= x[p] && y_min <= y[p] && y_max >= y[p] && z_min <= z[p] && z_max >= z[p]){ // to use for PBC
	  //cout << "adding " << p << " to " << n << endl;
	  unsigned long ID;
	  if (TrackParticleID) {ID= ParticleID[p];} else {ID=0;}
	  addP(ch, n, x[p], y[p], z[p], u[p], v[p], w[p], q[p], ID, vct);
	  break; // to make sure that a particle is added only to one msg
	} 
      }// end  for (int n=0; n< CG_numPBCMessages[ch]; n++)
    } // end for (int p=0; p<nop; p++){

    // post
    // the last one, to mark the end of the the meaningful part --> q=0.0 
    for (int m=0; m< CG_numPBCMessages[ch]; m++ ){
      int NUM= nopPCGMsg[ch][m];
      PCGMsg[ch][m][NUM].q= EXIT_VAL;

      //cout << "Grid " << numGrid <<" core " << vct->getCartesian_rank() << " on the local grid communicator is sending " << nopPCGMsg[ch][m] << " particles as BC to grid " << childNum << " core " << CG_Info[ch][m].RG_core << " (msg " << m << " of " <<CG_numPBCMessages[ch] <<")" << endl; 
      
    }

  } // end f (CG_numPBCMessages[ch]>0){ // this particular core has to send BC   
  
}

void Particles3Dcomm::addP(int ch, int n, double x, double y, double z, double u, double v, double w, double q, unsigned long ID, VirtualTopology3D* vct){


  int pos= nopPCGMsg[ch][n];
  PCGMsg[ch][n][pos].x=x;
  PCGMsg[ch][n][pos].y=y;
  PCGMsg[ch][n][pos].z=z;

  PCGMsg[ch][n][pos].u=u;
  PCGMsg[ch][n][pos].v=v;
  PCGMsg[ch][n][pos].w=w;

  PCGMsg[ch][n][pos].q=q;
  
  if (TrackParticleID)
    PCGMsg[ch][n][pos].ID=ID;
  
  nopPCGMsg[ch][n]= nopPCGMsg[ch][n]+1;

  //cout << "inside addP: pos: "<<pos <<"  PCGMsg[ch][n][pos].q " << PCGMsg[ch][n][pos].q << endl;
  if (nopPCGMsg[ch][n] == sizeCG_PBCMsg){
    if (AllowPMsgResize){
      resize_CG_MLMD_buffers(vct);
    }else{
      cout << "in addP, numGrid " << numGrid << " core " << vct->getCartesian_rank()  << " in the local grid communicator, you plan on passing too many particles as BC; " << endl << "ENABLE RESIZE --> AllowPMsgResize = 1 in the inputfile" << endl <<"Aborting now..." << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  return;

}

void Particles3Dcomm::ApplyPBC(VirtualTopology3D* vct, Grid* grid){

  int RR= vct->getCartesian_rank();
  /*if (RR==0){
    cout << "G" << numGrid <<"R" << RR << " INSIDE APPLYPBC" << endl;
    }*/
  //cout << "R" << RR <<" G" << numGrid << ", nop before applying BC: " << nop <<endl;

  for (int m=0; m< RG_numPBCMessages; m++){
    int count =0; 

    while (fabs(PRGMsg[m][count].q) > 1.1*EXIT_VAL){
      //cout << "m: " <<m << " count "<< count <<" PRGMsg[m][count].q: " << PRGMsg[m][count].q << endl;
      SplitPBC(vct, grid, PRGMsg[m][count]);
      count ++;
    }
    nopPRGMsg[m]= count;

    //cout <<"G"<<numGrid << "R" <<RR << " ns " << ns  << " m " << m << " nopPRGMsg[m] " << nopPRGMsg[m] << endl;
  } // end for (int m=0; m< RG_numPBCMessages; m++){

  //cout << "R" << RR <<" G" << numGrid << ", nop after applying BC: " << nop <<endl;
}

void Particles3Dcomm::MPI_Barrier_ParentChild(VirtualTopology3D* vct){

  int RS= vct->getSystemWide_rank();
  int RR= vct->getCartesian_rank();
  int TOT= vct->getNprocs();

  //cout <<"R" << RR <<"/" << TOT <<" G" <<numGrid << " before barrier CommToParent_P" << endl;

  if (CommToParent_P != MPI_COMM_NULL)  MPI_Barrier(CommToParent_P);
  /*if (RR==0)
    cout << "Grid " << numGrid << " after barrier CommToParent_P" << endl;
  
    cout <<"R" << RR <<"/" << TOT <<" G" <<numGrid << " after barrier CommToParent_P" << endl;*/

  for (int ch=0; ch< numChildren; ch++)
    if (CommToChild_P[ch] != MPI_COMM_NULL)
      MPI_Barrier(CommToChild_P[ch]);

  /*if (RR==0)
    cout << "Grid " << numGrid << " after barrier CommToChild_P[ch]" << endl;
    cout <<"R" <<RR <<"/" << TOT << " G" <<numGrid << " after barrier CommToParent_P" << endl;
  if (RS==0)
  cout << "Everybody at the end of particles3Dcomm::MPI_Barrier_ParentChild, ns " << ns << endl;*/
}

void Particles3Dcomm::CheckSentReceivedParticles(VirtualTopology3D* vct){

  int RR= vct->getCartesian_rank();

  // RG Side
  int RGPPerCore=0;
  int RGPGridWide=0;
  if (CommToParent_P != MPI_COMM_NULL){
    for (int m=0; m<RG_numPBCMessages; m++){
      RGPPerCore+= nopPRGMsg[m];
    }
    //cout << "G" << numGrid <<"R" << RR << " ns " << ns <<": RGPPerCore is "<< RGPPerCore <<endl;
    // MPI_Allreduce
    MPI_Allreduce(&RGPPerCore, &RGPGridWide, 1, MPI_INT, MPI_SUM, vct->getComm());

    if (RR==0){
      cout << "Grid "<< numGrid << " received " << RGPGridWide << " particles ns " << ns << " from grid " << vct->getParentGridNum()<< endl;
      }
  }
  // CG Side
  int * CGPPerCore= new int[numChildren];
  int * CGPGridWide= new int[numChildren];
  for (int ch=0; ch< numChildren; ch++){
    CGPPerCore[ch]=0;
    if (CommToChild_P[ch] != MPI_COMM_NULL){
      for (int m=0; m< CG_numPBCMessages[ch]; m++){
	CGPPerCore[ch] += nopPCGMsg[ch][m];
      } // end for (int m=0; m< CG_numPBCMessages[m]; m++)
      // MPI_Allreduce
      MPI_Allreduce(&(CGPPerCore[ch]), &(CGPGridWide[ch]), 1, MPI_INT, MPI_SUM, vct->getComm());
     
      if (RR==0){
	cout << "Grid "<< numGrid << " sent " << CGPGridWide[ch] << " particles ns " << ns << " to grid " << vct->getChildGridNum(ch)<< endl;
	}
    } // end if (CommToChild_P[ch] != MPI_COMM_NULL)
  } // end for (int ch=0; ch< numChildren; ch++)


  MPI_Status status;
  // compare the values in the two grids
  // RG sends info to CG
  if (CommToParent_P != MPI_COMM_NULL and RR==0){
    MPI_Send(&RGPGridWide, 1, MPI_INT, 0, 8, CommToParent_P);
    cout <<  " RGPGridWide: " << RGPGridWide  << endl;
  }
  // CG receives and check
  int *recvFromRG= new int[numChildren];
  for (int ch=0; ch < numChildren; ch++){
    if (CommToChild_P[ch] != MPI_COMM_NULL and RR==0){
      MPI_Recv(&(recvFromRG[ch]), 1, MPI_INT, XLEN*YLEN*ZLEN, 8, CommToChild_P[ch], &status);
      //cout << numGrid << " receives " << recvFromRG[ch] <<endl;
    }
  }

  for (int i=0; i<numChildren; i++){
    if (CommToChild_P[i] != MPI_COMM_NULL and RR==0){
      if (recvFromRG[i] != CGPGridWide[i]){
	cout << "FATAL ERROR: Grid " << numGrid << " there is a mismatch between # of particles sent/ received to child grid "<<vct->getChildGridNum(i)  << endl;
	cout <<"Grid " << numGrid <<" sent " << CGPGridWide[i] << " , child received " << recvFromRG[i] <<", abort..." << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
      else{
	//cout << "Grid " << numGrid << " and its child "<<vct->getChildGridNum(i) << " exchanged " << CGPGridWide[i] << " RG BC particles: same # of particles sent and received" << endl;
	
      }
    }
  }

  delete[]recvFromRG;
  delete[]CGPPerCore;
  delete[]CGPGridWide;

  // here, I check if particles in RG are the same # after receiving repopulated particles and after exchanging them (pre- and post-communicateRepopulatedParticles)

  int Tot_nop_AfterReceivePBC;
  int Tot_nop_EndCommunicate;
  int Tot_nop;

  MPI_Allreduce(&nop_AfterReceivePBC, &Tot_nop_AfterReceivePBC, 1, MPI_INT, MPI_SUM, vct->getComm());
  MPI_Allreduce(&nop, &Tot_nop, 1, MPI_INT, MPI_SUM, vct->getComm());
  MPI_Allreduce(&nop_EndCommunicate, &Tot_nop_EndCommunicate, 1, MPI_INT, MPI_SUM, vct->getComm());

  if (RR==XLEN*YLEN*ZLEN-1){
    if (Tot_nop_AfterReceivePBC != Tot_nop){
      cout << "Grid " << numGrid << " ns " << ns << " FATAL ERROR: I have changed the # of particles in communicateRepopulatedParticles, from " << Tot_nop_AfterReceivePBC << " to " << Tot_nop  << ", nop_EndCommunicate = " << Tot_nop_EndCommunicate << endl;
      cout << "Aborting ..." << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    } else{
      cout << "Grid " << numGrid << " check on nop pre-post communicateRepopulatedParticles passed "<< endl;
    }
  }
}
/** split each received particles **/
void Particles3Dcomm::SplitPBC(VirtualTopology3D * vct, Grid* grid, RepP_struct p){
  // to prevent the repopulation of particles which would try to accumulate outside the grid

  double PM= 0.0001;

  double xTmp;
  double yTmp;
  double zTmp;
  double qTmp;
  
  bool StX, StY, StZ;

  for (int i=0; i< ceil(RFx); i++){
    xTmp= p.x - DxP/2.0 + dx*(1./2. + i)- grid->getOx();
        
    // if outside the domain, no  point in continuing splitting
    if (xTmp < Coord_XLeft_Start+dx*PM or xTmp>Coord_XRight_End-dx*PM)
      continue;
    
    // am i inside a X PRA?
    StX= (xTmp> Coord_XLeft_Start+dx*PM and xTmp < Coord_XLeft_End) or (xTmp> Coord_XRight_Start and xTmp < Coord_XRight_End-dx*PM);
    /*if (! StX ) continue;*/

    for (int j=0; j< ceil(RFy); j++){
      yTmp= p.y - DyP/2.0 + dy*(1./2. + j)- grid->getOy();

      // if outside the domain, no  point in continuing splitting
      if (yTmp < Coord_YLeft_Start+dy*PM or yTmp>Coord_YRight_End-dy*PM)
	continue;
      
      // am i inside a Y PRA?
      StY= (yTmp> Coord_YLeft_Start+dy*PM and yTmp < Coord_YLeft_End) or (yTmp> Coord_YRight_Start and yTmp < Coord_YRight_End-dy*PM);
      /*if (! StY ) continue;*/

      for (int k=0; k< ceil(RFz); k++){
	zTmp= p.z - DzP/2.0 + dz*(1./2. + k)- grid->getOz();
	
	// if outside the domain, no  point in continuing splitting
	if (zTmp < Coord_ZLeft_Start+dz*PM or zTmp>Coord_ZRight_End-dz*PM)
	  continue;

	// am i inside a Z PRA?
	StZ= (zTmp> Coord_ZLeft_Start+dz*PM and zTmp < Coord_ZLeft_End) or (zTmp> Coord_ZRight_Start and zTmp < Coord_ZRight_End-dz*PM);
	
	//cout << " Coord_ZLeft_Start " << Coord_ZLeft_Start << " Coord_ZLeft_End " << Coord_ZLeft_End << " Coord_ZRight_Start " << Coord_ZRight_Start << " Coord_ZRight_End " << Coord_ZRight_End << endl;

	// if i am inside any PRA (if i arrived here i am inside the extended domain)
	bool Keep= false;
	// i should really separate by direction but i will do that later
	if (bcPfaceXleft== -1 ) { 
	  /* this is when i repopulate all the PRA
	     if you arrive here, you should be kep */
	  Keep= true;
	}
	else if (bcPfaceXleft== -2){ // here i keep only the particles that have entered in the last dt
	  Keep= RepopulatedParticleHasEnteredRG(xTmp, yTmp, zTmp, p.u, p.v, p.w);
	}

	if ((StX or StY or StZ) and Keep){
	  //cout << "Grid "<< numGrid <<" R" <<vct->getCartesian_rank() <<" Particle split and accepted " <<endl; 
	  
	  // if it arrives here, particle is to add to the PRA, not necessarily to this core
	  PRA_PAdded++;
	  
	  // add to core particle system
	  x[nop]= xTmp;
	  y[nop]= yTmp;
	  z[nop]= zTmp;
	  
	  u[nop]= p.u;
	  v[nop]= p.v;
	  w[nop]= p.w;
	  
	  q[nop]= p.q/ RFx/RFy/RFz;
	  
	  if (TrackParticleID)
	    ParticleID[nop]= p.ID;
	  
	  // end add to core particle system 
	  // check nop is not exceeded                                                                                              
	  nop++;
	  
	  if (nop > npmax) {
	    cout << "Grid " << numGrid <<": Number of particles in the domain " << nop << " and maxpart = " << npmax << endl;
	    cout << "Aborting ..." << endl;
	    MPI_Abort(MPI_COMM_WORLD, -1);
	    return;              // end the simulation because you dont have enough space on the array           
	  }
	}// end if (StX or StY or StZ)
      } // end for (int k=0; k< ceil(RFz); k++){
    } // end for (int j=0; j< ceil(RFy); j++){
  } // end for (int i=0; i< ceil(RFx); i++){

}

void Particles3Dcomm::CheckAfterInitWeightPBC(VirtualTopology3D * vct){
  // check if number of msgs sent and receive from the grids correspond
  int RR= vct->getCartesian_rank();

  // RG Side
  int RGMsgGW=0;
  if (CommToParent_P!= MPI_COMM_NULL){

    // MPI_Allreduce
    MPI_Allreduce(&RG_numPBCMessages, &RGMsgGW, 1, MPI_INT, MPI_SUM, vct->getComm());

    /*if (RR==0){
      cout << "Grid "<< numGrid << " must receive " << RGMsgGW  << " PBC msgs ns " << ns << " from grid " << vct->getParentGridNum()<< endl;
      }*/
  }
  // CG Side
  int * CGMsgGW= new int[numChildren];
  for (int ch=0; ch< numChildren; ch++){
    CGMsgGW[ch]=0;
    if (CommToChild_P[ch] != MPI_COMM_NULL){
     
      // MPI_Allreduce
      MPI_Allreduce(&(CG_numPBCMessages[ch] ), &(CGMsgGW[ch]), 1, MPI_INT, MPI_SUM, vct->getComm());
     
      /*if (RR==0){
	cout << "Grid "<< numGrid << " must send " << CGMsgGW[ch] << " PBC msgs ns " << ns << " to grid " << vct->getChildGridNum(ch)<< endl;
	}*/
    } // end if (CommToChild_P[ch] != MPI_COMM_NULL)
  } // end for (int ch=0; ch< numChildren; ch++)

  MPI_Status status;
  // compare the values in the two grids
  // RG sends info to CG
  if (CommToParent_P != MPI_COMM_NULL and RR==0){
    MPI_Send(&RGMsgGW, 1, MPI_INT, 0, 6, CommToParent_P);
    //cout << numGrid<< " RGMsgGW: " << RGMsgGW << endl;
  }
  // CG receives and check
  int *recvFromRG= new int[numChildren];
  for (int ch=0; ch < numChildren; ch++){
    if (CommToChild_P[ch] != MPI_COMM_NULL and RR==0){
      MPI_Recv(&(recvFromRG[ch]), 1, MPI_INT, XLEN*YLEN*ZLEN, 6, CommToChild_P[ch], &status);
      //cout << numGrid << " receives " << recvFromRG[ch] <<endl;
    }
  }

  for (int i=0; i<numChildren; i++){
    if (CommToChild_P[i] != MPI_COMM_NULL and RR==0){
      if (recvFromRG[i] != CGMsgGW[i]){
	cout << "FATAL ERROR: Grid " << numGrid << " there is a mismatch between PBC sent/ received to child grid "<<vct->getChildGridNum(i)  << endl;
	cout << numGrid <<" sends " << CGMsgGW[i] << " , child wants to receive " << recvFromRG[i] <<", abort..." << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
      else{
	//cout << numGrid << " and its child "<<vct->getChildGridNum(i) << " exchange " << recvFromRG[i] << " msgs" << endl;
      }
    }
  }

  delete[]CGMsgGW;
  delete[]recvFromRG;

}

void Particles3Dcomm::communicateRepopulatedParticles(Grid* grid, VirtualTopology3D * vct){
  if ( RG_numPBCMessages <1) return;

  int TAG=ns;
  int count;
  int rank_local= vct->getCartesian_rank();
  
  int HighestRank= XLEN*YLEN*ZLEN-1;

  bool RightCore;
  bool Id;
  long long np_current = nop_EndCommunicate, nplast = nop - 1;

  //cout <<"Grid " << numGrid <<" R " << rank_local << " nop_PreApplyPBC: " << " nop:" << nop << endl;
  
  int *Coord_P = new int[3];

  double xlen=Lx/ XLEN;
  double ylen=Ly/ YLEN;
  double zlen=Lz/ ZLEN;

  int DestCore;
  MPI_Status status;
  MPI_Request request;
  np_ToCoreH_CRP= 0;

  // this is set to 0 here because HighestRank starts updating it 
  //in the while (np_current < nplast+1)
  if (rank_local== HighestRank){
    for (int i=0; i<XLEN*YLEN*ZLEN; i++){
      H_num_CRP_nop[i]=0;
    }
    /*cout <<"Grid " << numGrid << ", Inside communicateRepopulatedParticles: the core that should participate are: " <<endl;
    for (int i=0; i<XLEN*YLEN*ZLEN; i++ ){
      cout << "R" <<i << ": " << H_CRP_cores[i] << endl;
      }*/
  }

  while (np_current < nplast+1){


    RightCore= (x[np_current] >= xStart_GC and x[np_current] <= xEnd_GC) and (y[np_current] >= yStart_GC and y[np_current] <= yEnd_GC) and (z[np_current] >= zStart_GC and z[np_current] <= zEnd_GC);

    if (! RightCore){
      // pack particle 
      //cout <<"Grid " << numGrid <<" R " << rank_local <<": particle " <<np_current << " will be moved" << endl;

      if (TrackParticleID)
	Id= ParticleID[np_current];
      else {Id=0;}
	
      Coord_P[0]= floor(x[np_current]/xlen);
      Coord_P[1]= floor(y[np_current]/ylen);
      Coord_P[2]= floor(z[np_current]/zlen);

      if (Coord_P[0]< 0 ) Coord_P[0]=0;
      if (Coord_P[0]> XLEN-1 ) Coord_P[0]=XLEN-1;
      if (Coord_P[1]<0 ) Coord_P[1]=0;
      if (Coord_P[1]> YLEN-1 ) Coord_P[1]=YLEN-1;
      if (Coord_P[2]< 0 ) Coord_P[2]=0;
      if (Coord_P[2]> ZLEN-1 ) Coord_P[2]=ZLEN-1;
      
      //Determines process rank in communicator given Cartesian location
      MPI_Cart_rank(vct->getComm(), Coord_P, &DestCore);

      if (rank_local != HighestRank){
	// check on vector size is done inside
	// np_ToCoreH_CRP is incremented inside
	addP_CRP(CRP_ToCoreH, &np_ToCoreH_CRP, x[np_current], y[np_current], z[np_current], u[np_current], v[np_current], w[np_current], q[np_current], Id, DestCore, vct);
      
      } // end if (rank_local != HighestRank){
      else{ // you are HighestRank, start packing it in the vector for the appropriate dest
	// check the destination is allowed
	if (H_CRP_cores[DestCore]!=1){
	  //cout << "L1: Grid "<<numGrid << ": in communicateRepopulatedParticles, HighestRank wants to send a msg to a forbidden destination... " << endl << "Check the inconsistency, aborting now ..." << endl;
	  cout << "L1: Grid "<<numGrid << ": in communicateRepopulatedParticles, HighestRank wants to send a msg to a forbidden destination... " << endl << "Check the inconsistency, now deleting particle and continuing ..." << endl;
	  
	  /*
	  cout <<"Cores to send msgs to:" << endl;
	  for (int j=0; j< XLEN*YLEN*ZLEN; j++)
	    cout <<" Core " <<j <<": " << H_CRP_cores[j] << endl;

	  cout << "But i want to message " << DestCore << endl;

	  if (DestCore == rank_local){
	    cout <<"WTF I am trying to msg myself??" << endl;
	    cout << "x[np_current]= " << x[np_current] << ", xStart_GC=" << xStart_GC <<", xEnd_GC=" <<xEnd_GC <<", -dx=" << -dx << ", Lx+dx= " << Lx+dx << endl;
	    cout << "y[np_current]= " << y[np_current] << ", yStart_GC=" << yStart_GC <<", yEnd_GC=" <<yEnd_GC <<", -dy=" << -dy << ", Ly+dy= " << Ly+dy << endl;
	    cout << "z[np_current]= " << z[np_current] << ", zStart_GC=" << zStart_GC <<", zEnd_GC=" <<zEnd_GC <<", -dz=" << -dz << ", Lz+dz= " << Lz+dz << endl;
	    cout << "Coordinates:" << endl;
	    cout << "X dir: " << vct->getCoordinates(0) << "/ " <<XLEN << endl;
	    cout << "Y dir: " << vct->getCoordinates(1) << "/ " <<YLEN << endl;
	    cout << "Z dir: " << vct->getCoordinates(2) << "/ " <<ZLEN << endl;
	    }
	  //MPI_Abort(MPI_COMM_WORLD, -1);
	  */
	  del_pack(np_current,&nplast);
	  continue;
	  
	}
	
	// I pack it in the vector from HighestRank to the appropriate core
	addP_CRP(H_CRP_Msg[DestCore], &(H_num_CRP_nop[DestCore]), x[np_current], y[np_current], z[np_current], u[np_current], v[np_current], w[np_current], q[np_current], Id, DestCore, vct);
      } // end else if (rank_local != HighestRank){
      
      del_pack(np_current,&nplast);
      
    }
    else{ // particle belongs here
      //cout <<"Grid " << numGrid <<" R " << rank_local <<": particle " <<np_current << " belongs here" << endl;
      np_current ++;
    }

  } // end while (np_current < nplast+1){

  nop= nplast +1;

  delete[]Coord_P;

  // intermediate check
  /*int nopOnCore=nop;
  if (rank_local == HighestRank){
    for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      nopOnCore+= H_num_CRP_nop[i];
    }
  }else{
    nopOnCore+=np_ToCoreH_CRP ;
  }

  int TotnopOnCore;
  int Totnop_AfterReceivePBC;
  MPI_Allreduce(&nopOnCore, &TotnopOnCore, 1, MPI_INT, MPI_SUM, vct->getComm());
  MPI_Allreduce(&nop_AfterReceivePBC, &Totnop_AfterReceivePBC, 1, MPI_INT, MPI_SUM, vct->getComm());
  
  if (rank_local == HighestRank){
    cout <<"Grid " << numGrid <<" Check in communicateRepopulatedParticles, after packing for send/ rcv to HR but before the actual send/rvc: " <<endl;
    if (Totnop_AfterReceivePBC== TotnopOnCore){
      cout << "Grid " <<numGrid <<": ns " << ns <<" intermediate test passed" << endl;
    }
    else{
      cout << "Grid " <<numGrid <<": ns " << ns <<" intermediate test NOT passed: TotnopOnCore: " << TotnopOnCore <<", Totnop_AfterReceivePBC: " << Totnop_AfterReceivePBC <<", Aborting..."<< endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  */

  // all cores which participate in PBC (i checked at the beginning on RG_numPBCMessages)
  // send their misplaced particles to Core HighestRank
  
  /*MPI_Barrier(vct->getComm());
  if (rank_local==HighestRank )
  cout << "Grid " << numGrid << " SONO QUI 0, grid barrier" << endl;*/

  /* -- exchange - resize_CRPbuffers_BSTH -- */
  if (AllowPMsgResize){ // if, from inputfile, the resize is an option, do it
    int New_size_CRP;
    // on the PBC communicators
    MPI_Allreduce(&size_CRP, &New_size_CRP, 1, MPI_INT, MPI_MAX, COMM_RG_PBCSubset_P);
    if (New_size_CRP> size_CRP) 
      resize_CRP_buffers(vct, New_size_CRP);
  }
  /* -- exchange - resize_CRPbuffers_BSTH -- */

  if (rank_local< HighestRank){
    MPI_Isend(CRP_ToCoreH, np_ToCoreH_CRP, MPI_CRP_struct, HighestRank, TAG, vct->getComm(), &request);  
    MPI_Wait(&request, &status);
    
    //cout << "CRP: Grid " <<numGrid << " R " <<rank_local << " of " << XLEN*YLEN*ZLEN  <<" sends " << np_ToCoreH_CRP << " repop particles to core " <<HighestRank << " ns "<< ns <<endl;
  }  

  if (rank_local == HighestRank){
    // receive from the other cores involved in communicateRepopulatedParticles
    // the msg with the particles to dispatch
    
    // num_H_CRP_cores has been initialised in initWeightPBC, it's the # of msg that HighestRank should expect
    // valid value only in HighestRank
    int count, src;
    for (int m=0; m< num_H_CRP_cores; m++){ 
      MPI_Recv(H_CRP_General, size_CRP, MPI_CRP_struct, MPI_ANY_SOURCE, TAG, vct->getComm(), &status);
      MPI_Get_count(&status, MPI_CRP_struct, &count);
      
      src= status.MPI_SOURCE;
      
      //cout <<"CRP: Grid "<<numGrid <<" R " << rank_local <<" has received " << count <<" repop particles from core " << src << " ns " <<ns << endl;
      
      // I am checking if this core is expected to participate in communicateRepopulatedParticles
      if (H_CRP_cores[src]!=1){
	cout << "L2: Grid "<<numGrid << ": in communicateRepopulatedParticles, HighestRank is receiving a msg from a forbidden source... " << endl << "Check the inconsistency, aborting now ..." << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
      
      // process the msg
      for (int mm=0; mm< count; mm++){
	
	// check if the destination (WhereToSend), as calculated from the originating core, is allowed
	int WTS= H_CRP_General[mm].Destination;
	
	if (WTS== rank_local){ // this core is the destination, unpack and go to next particle
	  unpack_CRP(H_CRP_General[mm], vct);
	}
	else { // you have to send 
	  
	  if (H_CRP_cores[WTS]!=1){
	    cout << "Grid "<<numGrid << ": in communicateRepopulatedParticles, HighestRank wants to send a msg to a forbidden destination... " << endl << "Check the inconsistency, aborting now ..." << endl;
	    MPI_Abort(MPI_COMM_WORLD, -1);
	  }
	  
	  // if so, assign to the right line in H_CRP_Msg and
	  unsigned long ID;
	  if (TrackParticleID) ID= H_CRP_General[mm].ID;
	  
	  addP_CRP(H_CRP_Msg[WTS], &(H_num_CRP_nop[WTS]), H_CRP_General[mm].x, H_CRP_General[mm].y, H_CRP_General[mm].z, H_CRP_General[mm].u, H_CRP_General[mm].v, H_CRP_General[mm].w, H_CRP_General[mm].q, ID, WTS, vct);
	  
	} // end if (WTS== rank_local){ 
      }// end for (int mm=0; mm< count; mm++) 
    } // end for (int m=0; m< num_H_CRP_cores; m++){ 
  } // end if (rank_local == HighestRank)
  
  /* -- exchange - resize_CRPbuffers_BSFH -- */
  if (AllowPMsgResize){ // if, from inputfile, the resize is an option, do it
    int New_size_CRP;
    // on the PBC communicators
    MPI_Allreduce(&size_CRP, &New_size_CRP, 1, MPI_INT, MPI_MAX, COMM_RG_PBCSubset_P);
    if (New_size_CRP> size_CRP) 
      resize_CRP_buffers(vct, New_size_CRP);
  }
  /* -- exchange - resize_CRPbuffers_BSFH -- */
    
  if (rank_local == HighestRank){
    // now, HighestRank sends directly to the right core
    int sent=0;
    for (int m=0; m< HighestRank; m++){  
      if (H_CRP_cores[m]!= 1) continue;
      
      MPI_Isend(H_CRP_Msg[m], H_num_CRP_nop[m], MPI_CRP_struct, m, ns, vct->getComm(), &request);
      MPI_Wait(&request, &status);
      //cout << "Grid "<< numGrid << " R "<< rank_local <<" has sent " << H_num_CRP_nop[m] << " particles to core " <<m << " ns "<<ns << endl;

      sent++;
    } // end for (int m=0; m< HighestRank; m++){
    if (sent != num_H_CRP_cores){
      cout << "Grid "<<numGrid << ": FATAL ERROR in communicateRepopulatedParticles, aborting ..." << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  } // end if (rank_local == HighestRank) 

  // here the next intermediate check, check nop for everybody Ã+ H_num_CRP_nopi] for HighestRank

  /*MPI_Barrier(vct->getComm());
  if (rank_local== HighestRank)
    cout << "Grid " << numGrid << " SONO QUI 2, grid barrier" << endl;

  nopOnCore=nop;
  if (rank_local == HighestRank){
    for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      nopOnCore+= H_num_CRP_nop[i];
    }
    }

  MPI_Allreduce(&nopOnCore, &TotnopOnCore, 1, MPI_INT, MPI_SUM, vct->getComm());
  MPI_Allreduce(&nop_AfterReceivePBC, &Totnop_AfterReceivePBC, 1, MPI_INT, MPI_SUM, vct->getComm());
  
  if (rank_local == HighestRank){
    cout <<"Grid " << numGrid <<" Second check in communicateRepopulatedParticles, after packing for send/ rcv to HR but before the actual send/rvc: " <<endl;
    if (Totnop_AfterReceivePBC== TotnopOnCore){
      cout << "Grid " <<numGrid <<": ns " << ns <<" Second intermediate test passed" << endl;
    }
    else{
      cout << "Grid " <<numGrid <<": ns " << ns <<" Second intermediate test NOT passed: TotnopOnCore: " << TotnopOnCore <<", Totnop_AfterReceivePBC: " << Totnop_AfterReceivePBC <<", Aborting..."<< endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }
  */
  // all the other cores now receive and unpack
  NH_num_CRP_nop=0;
  if (true){
  if (rank_local < HighestRank){
    // reuse the same bugger used for the send to HighestRank
    MPI_Recv(CRP_ToCoreH, size_CRP, MPI_CRP_struct, HighestRank, ns, vct->getComm(), &status);
    MPI_Get_count(&status, MPI_CRP_struct, &count);

    //cout << "LR: Grid "<< numGrid <<" R "<<rank_local <<" has received " << count << " particles from core " <<HighestRank << " ns " << ns << endl;
    
    for (int i=0; i< count; i++){
      // process each msg, checks are inside 
      unpack_CRP(CRP_ToCoreH[i], vct);
    
    }
  } // end if (rank_local < HighestRank){
  }
  /*MPI_Barrier(vct->getComm());
    cout << "Grid " << numGrid << " SONO QUI 3, grid barrier" << endl;*/
  
  //cout << "Grid "<< numGrid <<" R " << rank_local << " has finished communicateRepopulatedParticles, ns " << ns << endl;
}

void Particles3Dcomm::unpack_CRP(CRP_struct p, VirtualTopology3D * vct){
  
  // check if p.Destination corresponds to this core
  if (p.Destination != vct->getCartesian_rank()){
    cout << "Grid " << numGrid <<": FATAL ERROR in unpack_CRP, aborting ..." <<endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  // check if the particle belongs to this core- use xStart_GC, etc because must particles are going to be in GC
  bool RightCore= (p.x >= xStart_GC and p.x <= xEnd_GC) and (p.y >= yStart_GC and p.y <= yEnd_GC) and (p.z >= zStart_GC and p.z <= zEnd_GC);

  if (RightCore== false){
    cout << "Grid " << numGrid <<": FATAL ERROR in unpack_CRP, aborting ..." <<endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  
  // now i can add the particle
  x[nop]= p.x;
  y[nop]= p.y;
  z[nop]= p.z;

  q[nop]= p.q;

  u[nop]= p.u;
  v[nop]= p.v;
  w[nop]= p.w;

  
  if (TrackParticleID)
    ParticleID[nop]= p.ID;

  nop++;

  if (nop > npmax) {
      cout << "Number of particles in the domain " << nop << " and maxpart = " << npmax << endl;
      cout << "Grid " << numGrid <<": FATAL ERROR in unpack_CRP, aborting ..." <<endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
      return ;              // end the simulation because you dont have enough space on the array
  }

}

void Particles3Dcomm::addP_CRP(CRP_struct * Vec, int *num, double x, double y, double z, double u, double v, double w, double q, unsigned long ID, int DestinationRank, VirtualTopology3D* vct){
 
  (Vec+(*num))->x=x;

  (Vec+(*num))->y=y;
  (Vec+(*num))->z=z;

  (Vec+(*num))->u=u;
  (Vec+(*num))->v=v;
  (Vec+(*num))->w=w;

  (Vec+(*num))->q=q;
  
  (Vec+(*num))->ID=ID;

  (Vec+(*num))->Destination= DestinationRank;
  
  (*num)++;
  // cout << "adding P: after adding *num= " << *num << endl;

  if (*num  == size_CRP){
    if (AllowPMsgResize){
      resize_CRP_buffers(vct);
    }else{
      cout << "in addP_CRP, numGrid " << numGrid << " core " << vct->getCartesian_rank()  << " in the local grid communicator, you plan on passing too many particles as BC; " << endl << "ENABLE RESIZE --> AllowPMsgResize = 1 in the inputfile" << endl <<"Aborting now..." << endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  return;
 
}

/* resize the buffers responsible for sending repopulated particle from the CG to the RG -- CG side*/
void Particles3Dcomm::resize_CG_MLMD_buffers(VirtualTopology3D * vct){
  /* I will multiply by 2 */
  int NEW_sizePBCMsg= sizeCG_PBCMsg*2;
  
  int RR= vct->getCartesian_rank();

  if (NEW_sizePBCMsg >= MAXsizePBCMsg){
    cout << "Grid " <<numGrid << " R " << RR <<": attempt to resize CG particle repopulation buffers failed because NEW_sizePBCMsg > allowed value, "<< MAXsizePBCMsg << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  }else{
    cout << "Grid " <<numGrid <<" R " << RR <<": resize_CG_MLMD_buffers from " <<sizeCG_PBCMsg << " to " << NEW_sizePBCMsg << " particles, max allowed " << MAXsizePBCMsg  << endl;
  }
  
  // the only buffer that i have to copy is PCGMsg
  RepP_struct ***tmp= newArr3(RepP_struct, numChildren, MaxNumMsg, sizeCG_PBCMsg);

  for (int ch=0; ch<numChildren; ch++)
    for (int m=0; m< CG_numPBCMessages[ch]; m++)
      for (int n=0; n< nopPCGMsg[ch][m]; n++)
	tmp[ch][m][n]= PCGMsg_ptr[ch][m][n];
      //memcpy(&(tmp[ch][m][0]), &(PCGMsg_ptr[ch][m][0]), nopPCGMsg[ch][m]*sizeof(RepP_struct));

  //cout << "inside resize_CG_MLMD_buffers, before deleting: PCGMsg: " <<PCGMsg << " PCGMsg_ptr " << PCGMsg_ptr<<endl;
  delArr3(PCGMsg, numChildren, MaxNumMsg);

  PCGMsg= newArr3(RepP_struct, numChildren, MaxNumMsg, NEW_sizePBCMsg);
  PCGMsg_ptr= PCGMsg;
  //cout << "inside resize_CG_MLMD_buffers, after the new: PCGMsg: " <<PCGMsg << " PCGMsg_ptr " << PCGMsg_ptr<<endl;

  for (int ch=0; ch< numChildren; ch++)
    for (int m=0; m< CG_numPBCMessages[ch]; m++)
      for (int n=0; n< nopPCGMsg[ch][m]; n++)
	PCGMsg[ch][m][n]= tmp[ch][m][n];
	//memcpy(&(PCGMsg[ch][m][0]), &(tmp[ch][m][0]), nopPCGMsg[ch][m]*sizeof(RepP_struct));
  
      
  delArr3(tmp, numChildren, MaxNumMsg);
  
  sizeCG_PBCMsg= NEW_sizePBCMsg;
}
/* resize the buffers responsible for sending repopulated particles from the CG to the RG -- RG side                                                  
   the difference: I do not need to preserve info, only the resize buffers */
void Particles3Dcomm::resize_RG_MLMD_buffers(int NEW_sizePBCMsg){

  cout << "Grid " <<numGrid <<": resize_RG_MLMD_buffers from " <<sizeRG_PBCMsg << " to " << NEW_sizePBCMsg << " particles " << endl;

  sizeRG_PBCMsg= NEW_sizePBCMsg;

  // PRGMsg
  delArr2(PRGMsg, RG_numPBCMessages);
  PRGMsg= newArr2(RepP_struct, RG_numPBCMessages, sizeRG_PBCMsg);
  
  // PRGMsg_General
  delete[]PRGMsg_General;
  PRGMsg_General= new RepP_struct[sizeRG_PBCMsg];
}

void Particles3Dcomm::resize_CRP_buffers(VirtualTopology3D * vct){
  /** I will now multiply by 2 **/
  int NEW_size_CRP= size_CRP*2;

  int RR= vct->getCartesian_rank();
  int HighestRank= XLEN*YLEN*ZLEN-1;

  if (NEW_size_CRP > MAX_np_CRP){
    cout << "Grid " <<numGrid << " R " << RR <<": attempt to resize CRP particle repopulation buffers failed because NEW_size_CRPg > allowed value, "<< MAX_np_CRP << endl;
    MPI_Abort(MPI_COMM_WORLD, -1);
  } else{
    cout << "Grid " <<numGrid <<" R " << RR <<": resize_CRP_buffers from " << size_CRP << " to " << NEW_size_CRP << " particles, max allowed " << MAX_np_CRP  << endl;
  }
  
  if (RR!= HighestRank){ 
    /* if called by addP_CRP with RR!= HighestRank, it is called while building msg to HighestRank
       resize & preserve CRP_ToCoreH */

    CRP_struct * tmp_CRP_ToCoreH= new CRP_struct[size_CRP];

    //memcpy(tmp_CRP_ToCoreH, CRP_ToCoreH, sizeof(CRP_struct)* np_ToCoreH_CRP);
    for (int i=0; i<np_ToCoreH_CRP; i++)
      tmp_CRP_ToCoreH[i]= CRP_ToCoreH_ptr[i];
    delete[] CRP_ToCoreH;

    CRP_ToCoreH= new CRP_struct[NEW_size_CRP];
    CRP_ToCoreH_ptr= CRP_ToCoreH;

    //memcpy(CRP_ToCoreH, tmp_CRP_ToCoreH, sizeof(CRP_struct)* np_ToCoreH_CRP);
    for (int i=0; i< np_ToCoreH_CRP; i++)
      CRP_ToCoreH[i]= tmp_CRP_ToCoreH[i];
    delete[] tmp_CRP_ToCoreH;


  } else{ // here, RR== HighestRank
    /* if called by addP_CRP with RR== HighestRank, it is called when HighestRank is building msg
       to everbybody else -
       resize H_CRP_General; resize & preserve H_CRP_Msg*/

    /* -- H_CRP_General -- */
    //delete[] H_CRP_General;
    //H_CRP_General= new CRP_struct[NEW_size_CRP];
    //H_CRP_General_ptr= H_CRP_General;

    CRP_struct *H_CRP_General_tmp= new CRP_struct[size_CRP];

    for (int i=0; i< size_CRP; i++)
      H_CRP_General_tmp[i]= H_CRP_General_ptr[i];
    
    delete[] H_CRP_General;
    
    H_CRP_General= new CRP_struct[NEW_size_CRP];
    H_CRP_General_ptr= H_CRP_General;
    
    for (int i=0; i< size_CRP; i++)
      H_CRP_General[i]= H_CRP_General_tmp[i];
    delete[]H_CRP_General_tmp;

    /* -- H_CRP_Msg -- */
    CRP_struct ** H_CRP_Msg_tmp= newArr2(CRP_struct, XLEN*YLEN*ZLEN, size_CRP);
    
    for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      //memcpy(&(H_CRP_Msg_tmp[i][0]), &(H_CRP_Msg[i][0]), sizeof(CRP_struct)* H_num_CRP_nop[i]);
      for (int j=0; j< H_num_CRP_nop[i]; j++)
	H_CRP_Msg_tmp[i][j]= H_CRP_Msg_ptr[i][j];
    }

    delArr2(H_CRP_Msg, XLEN*YLEN*ZLEN);
    H_CRP_Msg= newArr2(CRP_struct, XLEN*YLEN*ZLEN, NEW_size_CRP);
    H_CRP_Msg_ptr= H_CRP_Msg;
    
    for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      //memcpy(&(H_CRP_Msg[i][0]), &(H_CRP_Msg_tmp[i][0]), sizeof(CRP_struct)* H_num_CRP_nop[i]);
      for (int j=0; j<  H_num_CRP_nop[i]; j++)
	H_CRP_Msg[i][j]= H_CRP_Msg_tmp[i][j];
    }
    delArr2(H_CRP_Msg_tmp, XLEN*YLEN*ZLEN);

  }

  size_CRP= NEW_size_CRP;
}

void Particles3Dcomm::resize_CRP_buffers(VirtualTopology3D * vct, int NewSize){

  if (size_CRP >= NewSize) return;

  cout << "Grid " <<numGrid <<": resize_CRP_buffers_BSTH from " << size_CRP << " to " << NewSize << " particles " << endl;

  int RR= vct->getCartesian_rank();
  int HighestRank= XLEN*YLEN*ZLEN-1;

  if (RR==HighestRank){
    /* resize - without preserving - H_CRP_General 
     resize & preserve H_CRP_Msg */
    
    /* -- H_CRP_General -- */
    //delete[] H_CRP_General;
    //H_CRP_General= new CRP_struct[NewSize];
    //H_CRP_General_ptr= H_CRP_General;

    CRP_struct *H_CRP_General_tmp= new CRP_struct[size_CRP];

    for (int i=0; i< size_CRP; i++)
      H_CRP_General_tmp[i]= H_CRP_General_ptr[i];
    
    delete[] H_CRP_General;
    
    H_CRP_General= new CRP_struct[NewSize];
    H_CRP_General_ptr= H_CRP_General;
    
    for (int i=0; i< size_CRP; i++)
      H_CRP_General[i]= H_CRP_General_tmp[i];
    delete[]H_CRP_General_tmp;

    /* -- H_CRP_Msg -- */
    CRP_struct ** H_CRP_Msg_tmp= newArr2(CRP_struct, XLEN*YLEN*ZLEN, size_CRP);
    
    for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      //memcpy(&(H_CRP_Msg_tmp[i][0]), &(H_CRP_Msg[i][0]), sizeof(CRP_struct)* H_num_CRP_nop[i]);
      for (int j=0; j< H_num_CRP_nop[i]; j++)
        H_CRP_Msg_tmp[i][j]= H_CRP_Msg_ptr[i][j];
    }

    delArr2(H_CRP_Msg, XLEN*YLEN*ZLEN);
    H_CRP_Msg= newArr2(CRP_struct, XLEN*YLEN*ZLEN, NewSize);
    H_CRP_Msg_ptr= H_CRP_Msg;
    
    for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      for (int j=0; j< H_num_CRP_nop[i]; j++)
        H_CRP_Msg[i][j]= H_CRP_Msg_tmp[i][j];
      //memcpy(&(H_CRP_Msg[i][0]), &(H_CRP_Msg_tmp[i][0]), sizeof(CRP_struct)* H_num_CRP_nop[i]);
    }
    delArr2(H_CRP_Msg_tmp, XLEN*YLEN*ZLEN);

  } else{ 
    /* not HighestRank: resize & preserve CRP_ToCoreH */

    CRP_struct * tmp_CRP_ToCoreH= new CRP_struct[size_CRP];

    //memcpy(tmp_CRP_ToCoreH, CRP_ToCoreH, sizeof(CRP_struct)* np_ToCoreH_CRP);
    for (int i=0; i<np_ToCoreH_CRP; i++)
      tmp_CRP_ToCoreH[i]= CRP_ToCoreH_ptr[i];
    delete[] CRP_ToCoreH;

    CRP_ToCoreH= new CRP_struct[NewSize];
    CRP_ToCoreH_ptr= CRP_ToCoreH;

    //memcpy(CRP_ToCoreH, tmp_CRP_ToCoreH, sizeof(CRP_struct)* np_ToCoreH_CRP);
    for (int i=0; i< np_ToCoreH_CRP; i++)
      CRP_ToCoreH[i]= tmp_CRP_ToCoreH[i];
    delete[] tmp_CRP_ToCoreH;

  }
  
  size_CRP= NewSize;
}


bool Particles3Dcomm::RepopulatedParticleHasEnteredRG(double xTmp, double yTmp, double zTmp, double u, double v, double w){
  // this is for repopulation method -2
  // i want to check if the repopulated particle has entered during the last dt
  // i return a false if the particle was already inside RG last dt
  // otherwise true

  double xB= xTmp- u*dt;
  double yB= yTmp- v*dt;
  double zB= zTmp- w*dt;

  if ((xB > Coord_XLeft_Start and xB < Coord_XRight_End) and (yB > Coord_YLeft_Start and yB < Coord_YRight_End) and (zB > Coord_ZLeft_Start and zB < Coord_ZRight_End) )
    return false;

  return true;
  
}
