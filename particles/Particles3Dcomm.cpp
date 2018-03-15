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
#include "MyClock.h"
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

extern MyClock *clocks;

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
  if (!FluidLikeRep){
    // RG side
    delArr2(PRGMsg, RG_numPBCMessages);
    delete[]nopPRGMsg;
    delete[]PRGMsgArrived;
    delete[]PRGMsg_General;
    // CG side
    delArr3(PCGMsg, numChildren, MaxNumMsg);
    delArr2(nopPCGMsg, numChildren);
  }
  if (RG_numPBCMessages and !FluidLikeRep){
    // RG side
    delete CRP_ToCoreH;  
    delArr2(H_CRP_Msg, XLEN*YLEN*ZLEN);
    delete[]H_num_CRP_nop;
    delete[]H_CRP_General;
    delete[]H_CRP_cores;
  }

  if (RG_numPBCMessages and FluidLikeRep){
    // RG side
    delArr2(rho_FBC, RG_numPBCMessages);
    delArr2(Jx_FBC, RG_numPBCMessages);
    delArr2(Jy_FBC, RG_numPBCMessages);
    delArr2(Jz_FBC, RG_numPBCMessages);
    delArr2(Pxx_FBC, RG_numPBCMessages);
    delArr2(Pxy_FBC, RG_numPBCMessages);
    delArr2(Pxz_FBC, RG_numPBCMessages);
    delArr2(Pyy_FBC, RG_numPBCMessages);
    delArr2(Pyz_FBC, RG_numPBCMessages);
    delArr2(Pzz_FBC, RG_numPBCMessages);

    delArr3(RHOP, nxc, nyc);
    delArr3(UP, nxc, nyc);
    delArr3(VP, nxc, nyc);
    delArr3(WP, nxc, nyc);
    delArr3(UTHP, nxc, nyc);
    delArr3(VTHP, nxc, nyc);
    delArr3(WTHP, nxc, nyc);

    delete[]RGFluidMsg;
    delArr3(REPO, nxc, nyc);
  }
  // put another condition so this happens only on RG
  if (FluidLikeRep){
    delete[]CGFluidMsg;
  }
  delArr3(Ex, nxn, nyn);
  delArr3(Ey, nxn, nyn);
  delArr3(Ez, nxn, nyn);

  delArr3(Bx, nxn, nyn);
  delArr3(By, nxn, nyn);
  delArr3(Bz, nxn, nyn);

  delArr3(Bx_ext, nxn, nyn);
  delArr3(By_ext, nxn, nyn);
  delArr3(Bz_ext, nxn, nyn);
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
  //cout << "numGrid: " << numGrid << ", RFx: "<< RFx <<", RFy: "<< RFy << ", RFz: " << RFz <<endl;
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

  nxc = grid->getNXC();
  nyc = grid->getNYC();
  nzc = grid->getNZC();
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

  dx_Ch=new double[numChildren];
  dy_Ch=new double[numChildren];
  dz_Ch=new double[numChildren];

  Ox_Ch=new double[numChildren];
  Oy_Ch=new double[numChildren];
  Oz_Ch=new double[numChildren];

  Lx_Ch=new double[numChildren];
  Ly_Ch=new double[numChildren];
  Lz_Ch=new double[numChildren];

  for (int i=0; i< numChildren; i++){
    // c: number of the child i in the SW  hierarchy
    c= vct->getChildGridNum(i);

    dx_Ch[i]= col->getDx_mlmd(c) ;
    dy_Ch[i]= col->getDy_mlmd(c);
    dz_Ch[i]= col->getDz_mlmd(c);

    // getOx_P(c) is the origin of grid c in the parent's coords
    // hence, these coordinates
    Ox_Ch[i]= col->getOx_P(c) ;
    Oy_Ch[i]= col->getOy_P(c);
    Oz_Ch[i]= col->getOz_P(c);

    Lx_Ch[i]= col->getLx_mlmd(c) ;
    Ly_Ch[i]= col->getLy_mlmd(c);
    Lz_Ch[i]= col->getLz_mlmd(c);
  }

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
  int MaxGridPer= vct->getMaxGridPer();
  MAX_RG_numPBCMessages= (vct->getMaxRF1()+8)* vct->getMaxRF1()* vct->getMaxRF1() + 8;
  MAX_RG_numPBCMessages_LevelWide= MAX_RG_numPBCMessages * MaxGridPer; 

  FluidLikeRep= col->getFluidLikeRep();
  numFBC= 10;

  // wether to allow resize of repopulated particle buffers
  // NB: if the repopulation is fluid, AllowPMsgResize is set to false 
  // otherwise, false in the inputfile means that resizing is done only the first three cycles, and then the values is kept
  // with true, the resizing is done every cycle
  // in both cases, i need the initial infrastructure, so i set to true here
  // the value is modified in UpdateAllowPMsgResize called in UpdateCycleInfo

  if (FluidLikeRep){
    AllowPMsgResize= false;
  } else {
    AllowPMsgResize= true;
  }


  // sizes of the PCGMsg values set here (but allocated only if needed) 
  // to have the send/ receive vectors with the same size, I cook up a number based 
  //   on the coarsest grid 
  // * start with a relatively low number, so it does not run out of memory *
  // if you are doing kinetic repopulation, communication between the grids and in the RG will size the buffers to a decent value
  int MaxNopMLMD= ceil(npcelx*npcely*npcelz *( col->getNxc_mlmd(0)/ col->getXLEN_mlmd(0) )*( col->getNyc_mlmd(0)/ col->getYLEN_mlmd(0) )*(col->getNzc_mlmd(0)/ col->getZLEN_mlmd(0)) *0.05 ); // *4*4 ; 

  sizeCG_PBCMsg = MaxNopMLMD; // CG side
  sizeRG_PBCMsg = MaxNopMLMD; // RG side
  MAXsizePBCMsg = 200*MaxNopMLMD; //20*MaxNopMLMD;


  // end sizes of the PCGMsg values set here 

  // here, to be able to used RGPBC_struct as an MPI_Datatype
  //MPI_RGPBC_struct_commit(); -- after unifying RGBC_struct and RGPBC_stuct in Grid
  MPI_RGPBC_struct_commit(); 

  // here, to be able to use the RepP_struct as an MPI_Datatype
  MPI_RepP_struct_commit();
  /* commit the structure for repopulated particle exchange within the RG */
  MPI_CRP_struct_commit();

  // PRA related initialisation
  // the # of PRA cells is a tmp variable
  //   the 'visible' variables are the index at which PRA starts/ ends 
  
  PRACells= col->getPRA();  

  /*cout << "I AM CHANGING PRA NUMBER FOR IONS" << endl;
    if (qom>0) PRACells=2;*/
  //cout << "Species " << ns <<" qom " << qom << "PRACells= " << PRACells << endl;

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


  if (vct->getPERIODICX_P()){  
    Coord_YLeft_Start=0;
    Coord_XLeft_End=0;
    Coord_XRight_Start=Lx;
    Coord_XRight_End=Lx;
  }

  if (vct->getPERIODICY_P()){  
    Coord_YLeft_Start=0;
    Coord_YLeft_End=0;
    Coord_YRight_Start=Ly;
    Coord_YRight_End=Ly;
  }

  if (vct->getPERIODICZ_P()){  
    Coord_ZLeft_Start=0;
    Coord_ZLeft_End=0;
    Coord_ZRight_Start=Lz;
    Coord_ZRight_End=Lz;
  }

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

  // the number will be sized automatically (see AllowPMsgResize)
  size_CRP=  ceil(nop*0.05); // it just needs to survive the init, then it will be resized
  


  // I am instantiating here local copies of Ex, Ey, Ezm Bx, By, Bz, Bx_ext, By_ext, Bz_ext
  // - used in the mover and only in one repopulation method
  Ex= newArr3(double, nxn, nyn, nzn);
  Ey= newArr3(double, nxn, nyn, nzn);
  Ez= newArr3(double, nxn, nyn, nzn);

  Bx= newArr3(double, nxn, nyn, nzn);
  By= newArr3(double, nxn, nyn, nzn);
  Bz= newArr3(double, nxn, nyn, nzn);

  Bx_ext= newArr3(double, nxn, nyn, nzn);
  By_ext= newArr3(double, nxn, nyn, nzn);
  Bz_ext= newArr3(double, nxn, nyn, nzn);

  // to output particle-related info
  string SaveDirName=col->getSaveDirName();
  stringstream num_grid_STR; num_grid_STR  << numGrid;
  stringstream num_sp_STR; num_sp_STR  << ns;
  stringstream num_rank_STR; num_rank_STR  << vct->getCartesian_rank();
  parInfo = SaveDirName + "/parInfo_G" + num_grid_STR.str() + "_sp" + num_sp_STR.str() + ".txt";
  if (vct->getCartesian_rank() == 0) {
    ofstream my_file(parInfo.c_str());
    my_file.close();
  }

  saveRepParFile= false;
  // to save repopulated particles
  if (saveRepParFile and numGrid >0){
    RepopulatedPar = SaveDirName + "/RepParInfo_G" + num_grid_STR.str() + "_sp" + num_sp_STR.str() + "_c" + num_rank_STR.str() + ".txt";

    ofstream my_file(RepopulatedPar.c_str());
    my_file.close();    
  }
  DiagnosticsOutputCycle= col->getDiagnosticsOutputCycle();

  TEST_FLUID_BC= false;

  // with this, the first d.f. of repopulated particles is = to the total distribution function
  nop_BeforeReceivePBC=0; 
  
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

  double qmax=0;
  double qmin=1000;

  int nop_InsideActive=0;
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

      /** some simple diagnostics **/
      if (fabs(q[i]) > qmax) qmax= fabs(q[i]);
      if (fabs(q[i]) < qmin) qmin= fabs(q[i]);
      if ((q[i]/ qom ) < 0) {
	cout <<"FATAL MISTAKE: IN INTERPP2G, AT LEAST ONE PARTICLE HAS CHARGE WITH SIGN != QOM: Q[I]= " <<q[i] <<", QOM= "<< qom <<". ABORTING NOW..." <<endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
	return;  
      }
      if (x[i]>= 0 and x[i]<=Lx and y[i]>= 0 and y[i]<=Ly and z[i]>= 0 and z[i]<=Lz){
	nop_InsideActive++;
      }
      /** some simple diagnostics **/

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
      EMf->addRho(weight, ix, iy, iz, ns);
      //EMf->addRho(weight, ix, iy, iz, ns, vct, x[i], y[i], z[i]);
      //EMf->addRho(weight, ix, iy, iz, ns, nxn, nyn, nzn);
      // add current density - X
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = u[i] * weight[ii][jj][kk];
      EMf->addJx(temp, ix, iy, iz, ns);
      //EMf->addJx(temp, ix, iy, iz, ns, nxn, nyn, nzn);
      // add current density - Y
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = v[i] * weight[ii][jj][kk];
      EMf->addJy(temp, ix, iy, iz, ns);
      //EMf->addJy(temp, ix, iy, iz, ns, nxn, nyn, nzn);
      // add current density - Z
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = w[i] * weight[ii][jj][kk];
      EMf->addJz(temp, ix, iy, iz, ns);
      //EMf->addJz(temp, ix, iy, iz, ns, nxn, nyn, nzn);  
      // Pxx - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = u[i] * u[i] * weight[ii][jj][kk];
      EMf->addPxx(temp, ix, iy, iz, ns);
      //EMf->addPxx(temp, ix, iy, iz, ns, nxn, nyn, nzn);
      // Pxy - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = u[i] * v[i] * weight[ii][jj][kk];
      EMf->addPxy(temp, ix, iy, iz, ns);
      //EMf->addPxy(temp, ix, iy, iz, ns, nxn, nyn, nzn);
      // Pxz - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = u[i] * w[i] * weight[ii][jj][kk];
      EMf->addPxz(temp, ix, iy, iz, ns);
      //EMf->addPxz(temp, ix, iy, iz, ns, nxn, nyn, nzn);
      // Pyy - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = v[i] * v[i] * weight[ii][jj][kk];
      EMf->addPyy(temp, ix, iy, iz, ns);
      //EMf->addPyy(temp, ix, iy, iz, ns, nxn, nyn, nzn);
      // Pyz - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = v[i] * w[i] * weight[ii][jj][kk];
      EMf->addPyz(temp, ix, iy, iz, ns);
      //EMf->addPyz(temp, ix, iy, iz, ns, nxn, nyn, nzn);
      // Pzz - add pressure tensor
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
            temp[ii][jj][kk] = w[i] * w[i] * weight[ii][jj][kk];
      EMf->addPzz(temp, ix, iy, iz, ns);
      //EMf->addPzz(temp, ix, iy, iz, ns, nxn, nyn, nzn);  
    }
    // change this to allow more parallelization after implementing array class
    //#pragma omp critical
    //EMf->addToSpeciesMoments(speciesMoments,ns);
  }
  // communicate contribution from ghost cells 

  EMf->communicateGhostP2G(ns, bcPfaceXright, bcPfaceXleft, bcPfaceYright, bcPfaceYleft, bcPfaceZright, bcPfaceZleft, vct, grid);

  //cout << "At the end of interpP2G, nop= " << nop << " in grid " << numGrid << " core "<< vct->getCartesian_rank() << endl;

#ifdef __PROFILING__

  double QMAX ; double QMIN; int nop_InsideActive_TOT;
  MPI_Allreduce(&qmin, &QMIN, 1, MPI_DOUBLE, MPI_MIN, vct->getCommGrid());
  MPI_Allreduce(&qmax, &QMAX, 1, MPI_DOUBLE, MPI_MAX, vct->getCommGrid());
  MPI_Allreduce(&nop_InsideActive, &nop_InsideActive_TOT, 1, MPI_INT, MPI_SUM, vct->getCommGrid());
  if (vct->getCartesian_rank()==0){

    ofstream my_file(parInfo.c_str(), fstream::app);
    my_file << endl << QMIN <<" " << QMAX <<" " << nop_InsideActive_TOT <<" ";
    cout << "QMIN: " << QMIN << " QMAX: " <<QMAX << endl;
    my_file.close();
    }
 #endif
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
    xMax=Lx; yMax=Ly; zMax=Lz;

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
    if (ptVCT->getPERIODICX_P()) { xMin=-Lx; xMax=2*Lx; } // so periodic particles will stay in the system and be communicate
    else if (bcPfaceXleft <0 and CommToParent_P!= MPI_COMM_NULL) { xMin= Coord_XLeft_Start; xMax= Coord_XRight_End; } // mlmd
    else {xMin= 0; xMax= Lx;} // non periodic, non mlmd

    if (ptVCT->getPERIODICY_P()) { yMin=-Ly; yMax=2*Ly; }
    else if (bcPfaceYleft <0 and CommToParent_P!= MPI_COMM_NULL) { yMin= Coord_YLeft_Start; yMax= Coord_YRight_End; }
    else {yMin= 0; yMax= Ly;}

    if (ptVCT->getPERIODICZ_P()) { zMin=-Lz; zMax=2*Lz; }
    else if (bcPfaceZleft <0 and CommToParent_P!= MPI_COMM_NULL) { zMin= Coord_ZLeft_Start; zMax= Coord_ZRight_End; }
    else {zMin= 0; zMax= Lz;}

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
  
  if (! (CommToParent_P!= MPI_COMM_NULL and (bcPfaceXleft <0 or bcPfaceXright <0 or bcPfaceYleft <0 or bcPfaceYright <0 or bcPfaceZleft <0 or bcPfaceZright <0 ))) return 0;
  
  MPI_Status status;

  double xMin, yMin, zMin;
  double xMax, yMax, zMax;

  long long np_current = 0, nplast = nop - 1;
  while (np_current < nplast+1){
     
    xMin=0; yMin=0; zMin=0;
    xMax=Lx; yMax=Ly; zMax=Lz;
    
    if (ptVCT->getPERIODICX_P()) { xMin=-Lx; xMax=2*Lx; } // so periodic particles will stay in the system and be communicated    
    else if ((bcPfaceXleft == -1 or bcPfaceXleft == -3 or bcPfaceXleft == -4) and CommToParent_P!= MPI_COMM_NULL) {
      xMin= Coord_XLeft_End; xMax= Coord_XRight_Start;
    }
    else if (bcPfaceXleft == -2 and CommToParent_P!= MPI_COMM_NULL) {
      xMin= Coord_XLeft_Start; xMax= Coord_XRight_End;
    }
    else {xMin= 0; xMax= Lx;} // non periodic, non mlmd              

    if (ptVCT->getPERIODICY_P()) { yMin=-Ly; yMax=2*Ly; } // so periodic particles will stay in the system and be communicated    
    else if ((bcPfaceYleft == -1 or bcPfaceYleft == -3 or bcPfaceYleft == -4) and CommToParent_P!= MPI_COMM_NULL) {
      yMin= Coord_YLeft_End; yMax= Coord_YRight_Start;
    }
    else if (bcPfaceYleft == -2 and CommToParent_P!= MPI_COMM_NULL) {
      yMin= Coord_YLeft_Start; yMax= Coord_YRight_End;
    }
    else {yMin= 0; yMax= Ly;} // non periodic, non mlmd     

    if (ptVCT->getPERIODICZ_P() ) { zMin=-Lz; zMax=2*Lz; } // so periodic particles will stay in the system and be communicated    
    else if ((bcPfaceZleft == -1 or bcPfaceZleft == -3 or bcPfaceZleft == -4) and CommToParent_P!= MPI_COMM_NULL) {
      zMin= Coord_ZLeft_End; zMax= Coord_ZRight_Start;
    }
    else if (bcPfaceZleft == -2 and CommToParent_P!= MPI_COMM_NULL) {
      zMin= Coord_ZLeft_Start; zMax= Coord_ZRight_End;
    }
    else {zMin= 0; zMax= Lz;} // non periodic, non mlmd 
    

    if (x[np_current] < xMin or x[np_current]> xMax or y[np_current] < yMin or y[np_current]> yMax or z[np_current] < zMin or z[np_current]> zMax){
      // particle to delete
      del_pack(np_current,&nplast);

    }
    else {
      // particle ok
      // particle is still in the domain, procede with the next particle
      np_current++;
    }
    
  }
  nop = nplast + 1;  
  /** do nor touch this otherwise mess in communicateRepopulatedParticles **/
  nop_EndCommunicate= nop;

  // i am removing sync points
#ifdef __PROFILING__
  int ppc=nop; int TotalP=0;
  MPI_Allreduce(&ppc, &TotalP, 1, MPI_INT, MPI_SUM, ptVCT->getCommGrid());
  if (ptVCT->getCartesian_rank()==0){
    cout << "Grid " << numGrid << " ns " <<ns  <<": total number of particles AFTER communicateAfterMover (when particles in PRA removed): " << TotalP << endl;

    ofstream my_file(parInfo.c_str(), fstream::app);
    my_file << TotalP <<" " ;
    my_file.close();
    }

#endif
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
double Particles3Dcomm::getMaxVelocity(MPI_Comm Comm) { 
  double localVel = 0.0;
  double maxVel = 0.0;
  for (long long i = 0; i < nop; i++){
    /** cosider only particles in the active part of the grid **/
    if (x[i]< 0 or x[i]>Lx or y[i]< 0 or y[i]>Ly or z[i]< 0 or z[i]>Lz)
      continue;
    localVel = max(localVel, sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]));
  }
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
    /** cosider only particles in the active part of the grid **/
    if (x[i]< 0 or x[i]>Lx or y[i]< 0 or y[i]>Ly or z[i]< 0 or z[i]>Lz)
      continue;

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
/* distribution function only of the repopulated particles -
   i test also the ones in the Ghost Cells */
unsigned long *Particles3Dcomm::getVelocityDistribution_RepPart(int nBins, double maxVel, MPI_Comm Comm) {
  unsigned long *f = new unsigned long[nBins];
  for (int i = 0; i < nBins; i++)
    f[i] = 0;
  double Vel = 0.0;
  double dv = maxVel / nBins;
  int bin = 0;

  //cout << "repopulated distribution function: " << (nop-nop_BeforeReceivePBC) <<" particles vs " << nop << " in this core " << endl;
  //cout << "NB: the first d.f. of repopulated particles is actually the entire d.f. of the RG" <<endl;
  /* test only the particles added after repopulation */
  for (long long i = nop_BeforeReceivePBC; i < nop; i++) {
    /** cosider only particles in the active part of the grid **/

    if (nop_BeforeReceivePBC==0){ // this happens only the first time
      // manually exclude everybody from the active area
      if ( (x[i]> -dx+PRACells*dx and x[i]<Lx+dx-PRACells*dx and y[i]> -dy+PRACells*dy and y[i]<Ly+dy-PRACells*dy and z[i]> -dz+PRACells*dx and z[i]<Lz+dz-PRACells*dz))
	continue;
    }
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


/* distribution function only of the NON- repopulated particles -
   i test also the ones in the Ghost Cells */
unsigned long *Particles3Dcomm::getVelocityDistribution_NonRepPart(int nBins, double maxVel, MPI_Comm Comm) {
  unsigned long *f = new unsigned long[nBins];
  for (int i = 0; i < nBins; i++)
    f[i] = 0;
  double Vel = 0.0;
  double dv = maxVel / nBins;
  int bin = 0;

  //cout << "NON repopulated distribution function: " << (nop_BeforeReceivePBC) <<" particles vs " << nop << " in this core " << endl;
  /* test only the particles added after repopulation */

  int END;
  bool ExtraCheck= false;
  END= nop_BeforeReceivePBC;
  if (END==0){    END=nop; ExtraCheck= true;}
  // so that the first time it prints initial distribution 

  for (long long i = 0; i < END; i++) {
    /** cosider only particles in the non-repopulation part of the grid -
	this only for the initial distribution, the others come automatically **/
    if (ExtraCheck){
      if (! (x[i]> -dx+PRACells*dx and x[i]<Lx+dx-PRACells*dx and y[i]> -dy+PRACells*dy and y[i]<Ly+dy-PRACells*dy and z[i]> -dz+PRACells*dx and z[i]<Lz+dz-PRACells*dz))

	continue;
    }

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
  MaxSizeCGFluidBCMsg= 0;
  int rank_local= vct->getCartesian_rank();  // on the grid communicator
  int HighestRank= XLEN*YLEN*ZLEN-1;
  MPI_Status status;
  int TAG_CG_RG= ns; // i need to put it here to be visible from both CG and RG

  /* phase 1: as a child */
  if (CommToParent_P != MPI_COMM_NULL) { // meaning: you are a child AND you want to receive PBC
    
    RGPBC_Info = new RGPBC_struct[MAX_RG_numPBCMessages];

    RG_MaxFluidMsgSize=0;
    if (FluidLikeRep){
      initWeightFluidPBC_Phase1(grid, vct, RGPBC_Info, &RG_numPBCMessages, &RG_MaxFluidMsgSize);
    }
    else{ // this is the option to use with -1 and -2       

      // change buildPBC accordingly
      initWeightPBC_Phase1(grid, vct, RGPBC_Info, &RG_numPBCMessages);
      // CAREFUL: initWeightPBC_Phase1_New does not work with many COARSE GRID cores- fix before using
      //initWeightPBC_Phase1_New(grid, vct, RGPBC_Info, &RG_numPBCMessages); 
    }

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
    // (but only if needed) -- these are used for kinetic repopulation
    if (RG_numPBCMessages>0 and !FluidLikeRep){
      PRGMsg= newArr2(RepP_struct, RG_numPBCMessages, sizeRG_PBCMsg);
      nopPRGMsg= new int[RG_numPBCMessages];
      PRGMsgArrived= new bool [RG_numPBCMessages];
      PRGMsg_General= new RepP_struct[sizeRG_PBCMsg];

    }

    // these are vectors needed to exchange repopulated particles inside the RG
    // needed only for kinetic repopulation
    if (RG_numPBCMessages>0 and !FluidLikeRep){
      if (rank_local==HighestRank){ 
	H_CRP_General= new CRP_struct[size_CRP];
	H_CRP_General_ptr= H_CRP_General; // for resize
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
    } // end if (RG_numPBCMessages>0 and !FluidLikeRep): buffers for exchange of repopulated particles

    if   (RG_numPBCMessages>0 and FluidLikeRep){
      // these are the vectors needed for the fluid repopulation
      rho_FBC = newArr2(double, RG_numPBCMessages, RG_MaxFluidMsgSize); 
      Jx_FBC = newArr2(double, RG_numPBCMessages, RG_MaxFluidMsgSize);
      Jy_FBC = newArr2(double, RG_numPBCMessages, RG_MaxFluidMsgSize);
      Jz_FBC = newArr2(double, RG_numPBCMessages, RG_MaxFluidMsgSize);
      Pxx_FBC = newArr2(double, RG_numPBCMessages, RG_MaxFluidMsgSize);
      Pxy_FBC = newArr2(double, RG_numPBCMessages, RG_MaxFluidMsgSize);
      Pxz_FBC = newArr2(double, RG_numPBCMessages, RG_MaxFluidMsgSize);
      Pyy_FBC = newArr2(double, RG_numPBCMessages, RG_MaxFluidMsgSize);
      Pyz_FBC = newArr2(double, RG_numPBCMessages, RG_MaxFluidMsgSize);
      Pzz_FBC = newArr2(double, RG_numPBCMessages, RG_MaxFluidMsgSize);

      // this is the vector used in the receives
      RGFluidMsg= new double[RG_MaxFluidMsgSize*numFBC];

      // this can be done better, but the rationale is to avoid to repopulate mulitple times per cell
      RHOP= newArr3(double, nxc, nyc, nzc);
      
      UP= newArr3(double, nxc, nyc, nzc);
      VP= newArr3(double, nxc, nyc, nzc);
      WP= newArr3(double, nxc, nyc, nzc);

      UTHP= newArr3(double, nxc, nyc, nzc);
      VTHP= newArr3(double, nxc, nyc, nzc);
      WTHP= newArr3(double, nxc, nyc, nzc);

      // for test
      for (int i=0; i<nxc; i++)
	for (int j=0; j<nyc; j++)
	  for (int k=0; k<nzc; k++)
	    RHOP[i][j][k]= -100.0;

    } // end if (RG_numPBCMessages>0 and !FluidLikeRep)

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
      MPI_Send(RGPBC_Info, RG_numPBCMessages+1, MPI_RGBC_struct, HighestRank, TAG_1b1c, vct->getComm());
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
	MPI_Recv(buffer_rcv, MAX_RG_numPBCMessages, MPI_RGBC_struct, MPI_ANY_SOURCE, TAG_1b1c, vct->getComm(), &status);
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
	if (recv_ThisMsg>0 and !FluidLikeRep){ // CRP ops only with kinetic repopulation
	  H_CRP_cores[src]= 1;
	}
	
      } // end for (int c=0; c< HighestRank; c++){

      if (!FluidLikeRep){ // CRP ops only with kinetic repopulation 
	num_H_CRP_cores=0;
	for (int i=0; i<XLEN*YLEN*ZLEN; i++){
	  num_H_CRP_cores+= H_CRP_cores[i];
	}
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
	MPI_Send(&(RGPBC_Info_ToCGCore[cg][0]), RG_numPBCMessages_ToCGCore[cg]+1, MPI_RGBC_struct, cg, TAG_CG_RG, CommToParent_P);
	//cout << "R " << rankAsChild << " on the PC communicator has just sent a msg to core " << cg << endl; 
      }

      delete[]RGPBC_Info_LevelWide;
      
      delArr2(RGPBC_Info_ToCGCore, ParentCoreNum);
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
      MPI_Recv(CG_buffer, MAX_RG_numPBCMessages, MPI_RGBC_struct, ChildHighestRank, TAG_CG_RG, CommToChild_P[ch], &status);

      //cout << "Grid " << numGrid << ", local rank on PC comm " << PCrank << " has recevied from ch " << ch << " tag " << TAG_CG_RG << endl;
      
      // process update the msg structure
      while (CG_buffer[CG_numPBCMessages[ch]].RG_core!= -1){
	//cout << "INSIDE WHILE: ch is " << ch << " CG_numPBCMessages[ch] is " << CG_numPBCMessages[ch] << " CG_buffer[CG_numPBCMessages[ch]].RG_core is " <<CG_buffer[CG_numPBCMessages[ch]].RG_core << endl;
	CG_Info[ch][CG_numPBCMessages[ch]]= CG_buffer[CG_numPBCMessages[ch]];

	/* used only for fluid PBC */
	int MsgSize=  CG_Info[ch][CG_numPBCMessages[ch]].np_x * CG_Info[ch][CG_numPBCMessages[ch]].np_y * CG_Info[ch][CG_numPBCMessages[ch]].np_z;
	if (MsgSize> MaxSizeCGFluidBCMsg) MaxSizeCGFluidBCMsg= MsgSize;
	/* end used only for fluid PBC */

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
    
    // MaxNumMsg >0: otherwise, i am not sending PBC and i don't need to allocate
    // !FluidLikeRep: these are needed only for kinetic repopulation
    if (MaxNumMsg >0 and !FluidLikeRep){       
      PCGMsg= newArr3(RepP_struct, numChildren, MaxNumMsg, sizeCG_PBCMsg);
      PCGMsg_ptr= PCGMsg; // for the resize
      nopPCGMsg = newArr2(int, numChildren, MaxNumMsg);
      //cout << "Grid " << numGrid << " core " << vct->getCartesian_rank() << " nopPCGMsg has sizes " << numChildren << " x " <<MaxNumMsg << endl;  
    }

    if (MaxNumMsg >0 and FluidLikeRep){
      CGFluidMsg= new double[MaxSizeCGFluidBCMsg* numFBC];
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

  /* finally instantiate H_CRP_Msg, but only the entries which are needed */
  // Fluid repopulation does not need communicate after repopulate
  if (CommToParent_P != MPI_COMM_NULL and RG_numPBCMessages>0 and rank_local==HighestRank and !FluidLikeRep) {
    H_CRP_Msg= newArr2_PA(CRP_struct, XLEN*YLEN*ZLEN, size_CRP, H_CRP_cores);      
    H_CRP_Msg_ptr= H_CRP_Msg;
    /*for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      cout << i << " of " << XLEN*YLEN*ZLEN <<": H_CRP_cores[i]: "<< H_CRP_cores[i] << endl;
      }*/
    /*for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      if (H_CRP_cores[i]!=1)
	cout << "Trying to provoke segm fault: " << H_CRP_Msg[i][3].q << endl;
	}*/
      

  }
  /* end finally instantiate H_CRP_Msg */

  
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
  
  // NB: AllowPMsgResize is overriden to false if fluid repopulation
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
  
  
  
  TrimInfoVector(vct);

  // allocate grid-sized vectors, according to the repopulation method choses

  MPI_Barrier(vct->getComm());
  if (vct->getCartesian_rank()==0){
    cout << "Grid " << numGrid <<" ns " << ns << " has finished initWeightPBC"<< endl;
  }
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
    // with extra RG dx
    /*i_s= 0-1; i_e= nxn-1+1;
    j_s= 0-1; j_e= nyn-1+1;
    k_s= PRA_ZLeft_Start-1; k_e= PRA_ZLeft_End+1;*/
    // just RG extension - add +/- 0.5 while building the CG msg
    i_s= 0; i_e= nxn-1;
    j_s= 0; j_e= nyn-1;
    k_s= PRA_ZLeft_Start; k_e= PRA_ZLeft_End;

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
    /*i_s=0-1; i_e= nxn-1+1;
    j_s=0-1; j_e= nyn-1+1;
    k_s= PRA_ZRight_Start-1; k_e= PRA_ZRight_End +1;*/
    i_s=0; i_e= nxn-1;
    j_s=0; j_e= nyn-1;
    k_s= PRA_ZRight_Start; k_e= PRA_ZRight_End;
   

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
    /*i_s= PRA_XLeft_Start-1; i_e= PRA_XLeft_End +1;
    j_s= 0-1; j_e= nyn-1+1;
    k_s= 0-1; k_e= nzn-1+1;*/
    i_s= PRA_XLeft_Start; i_e= PRA_XLeft_End;
    j_s= 0; j_e= nyn-1;
    k_s= 0; k_e= nzn-1;

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
    /*i_s= PRA_XRight_Start-1; i_e= PRA_XRight_End+1;
    j_s= 0-1; j_e= nyn-1+1;
    k_s= 0-1; k_e= nzn-1+1;*/
    i_s= PRA_XRight_Start; i_e= PRA_XRight_End;
    j_s= 0; j_e= nyn-1;
    k_s= 0; k_e= nzn-1;

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
    /*i_s= 0-1; i_e= nxn-1+1;
    j_s= PRA_YLeft_Start-1; j_e= PRA_YLeft_End+1;
    k_s= 0-1; k_e= nzn-1+1;*/
    i_s= 0; i_e= nxn-1;
    j_s= PRA_YLeft_Start; j_e= PRA_YLeft_End;
    k_s= 0; k_e= nzn-1;

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
    /*i_s= 0-1; i_e= nxn-1+1;
    j_s= PRA_YRight_Start-1; j_e= PRA_YRight_End+1;
    k_s= 0-1; k_e= nzn-1+1;*/
    i_s= 0; i_e= nxn-1;
    j_s= PRA_YRight_Start; j_e= PRA_YRight_End;
    k_s= 0; k_e= nzn-1;
    
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


/* Phase 1: RG cores build their side of the map for PBC,
differently form before, it also includes the DX/2, which should NOT be present anymore when building the CG msg-
this version should work also in the (unlikely) case that the extra DX/2 is outside the CG core in the msg */
void Particles3Dcomm::initWeightPBC_Phase1_New(Grid *grid, VirtualTopology3D * vct, RGPBC_struct *RGPBC_Info, int *RG_numPBCMessages){

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
    // just RG extension - add +/- 0.5 while building the CG msg
    /*i_s= 0; i_e= nxn-1;
    j_s= 0; j_e= nyn-1;
    k_s= PRA_ZLeft_Start; k_e= PRA_ZLeft_End;*/
    // RG + CG extension - remove that while building the CG msg
    i_s= 0-ceil(RFx/2); i_e= nxn-1+ceil(RFx/2);
    j_s= 0-ceil(RFy/2); j_e= nyn-1+ceil(RFy/2);
    k_s= PRA_ZLeft_Start-ceil(RFz/2); k_e= PRA_ZLeft_End+ceil(RFz/2);

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
    /*i_s=0; i_e= nxn-1;
    j_s=0; j_e= nyn-1;
    k_s= PRA_ZRight_Start; k_e= PRA_ZRight_End;*/
    i_s=0-ceil(RFx/2); i_e= nxn-1+ceil(RFx/2);
    j_s=0-ceil(RFy/2); j_e= nyn-1+ceil(RFy/2);
    k_s= PRA_ZRight_Start-ceil(RFz/2); k_e= PRA_ZRight_End+ceil(RFz/2);
   

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
    /*i_s= PRA_XLeft_Start; i_e= PRA_XLeft_End;
    j_s= 0; j_e= nyn-1;
    k_s= 0; k_e= nzn-1;*/
    i_s= PRA_XLeft_Start-ceil(RFx/2); i_e= PRA_XLeft_End+ceil(RFx/2);
    j_s= 0-ceil(RFy/2); j_e= nyn-1+ceil(RFy/2);
    k_s= 0-ceil(RFz/2); k_e= nzn-1+ceil(RFz/2);

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
    /*i_s= PRA_XRight_Start; i_e= PRA_XRight_End;
    j_s= 0; j_e= nyn-1;
    k_s= 0; k_e= nzn-1;*/
    i_s= PRA_XRight_Start-ceil(RFx/2); i_e= PRA_XRight_End+ceil(RFx/2);
    j_s= 0-ceil(RFy/2); j_e= nyn-1+ceil(RFy/2);
    k_s= 0-ceil(RFz/2); k_e= nzn-1+ceil(RFz/2);

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
    /*i_s= 0; i_e= nxn-1;
    j_s= PRA_YLeft_Start; j_e= PRA_YLeft_End;
    k_s= 0; k_e= nzn-1;*/
    i_s= 0-ceil(RFx/2); i_e= nxn-1+ceil(RFx/2);
    j_s= PRA_YLeft_Start-ceil(RFy/2); j_e= PRA_YLeft_End+ceil(RFy/2);
    k_s= 0-ceil(RFz/2); k_e= nzn-1+ceil(RFz/2);

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
    /*i_s= 0; i_e= nxn-1;
    j_s= PRA_YRight_Start; j_e= PRA_YRight_End;
    k_s= 0; k_e= nzn-1;*/
    i_s= 0-ceil(RFx/2); i_e= nxn-1+ceil(RFx/2);
    j_s= PRA_YRight_Start-ceil(RFy/2); j_e= PRA_YRight_End+ceil(RFy/2);
    k_s= 0-ceil(RFz/2); k_e= nzn-1+ceil(RFz/2);
    
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
	    MPI_Isend(&MaxGrid_sizeCG_PBCMsg, 1, MPI_INT, RcvList[ch][i], 200+ns, CommToChild_P[ch], &request);
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
void Particles3Dcomm::ReceivePBC_NoApply(Grid* grid, VirtualTopology3D * vct, int cycle){

  if (saveRepParFile and numGrid >0 and (cycle % DiagnosticsOutputCycle == 1)){
    ofstream my_file(RepopulatedPar.c_str(), fstream::app);
    my_file <<endl<< "Cycle "<< cycle <<endl;
    my_file.close();

  }
  
  // to be mover in: Before Applying 
  //nop_BeforeReceivePBC=nop;
  PRA_PAdded=0;

  int RR= vct->getCartesian_rank();

  // if not, exits
  if (CommToParent_P!= MPI_COMM_NULL and RG_numPBCMessages >0){

    MPI_Status status;

    if (AllowPMsgResize){ // to do before anybody has started receiving, so i don't have to copy info
      int NEW_sizePBCMsg;
      // as it is set now, each RG core receives a msg, so just do a rcv 
      MPI_Recv(&NEW_sizePBCMsg, 1, MPI_INT, MPI_ANY_SOURCE, 200+ns, CommToParent_P, &status);

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


  }
}


void Particles3Dcomm::ApplyPBC_NoReceive_NoCommunicate(Grid* grid, VirtualTopology3D * vct, int cycle){
  

  nop_BeforeReceivePBC=nop;

  int RR= vct->getCartesian_rank();

  // if not, exits
  if (CommToParent_P!= MPI_COMM_NULL and RG_numPBCMessages >0){
    
    // here, PBC have been successfully received
    // now I have to apply them (split the particles )
    ApplyPBC(vct, grid, cycle);


  }
  
  nop_AfterReceivePBC= nop;

  //cout << "Grid " << numGrid <<" R " << vct->getCartesian_rank() << " ns " << ns << " added " << PRA_PAdded << " particles, deleted " << PRA_deleted << endl;

  // here, check how many PBC particles have been added at grid level
  if (false){
    int ToTPRA_PAdded;
    int TotPRA_deleted;
    int TotNOP;
    MPI_Allreduce(&PRA_PAdded, &ToTPRA_PAdded, 1, MPI_INT, MPI_SUM, vct->getComm());
    MPI_Allreduce(&PRA_deleted, &TotPRA_deleted, 1, MPI_INT, MPI_SUM, vct->getComm());
    MPI_Allreduce(&nop, &TotNOP, 1, MPI_INT, MPI_SUM, vct->getComm());
    if (RR== XLEN*YLEN*ZLEN-1){
      cout << "Grid " << numGrid << " ns " << ns << " added " << ToTPRA_PAdded << " particles, deleted " << TotPRA_deleted << " (this before communicateRepopulatedParticles)" << ", total nop: " << TotNOP << endl;
    }
  } // end debug stuff

  

}


void Particles3Dcomm::communicateRepopulatedParticles_Wrap(Grid* grid, VirtualTopology3D * vct, int cycle){

  if (CommToParent_P!= MPI_COMM_NULL and XLEN*YLEN*ZLEN>1) {
    communicateRepopulatedParticles(grid, vct);
  }

}

/* RG receives PBC msg and acts accordingly */
void Particles3Dcomm::ReceivePBC(Grid* grid, VirtualTopology3D * vct, int cycle){

  ReceivePBC_NoApply(grid, vct, cycle);
  
  ApplyPBC_NoReceive_NoCommunicate(grid, vct, cycle);

  communicateRepopulatedParticles_Wrap(grid, vct, cycle);

}

void Particles3Dcomm::buildPBCMsg(Grid* grid, VirtualTopology3D * vct, int ch){
  
  /* returns the number - not the level or the order in the children vector - of the child grid n */
  int childNum= vct->getChildGridNum(ch);
  // This is a fundamental mistake! i need to pass 0.5 DX more
  // these are the dx, dy, dz of the child
  //double dx= grid->getDx_mlmd(childNum);
  //double dy= grid->getDy_mlmd(childNum);
  //double dz= grid->getDz_mlmd(childNum);

  double x_min, x_max, y_min, y_max, z_min, z_max;
  /*double EXTRAL=0; // 0 should be the "good" value; put extra for experimenting
  EXTRAL=-0.5; // like this i repopulate all with RF=8
  double EXTRAR=0;
  EXTRAR=-0.5;
  //EXTRA=0.5;
  cout << "EXTRAL: " << EXTRAL << " EXTRAR: " << EXTRAR << endl;*/

  if (CG_numPBCMessages[ch]>0){ // this particular core has to send BC

    // pre
    for (int m=0; m< CG_numPBCMessages[ch]; m++){
      nopPCGMsg[ch][m]= 0;
    } 

    // core
    for (int p=0; p<nop; p++){
      // here, i have to check all the msgs in else if, to make sure that I am sending 
      // this particular particle only to one core
      
      // first, coarse check to avoid spending too much time here
      if (x[p]< Ox_Ch[ch]- dx or x[p]> Ox_Ch[ch]+Lx_Ch[ch] + dx or y[p]< Oy_Ch[ch]- dy or y[p]> Oy_Ch[ch]+Ly_Ch[ch] + dy or z[p]< Oz_Ch[ch]- dz or z[p]> Oz_Ch[ch]+Lz_Ch[ch] + dz)
	continue;
      

      for (int n=0; n< CG_numPBCMessages[ch]; n++){ 


	// NB: dx*(CG_Info[ch][n].np_x) gives you the extra DX; then i remove 1/2

	// use this with initWeightPBC_Phase1
	x_min= CG_Info[ch][n].CG_x_first -0.5*dx;
	x_max= CG_Info[ch][n].CG_x_first+ dx*(CG_Info[ch][n].np_x) -0.5*dx;
	y_min= CG_Info[ch][n].CG_y_first -0.5*dy;
	y_max= CG_Info[ch][n].CG_y_first+ dy*(CG_Info[ch][n].np_y) -0.5*dy;
	z_min= CG_Info[ch][n].CG_z_first -0.5*dz;
	z_max= CG_Info[ch][n].CG_z_first+ dz*(CG_Info[ch][n].np_z) -0.5*dz;
	
	// with this, there is a risk if x_min is < xstart, x_max > xend
	// put a warning and exit in case
	// this is needed only with initWeightPBC_Phase1, not with _New (when it works)
	if (x_min < xstart and (Ox_Ch[ch]> xstart and Ox_Ch[ch]-0.5*dx < xstart ) ) {
	  cout << "Move child " << ch <<" of grid " << numGrid << " to the right in x"<<endl;
	  cout << "Aborting now" <<endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);}
	
	if (x_max > xend and (Ox_Ch[ch] + Lx_Ch[ch] < xend and Ox_Ch[ch] + Lx_Ch[ch]+ 0.5*dx > xend)){
	  cout << "Move child " << ch <<" of grid " << numGrid << " to the left in x"<<endl;
	  cout << "Aborting now" <<endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);}

	if (y_min < ystart and (Oy_Ch[ch]> ystart and Oy_Ch[ch]-0.5*dy < ystart ) ) {
	  cout << "Move child " << ch <<" of grid " << numGrid << " to the right in y"<<endl;
	  cout << "Aborting now" <<endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);}
	
	if (y_max > yend and (Oy_Ch[ch] + Ly_Ch[ch] < yend and Oy_Ch[ch] + Ly_Ch[ch]+ 0.5*dy > yend)){
	  cout << "Move child " << ch <<" of grid " << numGrid << " to the left in y"<<endl;
	  cout << "Aborting now" <<endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);}

	if (z_min < zstart and (Oz_Ch[ch]> zstart and Oz_Ch[ch]-0.5*dz < zstart ) ) {
	  cout << "Move child " << ch <<" of grid " << numGrid << " to the right in z"<<endl;
	  cout << "Aborting now" <<endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);}
	
	if (z_max > zend and (Oz_Ch[ch] + Lz_Ch[ch] < zend and Oz_Ch[ch] + Lz_Ch[ch]+ 0.5*dz > zend)){
	  cout << "Move child " << ch <<" of grid " << numGrid << " to the left in z"<<endl;
	  cout << "z_max " << z_max << " Lz_Ch[ch]+ 0.5*dz: " << Lz_Ch[ch]+ 0.5*dz << endl;
	  cout << "Aborting now" <<endl;
	  MPI_Abort(MPI_COMM_WORLD, -1);}
	
	  

	
	// use this with initWeightPBC_Phase1_New: careful at the moment it does not work with multiple cores	
	/*x_min= CG_Info[ch][n].CG_x_first ;
	x_max= CG_Info[ch][n].CG_x_first+ dx*(CG_Info[ch][n].np_x) -dx;
	y_min= CG_Info[ch][n].CG_y_first ;
	y_max= CG_Info[ch][n].CG_y_first+ dy*(CG_Info[ch][n].np_y) -dy;
	z_min= CG_Info[ch][n].CG_z_first ;
	z_max= CG_Info[ch][n].CG_z_first+ dz*(CG_Info[ch][n].np_z) -dz;*/

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

void Particles3Dcomm::ApplyPBC(VirtualTopology3D* vct, Grid* grid, int cycle){

  int RR= vct->getCartesian_rank();
  /*if (RR==0){
    cout << "G" << numGrid <<"R" << RR << " INSIDE APPLYPBC" << endl;
    }*/
  //cout << "R" << RR <<" G" << numGrid << ", nop before applying BC: " << nop <<endl;

  for (int m=0; m< RG_numPBCMessages; m++){
    int count =0; 

    while (fabs(PRGMsg[m][count].q) > 1.1*EXIT_VAL){
      //cout << "m: " <<m << " count "<< count <<" PRGMsg[m][count].q: " << PRGMsg[m][count].q << endl;
      SplitPBC(vct, grid, PRGMsg[m][count], cycle);
      count ++;
    }
    nopPRGMsg[m]= count;

    //cout <<"G"<<numGrid << "R" <<RR << " ns " << ns  << " m " << m << " nopPRGMsg[m] " << nopPRGMsg[m] << Endl;
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
void Particles3Dcomm::SplitPBC(VirtualTopology3D * vct, Grid* grid, RepP_struct p, int cycle){
  // to prevent the repopulation of particles which would try to accumulate outside the grid
  double inv_dx= 1./dx;
  double inv_dy= 1./dy;
  double inv_dz= 1./dz;

  // there is no need to add an extra buffer here; keep it 0
  double PM= 0.0; 

  double xTmp;
  double yTmp;
  double zTmp;
  double qTmp;
  
  bool StX, StY, StZ;

  double Sign= -1;
  double *Perc = new double[4];
  Perc[0]= 0.05; 
  Perc[1]= 0.01;
  Perc[2]= 0.15;
  Perc[3]= 0.2;
  // with any of the two true, nothing seems to change
  bool KickVelocity= false;
  bool KickPosition= false;
  int Index;

  for (int i=0; i< ceil(RFx); i++){
    xTmp= p.x - DxP/2.0 + dx*(1./2. + i)- grid->getOx();
        
    // if outside the domain, no  point in continuing splitting
    if (xTmp < Coord_XLeft_Start-dx*PM or xTmp>Coord_XRight_End+dx*PM)
      continue;
    
    // am i inside a X PRA?
    StX= (xTmp> Coord_XLeft_Start-dx*PM and xTmp <= Coord_XLeft_End) or (xTmp>= Coord_XRight_Start and xTmp < Coord_XRight_End+dx*PM);
    //if (! StX ) continue; DO NOT DARE UNCOMMENTING THIS; I will lose particles which are not in the X PRA, but are in an other

    for (int j=0; j< ceil(RFy); j++){
      yTmp= p.y - DyP/2.0 + dy*(1./2. + j)- grid->getOy();

      // if outside the domain, no  point in continuing splitting
      if (yTmp < Coord_YLeft_Start-dy*PM or yTmp>Coord_YRight_End+dy*PM)
	continue;
      
      // am i inside a Y PRA?
      StY= (yTmp> Coord_YLeft_Start-dy*PM and yTmp <= Coord_YLeft_End) or (yTmp>= Coord_YRight_Start and yTmp < Coord_YRight_End+dy*PM);
      //if (! StY ) continue; DO NOT DARE UNCOMMENTING THIS; I will lose particles which are in another PRA but not in thos

      for (int k=0; k< ceil(RFz); k++){
	zTmp= p.z - DzP/2.0 + dz*(1./2. + k)- grid->getOz();
	
	// if outside the domain, no  point in continuing splitting
	if (zTmp < Coord_ZLeft_Start-dz*PM or zTmp>Coord_ZRight_End+dz*PM)
	  continue;

	// am i inside a Z PRA?
	StZ= (zTmp> Coord_ZLeft_Start-dz*PM and zTmp <= Coord_ZLeft_End) or (zTmp>= Coord_ZRight_Start and zTmp < Coord_ZRight_End+dz*PM);

	
	//if (! StZ ) continue; DO NOT DARE UNCOMMENTING THIS
	

	// if i am inside any PRA (if i arrived here i am inside the extended domain)
	bool Keep= false;
	// i should really separate by direction but i will do that later
	if (bcPfaceXleft== -1 or bcPfaceXleft== -4 ) { 
	  /* this is when i repopulate all the PRA
	     if you arrive here, you should be kept */
	  Keep= true;
	}
	else if (bcPfaceXleft== -2){ // here i keep only the particles that have entered in the last dt
	  Keep= RepopulatedParticleHasEnteredRG(grid, xTmp, yTmp, zTmp, p.u, p.v, p.w);
	}
	
	// this is an OR not an AND; notice i have already filtered out particles
	// outside the RG
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
	  if (KickVelocity){
	    if (Sign== -1){
	      Index= rand() % 4;
	    }
	    //cout << "old u: " << u[nop] <<"; " ;
	    u[nop]= u[nop]*(1+Sign*Perc[Index]);
	    v[nop]= v[nop]*(1+Sign*Perc[Index]);
	    w[nop]= w[nop]*(1+Sign*Perc[Index]);
	    Sign= Sign*(-1); // for the next round
	    //cout << "new u: " << u[nop] << endl;

	    /*if (u[nop]== 0){
	      cout << "u: " << u[nop] << ", v: "<< v[nop] << ", w: " << w[nop] << ", x: " <<x[nop] << ", y: " <<y[nop] << ", z: "<< z[nop] << endl;
	      cout << "Original p: u: " << p.u << ", v: " << p.v << ", w: "<< p.w <<endl;
	      }*/
	  }
	  if (KickPosition){
	    const double ixd = floor((x[nop] - xstart) * inv_dx);
	    const double iyd = floor((y[nop] - ystart) * inv_dy);
	    const double izd = floor((z[nop] - zstart) * inv_dz);
	    int ix = 2 + int (ixd)-1; // -1, not sure about this
	    int iy = 2 + int (iyd)-1;
	    int iz = 2 + int (izd)-1;
	    
	    double rX = ((double)rand() / (double)(RAND_MAX));
	    double rY = ((double)rand() / (double)(RAND_MAX));
	    double rZ = ((double)rand() / (double)(RAND_MAX));

	    x[nop]= grid->getXN_XT(ix, iy, iz)+ dx*rX;
	    y[nop]= grid->getYN_XT(ix, iy, iz)+ dy*rY;
	    z[nop]= grid->getZN_XT(ix, iy, iz)+ dz*rZ;

	    const double ixd2 = floor((x[nop] - xstart) * inv_dx);
	    const double iyd2 = floor((y[nop] - ystart) * inv_dy);
	    const double izd2 = floor((z[nop] - zstart) * inv_dz);

	    if (ixd != ixd2 or iyd != iyd2 or izd != izd2){
	      cout <<"ixd: " << ixd <<", ixd2: " << ixd2 << endl;
	      cout <<"iyd: " << iyd <<", iyd2: " << iyd2 << endl;
	      cout <<"izd: " << izd <<", izd2: " << izd2 << endl;
	    }
	    
	    /* ok, this i check and nobody ends up here
	    if (x[nop]< -dx or y[nop]< -dy or z[nop]< -dz){
	      cout << "The particle is outside the PRA " << endl;
	      }*/
	  }
                                          
	  nop++; 
	  //cout << "Particle accepted, x: "<<x[nop] << ", y: " << y[nop] <<", z: "<< z[nop] << "(Lz: " << Lz << ")" << endl;

	  /* save data regarding this particle */
	  if (saveRepParFile and (cycle % DiagnosticsOutputCycle == 1)){
	    ofstream my_file(RepopulatedPar.c_str(), fstream::app);
	    my_file << endl << x[nop] <<" "<< y[nop] <<" "<< z[nop] <<" "<< u[nop] <<" "<< v[nop] <<" "<< w[nop] <<" "<< q[nop] <<" "  <<" ";
	    my_file.close();

	  }
	  /* end save data regarding this particle */
	 	 
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

  delete[]Perc;

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
    //MPI_Abort(MPI_COMM_WORLD, -1);
    cout << "now i am continuing, but check the inconsistency..." << endl; return;
  }
  // check if the particle belongs to this core- use xStart_GC, etc because must particles are going to be in GC
  bool RightCore= (p.x >= xStart_GC and p.x <= xEnd_GC) and (p.y >= yStart_GC and p.y <= yEnd_GC) and (p.z >= zStart_GC and p.z <= zEnd_GC);

  if (RightCore== false){
    cout << "Grid " << numGrid <<": FATAL ERROR in unpack_CRP, aborting ..." <<endl;
    //MPI_Abort(MPI_COMM_WORLD, -1);
    cout << "now i am continuing, but check the inconsistency..." << endl; return;
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
    CRP_struct ** H_CRP_Msg_tmp= newArr2_PA(CRP_struct, XLEN*YLEN*ZLEN, size_CRP, H_CRP_cores); 
    
    for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      if (H_CRP_cores[i]!= 1) continue;
      for (int j=0; j< H_num_CRP_nop[i]; j++)
	H_CRP_Msg_tmp[i][j]= H_CRP_Msg_ptr[i][j];
    }

    delArr2(H_CRP_Msg, XLEN*YLEN*ZLEN);
    H_CRP_Msg= newArr2_PA(CRP_struct, XLEN*YLEN*ZLEN, NEW_size_CRP, H_CRP_cores); 
    H_CRP_Msg_ptr= H_CRP_Msg;
    
    for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      if (H_CRP_cores[i]!= 1) continue;
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
    CRP_struct ** H_CRP_Msg_tmp= newArr2_PA(CRP_struct, XLEN*YLEN*ZLEN, size_CRP, H_CRP_cores); 
    
    for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      if (H_CRP_cores[i]!= 1) continue;
      for (int j=0; j< H_num_CRP_nop[i]; j++)
	H_CRP_Msg_tmp[i][j]= H_CRP_Msg_ptr[i][j];
    }

    delArr2(H_CRP_Msg, XLEN*YLEN*ZLEN);
    H_CRP_Msg= newArr2_PA(CRP_struct, XLEN*YLEN*ZLEN, NewSize, H_CRP_cores); 
    H_CRP_Msg_ptr= H_CRP_Msg;
    
    for (int i=0; i< XLEN*YLEN*ZLEN; i++){
      if (H_CRP_cores[i]!= 1) continue;
      for (int j=0; j<  H_num_CRP_nop[i]; j++)
	H_CRP_Msg[i][j]= H_CRP_Msg_tmp[i][j];
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

/*bool Particles3Dcomm::RepopulatedParticleHasEnteredRG(Grid* grid, double xTmp, double yTmp, double zTmp, double u, double v, double w){
  // implementation based on the old mover //

  const double dto2 = .5 * dt, qomdt2 = qom * dto2 / c;
  const double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy, inv_dz = 1.0 / dz;
  // don't bother trying to push any particles simultaneously;                                                                          
  // MIC already does vectorization automatically, and trying                                                                           
  // to do it by hand only hurts performance.                                                                                           
  double xp = xTmp;
  double yp = yTmp;
  double zp = zTmp;
  double up = u;
  double vp = v;
  double wp = w;
  const double xptilde = xTmp;
  const double yptilde = yTmp;
  const double zptilde = zTmp;
  double uptilde;
  double vptilde;
  double wptilde;
  // calculate the average velocity iteratively
  for (int innter = 0; innter < NiterMover; innter++) {
    // interpolation G-->P                                                                                                            
    const double ixd = floor((xp - xstart) * inv_dx);
    const double iyd = floor((yp - ystart) * inv_dy);
    const double izd = floor((zp - zstart) * inv_dz);
    int ix = 2 + int (ixd);
    int iy = 2 + int (iyd);
    int iz = 2 + int (izd);
    if (ix < 1)
      ix = 1;
    if (iy < 1)
      iy = 1;
    if (iz < 1)
      iz = 1;
    if (ix > nxn - 1)
      ix = nxn - 1;
    if (iy > nyn - 1)
      iy = nyn - 1;
    if (iz > nzn - 1)
      iz = nzn - 1;
    
    double xi  [2];
    double eta [2];
    double zeta[2];
    
    xi  [0] = xp - grid->getXN(ix-1,iy  ,iz  );
    eta [0] = yp - grid->getYN(ix  ,iy-1,iz  );
    zeta[0] = zp - grid->getZN(ix  ,iy  ,iz-1);
    xi  [1] = grid->getXN(ix,iy,iz) - xp;
    eta [1] = grid->getYN(ix,iy,iz) - yp;
    zeta[1] = grid->getZN(ix,iy,iz) - zp;
    
    double Exl = 0.0;
    double Eyl = 0.0;
    double Ezl = 0.0;
    double Bxl = 0.0;
    double Byl = 0.0;
    double Bzl = 0.0;
    
    const double weight000 = xi[0] * eta[0] * zeta[0] * invVOL;
    const double weight001 = xi[0] * eta[0] * zeta[1] * invVOL;
    const double weight010 = xi[0] * eta[1] * zeta[0] * invVOL;
    const double weight011 = xi[0] * eta[1] * zeta[1] * invVOL;
    const double weight100 = xi[1] * eta[0] * zeta[0] * invVOL;
    const double weight101 = xi[1] * eta[0] * zeta[1] * invVOL;
    const double weight110 = xi[1] * eta[1] * zeta[0] * invVOL;
    const double weight111 = xi[1] * eta[1] * zeta[1] * invVOL;
    //                                                                                                                                
    Bxl += weight000 * (Bx[ix][iy][iz]             + Fext*Bx_ext[ix][iy][iz]);
    Bxl += weight001 * (Bx[ix][iy][iz - 1]         + Fext*Bx_ext[ix][iy][iz-1]);
    Bxl += weight010 * (Bx[ix][iy - 1][iz]         + Fext*Bx_ext[ix][iy-1][iz]);
    Bxl += weight011 * (Bx[ix][iy - 1][iz - 1]     + Fext*Bx_ext[ix][iy-1][iz-1]);
    Bxl += weight100 * (Bx[ix - 1][iy][iz]         + Fext*Bx_ext[ix-1][iy][iz]);
    Bxl += weight101 * (Bx[ix - 1][iy][iz - 1]     + Fext*Bx_ext[ix-1][iy][iz-1]);
    Bxl += weight110 * (Bx[ix - 1][iy - 1][iz]     + Fext*Bx_ext[ix-1][iy-1][iz]);
    Bxl += weight111 * (Bx[ix - 1][iy - 1][iz - 1] + Fext*Bx_ext[ix-1][iy-1][iz-1]);
    //                                                                                                                                
    Byl += weight000 * (By[ix][iy][iz]             + Fext*By_ext[ix][iy][iz]);
    Byl += weight001 * (By[ix][iy][iz - 1]         + Fext*By_ext[ix][iy][iz-1]);
    Byl += weight010 * (By[ix][iy - 1][iz]         + Fext*By_ext[ix][iy-1][iz]);
    Byl += weight011 * (By[ix][iy - 1][iz - 1]     + Fext*By_ext[ix][iy-1][iz-1]);
    Byl += weight100 * (By[ix - 1][iy][iz]         + Fext*By_ext[ix-1][iy][iz]);
    Byl += weight101 * (By[ix - 1][iy][iz - 1]     + Fext*By_ext[ix-1][iy][iz-1]);
    Byl += weight110 * (By[ix - 1][iy - 1][iz]     + Fext*By_ext[ix-1][iy-1][iz]);
    Byl += weight111 * (By[ix - 1][iy - 1][iz - 1] + Fext*By_ext[ix-1][iy-1][iz-1]);
    //                                                                                                                                
    Bzl += weight000 * (Bz[ix][iy][iz]             + Fext*Bz_ext[ix][iy][iz]);
    Bzl += weight001 * (Bz[ix][iy][iz - 1]         + Fext*Bz_ext[ix][iy][iz-1]);
    Bzl += weight010 * (Bz[ix][iy - 1][iz]         + Fext*Bz_ext[ix][iy-1][iz]);
    Bzl += weight011 * (Bz[ix][iy - 1][iz - 1]     + Fext*Bz_ext[ix][iy-1][iz-1]);
    Bzl += weight100 * (Bz[ix - 1][iy][iz]         + Fext*Bz_ext[ix-1][iy][iz]);
    Bzl += weight101 * (Bz[ix - 1][iy][iz - 1]     + Fext*Bz_ext[ix-1][iy][iz-1]);
    Bzl += weight110 * (Bz[ix - 1][iy - 1][iz]     + Fext*Bz_ext[ix-1][iy-1][iz]);
    Bzl += weight111 * (Bz[ix - 1][iy - 1][iz - 1] + Fext*Bz_ext[ix-1][iy-1][iz-1]);
    
    Exl += weight000 * Ex[ix][iy][iz];
    Exl += weight001 * Ex[ix][iy][iz - 1];
    Exl += weight010 * Ex[ix][iy - 1][iz];
    Exl += weight011 * Ex[ix][iy - 1][iz - 1];
    Exl += weight100 * Ex[ix - 1][iy][iz];
    Exl += weight101 * Ex[ix - 1][iy][iz - 1];
    Exl += weight110 * Ex[ix - 1][iy - 1][iz];
    Exl += weight111 * Ex[ix - 1][iy - 1][iz - 1];
    //                                                                                                                                
    Eyl += weight000 * Ey[ix][iy][iz];
    Eyl += weight001 * Ey[ix][iy][iz - 1];
    Eyl += weight010 * Ey[ix][iy - 1][iz];
    Eyl += weight011 * Ey[ix][iy - 1][iz - 1];
    Eyl += weight100 * Ey[ix - 1][iy][iz];
    Eyl += weight101 * Ey[ix - 1][iy][iz - 1];
    Eyl += weight110 * Ey[ix - 1][iy - 1][iz];
    Eyl += weight111 * Ey[ix - 1][iy - 1][iz - 1];
    
    Ezl += weight000 * Ez[ix][iy][iz];
    Ezl += weight001 * Ez[ix][iy][iz - 1];
    Ezl += weight010 * Ez[ix][iy - 1][iz];
    Ezl += weight011 * Ez[ix][iy - 1][iz - 1];
    Ezl += weight100 * Ez[ix - 1][iy][iz];
    Ezl += weight101 * Ez[ix - 1][iy][iz - 1];
    Ezl += weight110 * Ez[ix - 1][iy - 1][iz];
    Ezl += weight111 * Ez[ix - 1][iy - 1][iz - 1];
    
    // end interpolation                                                                                                              
    const double omdtsq = qomdt2 * qomdt2 * (Bxl * Bxl + Byl * Byl + Bzl * Bzl);
    const double denom = 1.0 / (1.0 + omdtsq);
    // solve the position equation                                                                                                    
    //const double ut = up + qomdt2 * Exl;
    //const double vt = vp + qomdt2 * Eyl;
    //const double wt = wp + qomdt2 * Ezl;

    const double ut = up - qomdt2 * Exl;
    const double vt = vp - qomdt2 * Eyl;
    const double wt = wp - qomdt2 * Ezl;
    const double udotb = ut * Bxl + vt * Byl + wt * Bzl;
    // solve the velocity equation                                                                                                    
    //uptilde = (ut + qomdt2 * (vt * Bzl - wt * Byl + qomdt2 * udotb * Bxl)) * denom;
    //vptilde = (vt + qomdt2 * (wt * Bxl - ut * Bzl + qomdt2 * udotb * Byl)) * denom;
    //wptilde = (wt + qomdt2 * (ut * Byl - vt * Bxl + qomdt2 * udotb * Bzl)) * denom;
    uptilde = (ut - qomdt2 * (vt * Bzl - wt * Byl + qomdt2 * udotb * Bxl)) * denom;
    vptilde = (vt - qomdt2 * (wt * Bxl - ut * Bzl + qomdt2 * udotb * Byl)) * denom;
    wptilde = (wt - qomdt2 * (ut * Byl - vt * Bxl + qomdt2 * udotb * Bzl)) * denom;
    // update position                                                                                                                
    //xp = xptilde + uptilde * dto2;
    //yp = yptilde + vptilde * dto2;
    //zp = zptilde + wptilde * dto2;
    xp = xptilde - uptilde * dto2;
    yp = yptilde - vptilde * dto2;
    zp = zptilde - wptilde * dto2;
  }                          // end of iteration                                                                                     
  // update the final position
  //xp = xptilde + uptilde * dt;
  //yp = yptilde + vptilde * dt;
  //zp = zptilde + wptilde * dt;
  xp = xptilde - uptilde * dt;
  yp = yptilde - vptilde * dt;
  zp = zptilde - wptilde * dt;
 

  // if the reconstructed past particle was inside, kill it
  if ((xp > Coord_XLeft_Start and xp < Coord_XRight_End) and (yp > Coord_YLeft_Start and yp < Coord_YRight_End) and (zp > Coord_ZLeft_Start and zp < Coord_ZRight_End) )
  return false;

  return true;
 

}*/


bool Particles3Dcomm::RepopulatedParticleHasEnteredRG(Grid* grid, double xTmp, double yTmp, double zTmp, double u, double v, double w){
   // implementation based on the new mover; tested vs the implementation based on the old mover; they are the same

  // this is for repopulation method -2
  // i want to check if the repopulated particle has entered during the last dt
  // i return a false if the particle was already inside RG last dt
  // otherwise true

  // this method is too simplestic and does not work
  //double xB= xTmp- u*dt;
  //double yB= yTmp- v*dt;
  //double zB= zTmp- w*dt;

  //if ((xB > Coord_XLeft_Start and xB < Coord_XRight_End) and (yB > Coord_YLeft_Start and yB < Coord_YRight_End) and (zB > Coord_ZLeft_Start and zB < Coord_ZRight_End) )
  //return false;

  // copied from mover BUT in the mover i have xn, vn and i want to go to xn+1 vn+1; here the opposite
  // so i have to put -'s (see calculations)
  double xp = xTmp;
  double yp = yTmp;
  double zp = zTmp;
  double up = u;
  double vp = v;
  double wp = w;
  const double xptilde = xTmp;
  const double yptilde = yTmp;
  const double zptilde = zTmp;
  double uptilde;
  double vptilde;
  double wptilde;

  double Exl = 0.0;
  double Eyl = 0.0;
  double Ezl = 0.0;
  double Bxl = 0.0;
  double Byl = 0.0;
  double Bzl = 0.0;
  int ix;
  int iy;
  int iz;

  double weights[2][2][2];

  // BEGIN OF SUBCYCLING LOOP

  get_weights(grid, xp, yp, zp, ix, iy, iz, weights);
  get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);

  const double B_mag      = sqrt(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
  double       dt_sub     = M_PI*c/(4*abs(qom)*B_mag);
  const int    sub_cycles = int(dt/dt_sub) + 1;

  dt_sub = dt/double(sub_cycles);
    
  const double dto2 = .5 * dt_sub, qomdt2 = qom * dto2 / c;
    
  // if (sub_cycles>1) cout << " >> sub_cycles = " << sub_cycles << endl;

  for (int cyc_cnt = 0; cyc_cnt < sub_cycles; cyc_cnt++) {

    // calculate the average velocity iteratively
    int nit = NiterMover;
    if (sub_cycles > 2*NiterMover) nit = 1;

    for (int innter = 0; innter < nit; innter++) {
      // interpolation G-->P

      get_weights(grid, xp, yp, zp, ix, iy, iz, weights);
      get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);
      get_El(weights, ix, iy, iz, Exl, Eyl, Ezl, Ex, Ey, Ez);

      // end interpolation
      const double omdtsq = qomdt2 * qomdt2 * (Bxl * Bxl + Byl * Byl + Bzl * Bzl);
      const double denom = 1.0 / (1.0 + omdtsq);
      // solve the position equation
      //const double ut = up + qomdt2 * Exl;
      //const double vt = vp + qomdt2 * Eyl;
      //const double wt = wp + qomdt2 * Ezl;
      // BACKWARDS
      const double ut = up - qomdt2 * Exl;
      const double vt = vp - qomdt2 * Eyl;
      const double wt = wp - qomdt2 * Ezl;

      const double udotb = ut * Bxl + vt * Byl + wt * Bzl;
      // solve the velocity equation 
      //uptilde = (ut + qomdt2 * (vt * Bzl - wt * Byl + qomdt2 * udotb * Bxl)) * denom;
      //vptilde = (vt + qomdt2 * (wt * Bxl - ut * Bzl + qomdt2 * udotb * Byl)) * denom;
      //wptilde = (wt + qomdt2 * (ut * Byl - vt * Bxl + qomdt2 * udotb * Bzl)) * denom;
      // BACKWARDS
      uptilde = (ut - qomdt2 * (vt * Bzl - wt * Byl + qomdt2 * udotb * Bxl)) * denom;
      vptilde = (vt - qomdt2 * (wt * Bxl - ut * Bzl + qomdt2 * udotb * Byl)) * denom;
      wptilde = (wt - qomdt2 * (ut * Byl - vt * Bxl + qomdt2 * udotb * Bzl)) * denom;
      // update position
      //xp = xptilde + uptilde * dto2;
      //yp = yptilde + vptilde * dto2;
      //zp = zptilde + wptilde * dto2;
      // BACKWARDS
      xp = xptilde - uptilde * dto2;
      yp = yptilde - vptilde * dto2;
      zp = zptilde - wptilde * dto2;
    }                           // end of iteration
    // update the final position and velocity
    // BACKWARDS -> i don't need velocity
    //up = 2.0 * uptilde - u[rest];
    //vp = 2.0 * vptilde - v[rest];
    //wp = 2.0 * wptilde - w[rest];
    //xp = xptilde + uptilde * dt_sub;
    //yp = yptilde + vptilde * dt_sub;
    //zp = zptilde + wptilde * dt_sub;
    // BACKWARDS
    xp = xptilde - uptilde * dt_sub;
    yp = yptilde - vptilde * dt_sub;
    zp = zptilde - wptilde * dt_sub;
  } // END  OF SUBCYCLING LOOP
  // end copied from mover

  // if the reconstructed past particle was inside, kill it
  if ((xp > Coord_XLeft_Start and xp < Coord_XRight_End) and (yp > Coord_YLeft_Start and yp < Coord_YRight_End) and (zp > Coord_ZLeft_Start and zp < Coord_ZRight_End) )
  return false;

  return true;
  
}


void Particles3Dcomm::get_weights(Grid * grid, double xp, double yp, double zp, int& ix, int& iy, int& iz, double weights[2][2][2]){

  const double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy, inv_dz = 1.0 / dz;
  const double ixd = floor((xp - xstart) * inv_dx);
  const double iyd = floor((yp - ystart) * inv_dy);
  const double izd = floor((zp - zstart) * inv_dz);

  ix = 2 + int (ixd);
  iy = 2 + int (iyd);
  iz = 2 + int (izd);

  if (ix < 1)
    ix = 1;
  if (iy < 1)
    iy = 1;
  if (iz < 1)
    iz = 1;
  if (ix > nxn - 1)
    ix = nxn - 1;
  if (iy > nyn - 1)
    iy = nyn - 1;
  if (iz > nzn - 1)
    iz = nzn - 1;

  double xi  [2];
  double eta [2];
  double zeta[2];

  xi  [0] = xp - grid->getXN(ix-1,iy  ,iz  );
  eta [0] = yp - grid->getYN(ix  ,iy-1,iz  );
  zeta[0] = zp - grid->getZN(ix  ,iy  ,iz-1);
  xi  [1] = grid->getXN(ix,iy,iz) - xp;
  eta [1] = grid->getYN(ix,iy,iz) - yp;
  zeta[1] = grid->getZN(ix,iy,iz) - zp;

  for (int ii = 0; ii < 2; ii++)
    for (int jj = 0; jj < 2; jj++)
      for (int kk = 0; kk < 2; kk++)
        weights[ii][jj][kk] = xi[ii] * eta[jj] * zeta[kk] * invVOL;
}


void Particles3Dcomm::get_Bl(const double weights[2][2][2], int ix, int iy, int iz, double& Bxl, double& Byl, double& Bzl, double*** Bx, double*** By, double*** Bz, double*** Bx_ext, double*** By_ext, double*** Bz_ext, double Fext){

  Bxl = 0.0;
  Byl = 0.0;
  Bzl = 0.0;

  int l = 0;
  for (int i=0; i<=1; i++)
    for (int j=0; j<=1; j++)
      for (int k=0; k<=1; k++) {
        Bxl += weights[i][j][k] * (Bx[ix-i][iy-j][iz-k] + Fext*Bx_ext[ix-i][iy-j][iz-k]);
        Byl += weights[i][j][k] * (By[ix-i][iy-j][iz-k] + Fext*By_ext[ix-i][iy-j][iz-k]);
        Bzl += weights[i][j][k] * (Bz[ix-i][iy-j][iz-k] + Fext*Bz_ext[ix-i][iy-j][iz-k]);
        l = l + 1;
      }
}

void Particles3Dcomm::get_El(const double weights[2][2][2], int ix, int iy, int iz, double& Exl, double& Eyl, double& Ezl, double*** Ex, double*** Ey, double*** Ez){

  Exl = 0.0;
  Eyl = 0.0;
  Ezl = 0.0;

  int l = 0;
  for (int i=0; i<=1; i++)
    for (int j=0; j<=1; j++)
      for (int k=0; k<=1; k++) {
        Exl += weights[i][j][k] * Ex[ix-i][iy-j][iz-k];
        Eyl += weights[i][j][k] * Ey[ix-i][iy-j][iz-k];
        Ezl += weights[i][j][k] * Ez[ix-i][iy-j][iz-k];
        l = l + 1;
      }

}

void Particles3Dcomm::TrimInfoVector(VirtualTopology3D *vct){
  
  int dim;
  MPI_Comm CommToParent= vct->getCommToParent_P(ns);

  // as a child
  if (CommToParent != MPI_COMM_NULL){
    dim= RG_numPBCMessages+1;

    RGPBC_struct* tmp= new RGPBC_struct[dim];

    for (int i=0; i< RG_numPBCMessages; i++){
      tmp[i]= RGPBC_Info[i];
    }
    delete[]RGPBC_Info;
    RGPBC_Info= new RGPBC_struct[dim];

    for (int i=0; i< RG_numPBCMessages; i++){
      RGPBC_Info[i]= tmp[i];
    }
    delete[]tmp;

  } // end if (CommToParent != MPI_COMM_NULL){
  if (numChildren>0){
    dim=0;
    for (int ch=0; ch< numChildren; ch++){
      if (CG_numPBCMessages[ch] >dim) {dim= CG_numPBCMessages[ch];}
    }
    dim=dim+1;
    RGPBC_struct **tmp2= newArr2(RGPBC_struct, numChildren, dim);
    for (int ch=0; ch< numChildren; ch++){
      for(int d=0; d<dim; d++){
        tmp2[ch][d]= CG_Info[ch][d];
      }
    }
    delArr2(CG_Info, numChildren);

    CG_Info= newArr2(RGPBC_struct, numChildren, dim);
    for (int ch=0; ch< numChildren; ch++){
      for(int d=0; d<dim; d++){
        CG_Info[ch][d]= tmp2[ch][d];
      }
    }

    delArr2(tmp2, numChildren);
  } // end if (numChildren>0){
}

void Particles3Dcomm::MPI_RGPBC_struct_commit(){

  RGBC_struct *a;
  MPI_Datatype type[13]={MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_INT};
  int blocklen[13]={1,1,1,1,1,1,1,1,1,1,1,1,1};
  // displacement in bytes                                                                                                                                                     
  MPI_Aint disp[13];

  disp[0]= (MPI_Aint) &(a->ix_first) - (MPI_Aint)a ;
  disp[1]= (MPI_Aint) &(a->iy_first) - (MPI_Aint)a ;
  disp[2]= (MPI_Aint) &(a->iz_first) - (MPI_Aint)a ;
  // BCside                                                                                                                                                                    
  disp[3]= (MPI_Aint) &(a->BCside) - (MPI_Aint)a ;
  // np_*                                                                                                                                                                      
  disp[4]= (MPI_Aint) &(a->np_x) - (MPI_Aint)a ;
  disp[5]= (MPI_Aint) &(a->np_y) - (MPI_Aint)a ;
  disp[6]= (MPI_Aint) &(a->np_z) - (MPI_Aint)a ;
  // CG_*_first                                                                                                                                                                
  disp[7]= (MPI_Aint) &(a->CG_x_first) - (MPI_Aint)a ;
  disp[8]= (MPI_Aint) &(a->CG_y_first) - (MPI_Aint)a ;
  disp[9]= (MPI_Aint) &(a->CG_z_first) - (MPI_Aint)a ;
  // the cores                                                                                                                                                                 
  disp[10]= (MPI_Aint) &(a->CG_core) - (MPI_Aint)a ;
  disp[11]= (MPI_Aint) &(a->RG_core) - (MPI_Aint)a ;
  // the msg id                                                                                                                                                                
  disp[12]= (MPI_Aint) &(a->MsgID) - (MPI_Aint)a ;

  MPI_Type_create_struct(13, blocklen, disp, type, &MPI_RGBC_struct);
  MPI_Type_commit(&MPI_RGBC_struct);

}

void Particles3Dcomm::initWeightFluidPBC_Phase1(Grid *grid, VirtualTopology3D *vct, RGBC_struct *RGBC_Info, int *numMsg, int *MaxSizeMsg){
  int DNS=1;
  bool DIR_0= true; // L/R  
  bool DIR_1= true; // B/F
  bool DIR_2= true; // top/ bottom  
  int XLEN= vct->getXLEN();
  int YLEN= vct->getYLEN();
  int ZLEN= vct->getZLEN();

  int rank_CTP= vct->getRank_CommToParent();
  int rank_G= vct->getSystemWide_rank();
  int rank_local= vct->getCartesian_rank();

  // NB: _s and _e are included!!!, so <= in the for        
  int i_s, i_e;
  int j_s, j_e;
  int k_s, k_e;

  string FACE;

  char DIR;

  REPO= newArr3(bool, nxc, nyc, nzc);
  for (int i=0; i<nxc; i++)
    for (int j=0; j<nyc; j++)
      for (int k=0; k<nzc; k++)
	REPO[i][j][k]= false;



  int CC=PRACells-1;
  // this is the bottom face                                                                                     
  if (vct->getCoordinates(2)==0 && vct->getZleft_neighbor_P()== MPI_PROC_NULL && DIR_2){
    
    i_s=0; i_e= nxc-1;
    j_s=0; j_e= nyc-1;
    k_s=0; k_e=CC;
        
    DIR= 'B';
    grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
    
    // porcheria per evitare di spezzare in troppi messaggi
    for (int i=i_s; i<= i_e; i++ )
      for (int j=j_s; j<= j_e; j++)
	for (int k= k_s; k<= k_e; k++)
	  REPO[i][j][k]= true;


  } // end bottom face                                                                                           
  
  // this is the top face                                                                                        
  if (vct->getCoordinates(2) ==ZLEN-1 && vct->getZright_neighbor_P() == MPI_PROC_NULL && DIR_2){
    
    i_s=0; i_e= nxc-1;
    j_s=0; j_e= nyc-1;
    k_e= nzc-1; k_s= nzc-1-CC;
    
    DIR= 'T';
    grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
    
    // porcheria per evitare di spezzare in troppi messaggi
    for (int i=i_s; i<= i_e; i++ )
      for (int j=j_s; j<= j_e; j++)
	for (int k= k_s; k<= k_e; k++)
	  REPO[i][j][k]= true;

  } // end top face 

  // this is the left face
  if (vct->getCoordinates(0) ==0  && vct->getXleft_neighbor_P() == MPI_PROC_NULL && DIR_0){
    
    j_s=0; j_e= nyc-1;
    k_s=0; k_e= nzc-1;
    i_s=0; i_e= CC;
    
    DIR= 'L';
    grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
    
    // porcheria per evitare di spezzare in troppi messaggi
    for (int i=i_s; i<= i_e; i++ )
      for (int j=j_s; j<= j_e; j++)
	for (int k= k_s; k<= k_e; k++)
	  REPO[i][j][k]= true;
    
  } // end left face                                                                                             
  
  // this is the right face                                                                                      
  if (vct->getCoordinates(0) ==XLEN-1 && vct->getXright_neighbor_P() == MPI_PROC_NULL && DIR_0){
    
    j_s=0; j_e= nyc-1;
    k_s=0; k_e= nzc-1;
    i_e=nxc-1; i_s= nxc-1-CC;
    
    DIR= 'R';
    grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
    
    // porcheria per evitare di spezzare in troppi messaggi
    for (int i=i_s; i<= i_e; i++ )
      for (int j=j_s; j<= j_e; j++)
	for (int k= k_s; k<= k_e; k++)
	  REPO[i][j][k]= true;
    
  } // end right face                                                                                            
  
  // this is the front face                                                                                      
  if (vct->getCoordinates(1) ==0 && vct->getYleft_neighbor_P() == MPI_PROC_NULL && DIR_1){
    
    i_s=0; i_e= nxc-1;
    k_s=0; k_e= nzc-1;
    j_s=0; j_e= CC;
    
    DIR= 'F';
    grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
    
    // porcheria per evitare di spezzare in troppi messaggi
    for (int i=i_s; i<= i_e; i++ )
      for (int j=j_s; j<= j_e; j++)
	for (int k= k_s; k<= k_e; k++)
	  REPO[i][j][k]= true;

  } // end front face 
  // this is the back face                                                                                       
  if (vct->getCoordinates(1) == YLEN-1 && vct->getYright_neighbor_P() == MPI_PROC_NULL && DIR_1){
    
    i_s=0; i_e= nxc-1;
    k_s=0; k_e= nzc-1;
    j_e=nyc-1; j_s= nyc-1-CC;
    
    DIR= 'b';
    grid->Explore3DAndCommit_Centers(i_s, i_e, j_s, j_e, k_s, k_e, RGBC_Info, numMsg, MaxSizeMsg, vct, DIR);
    
    // porcheria per evitare di spezzare in troppi messaggi
    for (int i=i_s; i<= i_e; i++ )
      for (int j=j_s; j<= j_e; j++)
	for (int k= k_s; k<= k_e; k++)
	  REPO[i][j][k]= true;
    
  } // end back face                                                                                             
  
  RGBC_Info[*numMsg].RG_core= -1;
  RGBC_Info[*numMsg].CG_core= -1;
  
  
}




void Particles3Dcomm::SendFluidPBC(Grid* grid, VirtualTopology3D * vct, Field * EMf){

  if (numChildren==0) return;

  if (TEST_FLUID_BC) return;

  rho= newArr3(double, nxn, nyn, nzn);
  Jx= newArr3(double, nxn, nyn, nzn);
  Jy= newArr3(double, nxn, nyn, nzn);
  Jz= newArr3(double, nxn, nyn, nzn);
  pxx= newArr3(double, nxn, nyn, nzn);
  pxy= newArr3(double, nxn, nyn, nzn);
  pxz= newArr3(double, nxn, nyn, nzn);
  pyy= newArr3(double, nxn, nyn, nzn);
  pyz= newArr3(double, nxn, nyn, nzn);
  pzz= newArr3(double, nxn, nyn, nzn);

  EMf->copyMoments(grid, vct, rho, Jx, Jy, Jz, pxx, pxy, pxz, pyy, pyz, pzz, ns);

  int dest;
  int tag;
  MPI_Request request;
  MPI_Status status;

  for (int ch=0; ch < numChildren; ch++){
    // this child wants PBC and this core participates
    if (CommToChild_P[ch] != MPI_COMM_NULL and CG_numPBCMessages[ch]>0){ 

      for (int m=0; m< CG_numPBCMessages[ch]; m++){

	int Size= CG_Info[ch][m].np_x * CG_Info[ch][m].np_y * CG_Info[ch][m].np_z;
	// build the fluid BC msg
	buildFluidBCMsg(vct, grid, EMf, ch, CG_Info[ch][m], Size, CGFluidMsg );
	
	// send the fluid BC msg
	dest= CG_Info[ch][m].RG_core;
	tag= CG_Info[ch][m].MsgID;

	MPI_Isend(CGFluidMsg, Size*numFBC, MPI_DOUBLE, dest, tag, CommToChild_P[ch], &request);
	MPI_Wait(&request, &status);


      }// end for (int m=0; m< CG_numPBCMessages[ch]; m++){
      
    } // end if (CommToChild_P[ch] != MPI_COMM_NULL and CG_numPBCMessages[ch]>0){
    
  }// end for (int ch=0; ch < numChildren; ch++){

  // deletes
  delArr3(rho, nxn, nyn);
  delArr3(Jx, nxn, nyn);
  delArr3(Jy, nxn, nyn);
  delArr3(Jz, nxn, nyn);
  delArr3(pxx, nxn, nyn);
  delArr3(pxy, nxn, nyn);
  delArr3(pxz, nxn, nyn);
  delArr3(pyy, nxn, nyn);
  delArr3(pyz, nxn, nyn);
  delArr3(pzz, nxn, nyn);
}


void Particles3Dcomm::buildFluidBCMsg(VirtualTopology3D *vct, Grid * grid, Field * EMf, int ch, RGBC_struct RGInfo, int Size, double *Msg ){
  
  // NB: if i send more stuff, do all the interpolation here
  // i use the three indexes, knowing that one will be zero

  // position of the initial RG point in the CG grid
  // points are displacement over this

  //cout << "Building msg starting from " << RGInfo.CG_x_first <<", " << RGInfo.CG_y_first << ", " <<RGInfo.CG_z_first <<" size "<< RGInfo.np_x << ", "  << RGInfo.np_y <<", " << RGInfo.np_z << endl;

  double x0= RGInfo.CG_x_first;
  double y0= RGInfo.CG_y_first;
  double z0= RGInfo.CG_z_first;

  double inv_dx= 1./dx;
  double inv_dy= 1./dy;
  double inv_dz= 1./dz;

  double xp, yp, zp;
  int count =0;

  //cout << "inside building msg RGInfo.np_x " << RGInfo.np_x << " RGInfo.np_y " << RGInfo.np_y << " RGInfo.np_z " << RGInfo.np_z <<endl;

  //cout <<"building msg, starts @ " << x0 << "; " << y0 << "; " << z0 << endl;

  for (int i= 0; i<RGInfo.np_x; i++){
    for (int j=0; j<RGInfo.np_y; j++){
      for (int k=0; k<RGInfo.np_z; k++){
	
	// ok, now each RG point is treated as a particle here
	// and i need the E at the point position

	xp= x0 + i*dx_Ch[ch];
	yp= y0 + j*dy_Ch[ch];
	zp= z0 + k*dz_Ch[ch]; 
	
	//cout << "building msg, point " << xp << "; " <<yp <<"; " << zp << endl;

	// this is copied from the mover
	const double ixd = floor((xp - xstart) * inv_dx);
	const double iyd = floor((yp - ystart) * inv_dy);
	const double izd = floor((zp - zstart) * inv_dz);
	int ix = 2 + int (ixd);
	int iy = 2 + int (iyd);
	int iz = 2 + int (izd);
	if (ix < 1)
	  ix = 1;
	if (iy < 1)
	  iy = 1;
	if (iz < 1)
	  iz = 1;
	if (ix > nxn - 1)
	  ix = nxn - 1;
	if (iy > nyn - 1)
	  iy = nyn - 1;
	if (iz > nzn - 1)
	  iz = nzn - 1;

	double xi  [2];
	double eta [2];
	double zeta[2];

	xi  [0] = xp - grid->getXN(ix-1,iy  ,iz  );
	eta [0] = yp - grid->getYN(ix  ,iy-1,iz  );
	zeta[0] = zp - grid->getZN(ix  ,iy  ,iz-1);
	xi  [1] = grid->getXN(ix,iy,iz) - xp;
	eta [1] = grid->getYN(ix,iy,iz) - yp;
	zeta[1] = grid->getZN(ix,iy,iz) - zp;

	double rhol = 0.0;
	double Jxl = 0.0;
	double Jyl = 0.0;
	double Jzl = 0.0;
	double pxxl = 0.0;
	double pxyl = 0.0;
	double pxzl = 0.0;
	double pyyl = 0.0;
	double pyzl = 0.0;
	double pzzl = 0.0;

	const double weight000 = xi[0] * eta[0] * zeta[0] * invVOL;
	const double weight001 = xi[0] * eta[0] * zeta[1] * invVOL;
	const double weight010 = xi[0] * eta[1] * zeta[0] * invVOL;
	const double weight011 = xi[0] * eta[1] * zeta[1] * invVOL;
	const double weight100 = xi[1] * eta[0] * zeta[0] * invVOL;
	const double weight101 = xi[1] * eta[0] * zeta[1] * invVOL;
	const double weight110 = xi[1] * eta[1] * zeta[0] * invVOL;
	const double weight111 = xi[1] * eta[1] * zeta[1] * invVOL;
	
	rhol += weight000 * (rho[ix][iy][iz]             );
        rhol += weight001 * (rho[ix][iy][iz - 1]         );
        rhol += weight010 * (rho[ix][iy - 1][iz]         );
        rhol += weight011 * (rho[ix][iy - 1][iz - 1]     );
        rhol += weight100 * (rho[ix - 1][iy][iz]         );
        rhol += weight101 * (rho[ix - 1][iy][iz - 1]     );
        rhol += weight110 * (rho[ix - 1][iy - 1][iz]     );
        rhol += weight111 * (rho[ix - 1][iy - 1][iz - 1] );

	Jxl += weight000 * (Jx[ix][iy][iz]             );
        Jxl += weight001 * (Jx[ix][iy][iz - 1]         );
        Jxl += weight010 * (Jx[ix][iy - 1][iz]         );
        Jxl += weight011 * (Jx[ix][iy - 1][iz - 1]     );
        Jxl += weight100 * (Jx[ix - 1][iy][iz]         );
        Jxl += weight101 * (Jx[ix - 1][iy][iz - 1]     );
        Jxl += weight110 * (Jx[ix - 1][iy - 1][iz]     );
        Jxl += weight111 * (Jx[ix - 1][iy - 1][iz - 1] );

	Jyl += weight000 * (Jy[ix][iy][iz]             );
        Jyl += weight001 * (Jy[ix][iy][iz - 1]         );
        Jyl += weight010 * (Jy[ix][iy - 1][iz]         );
        Jyl += weight011 * (Jy[ix][iy - 1][iz - 1]     );
        Jyl += weight100 * (Jy[ix - 1][iy][iz]         );
        Jyl += weight101 * (Jy[ix - 1][iy][iz - 1]     );
        Jyl += weight110 * (Jy[ix - 1][iy - 1][iz]     );
        Jyl += weight111 * (Jy[ix - 1][iy - 1][iz - 1] );

	Jzl += weight000 * (Jz[ix][iy][iz]             );
        Jzl += weight001 * (Jz[ix][iy][iz - 1]         );
        Jzl += weight010 * (Jz[ix][iy - 1][iz]         );
        Jzl += weight011 * (Jz[ix][iy - 1][iz - 1]     );
        Jzl += weight100 * (Jz[ix - 1][iy][iz]         );
        Jzl += weight101 * (Jz[ix - 1][iy][iz - 1]     );
        Jzl += weight110 * (Jz[ix - 1][iy - 1][iz]     );
        Jzl += weight111 * (Jz[ix - 1][iy - 1][iz - 1] );

	pxxl += weight000 * (pxx[ix][iy][iz]             );
        pxxl += weight001 * (pxx[ix][iy][iz - 1]         );
        pxxl += weight010 * (pxx[ix][iy - 1][iz]         );
        pxxl += weight011 * (pxx[ix][iy - 1][iz - 1]     );
        pxxl += weight100 * (pxx[ix - 1][iy][iz]         );
        pxxl += weight101 * (pxx[ix - 1][iy][iz - 1]     );
        pxxl += weight110 * (pxx[ix - 1][iy - 1][iz]     );
        pxxl += weight111 * (pxx[ix - 1][iy - 1][iz - 1] );

	pxyl += weight000 * (pxy[ix][iy][iz]             );
        pxyl += weight001 * (pxy[ix][iy][iz - 1]         );
        pxyl += weight010 * (pxy[ix][iy - 1][iz]         );
        pxyl += weight011 * (pxy[ix][iy - 1][iz - 1]     );
        pxyl += weight100 * (pxy[ix - 1][iy][iz]         );
        pxyl += weight101 * (pxy[ix - 1][iy][iz - 1]     );
        pxyl += weight110 * (pxy[ix - 1][iy - 1][iz]     );
        pxyl += weight111 * (pxy[ix - 1][iy - 1][iz - 1] );

	pxzl += weight000 * (pxz[ix][iy][iz]             );
        pxzl += weight001 * (pxz[ix][iy][iz - 1]         );
        pxzl += weight010 * (pxz[ix][iy - 1][iz]         );
        pxzl += weight011 * (pxz[ix][iy - 1][iz - 1]     );
        pxzl += weight100 * (pxz[ix - 1][iy][iz]         );
        pxzl += weight101 * (pxz[ix - 1][iy][iz - 1]     );
        pxzl += weight110 * (pxz[ix - 1][iy - 1][iz]     );
        pxzl += weight111 * (pxz[ix - 1][iy - 1][iz - 1] );

	pyyl += weight000 * (pyy[ix][iy][iz]             );
        pyyl += weight001 * (pyy[ix][iy][iz - 1]         );
        pyyl += weight010 * (pyy[ix][iy - 1][iz]         );
        pyyl += weight011 * (pyy[ix][iy - 1][iz - 1]     );
        pyyl += weight100 * (pyy[ix - 1][iy][iz]         );
        pyyl += weight101 * (pyy[ix - 1][iy][iz - 1]     );
        pyyl += weight110 * (pyy[ix - 1][iy - 1][iz]     );
        pyyl += weight111 * (pyy[ix - 1][iy - 1][iz - 1] );

	pyzl += weight000 * (pyz[ix][iy][iz]             );
        pyzl += weight001 * (pyz[ix][iy][iz - 1]         );
        pyzl += weight010 * (pyz[ix][iy - 1][iz]         );
        pyzl += weight011 * (pyz[ix][iy - 1][iz - 1]     );
        pyzl += weight100 * (pyz[ix - 1][iy][iz]         );
        pyzl += weight101 * (pyz[ix - 1][iy][iz - 1]     );
        pyzl += weight110 * (pyz[ix - 1][iy - 1][iz]     );
        pyzl += weight111 * (pyz[ix - 1][iy - 1][iz - 1] );

	pzzl += weight000 * (pzz[ix][iy][iz]             );
        pzzl += weight001 * (pzz[ix][iy][iz - 1]         );
        pzzl += weight010 * (pzz[ix][iy - 1][iz]         );
        pzzl += weight011 * (pzz[ix][iy - 1][iz - 1]     );
        pzzl += weight100 * (pzz[ix - 1][iy][iz]         );
        pzzl += weight101 * (pzz[ix - 1][iy][iz - 1]     );
        pzzl += weight110 * (pzz[ix - 1][iy - 1][iz]     );
        pzzl += weight111 * (pzz[ix - 1][iy - 1][iz - 1] );


	Msg[0*Size +count]= rhol;
	
	Msg[1*Size +count]= Jxl;
	Msg[2*Size +count]= Jyl;
	Msg[3*Size +count]= Jzl;

	Msg[4*Size +count]= pxxl;
	Msg[5*Size +count]= pxyl;
	Msg[6*Size +count]= pxzl;
	Msg[7*Size +count]= pyyl;
	Msg[8*Size +count]= pyzl;
	Msg[9*Size +count]= pzzl;

	count ++;

      }
    }
  }
  //cout <<"R" << vct->getSystemWide_rank() << ", count " << count << "for child " << ch << endl;
}

/** refined grid receives fluid BCs from the coarse grid **/
void Particles3Dcomm::ReceiveFluidBC(Grid *grid, VirtualTopology3D *vct){

  PRA_PAdded=0;
  nop_BeforeReceivePBC=nop;

  MPI_Comm CommToParent= vct->getCommToParent_P(ns);

  if (CommToParent == MPI_COMM_NULL) {return;}  // if you are not a refined grid, no need to be here

  MPI_Status status;
  bool found;
  int countExp;
  int count;

  if (TEST_FLUID_BC) return;
  
  bool Testing= true; // put to false during production
  
  
  for (int m=0; m< RG_numPBCMessages; m++){
    
    MPI_Recv(RGFluidMsg, RG_MaxFluidMsgSize*numFBC, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, CommToParent, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &count );
    
    found= false;
    // check with stored msgs
    for (int i=0; i< RG_numPBCMessages; i++){
      
      if (RGPBC_Info[i].MsgID == status.MPI_TAG and RGPBC_Info[i].CG_core == status.MPI_SOURCE){ // NB: the tag gives you the msg position in the BC vector
	
	found= true;
	
	countExp= RGPBC_Info[i].np_x * RGPBC_Info[i].np_y * RGPBC_Info[i].np_z  ;
	
	if (Testing){
	  if (countExp *numFBC  != count){
	    cout << "R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent_P(ns) <<" : msg recv from core " << status.MPI_SOURCE <<" with tag " << status.MPI_TAG << " but fluid PBC msg size does not check: received: " << count << " expected " << countExp*numFBC << ", aborting ..." << endl;
	    MPI_Abort(MPI_COMM_WORLD, -1);
	  }
	}
	
	// put BC here
	// the tag gives you the position in the BC vector
	
	// rho
	memcpy(rho_FBC[status.MPI_TAG], RGFluidMsg, sizeof(double)*countExp );
	// J
	memcpy(Jx_FBC[status.MPI_TAG], RGFluidMsg +   countExp, sizeof(double)*countExp );
	memcpy(Jy_FBC[status.MPI_TAG], RGFluidMsg + 2*countExp, sizeof(double)*countExp );
	memcpy(Jz_FBC[status.MPI_TAG], RGFluidMsg + 3*countExp, sizeof(double)*countExp );
	// P
	memcpy(Pxx_FBC[status.MPI_TAG], RGFluidMsg + 4*countExp, sizeof(double)*countExp );
	memcpy(Pxy_FBC[status.MPI_TAG], RGFluidMsg + 5*countExp, sizeof(double)*countExp );
	memcpy(Pxz_FBC[status.MPI_TAG], RGFluidMsg + 6*countExp, sizeof(double)*countExp );
	memcpy(Pyy_FBC[status.MPI_TAG], RGFluidMsg + 7*countExp, sizeof(double)*countExp );
	memcpy(Pyz_FBC[status.MPI_TAG], RGFluidMsg + 8*countExp, sizeof(double)*countExp );
	memcpy(Pzz_FBC[status.MPI_TAG], RGFluidMsg + 9*countExp, sizeof(double)*countExp );
	break;
	
      } // end if (RGBC_Info_Active[i].MsgID == status.MPI_TAG )
      
    } //for (int i=0; i< RG_numBCMessages_Active; i++){ // here end search of the msg
    
    if (Testing){
      if (found == false){
	cout <<"R" << vct->getSystemWide_rank() << " numGrid " << numGrid <<" PC rank " << vct->getRank_CommToParent() << " I have received an active msg I cannot match with my record, aborting..." << endl;
	MPI_Abort(MPI_COMM_WORLD, -1);
      }
    }
    
  } // here is when i receive it 
  
  
}

void Particles3Dcomm::ApplyFluidPBC(Grid *grid, VirtualTopology3D *vct, Field * EMf){
  
  /* if true, use initial profile for repopulation -
     - do not use for production */
  


  MPI_Comm CommToParent= vct->getCommToParent_P(ns);
  
  if (CommToParent == MPI_COMM_NULL or FluidLikeRep == false or RG_numPBCMessages <1 )
    return;

  if (!TEST_FLUID_BC) {
    
    for (int m=0; m< RG_numPBCMessages; m++ )  {
      
      int II= RGPBC_Info[m].np_x;
      int JJ= RGPBC_Info[m].np_y;
      int KK= RGPBC_Info[m].np_z;
      
      int ix_f= RGPBC_Info[m].ix_first;
      int iy_f= RGPBC_Info[m].iy_first;
      int iz_f= RGPBC_Info[m].iz_first;
      
      int countMsg=0;
      
      for (int i= 0; i<II; i++)
	for (int j=0; j<JJ; j++)
	  for (int k=0; k<KK; k++){
	    
	    RHOP[ix_f + i][iy_f + j][iz_f + k]= rho_FBC[m][countMsg];
	    
	    UP[ix_f + i][iy_f + j][iz_f + k]= Jx_FBC[m][countMsg]/ rho_FBC[m][countMsg];
	    VP[ix_f + i][iy_f + j][iz_f + k]= Jy_FBC[m][countMsg]/ rho_FBC[m][countMsg];
	    WP[ix_f + i][iy_f + j][iz_f + k]= Jz_FBC[m][countMsg]/ rho_FBC[m][countMsg];
	    
	    double PTH_X= (Pxx_FBC[m][countMsg]- Jx_FBC[m][countMsg]*Jx_FBC[m][countMsg]/ rho_FBC[m][countMsg]  )/qom;
	    double PTH_Y= (Pyy_FBC[m][countMsg]- Jy_FBC[m][countMsg]*Jy_FBC[m][countMsg]/ rho_FBC[m][countMsg]  )/qom;
	    double PTH_Z= (Pzz_FBC[m][countMsg]- Jz_FBC[m][countMsg]*Jz_FBC[m][countMsg]/ rho_FBC[m][countMsg]  )/qom;
	    
	    UTHP[ix_f + i][iy_f + j][iz_f + k]= sqrt(abs(PTH_X/ rho_FBC[m][countMsg] *qom));
	    VTHP[ix_f + i][iy_f + j][iz_f + k]= sqrt(abs(PTH_Y/ rho_FBC[m][countMsg] *qom));
	    WTHP[ix_f + i][iy_f + j][iz_f + k]= sqrt(abs(PTH_Z/ rho_FBC[m][countMsg] *qom));
	    
	    countMsg++;
	  } // end cycle on i, j, k
    } // end  for (int m=0; m< RG_numPBCMessages; m++ ) 
  } else { // if TEST_FLUID_BC
    // if (TEST_FLUID_BC  and RG_numPBCMessages>0 ){
    for (int i=0; i< nxc; i++)
      for (int j=0; j< nyc; j++)
	for (int k=0; k< nzc; k++){
	  
	  UTHP[i][j][k]= uth;
	  VTHP[i][j][k]= vth;
	  WTHP[i][j][k]= wth;
	  
	  UP[i][j][k]= u0;
	  VP[i][j][k]= v0;
	  WP[i][j][k]= w0;
	  
	  RHOP[i][j][k]= EMf->getRHOINIT(ns, i, j, k);
	  
	}
  } // end if (TEST_FLUID_BC){
  
  // now i have the data to regenerate
  
  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);
  
  double harvest;
  double prob, theta, sign;
  long long counter = nop;

  int i_s=1, i_e= grid->getNXC() - 1;
  int j_s=1, j_e= grid->getNYC() - 1;
  int k_s=1, k_e= grid->getNZC() - 1;

  if ( vct->getXleft_neighbor_P() == MPI_PROC_NULL) i_s=0;
  if ( vct->getXright_neighbor_P() == MPI_PROC_NULL) i_e=grid->getNXC();

  if ( vct->getYleft_neighbor_P() == MPI_PROC_NULL) j_s=0;
  if ( vct->getYright_neighbor_P() == MPI_PROC_NULL) j_e=grid->getNYC();

  if ( vct->getZleft_neighbor_P() == MPI_PROC_NULL) k_s=0;
  if ( vct->getZright_neighbor_P() == MPI_PROC_NULL) k_e=grid->getNZC();

  for (int i = i_s; i < i_e; i++)
    for (int j = j_s; j < j_e; j++)
      for (int k = k_s; k < k_e; k++){

	if (REPO[i][j][k]== false) continue; // repopulate only the flagged cells

        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {

	      if (RHOP[i][j][k]== -100.0){
		cout << "I have messed up in fluid repopulation aborting " << endl;
		MPI_Abort(MPI_COMM_WORLD, -1);
		return ;     
	      }
	      
	      // init with random position in the cell
	      double rX = ((double)rand() / (double)(RAND_MAX));
	      double rY = ((double)rand() / (double)(RAND_MAX));
	      double rZ = ((double)rand() / (double)(RAND_MAX));

	      x[counter]= grid->getXN(i, j, k)+ dx*rX;
	      y[counter]= grid->getYN(i, j, k)+ dy*rY;
	      z[counter]= grid->getZN(i, j, k)+ dz*rZ;
	            
              // q = charge
              q[counter] = (qom / fabs(qom)) * (fabs(RHOP[i][j][k]) / npcel) * (1.0 / grid->getInvVOL());
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              u[counter] = UP[i][j][k] + UTHP[i][j][k] * prob * cos(theta);
              // v
              v[counter] = VP[i][j][k] + VTHP[i][j][k] * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              w[counter] = WP[i][j][k] + WTHP[i][j][k] * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];


              counter++;
	      PRA_PAdded++;
            }
      }
  // to set number of particles in allocate, keeping in mind that in mlmd there can be particles in the GC
  nop= counter;

}

void Particles3Dcomm::UpdateAllowPMsgResize(Collective * col, int cycle){

  if (FluidLikeRep) { AllowPMsgResize= false; return;}
 
  // so buffers are sized to a decent value
  if (cycle <3) {AllowPMsgResize= true; return;}

  if (cycle >=3){
    if (col->getAllowPMsgResize()) AllowPMsgResize= true;
    else
      AllowPMsgResize= false;
  }
  return;

}

int Particles3Dcomm::communicate_DepopulatePRA(VirtualTopology3D * ptVCT) {
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

  bool NoBCPart= false;
  if ( CommToParent_P!= MPI_COMM_NULL and bcPfaceXleft <0) NoBCPart=true;

  while (np_current < nplast+1){
     
    xMin=0; yMin=0; zMin=0;
    xMax=Lx; yMax=Ly; zMax=Lz;

    // i do not need to apply BC if repopulating
    if (NoBCPart== false ){

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
    }

    if (ptVCT->getPERIODICX_P()) { xMin=-Lx; xMax=2*Lx; } // so periodic particles will stay in the system and be communicated    
    else if ((bcPfaceXleft == -1 or bcPfaceXleft == -3 or bcPfaceXleft == -4) and CommToParent_P!= MPI_COMM_NULL) {
      xMin= Coord_XLeft_End; xMax= Coord_XRight_Start;
    }
    else if (bcPfaceXleft == -2 and CommToParent_P!= MPI_COMM_NULL) {
      xMin= Coord_XLeft_Start; xMax= Coord_XRight_End;
    }
    else {xMin= 0; xMax= Lx;} // non periodic, non mlmd              

    if (ptVCT->getPERIODICY_P()) { yMin=-Ly; yMax=2*Ly; } // so periodic particles will stay in the system and be communicated    
    else if ((bcPfaceYleft == -1 or bcPfaceYleft == -3 or bcPfaceYleft == -4) and CommToParent_P!= MPI_COMM_NULL) {
      yMin= Coord_YLeft_End; yMax= Coord_YRight_Start;
    }
    else if (bcPfaceYleft == -2 and CommToParent_P!= MPI_COMM_NULL) {
      yMin= Coord_YLeft_Start; yMax= Coord_YRight_End;
    }
    else {yMin= 0; yMax= Ly;} // non periodic, non mlmd     

    if (ptVCT->getPERIODICZ_P() ) { zMin=-Lz; zMax=2*Lz; } // so periodic particles will stay in the system and be communicated    
    else if ((bcPfaceZleft == -1 or bcPfaceZleft == -3 or bcPfaceZleft == -4) and CommToParent_P!= MPI_COMM_NULL) {
      zMin= Coord_ZLeft_End; zMax= Coord_ZRight_Start;
    }
    else if (bcPfaceZleft == -2 and CommToParent_P!= MPI_COMM_NULL) {
      zMin= Coord_ZLeft_Start; zMax= Coord_ZRight_End;
    }
    else {zMin= 0; zMax= Lz;} // non periodic, non mlmd 
    

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
  /** do nor touch this otherwise mess in communicateRepopulatedParticles **/
  nop_EndCommunicate= nop;

    // i am removing sync points
#ifdef __PROFILING__
  int ppc=nop; int TotalP=0;
  MPI_Allreduce(&ppc, &TotalP, 1, MPI_INT, MPI_SUM, ptVCT->getCommGrid());
  if (ptVCT->getCartesian_rank()==0){
    cout << "Grid " << numGrid << " ns " <<ns  <<": total number of particles AFTER communicateAfterMover (when particles in PRA removed): " << TotalP << endl;

    ofstream my_file(parInfo.c_str(), fstream::app);
    my_file << TotalP <<" " ;
    my_file.close();
    }

#endif

  return(0);

}
