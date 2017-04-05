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
  delete[]RGPBC_Info;
  delArr3(PCGMsg, numChildren, MaxNumMsg);
  delArr2(nopPCMsg, numChildren);
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

  // This if is necessary to restart with H5hut-io
  if (initnpmax==0){
    /*! pre-mlmd
    long ncproc = int(col->getNxc()/col->getXLEN()) *
                  int(col->getNyc()/col->getYLEN()) *
                  int(col->getNzc()/col->getZLEN()); */
    
    long ncproc = int(col->getNxc_mlmd(numGrid)/vct->getXLEN()) *
      int(col->getNyc_mlmd(numGrid)/vct->getYLEN()) *
      int(col->getNzc_mlmd(numGrid)/vct->getZLEN());

    nop   = ncproc * npcel;
    npmax = nop * col->getNpMaxNpRatio();
  }
  else {
    npmax = initnpmax*col->getNpMaxNpRatio();
    nop   = initnpmax;
  }

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
  
  int MaxGridCoreN= vct->getMaxGridCoreN();
  MAX_RG_numPBCMessages= (int) (MaxGridCoreN*6+1);
  MAX_RG_numPBCMessages_LevelWide= MAX_RG_numPBCMessages*4;

  numChildren= vct->getNumChildren();
  // here, to be able to used RGPBC_struct as an MPI_Datatype
  MPI_RGPBC_struct_commit();
  // here, to be able to use the RepP_struct as an MPI_Datatype
  MPI_RepP_struct_commit();


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
      EMf->addRho(weight, ix, iy, iz, ns);
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

  EMf->communicateGhostP2G(ns, 0, 0, 0, 0, vct);
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

  while (np_current < nplast+1){


    // BC on particles
    if (x[np_current] < 0 && ptVCT->getXleft_neighbor_P() == MPI_PROC_NULL)
      BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft);
    else if (x[np_current] > Lx && ptVCT->getXright_neighbor_P() == MPI_PROC_NULL)
      BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft); 
    if (y[np_current] < 0 && ptVCT->getYleft_neighbor_P() == MPI_PROC_NULL)  // check it here
      BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft);
    else if (y[np_current] > Ly && ptVCT->getYright_neighbor_P() == MPI_PROC_NULL) //check it here
      BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft); 
    if (z[np_current] < 0 && ptVCT->getZleft_neighbor_P() == MPI_PROC_NULL)  // check it here
      BCpart(&z[np_current],&w[np_current],&u[np_current],&v[np_current],Lz,wth,uth,vth,bcPfaceZright,bcPfaceZleft);
    else if (z[np_current] > Lz && ptVCT->getZright_neighbor_P() == MPI_PROC_NULL) //check it here
      BCpart(&z[np_current],&w[np_current],&u[np_current],&v[np_current],Lz,wth,uth,vth,bcPfaceZright,bcPfaceZleft);

    // if the particle exits, apply the boundary conditions add the particle to communication buffer
    if (x[np_current] < xstart || x[np_current] >xend){
      // communicate if they don't belong to the domain
      if (x[np_current] < xstart && ptVCT->getXleft_neighbor_P() != MPI_PROC_NULL){
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
      else if (x[np_current] < xstart && ptVCT->getXleft_neighbor_P() == MPI_PROC_NULL){
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
      else if (x[np_current] > xend && ptVCT->getXright_neighbor_P() == MPI_PROC_NULL){
        del_pack(np_current,&nplast);
        npExitXright++;
      }

    } else  if (y[np_current] < ystart || y[np_current] >yend){
      // communicate if they don't belong to the domain
      if (y[np_current] < ystart && ptVCT->getYleft_neighbor_P() != MPI_PROC_NULL){
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
      else if (y[np_current] < ystart && ptVCT->getYleft_neighbor_P() == MPI_PROC_NULL){
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
      else if (y[np_current] > yend && ptVCT->getYright_neighbor_P() == MPI_PROC_NULL){
        // delete the particle and pack the particle array, the value of nplast changes
        del_pack(np_current,&nplast);
        npExitYright++;
      }
    } else  if (z[np_current] < zstart || z[np_current] >zend){
      // communicate if they don't belong to the domain
      if (z[np_current] < zstart && ptVCT->getZleft_neighbor_P() != MPI_PROC_NULL){
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
      else if (z[np_current] < zstart && ptVCT->getZleft_neighbor_P() == MPI_PROC_NULL){
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
      }
      else if (z[np_current] > zend && ptVCT->getZright_neighbor_P() == MPI_PROC_NULL){
        del_pack(np_current,&nplast);
        npExitZright++;
      }
    }  else {
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
/** put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
void Particles3Dcomm::bufferXleft(double *b_, long long np_current, VirtualTopology3D * vct) {
  if (x[np_current] < 0)
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
  if (x[np_current] > Lx)
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
  if (y[np_current] < 0)
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
  if (y[np_current] > Ly)
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
  if (z[np_current] < 0)
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
  if (z[np_current] > Lz)
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
    if (x[nop] < xstart || x[nop] > xend || y[nop] < ystart || y[nop] > yend || z[nop] < zstart || z[nop] > zend)
      rightDomain++;            // the particle is not in the domain
    nop++;
    if (nop > npmax) {
      cout << "Number of particles in the domain " << nop << " and maxpart = " << npmax << endl;
      MPI_Abort(Comm, -1);
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

  /* sizes of the PCGMsg values set here */
  sizePCGMsg = nop;
  MAXsizePCMsg = nop;
  /* end sizes of the PCGMsg values set here */
  
  RG_numPBCMessages= 0;
  int rank_local= vct->getCartesian_rank();  // on the grid communicator
  int HighestRank= vct->getXLEN()*vct->getYLEN()*vct->getZLEN()-1;
  MPI_Comm CommToParent_P= vct->getCommToParent_P();
  MPI_Status status;
  int TAG_CG_RG= ns; // i need to put it here to be visible from both CG and RG

  /* phase 1: as a child */
  if (CommToParent_P != MPI_COMM_NULL) { // meaning: you are a child AND you want to receive PBC
    
    RGPBC_Info = new RGPBC_struct[MAX_RG_numPBCMessages];
    
    // this number is at the moment arbitrarily decided
    PRA_XLeft = 2;
    PRA_XRight = 2;
    PRA_YLeft = 2;
    PRA_YRight = 2;
    PRA_ZLeft = 2;
    PRA_ZRight = 2;

    initWeightPBC_Phase1(grid, vct, RGPBC_Info, &RG_numPBCMessages);
    
    

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
	abort();
      }
      
    }

    // here i should decide wether i want to instantiate vectors to receive particles, or create them on the spot; also, vectors to store particle while processing them

  } // end if (CommToParent_P != MPI_COMM_NULL) { 


  if (CommToParent_P != MPI_COMM_NULL) { // meaning: you are a child AND you want to receive PBC
    int TAG_1b1c= ns;
    // phase 1b: all RG cores (but itself) send their msgs to highest ranking core 
    if (rank_local < HighestRank){
      /* send one message more; the last message has -1 in the RG_core   
	 to signal end of 'valid' messages */
      int dest= vct->getXLEN()*vct->getYLEN()*vct->getZLEN()-1;
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
      
      // i am receiveing XLEN*YLEN*ZLEN-1 --> HighestRank msg
      for (int c=0; c< HighestRank; c++){
	int recv=0;
	MPI_Recv(buffer_rcv, MAX_RG_numPBCMessages, MPI_RGPBC_struct, MPI_ANY_SOURCE, TAG_1b1c, vct->getComm(), &status);
	
	while(buffer_rcv[recv].RG_core != -1){
	  RGPBC_Info_LevelWide[RG_numPBCMessages_LevelWide]= buffer_rcv[recv];
	  RG_numPBCMessages_LevelWide++;

	  recv++;

	  if (RG_numPBCMessages_LevelWide== MAX_RG_numPBCMessages_LevelWide){
	    cout << "initWeightPBC: MAX_RG_numPBCMessages_LevelWide is " << MAX_RG_numPBCMessages_LevelWide <<", but you need more..."<< endl << "Aborting...";
	    abort();
	  }
	  
	} // end while
	
      } // end for (int c=0; c< HighestRank; c++){

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
	  abort();
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

      cout << "ParentCoreNum= "<< ParentCoreNum<< endl;
      for (int cg=0; cg< ParentCoreNum; cg++){
	MPI_Send(&(RGPBC_Info_ToCGCore[cg][0]), RG_numPBCMessages_ToCGCore[cg]+1, MPI_RGPBC_struct, cg, TAG_CG_RG, CommToParent_P);
	cout << "R " << rankAsChild << " on thr PC communicator has just sent a msg to core " << cg << endl; 
      }

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
    if (vct->getCommToChild_P(ch) != MPI_COMM_NULL){ // it may happen that a particular child does not want PC
      
      int ChildHighestRank;
      MPI_Comm_size(vct->getCommToChild_P(ch), &ChildHighestRank); ChildHighestRank--;
      int PCrank= vct->getRank_CommToChildren_P(ch);
      
      //cout << "Grid " << numGrid << ", local rank on PC comm " << PCrank << " is trying to receive from ch " << ch << endl;

      // receive one msg from HighestRank on the PC comm

      // to put ChildHighestRank as source is just a precaution, MPI_ANY_Source should work
      MPI_Recv(CG_buffer, MAX_RG_numPBCMessages, MPI_RGPBC_struct, ChildHighestRank, TAG_CG_RG, vct->getCommToChild_P(ch), &status);

      //cout << "Grid " << numGrid << ", local rank on PC comm " << PCrank << " has recevied from ch " << ch << " tag " << TAG_CG_RG << endl;

      
      // process update the msg structure
      while (CG_buffer[CG_numPBCMessages[ch]].RG_core!= -1){
	//cout << "INSIDE WHILE: ch is " << ch << " CG_numPBCMessages[ch] is " << CG_numPBCMessages[ch] << " CG_buffer[CG_numPBCMessages[ch]].RG_core is " <<CG_buffer[CG_numPBCMessages[ch]].RG_core << endl;
	CG_Info[ch][CG_numPBCMessages[ch]]= CG_buffer[CG_numPBCMessages[ch]];
	CG_numPBCMessages[ch]++;
      }

    } // end if (vct->getCommToChild_P()!= MPI_COMM_NULL)
    
    // allocate CG vectors

    delete[]CG_buffer;
    
  } // end for (int ch=0; ch< numChildren; ch++)

  if (numChildren>0){
    MaxNumMsg=0;
    for (int i=0; i< numChildren; i++){
      if (CG_numPBCMessages[i]> MaxNumMsg) MaxNumMsg= CG_numPBCMessages[i];
    }
    
    PCGMsg= newArr3(RepP_struct, numChildren, MaxNumMsg, sizePCGMsg);
    nopPCMsg = newArr2(int, numChildren, MaxNumMsg);
    cout << "Grid " << numGrid << " core " << vct->getCartesian_rank() << " nopPCMsg has sizes " << numChildren << " x " <<MaxNumMsg << endl;  
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

/* Phase 1: RG cores build their side of the map for PBC */
void Particles3Dcomm::initWeightPBC_Phase1(Grid *grid, VirtualTopology3D * vct, RGPBC_struct *RGPBC_Info, int *RG_numPBCMessages){

  // NB: when build the PBC msg on the CG side, I have to pay attention not to replicate the info; that was not a problem with fields, but it is now
  int XLEN= vct->getXLEN();
  int YLEN= vct->getYLEN();
  int ZLEN= vct->getZLEN();

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
    /*i_s= 0; i_e= nxn-1;
    j_s= 0; j_e= nyn-1;
    k_s= 0; k_e= PRA_ZLeft;*/
    // 
    i_s= 0-1; i_e= nxn-1+1;
    j_s= 0-1; j_e= nyn-1+1;
    k_s= 0-1; k_e= PRA_ZLeft+1;

    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    cout <<"FACE " << FACE << endl;
  } // end bottom face
  
  // this is the top face
  if (vct->getCoordinates(2) ==ZLEN-1 && vct->getZright_neighbor_P() == MPI_PROC_NULL && DIR_2){
    
    FACE= "TOP";
    /*i_s=0; i_e= nxn-1;
    j_s=0; j_e= nyn-1;
    k_s= nzn-1-PRA_ZRight; k_e= nzn-1;*/
    i_s=0-1; i_e= nxn-1+1;
    j_s=0-1; j_e= nyn-1+1;
    k_s= nzn-1-PRA_ZRight-1; k_e= nzn-1+1;

    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    cout<<"FACE " << FACE << endl;
  } // end top face

  // this is the left face
  if (vct->getCoordinates(0) ==0  && vct->getXleft_neighbor_P() == MPI_PROC_NULL && DIR_0){

    FACE= "LEFT";
    /*i_s= 0; i_e= PRA_XLeft;
    j_s= 0; j_e= nyn-1;
    k_s= 0; k_e= nzn-1;*/
    i_s= 0-1; i_e= PRA_XLeft+1;
    j_s= 0-1; j_e= nyn-1+1;
    k_s= 0-1; k_e= nzn-1+1;

    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    cout<<"FACE " << FACE << endl;
  } // end left face

  // this is the right face
  if (vct->getCoordinates(0) ==XLEN-1 && vct->getXright_neighbor_P() == MPI_PROC_NULL && DIR_0){

    FACE= "RIGHT";
    /*i_s= nxn-1-PRA_XRight; i_e= nxn-1;
    j_s= 0; j_e= nyn-1;
    k_s= 0; k_e= nzn-1;*/
    i_s= nxn-1-PRA_XRight-1; i_e= nxn-1+1;
    j_s= 0-1; j_e= nyn-1+1;
    k_s= 0-1; k_e= nzn-1+1;

    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    cout<<"FACE " << FACE << endl;
  } // end right face
  
  // this is the front face
  if (vct->getCoordinates(1) ==0 && vct->getYleft_neighbor_P() == MPI_PROC_NULL && DIR_1){ 

    FACE= "FRONT";
    /*i_s= 0; i_e= nxn-1;
    j_s= 0; j_e= PRA_YLeft;
    k_s= 0; k_e= nzn-1;*/
    i_s= 0-1; i_e= nxn-1+1;
    j_s= 0-1; j_e= PRA_YLeft+1;
    k_s= 0-1; k_e= nzn-1+1;

    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    cout<<"FACE " << FACE << endl;
  } // end front face

  // this is the back face
  if (vct->getCoordinates(1) == YLEN-1 && vct->getYright_neighbor_P() == MPI_PROC_NULL && DIR_1){

    FACE= "BACK";
    /*i_s= 0; i_e= nxn-1;
    j_s= nyn-1-PRA_YRight; j_e= nyn-1;
    k_s= 0; k_e= nzn-1;*/
    i_s= 0-1; i_e= nxn-1+1;
    j_s= nyn-1-PRA_YRight-1; j_e= nyn-1+1;
    k_s= 0-1; k_e= nzn-1+1;
    
    Explore3DAndCommit(grid, i_s, i_e, j_s, j_e, k_s, k_e, RG_numPBCMessages, &MAX_RG_numPBCMessages, vct);
    cout<<"FACE " << FACE << endl;
  } // end back face
  
  /* for further use, i need to set the RG_core field of the first unused slot to -1  
     but DO NOT MODIFY THE NUMBER OF MSGs;                             
     I will just send a +1 */

  cout << "R" <<SW_rank <<" RG_numPBCMessages= " <<*RG_numPBCMessages <<endl;
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

  /*  cout << "Inside explore and commit: " << endl;
  cout << "i_s= " <<i_s <<" -> " << CG_C_i_s << endl;
  cout << "i_e= " <<i_e <<" -> " << CG_C_i_e <<endl;

  cout << "j_s= " <<j_s <<" -> " << CG_C_j_s << endl;
  cout << "j_e= " <<j_e <<" -> " << CG_C_j_e <<endl;

  cout << "k_s= " <<k_s <<" -> " << CG_C_k_s << endl;
  cout << "k_e= " <<k_e <<" -> " << CG_C_k_e <<endl;*/

  // Z dir / Dir 1
  grid->RGBCExploreDirection(vct, FACE, 2, k_s, k_e, i_s, j_s, Dir1_SPXperC, Dir1_SPYperC, Dir1_SPZperC, Dir1_NPperC, Dir1_rank, &Dir1_Ncores, Dir1_IndexFirstPointperC);

  /*
  // debug
  cout << "G"<< numGrid <<"R" << rank_CTP <<": Z dir, " << Dir1_Ncores << " cores " << endl;
  for (int ii=0; ii< Dir1_Ncores; ii++){
    cout << "G"<< numGrid <<"R" << rank_CTP << " " << ii <<": " << Dir1_rank[ii] << " first coord: " <<Dir1_SPZperC[ii] << ", n points: " << Dir1_NPperC[ii] << " which means up to " << Dir1_SPZperC[ii]+(Dir1_NPperC[ii] -1)*dz << endl;
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

  for (int ch=0; ch < numChildren; ch++){
    if (vct->getCommToChild_P(ch) != MPI_COMM_NULL){ // this child wants PBC
      buildPBCMsg(grid, vct, ch);
    }
  }
  // now send; send an extra msg for the -1
  
  //...
  return;
}

void Particles3Dcomm::buildPBCMsg(Grid* grid, VirtualTopology3D * vct, int ch){
  
  /* returns the number - not the level or the order in the children vector - of the child grid n */
  int childNum= vct->getChildGridNum(ch);
  // these are the dx, dy, dz of the child
  double dx= grid->getDx_mlmd(childNum);
  double dy= grid->getDy_mlmd(childNum);
  double dz= grid->getDz_mlmd(childNum);

  double x_min, x_max, y_min, y_max, z_min, z_max;
 
  // set number of particles to 0
  for (int m=0; m< CG_numPBCMessages[ch]; m++){
    nopPCMsg[ch][m]= 0;
  }
  
  if (CG_numPBCMessages[ch]>0){ // this particular core has to send BC
    for (int p=0; p<nop; p++){
      // here, i have to check all the msgs in else if, to make sure that I am sending 
      // this particular particle only to one core
      
      for (int n=0; n< CG_numPBCMessages[ch]; n++){ 
	  
	x_min= CG_Info[ch][n].CG_x_first;
	x_max= CG_Info[ch][n].CG_x_first+ dx*(CG_Info[ch][n].np_x-1);
	y_min= CG_Info[ch][n].CG_y_first;
	y_max= CG_Info[ch][n].CG_y_first+ dy*(CG_Info[ch][n].np_y-1);
	z_min= CG_Info[ch][n].CG_z_first;
	z_max= CG_Info[ch][n].CG_z_first+ dz*(CG_Info[ch][n].np_z-1);
	
	/*cout << "QUA" << endl;
	cout <<"p " << p  << " x_min " << x_min << " x_max " << x_max << " x[p] " << x[p] << " xstart " << xstart << " xend "<< xend << endl;
	cout <<"p " << p  << " y_min " << y_min << " y_max " << y_max << " y[p] " << y[p] << " ystart " << ystart << " yend "<< yend << endl;
	cout <<"p " << p  << " z_min " << z_min << " z_max " << z_max << " z[p] " << z[p] << " zstart " << zstart << " zend "<< zend << endl;
	cout << "End QUA" << endl;*/

	if (x_min <= x[p] && x_max >= x[p] && y_min <= y[p] && y_max >= y[p] && z_min <= z[p] && z_max >= z[p]){ // to use for PBC
	  cout << "adding " << p << " to " << n << endl;
	  unsigned long ID;
	  if (TrackParticleID) ID= ParticleID[p]; else ID=0;
	  addP(PCGMsg[ch][n], &(nopPCMsg[ch][n]), x[p], y[p], z[p], u[p], v[p], w[p], q[p], ID, vct);
	  break; // to make sure that a particle is added only to one msg
	} 
	
	
      }// end  for (int n=0; n< CG_numPBCMessages[ch]; n++)
      
      
    } // end for (int p=0; p<nop; p++){
  } // end f (CG_numPBCMessages[ch]>0){ // this particular core has to send BC   
  
  /*
    // the last one, to mark the end of the the meaningful part --> q=0.0
  for (int m=0; m< CG_numPBCMessages[ch]; m++ ){
    //addP(PCGMsg[ch][m], &(nopPCMsg[ch][m]), -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, 0.0, 0, vct);
    PCGMsg[ch][m][nopPCMsg[ch][m]].q= 0.0;
  }
  
  for (int m=0; m< CG_numPBCMessages[ch]; m++){
    cout << "Grid " << numGrid <<" core " << vct->getCartesian_rank() << " on the local grid communicator is sending " << nopPCMsg[ch][m] << " particles as BC to grid " << childNum << " core " << CG_Info[ch][m].RG_core << " (msg " << m << " of " <<CG_numPBCMessages[ch] <<")" << endl; 
  }
  */
  
  
}

void Particles3Dcomm::addP(RepP_struct * Vec, int *num, double x, double y, double z, double u, double v, double w, double q, unsigned long ID, VirtualTopology3D* vct){

  //  cout << "adding P: before adding *num= " << *num << endl; 
 
  (Vec+(*num))->x=x;

  (Vec+(*num))->y=y;
  (Vec+(*num))->z=z;

  (Vec+(*num))->u=u;
  (Vec+(*num))->v=v;
  (Vec+(*num))->w=w;

  (Vec+(*num))->q=q;
  
  (Vec+(*num))->ID=ID;
  
  (*num)++;
  // cout << "adding P: after adding *num= " << *num << endl;

  if (*num > sizePCGMsg){
    // TO DO:RESIZE
    cout << "in addP, numGrid " << numGrid << " core " << vct->getCartesian_rank()  << " in the local grid communicator, you plan on passing too many particles as BC; IMPLEMENT RESIZE; aborting now..." << endl;
    abort();
  }
  return;

}


