/*******************************************************************************************
  Particles3Dcomm.cpp  -  Class for communication of particles of the same species in 3D
  -------------------
developers: Stefano Markidis, Giovanni Lapenta.
 ********************************************************************************************/
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <math.h>
#include "VirtualTopology3D.h"
#include "VCtopology3D.h"
#include "Collective.h"
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
#include <climits>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define INVALID_PARTICLE   -4.0e32
/**
 * 
 * Class for communication of particles of the same species in 3D
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
}
/** constructors fo a single species*/
void Particles3Dcomm::allocate(int species, long long initnpmax, Collective * col, VirtualTopology3D * vct, Grid * grid) {
  // info from collectiveIO
  ns = species;
  npcel = col->getNpcel(species);
  npcelx = col->getNpcelx(species);
  npcely = col->getNpcely(species);
  npcelz = col->getNpcelz(species);

  // This if is necessary to restart with H5hut-io
  if (initnpmax==0) {
    long ncproc = int(col->getNxc()/col->getXLEN()) *
                  int(col->getNyc()/col->getYLEN()) *
                  int(col->getNzc()/col->getZLEN());
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
  Lx = col->getLx();
  Ly = col->getLy();
  Lz = col->getLz();
  dx = grid->getDX();
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

  dx = grid->getDX();
  dy = grid->getDY();
  dz = grid->getDZ();

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

  x_degenerated = (vct->getXleft_neighbor_P() == vct->getCartesian_rank());
  y_degenerated = (vct->getYleft_neighbor_P() == vct->getCartesian_rank());
  z_degenerated = (vct->getZleft_neighbor_P() == vct->getCartesian_rank());

  x_leftmost = x_degenerated || (vct->getXleft_neighbor_P() > vct->getCartesian_rank());
  y_leftmost = y_degenerated || (vct->getYleft_neighbor_P() > vct->getCartesian_rank());
  z_leftmost = z_degenerated || (vct->getZleft_neighbor_P() > vct->getCartesian_rank());
  x_rightmost = x_degenerated || (vct->getXright_neighbor_P() < vct->getCartesian_rank());
  y_rightmost = y_degenerated || (vct->getYright_neighbor_P() < vct->getCartesian_rank());
  z_rightmost = z_degenerated || (vct->getZright_neighbor_P() < vct->getCartesian_rank());

  no_x_left = x_leftmost || (vct->getXleft_neighbor_P() == MPI_PROC_NULL);
  no_y_left = y_leftmost || (vct->getYleft_neighbor_P() == MPI_PROC_NULL);
  no_z_left = z_leftmost || (vct->getZleft_neighbor_P() == MPI_PROC_NULL);
  no_x_right = x_rightmost || (vct->getXright_neighbor_P() == MPI_PROC_NULL);
  no_y_right = y_rightmost || (vct->getYright_neighbor_P() == MPI_PROC_NULL);
  no_z_right = z_rightmost || (vct->getZright_neighbor_P() == MPI_PROC_NULL);

  x_mirror = !x_degenerated && ((bcPfaceXleft == 1) || (bcPfaceXright == 1));
  y_mirror = !y_degenerated && ((bcPfaceYleft == 1) || (bcPfaceYright == 1));
  z_mirror = !z_degenerated && ((bcPfaceZleft == 1) || (bcPfaceZright == 1));
  x_reemission = !x_degenerated && ((bcPfaceXleft == 2) || (bcPfaceXright == 2));
  y_reemission = !y_degenerated && ((bcPfaceYleft == 2) || (bcPfaceYright == 2));
  z_reemission = !z_degenerated && ((bcPfaceZleft == 2) || (bcPfaceZright == 2));

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
  buffer_size_x = (long long) (0.025 * nop) + 10; // assume 2.5% of the particles in the processors are leaving
  buffer_size_y = buffer_size_x;
  buffer_size_z = buffer_size_x;

  b_X_RIGHT.resize (buffer_size_x*nVar, INVALID_PARTICLE);
  b_X_LEFT.resize (buffer_size_x*nVar, INVALID_PARTICLE);
  b_Y_RIGHT.resize (buffer_size_y*nVar, INVALID_PARTICLE);
  b_Y_LEFT.resize (buffer_size_y*nVar, INVALID_PARTICLE);
  b_Z_RIGHT.resize (buffer_size_z*nVar, INVALID_PARTICLE);
  b_Z_LEFT.resize (buffer_size_z*nVar, INVALID_PARTICLE);

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
void c_vDist::init(int ispec, double vX, double vY, double vZ, int bi, int bj, int bk, double vR, double vFact, Collective * col, Grid * grid) {
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

void c_vDist::add(double x, double y, double z, double u, double v, double w) {

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
void Particles3Dcomm::Add_vDist3D() {

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
void Particles3Dcomm::Write_vDist3D(string SaveDirName) {

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

/** communicate buffers or apply boundary conditions for particles: */
void Particles3Dcomm::communicate(VirtualTopology3D * vct) {
  // number of particles in the wrong domain
  long long comm_x, comm_y, comm_z;

  npExitXright = 0L, npExitXleft = 0L, npExitYright = 0L, npExitYleft = 0L, npExitZright = 0L, npExitZleft = 0L;
  long long np_current = 0L;
  int result;

  while (np_current < nop) {

    // check for boundary conditions
    result = bc_apply (&x[np_current], &y[np_current], &z[np_current], &u[np_current], &v[np_current], &w[np_current]);
    if (result == -1) { del_pack(np_current); continue; }

    // particle leaving the domain need to be communicated in a large enough buffer
    if (x_out_left) {
      npExitXleft++;
      if (x_leftmost) x[np_current] += Lx;
      if (npExitXleft >= buffer_size_x) resize_buffers(b_X_LEFT, b_X_RIGHT, buffer_size_x, npExitXleft);
      buffer_leaving(b_X_LEFT, (npExitXleft-1)*nVar, np_current);
    }
    else if (x_out_right) {
      npExitXright++;
      if (x_rightmost) x[np_current] -= Lx;
      if (npExitXright >= buffer_size_x) resize_buffers(b_X_LEFT, b_X_RIGHT, buffer_size_x, npExitXright);
      buffer_leaving(b_X_RIGHT, (npExitXright-1)*nVar, np_current);
    }
    else if (y_out_left) {
      npExitYleft++;
      if (y_leftmost) y[np_current] += Ly;
      if (npExitYleft >= buffer_size_y) resize_buffers(b_Y_LEFT, b_Y_RIGHT, buffer_size_y, npExitYleft);
      buffer_leaving(b_Y_LEFT, (npExitYleft-1)*nVar, np_current);
    }
    else if (y_out_right) {
      npExitYright++;
      if (y_rightmost) y[np_current] -= Ly;
      if (npExitYright >= buffer_size_y) resize_buffers(b_Y_LEFT, b_Y_RIGHT, buffer_size_y, npExitYright);
      buffer_leaving(b_Y_RIGHT, (npExitYright-1)*nVar, np_current);
    }
    else if (z_out_left) {
      npExitZleft++;
      if (z_leftmost) z[np_current] += Lz;
      if (npExitZleft >= buffer_size_z) resize_buffers(b_Z_LEFT, b_Z_RIGHT, buffer_size_z, npExitZleft);
      buffer_leaving(b_Z_LEFT, (npExitZleft-1)*nVar, np_current);
    }
    else if (z_out_right) {
      npExitZright++;
      if (z_rightmost) z[np_current] -= Lz;
      if (npExitZright >= buffer_size_z) resize_buffers(b_Z_LEFT, b_Z_RIGHT, buffer_size_z, npExitZright);
      buffer_leaving(b_Z_RIGHT, (npExitZright-1)*nVar, np_current);
    }
    else {
      // particle is still in the domain, proceed with the next particle
      np_current++;
    }
  }

  // put end markers into buffers
  b_X_LEFT[npExitXleft * nVar] = INVALID_PARTICLE;
  b_Y_LEFT[npExitYleft * nVar] = INVALID_PARTICLE;
  b_Z_LEFT[npExitZleft * nVar] = INVALID_PARTICLE;
  b_X_RIGHT[npExitXright * nVar] = INVALID_PARTICLE;
  b_Y_RIGHT[npExitYright * nVar] = INVALID_PARTICLE;
  b_Z_RIGHT[npExitZright * nVar] = INVALID_PARTICLE;

  // local maximum number of particles communicated in each direction
  comm_x = max(npExitXleft, npExitXright);
  comm_y = max(npExitYleft, npExitYright);
  comm_z = max(npExitZleft, npExitZright);

  iterate_communication (b_X_LEFT, b_X_RIGHT, b_Y_LEFT, b_Y_RIGHT, b_Z_LEFT, b_Z_RIGHT, comm_x, comm_y, comm_z, buffer_size_x, buffer_size_y, buffer_size_z, vct);
}

/** iterate communication of buffers and unbuffer received buffers to local domain: */
int Particles3Dcomm::iterate_communication(std::vector<double>& bxl, std::vector<double>& bxr, std::vector<double>& byl, std::vector<double>& byr, std::vector<double>& bzl, std::vector<double>& bzr, long long& num_x, long long& num_y, long long& num_z, long long& size_x, long long& size_y, long long& size_z, VirtualTopology3D * vct, int add_size) {

  // number of particles still in the wrong domain
  long long wrong_x = 0L, wrong_y = 0L, wrong_z = 0L;
  // vector of particles in the wrong domain
  std::vector<double> wxl, wxr, wyl, wyr, wzl, wzr;
  // variable for memory availability of space for new particles
  long long avail;

  // global maximum number of particles communicated in each direction
  num_x = globalMaximum(num_x);
  num_y = globalMaximum(num_y);
  num_z = globalMaximum(num_z);

  if (num_x + num_y + num_z > 0L) {

    // resize buffers, if necessary
    if (num_x+add_size > size_x) {
      if (vct->getCartesian_rank() == 0)
        cout << "resizing X-buffer: " << size_x << " => " << num_x << " particles" << endl;
      resize_buffers(bxl, bxr, size_x, num_x, false);
    }
    if (num_y+add_size > size_y) {
      if (vct->getCartesian_rank() == 0)
        cout << "resizing Y-buffer: " << size_y << " => " << num_y << " particles" << endl;
      resize_buffers(byl, byr, size_y, num_y, false);
    }
    if (num_z+add_size > size_z) {
      if (vct->getCartesian_rank() == 0)
        cout << "resizing Z-buffer: " << size_z << " => " << num_z << " particles" << endl;
      resize_buffers(bzl, bzr, size_z, num_z, false);
    }

    // communicate in the X direction
    if (num_x > 0L) {
      if (num_x*nVar > INT_MAX) {
        cout << "ERROR: X-buffer too large: " << num_x*nVar << " > " << INT_MAX << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
        return (-1);
      }
      communicateParticlesDIR((int) num_x*nVar, vct->getCartesian_rank(), vct->getXright_neighbor_P(), vct->getXleft_neighbor_P(), 0, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), bxr.data(), bxl.data());
    }

    // communicate in the Y direction
    if (num_y > 0L) {
      if (num_y*nVar > INT_MAX) {
        cout << "ERROR: Y-buffer too large: " << num_y*nVar << " > " << INT_MAX << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
        return (-1);
      }
      communicateParticlesDIR((int) num_y*nVar, vct->getCartesian_rank(), vct->getYright_neighbor_P(), vct->getYleft_neighbor_P(), 1, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), byr.data(), byl.data());
    }

    // communicate in the Z direction
    if (num_z > 0L) {
      if (num_z*nVar > INT_MAX) {
        cout << "ERROR: Z-buffer too large: " << num_z*nVar << " > " << INT_MAX << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
        return (-1);
      }
      communicateParticlesDIR((int) num_z*nVar, vct->getCartesian_rank(), vct->getZright_neighbor_P(), vct->getZleft_neighbor_P(), 2, vct->getXLEN(), vct->getYLEN(), vct->getZLEN(), bzr.data(), bzl.data());
    }

    // put received particles in the local domain
    wxl.reserve(100*nVar);
    wxr.reserve(100*nVar);
    wyl.reserve(100*nVar);
    wyr.reserve(100*nVar);
    wzl.reserve(100*nVar);
    wzr.reserve(100*nVar);
    avail  = unbuffer(bxr, wxl, wxr, wyl, wyr, wzl, wzr, wrong_x, wrong_y, wrong_z);
    avail += unbuffer(bxl, wxl, wxr, wyl, wyr, wzl, wzr, wrong_x, wrong_y, wrong_z);
    avail += unbuffer(byr, wxl, wxr, wyl, wyr, wzl, wzr, wrong_x, wrong_y, wrong_z);
    avail += unbuffer(byl, wxl, wxr, wyl, wyr, wzl, wzr, wrong_x, wrong_y, wrong_z);
    avail += unbuffer(bzr, wxl, wxr, wyl, wyr, wzl, wzr, wrong_x, wrong_y, wrong_z);
    avail += unbuffer(bzl, wxl, wxr, wyl, wyr, wzl, wzr, wrong_x, wrong_y, wrong_z);

    // if any of these numbers is negative there is not enough space to store the incoming particles
    avail = globalSum(avail);
    if (avail < 0) return (-1); // save data and stop simulation

    // recursive call to treat and more particles that are in the wrong domain
    iterate_communication (wxl, wxr, wyl, wyr, wzl, wzr, wrong_x, wrong_y, wrong_z, wrong_x, wrong_y, wrong_z, vct, 0);
  }
  return (0);
}

/** apply boundary conditions */
inline int Particles3Dcomm::bc_apply(double* x, double* y, double* z, double* u, double* v, double* w) {

    x_out_left = (*x < xstart);
    x_out_right = (*x > xend);
    y_out_left = (*y < ystart);
    y_out_right = (*y > yend);
    z_out_left = (*z < zstart);
    z_out_right = (*z > zend);

    if (no_x_left && x_out_left) {
      if (x_degenerated) BCpart_left_degenerated(x,Lx);
      else if (x_mirror) BCpart_left_mirror(x,u,Lx);
      else if (x_reemission) BCpart_left_reemission(x,u,v,w,Lx,uth,vth,wth);
      else { return (-1); }
      x_out_left = false;
    }
    else if (no_x_right && x_out_right) {
      if (x_degenerated) BCpart_right_degenerated(x,Lx);
      else if (x_mirror) BCpart_right_mirror(x,u,Lx);
      else if (x_reemission) BCpart_right_reemission(x,u,v,w,Lx,uth,vth,wth);
      else { return (-1); }
      x_out_right = false;
    }
    if (no_y_left && y_out_left) {
      if (y_degenerated) BCpart_left_degenerated(y,Ly);
      else if (y_mirror) BCpart_left_mirror(y,v,Ly);
      else if (y_reemission) BCpart_left_reemission(y,v,u,w,Ly,vth,uth,wth);
      else { return (-1); }
      y_out_left = false;
    }
    else if (no_y_right && y_out_right) {
      if (y_degenerated) BCpart_right_degenerated(y,Ly);
      else if (y_mirror) BCpart_right_mirror(y,v,Ly);
      else if (y_reemission) BCpart_right_reemission(y,v,u,w,Ly,vth,uth,wth);
      else { return (-1); }
      y_out_right = false;
    }
    if (no_z_left && z_out_left) {
      if (z_degenerated) BCpart_left_degenerated(z,Lz);
      else if (z_mirror) BCpart_left_mirror(z,w,Lz);
      else if (z_reemission) BCpart_left_reemission(z,w,u,v,Lz,wth,uth,vth);
      else { return (-1); }
      z_out_left = false;
    }
    else if (no_z_right && z_out_right) {
      if (z_degenerated) BCpart_right_degenerated(z,Lz);
      else if (z_mirror) BCpart_right_mirror(z,w,Lz);
      else if (z_reemission) BCpart_right_reemission(z,w,u,v,Lz,wth,uth,vth);
      else { return (-1); }
      z_out_right = false;
    }

  return (0);
}

/** resize a buffer */
void Particles3Dcomm::resize_buffers(std::vector<double>& b_left, std::vector<double>& b_right, long long& size, long long request_size, bool extend) {
  double *temp = NULL;
  long long old_size = size * nVar;
  size = request_size + 1;
  if (extend) size += ((long long) (request_size*0.1 + 0.025*nop)) + 100;
  long long new_size = size * nVar;

  // resize and initialize
  b_left.resize(new_size, INVALID_PARTICLE);
  b_right.resize(new_size, INVALID_PARTICLE);
}

/** put a leaving particle to the communication buffer */
inline void Particles3Dcomm::buffer_leaving(std::vector<double>& buffer, long long pos, long long& np_current) {
  buffer[pos] = x[np_current];
  buffer[pos+1] = y[np_current];
  buffer[pos+2] = z[np_current];
  buffer[pos+3] = u[np_current];
  buffer[pos+4] = v[np_current];
  buffer[pos+5] = w[np_current];
  buffer[pos+6] = q[np_current];
  if (TrackParticleID) buffer[pos+7] = ParticleID[np_current];
  del_pack(np_current);
}

/** Unbuffer the last communication */
int Particles3Dcomm::unbuffer(std::vector<double>& buffer, std::vector<double>& wxl, std::vector<double>& wxr, std::vector<double>& wyl, std::vector<double>& wyr, std::vector<double>& wzl, std::vector<double>& wzr, long long& wrong_x, long long& wrong_y, long long& wrong_z) {
  double *start;
  long long pos = 0L, size;
  int result;
  // put the new particles at the end of the array, and update the number of particles
  start = buffer.data();
  size = buffer.size();
  while ((*start != INVALID_PARTICLE) && (pos < size)) {

    result = bc_apply (start, start+1, start+2, start+3, start+4, start+5);
    if (result == 0) {
      // these particles need further communication
      if (x_out_left) {
        rebuffer(start, wxl, wrong_x);
      }
      else if (x_out_right) {
        rebuffer(start, wxr, wrong_x);
      }
      else if (y_out_left) {
        rebuffer(start, wyl, wrong_y);
      }
      else if (y_out_right) {
        rebuffer(start, wyr, wrong_y);
      }
      else if (z_out_left) {
        rebuffer(start, wzl, wrong_z);
      }
      else if (z_out_right) {
        rebuffer(start, wzr, wrong_z);
      }
      else {
        x[nop] = start[0];
        y[nop] = start[1];
        z[nop] = start[2];
        u[nop] = start[3];
        v[nop] = start[4];
        w[nop] = start[5];
        q[nop] = start[6];
        if (TrackParticleID) ParticleID[nop] = (unsigned long) start[7];
        nop++;
        if (nop > npmax) {
          cout << "Number of particles in the domain " << nop << " and maxpart = " << npmax << endl;
          MPI_Abort(MPI_COMM_WORLD, -1);
          return (-1);              // end the simulation because you dont have enough space on the array
        }
      }
    }
    start += nVar;
    pos += nVar;
  }
  buffer[0] = INVALID_PARTICLE;
  return (0);                   // everything was fine
}

/** This unbuffer the last communication */
inline void Particles3Dcomm::rebuffer(double *start, std::vector<double>& buffer, long long& wrong) {
  buffer.push_back (start[0]);
  buffer.push_back (start[1]);
  buffer.push_back (start[2]);
  buffer.push_back (start[3]);
  buffer.push_back (start[4]);
  buffer.push_back (start[5]);
  buffer.push_back (start[6]);
  if (TrackParticleID) buffer.push_back (start[7]);
  wrong++;
}

/** Delete the a particle from the array and pack the array,
 * update the number of particles that are leaving.
 * For deleting the particle from the array take the last particle and
 * put it in the position of the particle you want to delete.
 * @param np = the index of the particle that must be deleted
 */
void Particles3Dcomm::del_pack(long long np_current) {
  nop--;
  x[np_current] = x[nop];
  y[np_current] = y[nop];
  z[np_current] = z[nop];
  u[np_current] = u[nop];
  v[np_current] = v[nop];
  w[np_current] = w[nop];
  q[np_current] = q[nop];
  if (TrackParticleID) ParticleID[np_current] = ParticleID[nop];
}

/** calculate the maximum number leaving from this domain */
long long Particles3Dcomm::maxNpExiting(long long *max_x, long long *max_y, long long *max_z) {
  *max_x = max(npExitXleft, npExitXright);
  *max_y = max(npExitYleft, npExitYright);
  *max_z = max(npExitZleft, npExitZright);
  long long max_xyz = max(*max_x, *max_y);
  max_xyz = max(max_xyz, *max_z);
  return (max_xyz);
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
double Particles3Dcomm::getKe() {
  double localKe = 0.0;
  double totalKe = 0.0;
  for (register long long i = 0; i < nop; i++)
    localKe += .5 * (q[i] / qom) * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
  MPI_Allreduce(&localKe, &totalKe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalKe);
}
/** return the total momentum */
double Particles3Dcomm::getP() {
  double localP = 0.0;
  double totalP = 0.0;
  for (register long long i = 0; i < nop; i++)
    localP += (q[i] / qom) * sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
  MPI_Allreduce(&localP, &totalP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalP);
}

/** return the highest kinetic energy */
double Particles3Dcomm::getMaxVelocity() {
  double localVel = 0.0;
  double maxVel = 0.0;
  for (long long i = 0; i < nop; i++)
    localVel = max(localVel, sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]));
  MPI_Allreduce(&localVel, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return (maxVel);
}


/** get energy spectrum */
unsigned long *Particles3Dcomm::getVelocityDistribution(int nBins, double maxVel) {
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

  MPI_Allreduce(MPI_IN_PLACE, f, nBins, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

  return f;
}


/** print particles info */
void Particles3Dcomm::Print(VirtualTopology3D * vct) const {
  cout << endl;
  cout << "Number of Particles: " << nop << endl;
  cout << "Subgrid (" << vct->getCoordinates(0) << "," << vct->getCoordinates(1) << "," << vct->getCoordinates(2) << ")" << endl;
  cout << "Xin = " << xstart << "; Xfin = " << xend << endl;
  cout << "Yin = " << ystart << "; Yfin = " << yend << endl;
  cout << "Zin = " << zstart << "; Zfin = " << zend << endl;
  cout << "Number of species = " << ns << endl;
  for (long long i = 0; i < nop; i++)
    cout << "Particles #" << i << " x=" << x[i] << " y=" << y[i] << " z=" << z[i] << " u=" << u[i] << " v=" << v[i] << " w=" << w[i] << endl;
  cout << endl;
}
/** print just the number of particles */
void Particles3Dcomm::PrintNp(VirtualTopology3D * vct)  const {
  cout << endl;
  cout << "Number of Particles of species " << ns << ": " << nop << endl;
  cout << "Subgrid (" << vct->getCoordinates(0) << "," << vct->getCoordinates(1) << "," << vct->getCoordinates(2) << ")" << endl;
  cout << endl;
}
