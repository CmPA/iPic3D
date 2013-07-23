/*******************************************************************************************
  Particles3Dcomm.cpp  -  Class for particles of the same species, in a 2D space and 3component velocity
  -------------------
developers: Stefano Markidis, Giovanni Lapenta.
 ********************************************************************************************/

#include <mpi.h>
#include <iostream>
#include <math.h>
#include <limits.h>
#include "asserts.h"
#include "VirtualTopology3D.h"
#include "VCtopology3D.h"
#include "CollectiveIO.h"
#include "Collective.h"
#include "ComParticles3D.h"
#include "Alloc.h"
#include "Basic.h"
#include "BcParticles.h"
#include "Grid.h"
#include "Grid3DCU.h"
#include "Field.h"
#include "Moments.h"
#include "MPIdata.h"
#include "ompdefs.h"

#include "Particles3Dcomm.h"

#include "hdf5.h"
#include <vector>
#include <complex>
#include "debug.h"

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
  // see allocate(int species, CollectiveIO* col, VirtualTopology3D* vct, Grid* grid)

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
}
/** constructors fo a single species*/
void Particles3Dcomm::allocate(int species, CollectiveIO * col, VirtualTopology3D * vct, Grid * grid) {
  // info from collectiveIO
  ns = species;
  npcel = col->getNpcel(species);
  npcelx = col->getNpcelx(species);
  npcely = col->getNpcely(species);
  npcelz = col->getNpcelz(species);
  nop = col->getNp(species) / (vct->getNprocs());
  np_tot = col->getNp(species);
  npmax = col->getNpMax(species) / (vct->getNprocs());
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
    dataset_id = H5Dopen(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
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
    dataset_id = H5Dopen(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, y);
    status = H5Dclose(dataset_id);

    // get z
    name_dataset = "/particles/species_" + species_name.str() + "/z/cycle_0";
    dataset_id = H5Dopen(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, z);
    status = H5Dclose(dataset_id);

    // get u
    name_dataset = "/particles/species_" + species_name.str() + "/u/cycle_0";
    dataset_id = H5Dopen(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, u);
    status = H5Dclose(dataset_id);
    // get v
    name_dataset = "/particles/species_" + species_name.str() + "/v/cycle_0";
    dataset_id = H5Dopen(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v);
    status = H5Dclose(dataset_id);
    // get w
    name_dataset = "/particles/species_" + species_name.str() + "/w/cycle_0";
    dataset_id = H5Dopen(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, w);
    status = H5Dclose(dataset_id);
    // get q
    name_dataset = "/particles/species_" + species_name.str() + "/q/cycle_0";
    dataset_id = H5Dopen(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, q);
    status = H5Dclose(dataset_id);
    // ID 
    if (TrackParticleID) {
      // herr_t (*old_func)(void*); // HDF 1.6
      H5E_auto_t old_func;      // HDF 1.8.8
      void *old_client_data;
      H5Eget_auto(H5E_DEFAULT, &old_func, &old_client_data);  // HDF 1.8.8
      /* Turn off error handling */
      // H5Eset_auto(NULL, NULL); // HDF 1.6
      H5Eset_auto(H5E_DEFAULT, 0, 0); // HDF 1.8
      name_dataset = "/particles/species_" + species_name.str() + "/ID/cycle_0";
      dataset_id = H5Dopen(file_id, name_dataset.c_str(), H5P_DEFAULT); // HDF 1.8.8

      // H5Eset_auto(old_func, old_client_data); // HDF 1.6
      H5Eset_auto(H5E_DEFAULT, old_func, old_client_data);
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

}


// A much faster version of this is at EMfields3D::sumMoments
//
void Particles3Dcomm::interpP2G(Field * EMf, Grid * grid, VirtualTopology3D * vct) {
  const double inv_dx = 1.0 / dx;
  const double inv_dy = 1.0 / dy;
  const double inv_dz = 1.0 / dz;
  const double nxn = grid->getNXN();
  const double nyn = grid->getNYN();
  const double nzn = grid->getNZN();
  assert_le(nop,INT_MAX); // else would need to use long long
  // to make memory use scale to a large number of threads we
  // could first apply an efficient parallel sorting algorithm
  // to the particles and then accumulate moments in smaller
  // subarrays.
  {
    assert_le(nop,INT_MAX); // else would need to use long long
    for (int i = 0; i < nop; i++)
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
  npExitingMax = reduceMaxNpExiting(npExitingMax);

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
    // message from XLEFT
    avail1 = unbuffer(b_X_RIGHT);
    avail2 = unbuffer(b_X_LEFT);
    avail3 = unbuffer(b_Y_RIGHT);
    avail4 = unbuffer(b_Y_LEFT);
    avail5 = unbuffer(b_Z_RIGHT);
    avail6 = unbuffer(b_Z_LEFT);
    // if one of these numbers is negative than there is not enough space for particles
    avail = avail1 + avail2 + avail3 + avail4 + avail5 + avail6;
    availALL = reduceNumberParticles(avail);
    if (availALL < 0)
      return (-1);              // too many particles coming, save data nad stop simulation
  }

  return (0);                   // everything was fine


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
int Particles3Dcomm::unbuffer(double *b_) {
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
  result = reduceNumberParticles(rightDomain);
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
/** return the number of particles */ long long Particles3Dcomm::getNOP()  const {
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
  unsigned long localN = 0;
  unsigned long totalN = 0;
  for (int i = 0; i < nBins; i++) {
    localN = f[i];
    MPI_Allreduce(&localN, &totalN, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    f[i] = totalN;
  }
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
