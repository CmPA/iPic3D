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
#include <algorithm> // for swap
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
#include "MPIdata.h"
#include "ompdefs.h"
#include "ipicmath.h"

#include "Particle.h"
#include "Particles3Dcomm.h"
#include "Parameters.h"

#include "hdf5.h"
#include <vector>
#include <complex>
#include "debug.h"
#include "TimeTasks.h"

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
Particles3Dcomm::Particles3Dcomm(){
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
  delete[]ParticleID;
  // AoS representation
  delete[]_pcls;
  // average position used in particle advance
  delete[]_xavg;
  delete[]_yavg;
  delete[]_zavg;
  // deallocate alternate storage
  delete[]_xtmp;
  delete[]_ytmp;
  delete[]_ztmp;
  delete[]_utmp;
  delete[]_vtmp;
  delete[]_wtmp;
  delete[]_qtmp;
  delete[]_ParticleIDtmp;
  delete[] _pclstmp;
  // extra xavg for sort
  delete[]_xavgtmp;
  delete[]_yavgtmp;
  delete[]_zavgtmp;
  // deallocate buffers
  delete[]b_X_RIGHT;
  delete[]b_X_LEFT;
  delete[]b_Y_RIGHT;
  delete[]b_Y_LEFT;
  delete[]b_Z_RIGHT;
  delete[]b_Z_LEFT;
  delete numpcls_in_bucket;
  delete numpcls_in_bucket_now;
  delete bucket_offset;
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
  // ensure that npmax is a multiple of AoS_PCLS_AT_A_TIME
  npmax = roundup_to_multiple(npmax,AoS_PCLS_AT_A_TIME);
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
  inv_dx = 1/dx;
  inv_dy = 1/dy;
  inv_dz = 1/dz;

  nxn = grid->getNXN();
  nyn = grid->getNYN();
  nzn = grid->getNZN();
  nxc = grid->getNXC();
  nyc = grid->getNYC();
  nzc = grid->getNZC();
  assert_eq(nxc,nxn-1);
  assert_eq(nyc,nyn-1);
  assert_eq(nzc,nzn-1);
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
  //
  // allocate arrays for sorting particles
  //
  numpcls_in_bucket = new array3_int(nxc,nyc,nzc);
  numpcls_in_bucket_now = new array3_int(nxc,nyc,nzc);
  bucket_offset = new array3_int(nxc,nyc,nzc);
  //num_threads = omp_get_max_threads();
  //numpcls_in_bucket_thr = (arr3_int*)malloc(sizeof(void*)*num_threads);
  //for(int i=0; i<num_threads; i++)
  //{
  //  numpcls_in_bucket_thr[i] = new array3_int(nxc,nyc,nzc);
  //}
  
  //
  // //////////////////////////////////////////////////////////////
  // ////////////// ALLOCATE ARRAYS /////////////////////////
  // //////////////////////////////////////////////////////////////
  //
  // SoA particle representation
  //
  // positions
  particleType = ParticleType::SoA;
  x = new double[npmax];
  y = new double[npmax];
  z = new double[npmax];
  // velocities
  u = new double[npmax];
  v = new double[npmax];
  w = new double[npmax];
  // charge
  q = new double[npmax];
  // average positions, used in iterative particle advance
  _xavg = 0;
  _yavg = 0;
  _zavg = 0;
  //if(Parameters::get_USING_XAVG())
  //{
  //  xavg = new double[npmax];
  //  yavg = new double[npmax];
  //  zavg = new double[npmax];
  //}
  _xtmp = 0;
  _ytmp = 0;
  _ztmp = 0;
  _utmp = 0;
  _vtmp = 0;
  _wtmp = 0;
  _qtmp = 0;
  _xavgtmp = 0;
  _yavgtmp = 0;
  _zavgtmp = 0;
  _pcls = 0;
  _pclstmp = 0;
  _ParticleIDtmp = 0;
  // accessors for data that should be allocated only if needed
  //
  if(Parameters::get_USING_XAVG())
  {
    _xavg=new double[npmax];
    _yavg=new double[npmax];
    _zavg=new double[npmax];
  }
  if(Parameters::get_SORTING_SOA())
  {
    _xtmp=new double[npmax];
    _ytmp=new double[npmax];
    _ztmp=new double[npmax];
    _utmp=new double[npmax];
    _vtmp=new double[npmax];
    _wtmp=new double[npmax];
    _qtmp=new double[npmax];
    _xavgtmp=new double[npmax];
    _yavgtmp=new double[npmax];
    _zavgtmp=new double[npmax];
    if(TrackParticleID)
    {
      _ParticleIDtmp = new long long[npmax];
    }
  }
  if(Parameters::get_USING_AOS())
  {
    assert_eq(sizeof(SpeciesParticle),64);
    //_pcls = AlignedAlloc(SpeciesParticle,npmax);
    _pcls = new SpeciesParticle[npmax];
    if(Parameters::get_SORTING_PARTICLES())
    {
      //_pclstmp = AlignedAlloc(SpeciesParticle,npmax);
      _pclstmp = new SpeciesParticle[npmax];
    }
    #ifdef __INTEL_COMPILER
      assert_eq(sizeof(SpeciesParticle),64);
      ALIGNED(_pcls);
    #endif
  }

  ParticleID = 0;
  // ID
  if (TrackParticleID) {
    ParticleID = new long long[npmax];
    //if(Parameters::get_SORTING_PARTICLES())
    //  _ParticleIDtmp = new long long[npmax];
    BirthRank[0] = vct->getCartesian_rank();
    if (vct->getNprocs() > 1)
      BirthRank[1] = (int) ceil(log10((double) (vct->getNprocs())));  // Number of digits needed for # of process in ID
    else
      BirthRank[1] = 1;
    if (BirthRank[1] + (int) ceil(log10((double) (npmax))) > 10 && BirthRank[0] == 0) {
      cerr << "Error: can't Track particles in Particles3Dcomm::allocate" << endl;
      cerr << "long long 'ParticleID' cannot store all the particles" << endl;
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
        for (int counter = 0; counter < nop; counter++)
          ParticleID[counter] = counter * (long long) pow(10.0, BirthRank[1]) + BirthRank[0];
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
  // assert_le(nop,(long long)INT_MAX); // else would need to use long long
  // to make memory use scale to a large number of threads we
  // could first apply an efficient parallel sorting algorithm
  // to the particles and then accumulate moments in smaller
  // subarrays.
  {
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
  int np_current = 0, nplast = nop - 1;

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
void Particles3Dcomm::bufferXleft(double *b_, int np_current, VirtualTopology3D * vct) {
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
void Particles3Dcomm::bufferXright(double *b_, int np_current, VirtualTopology3D * vct) {
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
inline void Particles3Dcomm::bufferYleft(double *b_, int np_current, VirtualTopology3D * vct) {
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
inline void Particles3Dcomm::bufferYright(double *b_, int np_current, VirtualTopology3D * vct) {
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
inline void Particles3Dcomm::bufferZleft(double *b_, int np_current, VirtualTopology3D * vct) {
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
inline void Particles3Dcomm::bufferZright(double *b_, int np_current, VirtualTopology3D * vct) {
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
  int np_current = 0;
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
      ParticleID[nop] = (long long) b_[nVar * np_current + 7];
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
void Particles3Dcomm::del_pack(int np_current, int *nplast) {
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

/** return the Kinetic energy */
double Particles3Dcomm::getKe() {
  double localKe = 0.0;
  double totalKe = 0.0;
  for (register int i = 0; i < nop; i++)
    localKe += .5 * (q[i] / qom) * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
  MPI_Allreduce(&localKe, &totalKe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalKe);
}
/** return the total momentum */
double Particles3Dcomm::getP() {
  double localP = 0.0;
  double totalP = 0.0;
  for (register int i = 0; i < nop; i++)
    localP += (q[i] / qom) * sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
  MPI_Allreduce(&localP, &totalP, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalP);
}

/** return the highest kinetic energy */
double Particles3Dcomm::getMaxVelocity() {
  double localVel = 0.0;
  double maxVel = 0.0;
  for (int i = 0; i < nop; i++)
    localVel = max(localVel, sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]));
  MPI_Allreduce(&localVel, &maxVel, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return (maxVel);
}


/** get energy spectrum */
long long *Particles3Dcomm::getVelocityDistribution(int nBins, double maxVel) {
  long long *f = new long long[nBins];
  for (int i = 0; i < nBins; i++)
    f[i] = 0;
  double Vel = 0.0;
  double dv = maxVel / nBins;
  int bin = 0;
  for (int i = 0; i < nop; i++) {
    Vel = sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
    bin = int (floor(Vel / dv));
    if (bin >= nBins)
      f[nBins - 1] += 1;
    else
      f[bin] += 1;
  }
  long long localN = 0;
  long long totalN = 0;
  for (int i = 0; i < nBins; i++) {
    localN = f[i];
    MPI_Allreduce(&localN, &totalN, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
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
  for (int i = 0; i < nop; i++)
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

/***** particle sorting routines *****/

void Particles3Dcomm::sort_particles_serial(Grid * grid, VirtualTopology3D * vct)
{
  switch(particleType)
  {
    case ParticleType::AoS:
      sort_particles_serial_AoS(grid,vct);
      break;
    case ParticleType::SoA:
      sort_particles_serial_SoA(grid,vct);
      break;
    default:
      unsupported_value_error(particleType);
  }
}

// need to sort and communicate particles after each iteration
void Particles3Dcomm::sort_particles_serial_AoS(
  Grid * grid, VirtualTopology3D * vct)
{
  SpeciesParticle* pcls = fetch_pcls();
  SpeciesParticle* pclstmp = fetch_pclstmp();
  {
    numpcls_in_bucket->setall(0);
    // iterate through particles and count where they will go
    for (int pidx = 0; pidx < nop; pidx++)
    {
      const SpeciesParticle& pcl = get_pcl(pidx);
      // get the cell indices of the particle
      int cx,cy,cz;
      get_safe_cell_for_pos(cx,cy,cz,pcl.get_x(),pcl.get_y(),pcl.get_z());

      // increment the number of particles in bucket of this particle
      (*numpcls_in_bucket)[cx][cy][cz]++;
    }

    // compute prefix sum to determine initial position
    // of each bucket (could parallelize this)
    //
    int accpcls=0;
    for(int cx=0;cx<nxc;cx++)
    for(int cy=0;cy<nyc;cy++)
    for(int cz=0;cz<nzc;cz++)
    {
      (*bucket_offset)[cx][cy][cz] = accpcls;
      accpcls += (*numpcls_in_bucket)[cx][cy][cz];
    }
    assert_eq(accpcls,nop);

    numpcls_in_bucket_now->setall(0);
    // put the particles where they are supposed to go
    for (int pidx = 0; pidx < nop; pidx++)
    {
      const SpeciesParticle& pcl = get_pcl(pidx);
      // get the cell indices of the particle
      int cx,cy,cz;
      get_safe_cell_for_pos(cx,cy,cz,pcl.get_x(),pcl.get_y(),pcl.get_z());

      // compute where the data should go
      const int numpcls_now = (*numpcls_in_bucket_now)[cx][cy][cz]++;
      const int outpidx = (*bucket_offset)[cx][cy][cz] + numpcls_now;
      assert_lt(outpidx, nop);
      assert_ge(outpidx, 0);
      assert_lt(pidx, nop);
      assert_ge(pidx, 0);

      // copy particle data to new location
      //
      pclstmp[outpidx] = pcl;
    }
    // swap the tmp particle memory with the official particle memory
    {
      swap(_pclstmp,_pcls);
    }

    // check if the particles were sorted incorrectly
    if(true)
    {
      for(int cx=0;cx<nxc;cx++)
      for(int cy=0;cy<nyc;cy++)
      for(int cz=0;cz<nzc;cz++)
      {
        assert_eq((*numpcls_in_bucket_now)[cx][cy][cz], (*numpcls_in_bucket)[cx][cy][cz]);
      }
    }
  }
}

// need to sort and communicate particles after each iteration
void Particles3Dcomm::sort_particles_serial_SoA(
  Grid * grid, VirtualTopology3D * vct)
{
  double * xtmp = fetch_xtmp();
  double * ytmp = fetch_ytmp();
  double * ztmp = fetch_ztmp();
  double * utmp = fetch_utmp();
  double * vtmp = fetch_vtmp();
  double * wtmp = fetch_wtmp();
  double * qtmp = fetch_qtmp();

  long long* ParticleIDtmp = 0;
  if (TrackParticleID)
  {
    assert(ParticleID);
    ParticleIDtmp = fetch_ParticleIDtmp();
    assert(fetch_ParticleIDtmp());
    assert(ParticleIDtmp);
  }

  // sort the particles
  {
    numpcls_in_bucket->setall(0);
    // iterate through particles and count where they will go
    for (int pidx = 0; pidx < nop; pidx++)
    {
      // get the cell indices of the particle
      //
      int cx,cy,cz;
      get_safe_cell_for_pos(cx,cy,cz,x[pidx],y[pidx],z[pidx]);
      //
      // is it better just to recompute this?
      //
      //xcell[pidx]=cx;
      //ycell[pidx]=cy;
      //zcell[pidx]=cz;

      // increment the number of particles in bucket of this particle
      (*numpcls_in_bucket)[cx][cy][cz]++;
    }

    // compute prefix sum to determine initial position
    // of each bucket (could parallelize this)
    //
    int accpcls=0;
    for(int cx=0;cx<nxc;cx++)
    for(int cy=0;cy<nyc;cy++)
    for(int cz=0;cz<nzc;cz++)
    {
      (*bucket_offset)[cx][cy][cz] = accpcls;
      accpcls += (*numpcls_in_bucket)[cx][cy][cz];
    }
    assert_eq(accpcls,nop);

    numpcls_in_bucket_now->setall(0);

    // put the particles where they are supposed to go
    for (int pidx = 0; pidx < nop; pidx++)
    {
      // get the cell indices of the particle
      //
      int cx,cy,cz;
      get_safe_cell_for_pos(cx,cy,cz,x[pidx],y[pidx],z[pidx]);
      //
      //cx = xcell[pidx];
      //cy = ycell[pidx];
      //cz = zcell[pidx];

      // compute where the data should go
      const int numpcls_now = (*numpcls_in_bucket_now)[cx][cy][cz]++;
      const int outpidx = (*bucket_offset)[cx][cy][cz] + numpcls_now;
      assert_lt(outpidx, nop);
      assert_ge(outpidx, 0);
      assert_lt(pidx, nop);
      assert_ge(pidx, 0);

      // copy particle data to new location
      //
      xtmp[outpidx] = x[pidx];
      ytmp[outpidx] = y[pidx];
      ztmp[outpidx] = z[pidx];
      utmp[outpidx] = u[pidx];
      vtmp[outpidx] = v[pidx];
      wtmp[outpidx] = w[pidx];
      qtmp[outpidx] = q[pidx];
      if (TrackParticleID)
      {
        ParticleIDtmp[outpidx] = ParticleID[pidx];
      }
    }
    // swap the tmp particle memory with the official particle memory
    {
      swap(_xtmp,x);
      swap(_ytmp,y);
      swap(_ztmp,z);
      swap(_utmp,u);
      swap(_vtmp,v);
      swap(_wtmp,w);
      swap(_qtmp,q);
      swap(_ParticleIDtmp,ParticleID);
    }

    // check that the number of bins was correct
    //
    if(true)
    {
      for(int cx=0;cx<nxc;cx++)
      for(int cy=0;cy<nyc;cy++)
      for(int cz=0;cz<nzc;cz++)
      {
        assert_eq((*numpcls_in_bucket_now)[cx][cy][cz], (*numpcls_in_bucket)[cx][cy][cz]);
      }
    }
    // confirm that the particles were sorted correctly
    if(false)
    {
      for(int cx=0;cx<nxc;cx++)
      for(int cy=0;cy<nyc;cy++)
      for(int cz=0;cz<nzc;cz++)
      {
        const int numpcls_in_cell = get_numpcls_in_bucket(cx,cy,cz);
        const int bucket_offset = get_bucket_offset(cx,cy,cz);
        const int bucket_end = bucket_offset+numpcls_in_cell;
        for(int pidx=bucket_offset; pidx<bucket_end; pidx++)
        {
          // confirm that particle is in correct cell
          {
            int cx_,cy_,cz_;
            get_safe_cell_for_pos(cx_,cy_,cz_,x[pidx],y[pidx],z[pidx]);
            if((cx_!=cx)
             ||(cy_!=cy)
             ||(cz_!=cz))
            {
              dprintf("\n\t cx =%d, cy =%d, cz =%d"
                      "\n\t cx_=%d, cy_=%d, cz_=%d"
                      "\n\t cxf=%f, cyf=%f, czf=%f",
                      cx,cy,cz,
                      cx_,cy_,cz_,
                      1.+(x[pidx]-xstart)*inv_dx,
                      1.+(y[pidx]-ystart)*inv_dy,
                      1.+(z[pidx]-zstart)*inv_dz);
            }
            assert_eq(cx_,cx);
            assert_eq(cy_,cy);
            assert_eq(cz_,cz);
          }
        }
      }
    }
  }
}

// need to sort and communicate particles after each iteration
void Particles3Dcomm::sort_particles_serial_SoA_by_xavg(
  Grid * grid, VirtualTopology3D * vct)
{
  double * xtmp = fetch_xtmp();
  double * ytmp = fetch_ytmp();
  double * ztmp = fetch_ztmp();
  double * utmp = fetch_utmp();
  double * vtmp = fetch_vtmp();
  double * wtmp = fetch_wtmp();
  double * qtmp = fetch_qtmp();
  double* xavg = fetch_xavg();
  double* yavg = fetch_yavg();
  double* zavg = fetch_zavg();
  double* xavgtmp = fetch_xavgtmp();
  double* yavgtmp = fetch_yavgtmp();
  double* zavgtmp = fetch_zavgtmp();

  long long* ParticleIDtmp = 0;
  if (TrackParticleID) ParticleIDtmp = fetch_ParticleIDtmp();

  // sort the particles
  {
    numpcls_in_bucket->setall(0);
    // iterate through particles and count where they will go
    for (int pidx = 0; pidx < nop; pidx++)
    {
      // get the cell indices of the particle
      //
      int cx,cy,cz;
      get_safe_cell_for_pos(cx,cy,cz,xavg[pidx],yavg[pidx],zavg[pidx]);
      //
      // is it better just to recompute this?
      //
      //xcell[pidx]=cx;
      //ycell[pidx]=cy;
      //zcell[pidx]=cz;

      // increment the number of particles in bucket of this particle
      (*numpcls_in_bucket)[cx][cy][cz]++;
    }

    // compute prefix sum to determine initial position
    // of each bucket (could parallelize this)
    //
    int accpcls=0;
    for(int cx=0;cx<nxc;cx++)
    for(int cy=0;cy<nyc;cy++)
    for(int cz=0;cz<nzc;cz++)
    {
      (*bucket_offset)[cx][cy][cz] = accpcls;
      accpcls += (*numpcls_in_bucket)[cx][cy][cz];
    }
    assert_eq(accpcls,nop);

    numpcls_in_bucket_now->setall(0);

    // put the particles where they are supposed to go
    for (int pidx = 0; pidx < nop; pidx++)
    {
      // get the cell indices of the particle
      //
      int cx,cy,cz;
      get_safe_cell_for_pos(cx,cy,cz,xavg[pidx],yavg[pidx],zavg[pidx]);
      //
      //cx = xcell[pidx];
      //cy = ycell[pidx];
      //cz = zcell[pidx];

      // compute where the data should go
      const int numpcls_now = (*numpcls_in_bucket_now)[cx][cy][cz]++;
      const int outpidx = (*bucket_offset)[cx][cy][cz] + numpcls_now;
      assert_lt(outpidx, nop);
      assert_ge(outpidx, 0);
      assert_lt(pidx, nop);
      assert_ge(pidx, 0);

      // copy particle data to new location
      //
      xtmp[outpidx] = x[pidx];
      ytmp[outpidx] = y[pidx];
      ztmp[outpidx] = z[pidx];
      utmp[outpidx] = u[pidx];
      vtmp[outpidx] = v[pidx];
      wtmp[outpidx] = w[pidx];
      qtmp[outpidx] = q[pidx];
      xavgtmp[outpidx] = xavg[pidx];
      yavgtmp[outpidx] = yavg[pidx];
      zavgtmp[outpidx] = zavg[pidx];
      if (TrackParticleID)
      {
        ParticleIDtmp[outpidx] = ParticleID[pidx];
      }
    }
    // swap the tmp particle memory with the official particle memory
    {
      swap(_xtmp,x);
      swap(_ytmp,y);
      swap(_ztmp,z);
      swap(_utmp,u);
      swap(_vtmp,v);
      swap(_wtmp,w);
      swap(_qtmp,q);
      swap(_ParticleIDtmp,ParticleID);
      swap(_xavgtmp,_xavg);
      swap(_yavgtmp,_yavg);
      swap(_zavgtmp,_zavg);
    }
    xavg = _xavg;
    yavg = _yavg;
    zavg = _zavg;

    // check that the number of bins was correct
    //
    if(true)
    {
      for(int cx=0;cx<nxc;cx++)
      for(int cy=0;cy<nyc;cy++)
      for(int cz=0;cz<nzc;cz++)
      {
        assert_eq((*numpcls_in_bucket_now)[cx][cy][cz], (*numpcls_in_bucket)[cx][cy][cz]);
      }
    }
    int serial_pidx=0;
    // confirm that the particles were sorted correctly
    for(int cx=0;cx<nxc;cx++)
    for(int cy=0;cy<nyc;cy++)
    for(int cz=0;cz<nzc;cz++)
    {
      const int numpcls_in_cell = get_numpcls_in_bucket(cx,cy,cz);
      const int bucket_offset = get_bucket_offset(cx,cy,cz);
      const int bucket_end = bucket_offset+numpcls_in_cell;
      for(int pidx=bucket_offset; pidx<bucket_end; pidx++)
      {
        // serial case: check that pidx is correct
        assert_eq(pidx,serial_pidx);
        serial_pidx++;
        // confirm that particle is in correct cell
        if(true)
        {
          int cx_,cy_,cz_;
          get_safe_cell_for_pos(cx_,cy_,cz_,xavg[pidx],yavg[pidx],zavg[pidx]);
          if((cx_!=cx)
           ||(cy_!=cy)
           ||(cz_!=cz))
          {
            dprint(cx)
            dprintf("\n\t cx =%d, cy =%d, cz =%d"
                    "\n\t cx_=%d, cy_=%d, cz_=%d"
                    "\n\t cxf=%f, cyf=%f, czf=%f",
                    cx,cy,cz,
                    cx_,cy_,cz_,
                    1.+(xavg[pidx]-xstart)*inv_dx,
                    1.+(yavg[pidx]-ystart)*inv_dy,
                    1.+(zavg[pidx]-zstart)*inv_dz);
            dprint(serial_pidx);
          }
          assert_eq(cx_,cx);
          assert_eq(cy_,cy);
          assert_eq(cz_,cz);
        }
      }
    }
  }
}

//void Particles3Dcomm::sort_particles_parallel(
//  double *xpos, double *ypos, double *zpos,
//  Grid * grid, VirtualTopology3D * vct)
//{
//  // should change this to first communicate particles so that
//  // they are in the correct process and all particles
//  // lie in this subdomain.
//
//  // count the number of particles to go in each bucket
//  numpcls_in_bucket.setall(0);
//  #pragma omp parallel
//  {
//    const int thread_num = omp_get_thread_num();
//    arr3_int numpcls_in_bucket_thr = fetch_numpcls_in_bucket_thr(thread_num);
//    numpcls_in_bucket_thr.setall(0);
//    // iterate through particles and count where they will go
//    #pragma omp for // nowait
//    for (int pidx = 0; pidx < nop; pidx++)
//    {
//      // get the cell indices of the particle
//      // (should change this to use xavg[pidx])
//      const pfloat xpos = xpos[pidx];
//      const pfloat ypos = ypos[pidx];
//      const pfloat zpos = zpos[pidx];
//      int cx,cy,cz;
//      get_safe_cell_for_pos(cx,cy,cz,xpos,ypos,zpos);
//
//      // need to allocate these
//      //
//      //xidx[pidx]=cx;
//      //yidx[pidx]=cy;
//      //zidx[pidx]=cz;
//
//      // increment the number of particles in bucket of this particle
//      numpcls_in_bucket_thr[cx][cy][cz]++;
//    }
//    // reduce the thread buckets into the main bucket
//    // #pragma omp critical (numpcls_in_bucket_reduction)
//    {
//      #pragma omp for collapse(2)
//      for(int cx=0;cx<nxc;cx++)
//      for(int cy=0;cy<nyc;cy++)
//      for(int th=0;th<num_threads;th++)
//      for(int cz=0;cz<nzc;cz++)
//      {
//        numpcls_in_bucket[cx][cy][cz]
//          += get_numpcls_in_bucket_thr(th)[cx][cy][cz];
//      }
//    }
//
//    // compute prefix sum to determine initial position
//    // of each bucket (could parallelize this)
//    //
//    int accpcls=0;
//    #pragma omp critical (bucket_offset_reduction)
//    for(int cx=0;cx<nxc;cx++)
//    for(int cy=0;cy<nyc;cy++)
//    for(int cz=0;cz<nzc;cz++)
//    {
//      bucket_offset[cx][cy][cz] = accpcls;
//      accpcls += numpcls_in_bucket[cx][cy][cz];
//    }
//
//    // cycle through the mesh cells mod 3
//    // (or mod(2*N+1), where N is number of mesh cells
//    // that a slow particle can move).
//    // This ensures that slow particles can be moved
//    // to their destinations without write conflicts
//    // among threads.  But what about cache contention?
//    //
//    for(int cxmod3=0; cxmod3<3; cxmod3++)
//    #pragma omp for collapse(2)
//    for(int cx=cxmod3; cx<nxc; cx+=3)
//    for(int cy=0; cy<nyc; cy++)
//    for(int cz=0; cz<nzc; cz++)
//    {
//      // put the slow particles where they are supposed to go and
//      // set aside the fast particles for separate processing.
//      // (to vectorize would need to sort separately in each
//      // dimension of space).
//      //
//      // problem: particles might have to be moved not because
//      // they are fast but because of an overall shift in the
//      // number of particles in a location, e.g. because of
//      // particles flowing in from a jet. Need a different
//      // approach, where memory is allocated for each cell.
//      _numpcls_in_bucket = numpcls_in_bucket[cx][cy][cz];
//      for(int pidx=bucket_offset[cx][cy][cz]; pidx<_numpcls_in_bucket; pidx++)
//      {
//        const int outcx = xidx[pidx];
//        const int outcy = yidx[pidx];
//        const int outcz = zidx[pidx];
//        const int cxlower = outcx <= 0 ? 0 : outcx-1;
//        const int cxupper = outcx >= (nxc-1) ? nxc-1 : outcx+1;
//        const int lowerindex = bucket_offset[cxlower][cylower][czlower];
//        const int upperoffset = bucket_offset[cxupper][cyupper][czupper];
//        const int upperindex = upperoffset + numpcls_in_bucket[outcx][outcy][outcz];
//        ...
//      }
//    }
//    // (1) put fast particles that must be moved more than one
//    // mesh cell at the end of the cell's list, and
//    // (2) put slow particles in the correct location
//
//    // count the number of particles that need to be moved
//    // more than one mesh cell and allocate a special buffer for them.
//    // (could change to count number of particles that need
//    // to move more than N mesh cells).
//    //
//    int numpcls_long_move_thr = 0;
//    #pragma omp for // nowait
//    for (int i = 0; i < nop; i++)
//    {
//      const int cx = xidx[pidx];
//      const int cy = yidx[pidx];
//      const int cz = zidx[pidx];
//
//      const int cxlower = cx <= 0 ? 0 : cx-1;
//      const int cxupper = cx >= (nxc-1) ? nxc-1 : cx+1;
//      const int lowerindex = bucket_offset[cxlower][cylower][czlower];
//      const int upperoffset = bucket_offset[cxupper][cyupper][czupper];
//      const int upperindex = upperoffset + numpcls_in_bucket[cx][cy][cz];
//      if(i < lowerindex || i > upperindex)
//      {
//        numpcls_long_move_thr++;
//      }
//    }
//  }
//}
//#endif

void Particles3Dcomm::copyParticlesToSoA()
{
  timeTasks_set_task(TimeTasks::TRANSPOSE_PCLS_TO_SOA);
  dprintf("copying to struct of arrays");
  SpeciesParticle const*const pcls = fetch_pcls();
  #pragma omp for
  for(int pidx=0; pidx<nop; pidx++)
  {
    const SpeciesParticle& pcl = pcls[pidx];
    if(ParticleID) ParticleID[pidx] = pcl.get_ID();
    x[pidx] = pcl.get_x(0);
    y[pidx] = pcl.get_x(1);
    z[pidx] = pcl.get_x(2);
    u[pidx] = pcl.get_u(0);
    v[pidx] = pcl.get_u(1);
    w[pidx] = pcl.get_u(2);
    q[pidx] = pcl.get_q();
  }
}

void Particles3Dcomm::copyParticlesToAoS()
{
  SpeciesParticle * pcls = fetch_pcls();
  timeTasks_set_task(TimeTasks::TRANSPOSE_PCLS_TO_AOS);
  dprintf("copying to array of structs");
  assert(pcls);
  if(TrackParticleID)
  {
    assert(ParticleID);
    #pragma omp for
    for(int pidx=0; pidx<nop; pidx++)
    {
      pcls[pidx].set( ParticleID[pidx],
        x[pidx],y[pidx],z[pidx],
        u[pidx],v[pidx],w[pidx],
        q[pidx]);
    }
  }
  else
  {
    #pragma omp for
    for(int pidx=0; pidx<nop; pidx++)
    {
      pcls[pidx].set( 0,
        x[pidx],y[pidx],z[pidx],
        u[pidx],v[pidx],w[pidx],
        q[pidx]);
    }
  }
}

void Particles3Dcomm::convertParticlesToAoS()
{
  if(particleType!=ParticleType::AoS)
  {
    assert_eq(particleType,ParticleType::SoA);
    copyParticlesToAoS();
    particleType = ParticleType::AoS;
  }
}

void Particles3Dcomm::convertParticlesToSoA()
{
  if(particleType != ParticleType::SoA)
  {
    assert_eq(particleType,ParticleType::AoS);
    copyParticlesToSoA();
    particleType = ParticleType::SoA;
  }
}

