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

  if (TrackSpecies){
    delete[]TrackSpBirthRank;
    delete[]TrackSpID;
  }
}
/** constructors fo a single species*/
void Particles3Dcomm::allocate(int species, long long initnpmax, Collective * col, VirtualTopology3D * vct, Grid * grid) {
  
  /* initialize random generator with different seed on different processor */
  seed = vct->getCartesian_rank()+ 2;  

  // info from collectiveIO
  ns = species;
  npcel = col->getNpcel(species);
  npcelx = col->getNpcelx(species);
  npcely = col->getNpcely(species);
  npcelz = col->getNpcelz(species);


  // these are used if you are restarting a non EB simulation from a EB simulation       
  EBRestart_RdivR0= col->getEBRestart_RdivR0();
  EBRestart_RelocPart= col->getEBRestart_RelocPart();

  // This if is necessary to restart with H5hut-io
  if (initnpmax==0){
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

  /* initial parameters for whistler */
  omega_r= col->getOmega_r();
  deltaB= col->getDeltaB();
  // I need the number of electron species to distribute current
  ElectronSpNumber=0;
  for (int ss=0; ss< col->getNs(); ss++)
    if (col->getQOM(ss)<0)
      ElectronSpNumber++;
  /* end initial parameters for whistler */




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
  npTracked = 0;
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
    ParticleID = new long long[npmax];
    partRank = new long long[npmax];
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

  /* NB: this tracking is different, because data not saved in hdf;
     saved in txt file; 
     it is expected that fewer particles will be tracked this way
     in an apposite species */
  TrackSpecies= col->getTrackSpecies(ns);
  TrackingSpOutputCycle= col->getTrackingOutputCycle();
  //cout <<"sp " << ns << " TrackSpecies: "<< TrackSpecies << ", TrackingSpOutputCycle: " << TrackingSpOutputCycle << endl;
  if (TrackSpecies){  
    TrackSpBirthRank = new unsigned long[npmax];
    TrackSpID        = new long long[npmax];
  }
  // BUFFERS
  // the buffer size should be decided depending on number of particles
  // the buffer size should be decided depending on number of particles
  if (TrackParticleID){
    nVar = 8;
  }
  else{
    nVar = 7;
  }
  
  if (TrackSpecies){
    nVar+=2;
    startTr= nVar-2;
  }

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

    // it checks internally if needed ()
    Do_EBRestart_RelocPart(vct);
  }

  // I copy these values here, even if the box is not expanding, for the tests
  UEB= col->getUEB_0(); 
  REB= col->getREB_0();  

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

#ifdef EB

	    /* the RESCALED moments have to be used in the big equation -
	       tested with electrostatic perturbation */
	    double R0= EMf->getREB_0_EB();
	    // this value is updated in UpdateEBVectors
	    double R= EMf->getR_EB();
	    /* rescaling the moments */
            weight[ii][jj][kk] = q[i] * xi[ii] * eta[jj] * zeta[kk] * invVOL/ (R*R/ R0/ R0);
	    if (ns ==0 && i==0 && vct->getCartesian_rank() ==0){
	      /*cout << "Species " << ns << " is rescaling moments (R/R0)^2 for EB, R= " << R << ", R0= "<< R0<<endl;;*/
	    }
#else
	    weight[ii][jj][kk] = q[i] * xi[ii] * eta[jj] * zeta[kk] * invVOL;
#endif
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

      // add energy flux density - X
      // careful: heat flux and pressure are defined in a non consistent way with regard to the 1/qom; for the pressure, it has to be multiplied in post-processing; heat flux is already OK; this implementation of HF is consistent with the one in the branch Coil
      for (int ii = 0; ii < 2; ii++)
	for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
	    temp[ii][jj][kk] = u[i] * 0.5/qom  * (u[i]*u[i] +v[i]*v[i]+w[i]*w[i]) * weight[ii][jj][kk];
      EMf->addEFx(temp, ix, iy, iz, ns);

      // add energy flux density - Y
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
	    temp[ii][jj][kk] = v[i] * 0.5/ qom * (u[i]*u[i] +v[i]*v[i]+w[i]*w[i]) * weight[ii][jj][kk];
      EMf->addEFy(temp, ix, iy, iz, ns);

      // add energy flux density - Z
      for (int ii = 0; ii < 2; ii++)
        for (int jj = 0; jj < 2; jj++)
          for (int kk = 0; kk < 2; kk++)
	    temp[ii][jj][kk] = w[i] * 0.5/qom * (u[i]*u[i] +v[i]*v[i]+w[i]*w[i]) * weight[ii][jj][kk];
      EMf->addEFz(temp, ix, iy, iz, ns);

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
  if (TrackParticleID){
    b_[npExitXleft * nVar + 7] = ParticleID[np_current];
    b_[npExitXleft * nVar + 8] = partRank[np_current];
  }
  if (TrackSpecies){
    b_[npExitXleft * nVar + startTr +0] = TrackSpBirthRank[np_current]; 
    b_[npExitXleft * nVar + startTr +1] = TrackSpID[np_current];
  }
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
  if (TrackParticleID){
    b_[npExitXright * nVar + 7] = ParticleID[np_current];
    b_[npExitXright * nVar + 8] = partRank[np_current];
  }
  if (TrackSpecies){
    b_[npExitXright * nVar + startTr +0] = TrackSpBirthRank[np_current];
    b_[npExitXright * nVar + startTr +1] = TrackSpID[np_current];
  }
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
  if (TrackParticleID){
    b_[npExitYleft * nVar + 7] = ParticleID[np_current];
    b_[npExitYleft * nVar + 8] = partRank[np_current];
  }
  if (TrackSpecies){
    b_[npExitYleft * nVar + startTr +0] = TrackSpBirthRank[np_current];
    b_[npExitYleft * nVar + startTr +1] = TrackSpID[np_current];
  }
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
  if (TrackParticleID){
    b_[npExitYright * nVar + 7] = ParticleID[np_current];
    b_[npExitYright * nVar + 8] = partRank[np_current];
  }
  if (TrackSpecies){
    b_[npExitYright * nVar + startTr +0] = TrackSpBirthRank[np_current];
    b_[npExitYright * nVar + startTr +1] = TrackSpID[np_current];
  }
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
  if (TrackParticleID){
    b_[npExitZleft * nVar + 7] = ParticleID[np_current];
    b_[npExitZleft * nVar + 8] = partRank[np_current];
  }
  if (TrackSpecies){
    b_[npExitZleft * nVar + startTr +0] = TrackSpBirthRank[np_current];
    b_[npExitZleft * nVar + startTr +1] = TrackSpID[np_current];
  }
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
  if (TrackParticleID){
    b_[npExitZright * nVar + 7] = ParticleID[np_current];
    b_[npExitZright * nVar + 8] = partRank[np_current];
  }
  if (TrackSpecies){
    b_[npExitZright * nVar + startTr +0] = TrackSpBirthRank[np_current];
    b_[npExitZright * nVar + startTr +1] = TrackSpID[np_current];
  }
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
    if (TrackParticleID){
      ParticleID[nop] = (unsigned long) b_[nVar * np_current + 7];
      partRank[nop] = (long) b_[nVar * np_current +8];
    }
    if (TrackSpecies){
      TrackSpBirthRank[nop] = (unsigned long) b_[nVar * np_current + startTr +0];
      TrackSpID[nop]        = (long long) b_[nVar * np_current + startTr +1];
    }
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
  if (TrackParticleID){
    ParticleID[np_current] = ParticleID[*nplast];
    partRank[np_current] = partRank[*nplast];
  }
  if (TrackSpecies){
    TrackSpBirthRank[np_current] = TrackSpBirthRank[*nplast];
    TrackSpID[np_current]        = TrackSpID[*nplast];
  }
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
long long *Particles3Dcomm::getParticleIDall()  const {
  return (ParticleID);
}
/**get proc of particle with label indexPart */
long long *Particles3Dcomm::getParticleRankall()  const {
  return (partRank);
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
long long Particles3Dcomm::getParticleID(long long indexPart)  const {
  return (ParticleID[indexPart]);
}
/**get rank ID of particle with label indexPart */
long long Particles3Dcomm::getParticleRank(long long indexPart)  const {
  return (partRank[indexPart]);
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
bool Particles3Dcomm::GetTrackSpecies(){
  return TrackSpecies;
}
int Particles3Dcomm::GetTrackingSpOutputCycle(){
  return TrackingSpOutputCycle;
}
void Particles3Dcomm::AssignParticlesID(VirtualTopology3D * vct){
  for (int i=0; i< nop; i++){
    TrackSpBirthRank[i] = vct->getCartesian_rank();
    TrackSpID[i]        = i;
  }
}

void Particles3Dcomm::WriteTracking(int cycle, VirtualTopology3D * vct, Collective * col, Field * EMf, Grid * grid){

  stringstream num_proc;
  stringstream num_sp;
  stringstream num_cyc;

  num_proc << vct->getCartesian_rank();
  num_sp << ns;
  num_cyc << cycle;

  cqTr = col->getSaveDirName() + "/Tracking_proc" + num_proc.str() + "_sp" + num_sp.str() + "_cyc" + num_cyc.str()+ ".txt";

  // I want to start writing also Ex, Ey, Ez, Bx, By, Bz @ particle position, so I have to do another interpolation
  const double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy, inv_dz = 1.0 / dz;
  
  double ***Ex = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEx());
  double ***Ey = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEy());
  double ***Ez = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEz());
  double ***Bx = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx());
  double ***By = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy());
  double ***Bz = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz());

  double ***Bx_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx_ext());
  double ***By_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy_ext());
  double ***Bz_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz_ext());

  double Fext = EMf->getFext();


  ofstream my_file(cqTr.c_str(), fstream::app);
  //my_file << cycle << endl;
  for (int i=0; i< nop; i++){

    /* interpolate fields at particle position*/
    // interpolation G-->P                                                                                            
    const double ixd = floor((x[i] - xstart) * inv_dx);
    const double iyd = floor((y[i] - ystart) * inv_dy);
    const double izd = floor((z[i] - zstart) * inv_dz);
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

    xi  [0] = x[i] - grid->getXN(ix-1,iy  ,iz  );
    eta [0] = y[i] - grid->getYN(ix  ,iy-1,iz  );
    zeta[0] = z[i] - grid->getZN(ix  ,iy  ,iz-1);
    xi  [1] = grid->getXN(ix,iy,iz) - x[i];
    eta [1] = grid->getYN(ix,iy,iz) - y[i];
    zeta[1] = grid->getZN(ix,iy,iz) - z[i];

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
    //                                                                                                                
    Ezl += weight000 * Ez[ix][iy][iz];
    Ezl += weight001 * Ez[ix][iy][iz - 1];
    Ezl += weight010 * Ez[ix][iy - 1][iz];
    Ezl += weight011 * Ez[ix][iy - 1][iz - 1];
    Ezl += weight100 * Ez[ix - 1][iy][iz];
    Ezl += weight101 * Ez[ix - 1][iy][iz - 1];
    Ezl += weight110 * Ez[ix - 1][iy - 1][iz];
    Ezl += weight111 * Ez[ix - 1][iy - 1][iz - 1];
    /* end interpolate fields at particle position;*/

    my_file << TrackSpBirthRank[i] <<" " << TrackSpID[i] <<" " << x[i]<<" " << y[i]<<" " << z[i]<<" " << q[i]<<" " << u[i]<<" " << v[i]<<" " << w[i] << " " << Exl << " " << Eyl << " " << Ezl << " " << Bxl << " " << Byl << " " << Bzl << endl;
  }
  cout << "Proc " << num_proc.str() << " wrote " << nop << " particles, sp " << num_sp << endl; 
  my_file.close();
}

/** this is used if you are restarting a non EB simulation from an EB simulation                                        
   the particle position in the transverse direction is scaled to keep into account volume transverse expansion;
   in the inputfile, scale L_trans= L_0 R/R_0**/
void Particles3Dcomm::Do_EBRestart_RelocPart(VirtualTopology3D * vct){
  // NB: this is done within a "if restart", so that has not to be checked

  if (!EBRestart_RelocPart)
    return;

  if (vct->getCartesian_rank() == 0 && ns == 0){
    cout << "RESTART: I AM RELOCATING PARTICLES TO ACCOUNT FOR VOLUME EXPANSION!!!!" << endl;
    cout << "(if you are not restarting a non EB sim, from an EB run, you should not see this)" << endl;
    cout << "Lx: " << Lx << ", Ly: "<< Ly << ", Lz: " << Lz << endl;
  }
  for (int i=0; i< nop; i++){
    y[i]= y[i]*EBRestart_RdivR0;
    z[i]= z[i]*EBRestart_RdivR0;
  }

}

/** return  <v^2*vx> */
double Particles3Dcomm::getqx() {
  double localqx = 0.0;
  double totalqx = 0.0;
  for (register long long i = 0; i < nop; i++)
    localqx += 0.5*(u[i]*u[i] + v[i]*v[i] + w[i]*w[i])*u[i]*q[i]/qom;
  MPI_Allreduce(&localqx, &totalqx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalqx);
}

/** return  <v^2*vy> */
double Particles3Dcomm::getqy() {
  double localqy = 0.0;
  double totalqy = 0.0;
  for (register long long i = 0; i < nop; i++)
    localqy += 0.5*(u[i]*u[i] + v[i]*v[i] + w[i]*w[i])*v[i]*q[i]/qom;
  MPI_Allreduce(&localqy, &totalqy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalqy);
}

/** return  <v^2*vz> */
double Particles3Dcomm::getqz() {
  double localqz = 0.0;
  double totalqz = 0.0;
  for (register long long i = 0; i < nop; i++)
    localqz += 0.5*(u[i]*u[i] + v[i]*v[i] + w[i]*w[i])*w[i]*q[i]/qom;
  MPI_Allreduce(&localqz, &totalqz, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (totalqz);
}
