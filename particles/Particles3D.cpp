/*******************************************************************************************
  Particles3D.cpp  -  Class for particles of the same species, in a 3D space and 3component velocity
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/

#include <mpi.h>
#include <iostream>
#include <math.h>
#include "../processtopology/VirtualTopology3D.h"
#include "../processtopology/VCtopology3D.h"
#include "../inputoutput/CollectiveIO.h"
#include "../inputoutput/Collective.h"
#include "../mathlib/Basic.h"
#include "../bc/BcParticles.h"
#include "../grids/Grid.h"
#include "../grids/Grid3DCU.h"
#include "../fields/Field.h"

#include "Particles3D.h"
#include "../utility/TimeTasks.h"


#include "hdf5.h"
#include <complex>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define MIN_VAL   1E-16
// particles processed together
#define P_SAME_TIME 2

/**
 * 
 * Class for particles of the same species
 * @date Fri Jun 4 2009
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */

/** constructor */
Particles3D::Particles3D() {
  // see allocate(int species, CollectiveIO* col, VirtualTopology3D* vct, Grid* grid)

}
/** deallocate particles */
Particles3D::~Particles3D() {
  delete[]x;
  delete[]y;
  delete[]z;
  delete[]u;
  delete[]v;
  delete[]w;
  delete[]q;
}

/** particles are uniformly distributed with zero velocity   */
void Particles3D::uniform_background(Grid * grid, Field * EMf) {
  long long counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
              u[counter] = 0.0;
              v[counter] = 0.0;
              w[counter] = 0.0;
              q[counter] = (qom / fabs(qom)) * (EMf->getRHOcs(i, j, k, ns) / npcel) * (1.0 / grid->getInvVOL());
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];
              counter++;
            }


  cout << "Velocity Maxwellian Distribution " << endl;
}
/** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
  <ul>
  <li> dim = 0 --> constant velocity on X direction </li>
  <li> dim = 1 --> constant velocity on Y direction </li>
  <li> dim = 2 --> constant velocity on Z direction </li>
  </ul>

*/
void Particles3D::constantVelocity(double vel, int dim, Grid * grid, Field * EMf) {
  switch (dim) {
    case 0:
      for (long long i = 0; i < nop; i++)
        u[i] = vel, v[i] = 0.0, w[i] = 0.0;
      break;
    case 1:
      for (register long long i = 0; i < nop; i++)
        u[i] = 0.0, v[i] = vel, w[i] = 0.0;
      break;
    case 2:
      for (register long long i = 0; i < nop; i++)
        u[i] = 0.0, v[i] = 0.0, w[i] = vel;
      break;

  }

}

/** alternative routine maxellian random velocity and uniform spatial distribution */
void Particles3D::alt_maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct) {


}

/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);

  double harvest;
  double prob, theta, sign;
	long long counter=0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);  // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
              // q = charge
              q[counter] = (qom / fabs(qom)) * (EMf->getRHOcs(i, j, k, ns) / npcel) * (1.0 / grid->getInvVOL());
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              u[counter] = u0 + uth * prob * cos(theta);
              // v
              v[counter] = v0 + vth * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              w[counter] = w0 + wth * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];


              counter++;
            }


}

/** Force Free initialization (JxB=0) for particles */
void Particles3D::force_free(Grid * grid, Field * EMf, VirtualTopology3D * vct) {


  double harvest, prob, theta;
  long long counter = 0;
  double shaperx, shapery, shaperz;
  double flvx = 1.0, flvy = 1.0, flvz = 1.0;


  /* initialize random generator */
  srand(vct->getCartesian_rank() + 1 + ns);
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              flvx = 1.0;
              flvy = 1.0;
              flvz = 1.0;
              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
              // q = charge
              q[counter] = (qom / fabs(qom)) * (EMf->getRHOcs(i, j, k, ns) / npcel) * (1.0 / invVOL);
              shaperx = tanh((y[counter] - Ly / 2) / delta) / cosh((y[counter] - Ly / 2) / delta) / delta;
              shaperz = 1.0 / (cosh((y[counter] - Ly / 2) / delta) * cosh((y[counter] - Ly / 2) / delta)) / delta;
              shapery = shapery;
              // new drift velocity to satisfy JxB=0
              flvx = u0 * flvx * shaperx;
              flvz = w0 * flvz * shaperz;
              flvy = v0 * flvy * shapery;
              u[counter] = c;
              v[counter] = c;
              w[counter] = c;
              while ((fabs(u[counter]) >= c) | (fabs(v[counter]) >= c) | (fabs(w[counter]) >= c)) {
                harvest = rand() / (double) RAND_MAX;
                prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
                harvest = rand() / (double) RAND_MAX;
                theta = 2.0 * M_PI * harvest;
                u[counter] = flvx + uth * prob * cos(theta);
                // v
                v[counter] = flvy + vth * prob * sin(theta);
                // w
                harvest = rand() / (double) RAND_MAX;
                prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
                harvest = rand() / (double) RAND_MAX;
                theta = 2.0 * M_PI * harvest;
                w[counter] = flvz + wth * prob * cos(theta);
              }
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

              counter++;
            }

}

/**Add a periodic perturbation in J exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
inline void Particles3D::AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0, Grid * grid) {

  // rescaling of amplitudes according to deltaBoB //
  double alpha;
  alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);
  jx_mod *= alpha;
  jy_mod *= alpha;
  jz_mod *= alpha;
  for (register long long i = 0; i < nop; i++) {
    u[i] += jx_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jx_phase);
    v[i] += jy_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jy_phase);
    w[i] += jz_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jz_phase);
  }
}


/** explicit mover */
void Particles3D::mover_explicit(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  // to be implemented

}
/** mover with a Predictor-Corrector scheme */
int Particles3D::mover_PC(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  double start_mover_PC = MPI_Wtime();
  doubleArr3 Ex(EMf->getEx());
  doubleArr3 Ey(EMf->getEy());
  doubleArr3 Ez(EMf->getEz());
  doubleArr3 Bx(EMf->getBx());
  doubleArr3 By(EMf->getBy());
  doubleArr3 Bz(EMf->getBz());
  doubleArr4 node_coordinate(grid->getN());
  const double dto2 = .5 * dt, qomdt2 = qom * dto2 / c;
  const double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy, inv_dz = 1.0 / dz;
  // don't bother trying to push any particles simultaneously;
  // MIC already does vectorization automatically, and trying
  // to do it by hand only hurts performance.
  #pragma omp parallel for
  #pragma simd // this just slows things down (why?)
  for (long long rest = 0; rest < nop; rest++) {
    // copy the particle
    double xp = x[rest];
    double yp = y[rest];
    double zp = z[rest];
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];
    const double xptilde = x[rest];
    const double yptilde = y[rest];
    const double zptilde = z[rest];
    double uptilde;
    double vptilde;
    double wptilde;
    // calculate the average velocity iteratively
    for (int innter = 0; innter < 1; innter++) {
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

      double xi[2];
      double eta[2];
      double zeta[2];
      xi[0]   = xp - node_coordinate.get(ix - 1, iy, iz, 0);
      eta[0]  = yp - node_coordinate.get(ix, iy - 1, iz, 1);
      zeta[0] = zp - node_coordinate.get(ix, iy, iz - 1, 2);
      xi[1]   = node_coordinate.get(ix,iy,iz,0) - xp;
      eta[1]  = node_coordinate.get(ix,iy,iz,1) - yp;
      zeta[1] = node_coordinate.get(ix,iy,iz,2) - zp;

      double Exl = 0.0; double Eyl = 0.0; double Ezl = 0.0; double Bxl = 0.0; double Byl = 0.0; double Bzl = 0.0;

      // MIC refuses to vectorize this ...
      //
      //double weight[2][2][2];
      //for (int ii = 0; ii < 2; ii++)
      //  for (int jj = 0; jj < 2; jj++)
      //    for (int kk = 0; kk < 2; kk++)
      //      weight[ii][jj][kk] = xi[ii] * eta[jj] * zeta[kk] * invVOL;
      //for (int ii = 0; ii < 2; ii++)
      //  for (int jj = 0; jj < 2; jj++)
      //    for (int kk = 0; kk < 2; kk++) {
      //      const double Exlp = weight[ii][jj][kk] * Ex.get(ix - ii, iy - jj, iz - kk);
      //      const double Eylp = weight[ii][jj][kk] * Ey.get(ix - ii, iy - jj, iz - kk);
      //      const double Ezlp = weight[ii][jj][kk] * Ez.get(ix - ii, iy - jj, iz - kk);
      //      const double Bxlp = weight[ii][jj][kk] * Bx.get(ix - ii, iy - jj, iz - kk);
      //      const double Bylp = weight[ii][jj][kk] * By.get(ix - ii, iy - jj, iz - kk);
      //      const double Bzlp = weight[ii][jj][kk] * Bz.get(ix - ii, iy - jj, iz - kk);
      //      Exl += Exlp;
      //      Eyl += Eylp;
      //      Ezl += Ezlp;
      //      Bxl += Bxlp;
      //      Byl += Bylp;
      //      Bzl += Bzlp;
      //    }

      // ... so we expand things out instead
      //
      const double weight000 = xi[0] * eta[0] * zeta[0] * invVOL;
      const double weight001 = xi[0] * eta[0] * zeta[1] * invVOL;
      const double weight010 = xi[0] * eta[1] * zeta[0] * invVOL;
      const double weight011 = xi[0] * eta[1] * zeta[1] * invVOL;
      const double weight100 = xi[1] * eta[0] * zeta[0] * invVOL;
      const double weight101 = xi[1] * eta[0] * zeta[1] * invVOL;
      const double weight110 = xi[1] * eta[1] * zeta[0] * invVOL;
      const double weight111 = xi[1] * eta[1] * zeta[1] * invVOL;
      //
      Bxl += weight000 * Bx.get(ix  , iy  , iz  );
      Bxl += weight001 * Bx.get(ix  , iy  , iz-1);
      Bxl += weight010 * Bx.get(ix  , iy-1, iz  );
      Bxl += weight011 * Bx.get(ix  , iy-1, iz-1);
      Bxl += weight100 * Bx.get(ix-1, iy  , iz  );
      Bxl += weight101 * Bx.get(ix-1, iy  , iz-1);
      Bxl += weight110 * Bx.get(ix-1, iy-1, iz  );
      Bxl += weight111 * Bx.get(ix-1, iy-1, iz-1);
      //
      Byl += weight000 * By.get(ix  , iy  , iz  );
      Byl += weight001 * By.get(ix  , iy  , iz-1);
      Byl += weight010 * By.get(ix  , iy-1, iz  );
      Byl += weight011 * By.get(ix  , iy-1, iz-1);
      Byl += weight100 * By.get(ix-1, iy  , iz  );
      Byl += weight101 * By.get(ix-1, iy  , iz-1);
      Byl += weight110 * By.get(ix-1, iy-1, iz  );
      Byl += weight111 * By.get(ix-1, iy-1, iz-1);
      //
      Bzl += weight000 * Bz.get(ix  , iy  , iz  );
      Bzl += weight001 * Bz.get(ix  , iy  , iz-1);
      Bzl += weight010 * Bz.get(ix  , iy-1, iz  );
      Bzl += weight011 * Bz.get(ix  , iy-1, iz-1);
      Bzl += weight100 * Bz.get(ix-1, iy  , iz  );
      Bzl += weight101 * Bz.get(ix-1, iy  , iz-1);
      Bzl += weight110 * Bz.get(ix-1, iy-1, iz  );
      Bzl += weight111 * Bz.get(ix-1, iy-1, iz-1);
      //
      Exl += weight000 * Ex.get(ix  , iy  , iz  );
      Exl += weight001 * Ex.get(ix  , iy  , iz-1);
      Exl += weight010 * Ex.get(ix  , iy-1, iz  );
      Exl += weight011 * Ex.get(ix  , iy-1, iz-1);
      Exl += weight100 * Ex.get(ix-1, iy  , iz  );
      Exl += weight101 * Ex.get(ix-1, iy  , iz-1);
      Exl += weight110 * Ex.get(ix-1, iy-1, iz  );
      Exl += weight111 * Ex.get(ix-1, iy-1, iz-1);
      //
      Eyl += weight000 * Ey.get(ix  , iy  , iz  );
      Eyl += weight001 * Ey.get(ix  , iy  , iz-1);
      Eyl += weight010 * Ey.get(ix  , iy-1, iz  );
      Eyl += weight011 * Ey.get(ix  , iy-1, iz-1);
      Eyl += weight100 * Ey.get(ix-1, iy  , iz  );
      Eyl += weight101 * Ey.get(ix-1, iy  , iz-1);
      Eyl += weight110 * Ey.get(ix-1, iy-1, iz  );
      Eyl += weight111 * Ey.get(ix-1, iy-1, iz-1);
      //
      Ezl += weight000 * Ez.get(ix  , iy  , iz  );
      Ezl += weight001 * Ez.get(ix  , iy  , iz-1);
      Ezl += weight010 * Ez.get(ix  , iy-1, iz  );
      Ezl += weight011 * Ez.get(ix  , iy-1, iz-1);
      Ezl += weight100 * Ez.get(ix-1, iy  , iz  );
      Ezl += weight101 * Ez.get(ix-1, iy  , iz-1);
      Ezl += weight110 * Ez.get(ix-1, iy-1, iz  );
      Ezl += weight111 * Ez.get(ix-1, iy-1, iz-1);

      // end interpolation
      const double omdtsq = qomdt2 * qomdt2 * (Bxl * Bxl + Byl * Byl + Bzl * Bzl);
      const double denom = 1.0 / (1.0 + omdtsq);
      // solve the position equation
      const double ut = up + qomdt2 * Exl;
      const double vt = vp + qomdt2 * Eyl;
      const double wt = wp + qomdt2 * Ezl;
      const double udotb = ut * Bxl + vt * Byl + wt * Bzl;
      // solve the velocity equation 
      uptilde = (ut + qomdt2 * (vt * Bzl - wt * Byl + qomdt2 * udotb * Bxl)) * denom;
      vptilde = (vt + qomdt2 * (wt * Bxl - ut * Bzl + qomdt2 * udotb * Byl)) * denom;
      wptilde = (wt + qomdt2 * (ut * Byl - vt * Bxl + qomdt2 * udotb * Bzl)) * denom;
      // update position
      xp = xptilde + uptilde * dto2;
      yp = yptilde + vptilde * dto2;
      zp = zptilde + wptilde * dto2;
    }                           // end of iteration
    // update the final position and velocity
    up = 2.0 * uptilde - u[rest];
    vp = 2.0 * vptilde - v[rest];
    wp = 2.0 * wptilde - w[rest];
    xp = xptilde + uptilde * dt;
    yp = yptilde + vptilde * dt;
    zp = zptilde + wptilde * dt;
    x[rest] = xp;
    y[rest] = yp;
    z[rest] = zp;
    u[rest] = up;
    v[rest] = vp;
    w[rest] = wp;
  }                             // END OF ALL THE PARTICLES

  // ********************//
  // COMMUNICATION 
  // *******************//
  timeTasks.start_communicate();
  const int avail = communicate(vct);
  if (avail < 0)
    return (-1);
  MPI_Barrier(MPI_COMM_WORLD);
  // communicate again if particles are not in the correct domain
  while (isMessagingDone(vct) > 0) {
    // COMMUNICATION
    const int avail = communicate(vct);
    if (avail < 0)
      return (-1);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  timeTasks.addto_communicate();
  return (0);                   // exit succcesfully (hopefully) 
}


/** relativistic mover with a Predictor-Corrector scheme */
int Particles3D::mover_relativistic(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  return (0);
}


/** interpolation Particle->Grid only for pressure tensor */
void Particles3D::interpP2G_onlyP(Field * EMf, Grid * grid, VirtualTopology3D * vct) {
  double weight[2][2][2];
  double temp[2][2][2];
  int ix, iy, iz, temp1, temp2, temp3;
  for (register long long i = 0; i < nop; i++) {
    ix = 2 + int (floor((x[i] - grid->getXstart()) / grid->getDX()));
    iy = 2 + int (floor((y[i] - grid->getYstart()) / grid->getDY()));
    iz = 2 + int (floor((z[i] - grid->getZstart()) / grid->getDZ()));
    calculateWeights(weight, x[i], y[i], z[i], ix, iy, iz, grid);
    scale(weight, q[i], 2, 2, 2);
    // Pxx
    eqValue(0.0, temp, 2, 2, 2);
    addscale(u[i] * u[i], temp, weight, 2, 2, 2);
    EMf->addPxx(temp, ix, iy, iz, ns);
    // Pxy
    eqValue(0.0, temp, 2, 2, 2);
    addscale(u[i] * v[i], temp, weight, 2, 2, 2);
    EMf->addPxy(temp, ix, iy, iz, ns);
    // Pxz
    eqValue(0.0, temp, 2, 2, 2);
    addscale(u[i] * w[i], temp, weight, 2, 2, 2);
    EMf->addPxz(temp, ix, iy, iz, ns);
    // Pyy
    eqValue(0.0, temp, 2, 2, 2);
    addscale(v[i] * v[i], temp, weight, 2, 2, 2);
    EMf->addPyy(temp, ix, iy, iz, ns);
    // Pyz
    eqValue(0.0, temp, 2, 2, 2);
    addscale(v[i] * w[i], temp, weight, 2, 2, 2);
    EMf->addPyz(temp, ix, iy, iz, ns);
    // Pzz
    eqValue(0.0, temp, 2, 2, 2);
    addscale(w[i] * w[i], temp, weight, 2, 2, 2);
    EMf->addPzz(temp, ix, iy, iz, ns);
  }
}
/** interpolation Particle->Grid only charge density, current */
void Particles3D::interpP2G_notP(Field * EMf, Grid * grid, VirtualTopology3D * vct) {
  double weight[2][2][2];
  double temp[2][2][2];
  int ix, iy, iz, temp2, temp1, temp3;
  for (register long long i = 0; i < nop; i++) {
    ix = 2 + int (floor((x[i] - grid->getXstart()) / grid->getDX()));
    iy = 2 + int (floor((y[i] - grid->getYstart()) / grid->getDY()));
    iz = 2 + int (floor((z[i] - grid->getZstart()) / grid->getDZ()));
    temp1 = (int) min(ix, nxn - 2);
    temp2 = (int) min(iy, nyn - 2);
    temp3 = (int) min(iz, nzn - 2);
    ix = (int) max(temp1, 2);
    iy = (int) max(temp2, 2);
    iz = (int) max(temp3, 2);
    calculateWeights(weight, x[i], y[i], z[i], ix, iy, iz, grid);
    scale(weight, q[i], 2, 2, 2);
    // rho
    EMf->addRho(weight, ix, iy, iz, ns);
    // Jx
    eqValue(0.0, temp, 2, 2, 2);
    addscale(u[i], temp, weight, 2, 2, 2);
    EMf->addJx(temp, ix, iy, iz, ns);
    // Jy
    eqValue(0.0, temp, 2, 2, 2);
    addscale(v[i], temp, weight, 2, 2, 2);
    EMf->addJy(temp, ix, iy, iz, ns);
    // Jz
    eqValue(0.0, temp, 2, 2, 2);
    addscale(w[i], temp, weight, 2, 2, 2);
    EMf->addJz(temp, ix, iy, iz, ns);

  }
  // communicate contribution from ghost cells 
  EMf->communicateGhostP2G(ns, 0, 0, 0, 0, vct);
}
/** apply a linear perturbation to particle distribution */
void Particles3D::linear_perturbation(double deltaBoB, double kx, double ky, double angle, double omega_r, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  double value1 = 0.0, value2 = 0.0, max_value = 0.0, min_value = 0.0, phi, n;
  long long counter = 0, total_generated = 0;
  bool rejected;
  double harvest, prob, theta;
  // rescaling of amplitudes according to deltaBoB //
  double alpha;
  double integral = 0.0;

  alpha = deltaBoB * sqrt(EMf->getBx(1, 1, 0) * EMf->getBx(1, 1, 0) + EMf->getBy(1, 1, 0) * EMf->getBy(1, 1, 0) + EMf->getBz(1, 1, 0) * EMf->getBz(1, 1, 0)) / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

  Ex_mod *= alpha;
  Ey_mod *= alpha;
  Ez_mod *= alpha;
  Bx_mod *= alpha;
  By_mod *= alpha;
  Bz_mod *= alpha;



  // find the maximum value of f=1+delta_f/f0
  for (register double vpar = -2 * uth; vpar <= 2 * uth; vpar += 0.0005)
    for (register double vperp = 1e-10; vperp <= 2 * vth; vperp += 0.0005)
      for (register double X = xstart; X <= xend; X += 2 * grid->getDX())
        for (register double Y = ystart; Y <= yend; Y += 2 * grid->getDY()) {
          value1 = 1 + delta_f(vpar, vperp, 0.0, X, Y, kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle, EMf) / f0(vpar, vperp);

          if (value1 > max_value)
            max_value = value1;


        }



  max_value *= 3.2;
  phi = 1.48409;
  n = 2.948687;                 // security factor...
  if (ns == 1) {
    max_value *= 3.0;
    phi = -1.65858;
    n = 2.917946;
  }                             // security factor...
  cout << "max-value=" << max_value << " min-value=" << min_value << endl;

  /* initialize random generator */
  srand(vct->getCartesian_rank() + 2);

  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int ii = 0; ii < npcelx + (int) (2 * n * (cos(2 * M_PI * 0.4125 * grid->getXN(i, j, 0) + phi))); ii++)
        for (int jj = 0; jj < npcely; jj++) {
          x[counter] = (ii + .5) * (dx / (npcelx + (int) (2 * n * (cos(2 * M_PI * 0.4125 * grid->getXN(i, j, 0) + phi))))) + grid->getXN(i, j, 0);
          y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, 0);
          q[counter] = (qom / fabs(qom)) * ((0.19635) / npcel) * (1.0 / invVOL);

          // apply rejection method in velocity space
          rejected = true;
          while (rejected) {
            total_generated++;
            harvest = rand() / (double) RAND_MAX;
            prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
            harvest = rand() / (double) RAND_MAX;
            theta = 2.0 * M_PI * harvest;
            // u
            u[counter] = u0 + uth * prob * cos(theta);
            // v
            v[counter] = v0 + vth * prob * sin(theta);
            // w
            harvest = rand() / (double) RAND_MAX;
            prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
            harvest = rand() / (double) RAND_MAX;
            theta = 2.0 * M_PI * harvest;
            w[counter] = w0 + wth * prob * cos(theta);

            // test: if rand < (1+delta_f/f0)/max_value --> accepted
            if (rand() / (double) RAND_MAX <= (1 + delta_f(u[counter], v[counter], w[counter], x[counter], y[counter], kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle, EMf) / f0(u[counter], sqrt(v[counter] * v[counter] + w[counter] * w[counter]))) / max_value)
              rejected = false;

          }
          if (TrackParticleID)
            ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];
          counter++;
        }
  nop = counter + 1;
  // if (vct->getCartesian_rank()==0)
  cout << "Rejection method: " << (counter + 1) / double (total_generated) * 100 << " % of particles are accepted for species " << ns << " counter=" << counter << endl;
}

/** Linear delta f for bi-maxwellian plasma */
double Particles3D::delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double theta, Field * EMf) {
  const complex < double >I(0.0, 1.0);
  const double vperp = sqrt(v * v + w * w);
  const double vpar = u;
  const double kpar = kx;
  double kperp;
  if (ky == 0.0)                // because this formula is not valid for exactly parallel
    kperp = 1e-9;
  else
    kperp = ky;
  const double om_c = qom / c * sqrt(EMf->getBx(1, 1, 0) * EMf->getBx(1, 1, 0) + EMf->getBy(1, 1, 0) * EMf->getBy(1, 1, 0)) / 2 / M_PI;
  const double phi = atan2(w, v);
  const double lambda = kperp * vperp / om_c;
  const complex < double >omega(omega_re, omega_i);

  const int lmax = 5;           // sum from -lmax to lmax

  double bessel_Jn_array[lmax + 2];
  double bessel_Jn_prime_array[lmax + 1];
  complex < double >a1[2 * lmax + 1], a2[2 * lmax + 1], a3[2 * lmax + 1];
  complex < double >factor, deltaf;

  // rotation of x,y
  double temp;
  temp = x;
  x = x * cos(theta) - y * sin(theta);
  y = temp * sin(theta) + y * cos(theta);


/** for compilation issues comment this part: PUT in the math stuff */
  // calc_bessel_Jn_seq(lambda, lmax, bessel_Jn_array, bessel_Jn_prime_array);
  factor = (kpar * vperp / omega * df0_dvpar(vpar, vperp) + (1.0 - (kpar * vpar / omega)) * df0_dvperp(vpar, vperp));
  for (register int l = -lmax; l < 0; l++) {  // negative index
    a1[l + lmax] = factor / lambda * pow(-1.0, -l) * bessel_Jn_array[-l];
    a1[l + lmax] *= (double) l;
    a2[l + lmax] = factor * I * 0.5 * pow(-1.0, -l) * (bessel_Jn_array[-l - 1] - bessel_Jn_array[-l + 1]);
    a3[l + lmax] = kperp / omega * (vpar * df0_dvperp(vpar, vperp) - vperp * df0_dvpar(vpar, vperp)) / lambda * pow(-1.0, -l) * bessel_Jn_array[-l];
    a3[l + lmax] *= (double) l;
    a3[l + lmax] += df0_dvpar(vpar, vperp) * pow(-1.0, -l) * bessel_Jn_array[-l];
  }

  for (register int l = 0; l < lmax + 1; l++) { // positive index
    a1[l + lmax] = factor / lambda * bessel_Jn_array[l];
    a1[l + lmax] *= (double) l;
    a2[l + lmax] = factor * I * bessel_Jn_prime_array[l];
    a3[l + lmax] = kperp / omega * (vpar * df0_dvperp(vpar, vperp) - vperp * df0_dvpar(vpar, vperp)) / lambda * bessel_Jn_array[l];
    a3[l + lmax] *= (double) l;
    a3[l + lmax] += df0_dvpar(vpar, vperp) * bessel_Jn_array[l];
  }

  deltaf = (0.0, 0.0);
  for (register int l = -lmax; l < lmax + 1; l++) {
    deltaf += (a3[l + lmax] * Ex_mod * exp(I * Ex_phase) + a1[l + lmax] * Ey_mod * exp(I * Ey_phase) + a2[l + lmax] * Ez_mod * exp(I * Ez_phase)) / (kpar * vpar + l * om_c - omega) * exp(-I * phi * (double) l);
  }
  deltaf *= I * qom * exp(I * lambda * sin(phi)) * exp(I * (2 * M_PI * kx * x + 2 * M_PI * ky * y));

  return (real(deltaf));
}

double Particles3D::df0_dvpar(double vpar, double vperp) {
  double result;
  result = -2 * (vpar - u0) / uth / uth * exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

double Particles3D::df0_dvperp(double vpar, double vperp) {
  double result;
  result = -2 * (vperp) / vth / vth * exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

double Particles3D::f0(double vpar, double vperp) {
  double result;
  result = exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

void Particles3D::RotatePlaneXY(double theta) {
  double temp, temp2;
  for (register long long s = 0; s < nop; s++) {
    temp = u[s];
    temp2 = v[s];
    u[s] = temp * cos(theta) + v[s] * sin(theta);
    v[s] = -temp * sin(theta) + temp2 * cos(theta);
  }
}
