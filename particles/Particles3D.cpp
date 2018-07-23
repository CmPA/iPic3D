/*******************************************************************************************
  Particles3D.cpp  -  Class for particles of the same species, in a 3D space and 3component velocity
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/


#include <iostream>
#include <math.h>

#include "VirtualTopology3D.h"
#include "VCtopology3D.h"
#include "Collective.h"
#include "Collective.h"
#include "Basic.h"
#include "BcParticles.h"
#include "Grid.h"
#include "Grid3DCU.h"
#include "MPIdata.h"
#include "TimeTasks.h"

#include "Particles3D.h"


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
  // see allocate(int species, Collective* col, VirtualTopology3D* vct, Grid* grid)

}
/** deallocate particles */
Particles3D::~Particles3D() {
  delete[]x;
  delete[]y;
  delete[]z;
  delete[]u;
  delete[]v;
  delete[]w;
  delete[]mxp;
  delete[]myp;
  delete[]mzp;
  delete[]q;
}

/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::maxwellian(Grid * grid, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2 + 10*ns);

  double harvest;
  double prob, theta, sign;
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
              // q = charge
              q[counter] = (qom / fabs(qom)) * (rhoINIT / npcel) * (1.0 / grid->getInvVOL());
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              mxp[counter] = u0 / sqrt(1.0 - u0*u0) + uth / sqrt(1.0 - uth*uth) * cos(theta);
              // v
              myp[counter] = v0 / sqrt(1.0 - v0*v0) + vth / sqrt(1.0 - vth*vth) * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              mzp[counter] = w0 / sqrt(1.0 - w0*w0) + wth / sqrt(1.0 - wth*wth) * cos(theta);

              double g = sqrt(1.0 + mxp[counter]*mxp[counter] + myp[counter]*myp[counter] + mzp[counter]*mzp[counter]);
              u[counter] = mxp[counter] / g;
              v[counter] = myp[counter] / g;
              w[counter] = mzp[counter] / g;

              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

              counter++;
            }
}

/** 1D Two-Stream instability and uniform spatial distribution */
void Particles3D::twostream1D(Grid * grid, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2 + 10*ns);

  double harvest;
  double prob, theta, sign;
  long long counter = 0;

  // Beam travelling right
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx/2; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + .5) * (dx / (npcelx/2)) + grid->getXN(i, j, k);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
              // q = charge
              q[counter] = (qom / fabs(qom)) * (rhoINIT / npcel) * (1.0 / grid->getInvVOL());
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              mxp[counter] = u0 / sqrt(1.0 - u0*u0) + uth / sqrt(1.0 - uth*uth) * cos(theta);
              // v
              myp[counter] = v0 / sqrt(1.0 - v0*v0) + vth / sqrt(1.0 - vth*vth) * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              mzp[counter] = w0 / sqrt(1.0 - w0*w0) + wth / sqrt(1.0 - wth*wth) * cos(theta);

              double g = sqrt(1.0 + mxp[counter]*mxp[counter] + myp[counter]*myp[counter] + mzp[counter]*mzp[counter]);
              u[counter] = mxp[counter] / g;
              v[counter] = myp[counter] / g;
              w[counter] = mzp[counter] / g;

              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

              counter++;
            }

  // Beam travelling left
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx/2; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + .5) * (dx / (npcelx/2)) + grid->getXN(i, j, k);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
              // q = charge
              q[counter] = (qom / fabs(qom)) * (rhoINIT / npcel) * (1.0 / grid->getInvVOL());
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              mxp[counter] = -u0 / sqrt(1.0 - u0*u0) + uth / sqrt(1.0 - uth*uth) * cos(theta);
              // v
              myp[counter] = -v0 / sqrt(1.0 - v0*v0) + vth / sqrt(1.0 - vth*vth) * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              mzp[counter] = -w0 / sqrt(1.0 - w0*w0) + wth / sqrt(1.0 - wth*wth) * cos(theta);

              double g = sqrt(1.0 + mxp[counter]*mxp[counter] + myp[counter]*myp[counter] + mzp[counter]*mzp[counter]);
              u[counter] = mxp[counter] / g;
              v[counter] = myp[counter] / g;
              w[counter] = mzp[counter] / g;

              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

              counter++;
            }
}

void Particles3D::get_weights(Grid * grid, double xp, double yp, double zp, int& ix, int& iy, int& iz, double weights[2][2][2]){

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
        weights[ii][jj][kk] = xi[ii] * eta[jj] * zeta[kk] * grid->getInvVOL();
}

void Particles3D::get_El(const double weights[2][2][2], int ix, int iy, int iz, double& Exl, double& Eyl, double& Ezl, double*** Ex, double*** Ey, double*** Ez){

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

void Particles3D::get_Bl(const double weights[2][2][2], int ix, int iy, int iz, double& Bxl, double& Byl, double& Bzl, double*** Bx, double*** By, double*** Bz, double*** Bx_ext, double*** By_ext, double*** Bz_ext, double Fext){

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

/** Relativistic EM mover with a Newton scheme: solve momentum equation */
/** Particles have been moved already by mover_relativistic_pos **/
/*
int Particles3D::mover_relativistic_mom_EM(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  double start_mover_rel = MPI_Wtime();
  double weights[2][2][2];
  // Get fields (at n+1/2 time level) for this processor
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

  for (long long rest = 0; rest < nop; rest++) {
    // copy the particle
    double xp = x[rest];
    double yp = y[rest];
    double zp = z[rest];
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];

    double Exl = 0.0;
    double Eyl = 0.0;
    double Ezl = 0.0;
    double Bxl = 0.0;
    double Byl = 0.0;
    double Bzl = 0.0;
    int ix;
    int iy;
    int iz;

    // Get fields (at n+1/2 time level) at the particle position
    get_weights(grid, xp, yp, zp, ix, iy, iz, weights);
    get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);
    get_El(weights, ix, iy, iz, Exl, Eyl, Ezl, Ex, Ey, Ez);

    // Precondition the newton variables
    double gp = 1. / sqrt(1. - up*up - vp*vp - wp*wp);
    double pxp = up*gp;
    double pyp = vp*gp;
    double pzp = wp*gp;
    double pxk = pxp + qom*dt*(Exl + vp*Bzl - wp*Byl);
    double pyk = pyp + qom*dt*(Eyl - up*Bzl - wp*Bxl);
    double pzk = pzp + qom*dt*(Ezl + up*Byl - vp*Bxl);
    double gk = sqrt(1. + pxk*pxk + pyk*pyk + pzk*pzk);

    double err = 1.;
    // Newton step: maximum number of iterations specified in input file
    for (int nk = 0; nk < NiterMover; nk++) {
      // Residuals
      double F1 = pxk - pxp - qom*dt*(Exl + (pyk+pyp)/(gk+gp)*Bzl - (pzk+pzp)/(gk+gp)*Byl);
      double F2 = pyk - pyp - qom*dt*(Eyl - (pxk+pxp)/(gk+gp)*Bzl + (pzk+pzp)/(gk+gp)*Bxl);
      double F3 = pzk - pzp - qom*dt*(Ezl + (pxk+pxp)/(gk+gp)*Byl - (pyk+pyp)/(gk+gp)*Bxl);

      // Jacobian
      double J11 = 1. - qom*dt*(-(pyk+pyp)/(gk+gp)*Bzl + (pzk+pzp)/(gk+gp)*Byl) * pxk/gk/(gk+gp);
      double J12 = - qom*dt*((-(pyk+pyp)/(gk+gp)*pyk/gk + 1.)*Bzl + (pzk+pzp)/(gk+gp)*pyk/gk*Byl) / (gk+gp);
      double J13 = - qom*dt*(-(pyk+pyp)/(gk+gp)*pzk/gk*Bzl - (-(pzk+pzp)/(gk+gp)*pzk/gk + 1.)*Byl) / (gk+gp);
      double J21 = - qom*dt*(-(-(pxk+pxp)/(gk+gp)*pxk/gk + 1.)*Bzl - (pzk+pzp)/(gk+gp)*pxk/gk*Bxl) /(gk+gp);
      double J22 = 1. - qom*dt*((pxk+pxp)/(gk+gp)*Bzl - (pzk+pzp)/(gk+gp)*Bxl) * pyk/gk / (gk+gp);
      double J23 = - qom*dt*((pxk+pxp)/(gk+gp)*pzk/gk*Bzl + (-(pzk+pzp)/(gk+gp)*pzk/gk + 1.)*Bxl) / (gk+gp);
      double J31 = - qom*dt*((-(pxk+pxp)/(gk+gp)*pxk/gk + 1.)*Byl + (pyk+pyp)/(gk+gp)*pxk/gk*Bxl) /(gk+gp);
      double J32 = - qom*dt*(-(pxk+pxp)/(gk+gp)*pyk/gk*Byl - (-(pyk+pyp)/(gk+gp)*pyk/gk + 1.)*Bxl) / (gk+gp);
      double J33 = 1. - qom*dt*(-(pxk+pxp)/(gk+gp)*Byl + (pyk+pyp)/(gk+gp)*Bxl) * pzk/gk / (gk+gp);

      // Inverse Jacobian
      double Det = J11*J22*J33 + J21*J32*J13 + J31*J12*J23 - J11*J32*J23 - J31*J22*J13 - J21*J12*J33;
      double iJ11 = (J22*J33 - J23*J32) / Det;
      double iJ12 = (J13*J32 - J12*J33) / Det;
      double iJ13 = (J12*J23 - J13*J22) / Det;
      double iJ21 = (J23*J31 - J21*J33) / Det;
      double iJ22 = (J11*J33 - J13*J31) / Det;
      double iJ23 = (J13*J21 - J11*J23) / Det;
      double iJ31 = (J21*J32 - J22*J31) / Det;
      double iJ32 = (J12*J31 - J11*J32) / Det;
      double iJ33 = (J11*J22 - J12*J21) / Det;
      
      // Compute increments dpk=-Jac^-1*Fk
      double dpxk = - (iJ11*F1 + iJ12*F2 + iJ13*F3);
      double dpyk = - (iJ21*F1 + iJ22*F2 + iJ23*F3);
      double dpzk = - (iJ31*F1 + iJ32*F2 + iJ33*F3);

      // Error
      err = sqrt(dpxk*dpxk + dpyk*dpyk + dpzk*dpzk);

      // Update iteration variables
      pxk += dpxk;
      pyk += dpyk;
      pzk += dpzk;
      gk = sqrt(1. + pxk*pxk + pyk*pyk + pzk*pzk);

      if (err < CGtol) break;
    } // END OF THE NEWTON ITERATION

    // Update the final position and velocity
    // Needs to distinguish between a final update (after the fields have converged) and a temporary one (during NK iteration)
    // E.g. with an input parameter "mode"
    // Also needs 3 additional variables mxp,myp,mzp for the particles
    mxp[rest] = pxk;
    myp[rest] = pyk;
    mzp[rest] = pzk;
  } // END OF ALL THE PARTICLES

  // timeTasks.addto_communicate();
  return (0);                   // exit succesfully (hopefully) 

}
*/

/** Relativistic ES mover: solve momentum equation */
/** Particles have been moved already by mover_relativistic_pos **/
int Particles3D::mover_relativistic_mom_ES(Grid * grid, VirtualTopology3D * vct, double*** Ex, double*** Ey, double*** Ez) {
  double start_mover_rel = MPI_Wtime();
  double weights[2][2][2];
//cout << "particle side E field" << endl;
//for (int i = 0; i < nxn; i++)
//    for (int j = 0; j < 1; j++)
//      for (int k = 0; k < 1; k++)
//        cout << Ex[i][j][k] << endl;
  for (long long rest = 0; rest < nop; rest++) {
    // copy the particle
    double xp = x[rest];
    double yp = y[rest];
    double zp = z[rest];
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];

    double Exl = 0.0;
    double Eyl = 0.0;
    double Ezl = 0.0;
    int ix;
    int iy;
    int iz;

    // Get fields (at n+1/2 time level) at the particle position
    get_weights(grid, xp, yp, zp, ix, iy, iz, weights);
    get_El(weights, ix, iy, iz, Exl, Eyl, Ezl, Ex, Ey, Ez);

    // Solve the momentum equation
    double gp = 1.0 / sqrt(1.0 - up*up - vp*vp - wp*wp);
    double pxp = up*gp;
    double pyp = vp*gp;
    double pzp = wp*gp;
    double pxk = pxp + qom*dt*Exl;
    double pyk = pyp + qom*dt*Eyl;
    double pzk = pzp + qom*dt*Ezl;
    // Update the momentum
    mxp[rest] = pxk;
    myp[rest] = pyk;
    mzp[rest] = pzk;
  } // END OF ALL THE PARTICLES

}

/** Relativistic EM mover: move particles */
int Particles3D::mover_relativistic_pos(Grid * grid, VirtualTopology3D * vct) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER species " << ns << endl;
  }
  //double start_mover_rel = MPI_Wtime();

  for (long long rest = 0; rest < nop; rest++) {
    // copy the particle
    double xp = x[rest];
    double yp = y[rest];
    double zp = z[rest];
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];

    // Move particles
    xp += up * dt/2.0;
    yp += vp * dt/2.0;
    zp += wp * dt/2.0;

    x[rest] = xp;
    y[rest] = yp;
    z[rest] = zp;
  }               // END OF ALL THE PARTICLES

  // ********************//
  // COMMUNICATION 
  // *******************//
  // timeTasks.start_communicate();
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
  // timeTasks.addto_communicate();
  return (0);                   // exit succesfully (hopefully) 

}

void Particles3D::mover_relativistic_vel(Grid*, VirtualTopology3D*)
{
  for (long long i=0; i<nop; i++) {
    double gk = sqrt(1.0 + mxp[i]*mxp[i] + myp[i]*myp[i] + mzp[i]*mzp[i]);
    u[i] = mxp[i] / gk;
    v[i] = myp[i] / gk;
    w[i] = mzp[i] / gk;
  }
}
