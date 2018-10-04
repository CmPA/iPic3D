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

/* !!!!!!!!!!! Mathematical functions !!!!!!!!!!! */

/* Modified Bessel functions of the first and second kind */
/*                      First kind                        */
inline double BesselI(double nu, double z)
{
  double tol = 1.e-14;
  double err = 1.;

  double add;
  double Inu = 0.;
  int fact1 = 1;

  int kk = 0;
  while (err > tol) {
    if (kk > 0) fact1 *= kk;
    add = pow(.25*z*z, kk) / fact1 / tgamma(nu+double(kk+1));
    Inu += add;

    err = std::abs(add);
    kk += 1;
  }
  Inu *= pow(.5*z, nu);

  return Inu;
}

/*                     Second kind                        */
inline double BesselKf(double nu, double z)
{
  double Inup = BesselI(nu, z);
  double Inum = BesselI(-nu, z);

  double Knuf = M_PI / 2. * (Inum - Inup) / sin(nu * M_PI);

  return Knuf;
}

/* Digamma function for (positive) integer argument: psi=d(log(Gamma(n)))/dn */
inline double psi(int m)
{
  double eulerg = 0.5772156649015328606065120;
  
  double psi = 0.;
  for (auto n=1; n<m; n++) {
    psi += 1. / double(n);
  }
  psi -= eulerg;
  
  return psi;
}

inline double BesselKi(int nu, double z)
{
  double tol = 1.e-14;
  double err = 1.;

  double add;
  double Knui = 0.;
  int fact1 = 1;
  int fact2 = 1;

  int kk = 0;
  while (err > tol) {
    if (kk > 0) fact1 *= kk;
    fact2 *= nu + kk;
    add = (psi(kk+1) + psi(kk+1+nu)) * pow(.25*z*z, kk) / fact1 / fact2;
    Knui += add;

    err = std::abs(add);
    kk += 1;
  }
  
  Knui *= pow(-1., nu) * .5*pow(.5*z, nu);
  double Inu = BesselI(double(nu), z);
  Knui += pow(-1., nu+1) * log(.5*z) * Inu;

  fact1 = 1;
  fact2 = 1;
  for (auto ii=0; ii<nu; ii++) {
    if (ii > 0) fact1 *= ii;
    if (ii != nu-1) fact2 *= nu - ii - 1;
    else fact2 = 1;
    Knui += .5 * pow(.5*z, -nu) * fact2 / fact1 * pow(-.25*z*z, ii);
  }

  return Knui;
}

inline double BesselK(double nu, double z)
{
  double Knu;
  if (ceilf(nu) != nu) Knu = BesselKf(nu, z); // Non-integer order
  else {
    if (nu < 0) {
      std::cout << " ERROR: Modified Bessel functions of the second kind not implemented for negative orders." << std::endl;
      abort();
    }
    else Knu = BesselKi((int) nu, z); // Integer order
  }

  return Knu;
}
/* !!!!!!!!!!! End of Mathematical functions !!!!!!!!!!! */

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
              mxp[counter] = u0 / sqrt(1.0 - u0*u0) + uth / sqrt(1.0 - uth*uth) * prob * cos(theta);
              // v
              myp[counter] = v0 / sqrt(1.0 - v0*v0) + vth / sqrt(1.0 - vth*vth) * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              mzp[counter] = w0 / sqrt(1.0 - w0*w0) + wth / sqrt(1.0 - wth*wth) * prob * cos(theta);

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
              mxp[counter] = u0 / sqrt(1.0 - u0*u0) + uth / sqrt(1.0 - uth*uth) * prob * cos(theta);
              // v
              myp[counter] = v0 / sqrt(1.0 - v0*v0) + vth / sqrt(1.0 - vth*vth) * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              mzp[counter] = w0 / sqrt(1.0 - w0*w0) + wth / sqrt(1.0 - wth*wth) * prob * cos(theta);

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
              mxp[counter] = -u0 / sqrt(1.0 - u0*u0) + uth / sqrt(1.0 - uth*uth) * prob * cos(theta);
              // v
              myp[counter] = -v0 / sqrt(1.0 - v0*v0) + vth / sqrt(1.0 - vth*vth) * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              mzp[counter] = -w0 / sqrt(1.0 - w0*w0) + wth / sqrt(1.0 - wth*wth) * prob * cos(theta);

              double g = sqrt(1.0 + mxp[counter]*mxp[counter] + myp[counter]*myp[counter] + mzp[counter]*mzp[counter]);
              u[counter] = mxp[counter] / g;
              v[counter] = myp[counter] / g;
              w[counter] = mzp[counter] / g;

              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

              counter++;
            }
}

/** Double periodic GEM distribution: Background (Maxwellian) + Drift (Maxwell-Juttner) */
void Particles3D::relgem2D(Grid * grid, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2 + 10*ns);

  double qbg = (qom / fabs(qom)) * (rhoINIT / npcel) * (1.0 / grid->getInvVOL()); // Charge of the background

  double harvest;
  double prob, theta, sign;
  long long counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        // Check if background
        if (w0 == 0) {         // Background: uniform distribution and Maxwellian velocity
          for (int ii = 0; ii < npcelx; ii++)
            for (int jj = 0; jj < npcely; jj++)
              for (int kk = 0; kk < npcelz; kk++) {
                x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
                y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
                z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
                // q = charge
                q[counter] = qbg;
                // u
                harvest = rand() / (double) RAND_MAX;
                prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
                harvest = rand() / (double) RAND_MAX;
                theta = 2.0 * M_PI * harvest;
                mxp[counter] = u0 / sqrt(1.0 - u0*u0) + uth / sqrt(1.0 - uth*uth) * prob * cos(theta);
                // v
                myp[counter] = v0 / sqrt(1.0 - v0*v0) + vth / sqrt(1.0 - vth*vth) * prob * sin(theta);
                // w
                harvest = rand() / (double) RAND_MAX;
                prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
                harvest = rand() / (double) RAND_MAX;
                theta = 2.0 * M_PI * harvest;
                mzp[counter] = w0 / sqrt(1.0 - w0*w0) + wth / sqrt(1.0 - wth*wth) * prob * cos(theta);

                double g = sqrt(1.0 + mxp[counter]*mxp[counter] + myp[counter]*myp[counter] + mzp[counter]*mzp[counter]);
                u[counter] = mxp[counter] / g;
                v[counter] = myp[counter] / g;
                w[counter] = mzp[counter] / g;

                if (TrackParticleID)
                  ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

                counter++;
              }
        }
        else {           // Current sheet
          // Check if rho is large enough here to place at least 1 particle
          const double eta = 3.0; // Overdensity parameter
          const double xB = grid->getXC(i, j, k) - .25 * Lx;
          const double xT = grid->getXC(i, j, k) - .75 * Lx;
          const double xBd = xB / delta;
          const double xTd = xT / delta;
          const double sech_xBd = 1. / cosh(xBd);
          const double sech_xTd = 1. / cosh(xTd);
          double rhoCS = eta * rhoINIT * sech_xBd * sech_xBd - eta * rhoINIT * sech_xTd * sech_xTd;
          int npcelCS = floor(double(npcel) * std::abs(rhoCS) / rhoINIT); // Particles in this cell
          // If yes, place them uniformly in the cell
          if (npcelCS >= 1) {
            int npcelyCS = (int) max(floor(sqrt(double(npcelCS) * dy/dx)), 1);
            int npcelxCS = (int) ceil(double(npcelCS) / double(npcelyCS));
            int npcelzCS = 1;
            for (int ii = 0; ii < npcelxCS; ii++)
              for (int jj = 0; jj < npcelyCS; jj++)
                for (int kk = 0; kk < npcelzCS; kk++) {
                  x[counter] = (ii + .5) * (dx / npcelxCS) + grid->getXN(i, j, k);
                  y[counter] = (jj + .5) * (dy / npcelyCS) + grid->getYN(i, j, k);
                  z[counter] = (kk + .5) * (dz / npcelzCS) + grid->getZN(i, j, k);
                  // q = charge
                  q[counter] = std::abs(rhoCS)/rhoCS; // Here use the charge as a proxy for the sign of the drift velocity
                  if (TrackParticleID)
                    ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];
                  
                  counter++;
                }
          }
        }    //  Done placing particles

  nop = counter; // Update number of particles

  // If CS particles were placed, initialise a Maxwell-Juttner velocity distribution
  // See Melzani et al. 2013
  if (w0 != 0 && counter > 0) {
    // CS parameters
    double muCS = 1.0 / wth;
    double g0CS = 1.0 / sqrt(1.0 - w0*w0);
    double besselkmu = BesselK(2., muCS);

    /* Marginal distribution function for this proc */
    int Npoints = min(counter,1000); // We want good statistics
    double zz [Npoints];
    double dzmax = 0.001;
    double tol = 1e-14;
    double z0 = -20.;
    double nmax = 10;

    /* Produce a table of values to draw pz from */
    for (auto ii=1; ii<Npoints; ii++) {
      double dz = dzmax;
      double z = z0;
      double J = 3.1553e-10;
      double t = double(ii) / double(Npoints);

      int nk = 0;
      while (std::abs(J-t) > tol) {
        nk += 1;
        double gz = sqrt(1. + z*z);
        double jZ = (1.+muCS*gz*g0CS)/(2.*muCS*g0CS*g0CS*g0CS*besselkmu)*exp(-muCS*g0CS*(gz-std::abs(w0)*z));

        J += dz*jZ;
        if (J < t) z += dz;
        else if (J > t) { J -= dz*jZ; dz /= 2.;}
      }

      zz[ii] = z; // Values stored in zz
    }

    /* Now initialise the particle velocity with the proper distribution */
    for (auto i=0; i<counter; i++) {

      // Drift component drawn from the values stored previously
      harvest = rand() / (double) RAND_MAX;
      int r = 1 + (int)(double(Npoints - 1) * harvest);
      double ppz = std::abs(qom)/qom * q[i] * zz[r];
      double gppz = sqrt(1. + ppz*ppz);
      
      // Other two components*: solve nonlinear equation with Newton
      // *not independent of the first one, which is the whole problem

      double az = muCS * g0CS * gppz;
      harvest = rand() / (double) RAND_MAX;
      double ww = harvest;
      // Apply Newton's method
      double err = 1.;
      int nk = 0;
      double uu = 1.;

      while (err > tol && nk <= nmax) {
        nk += 1;
        double F = (1.+az*uu) * exp(-az*uu) - ww * (1.+az) * exp(-az);
        double D = az * exp(-az*uu) - (1.+az*uu) * az * exp(-az*uu);
        double du = - F / D;
        uu += du;
        err = std::abs(du);
      }

      double rr = sqrt(uu*uu - 1.);
      harvest = rand() / (double) RAND_MAX;
      double th = 2. * M_PI * harvest;
      double ppx = rr * cos(th) * gppz;
      double ppy = rr * sin(th) * gppz;   
      double gpp = sqrt(ppx*ppx + ppy*ppy + ppz*ppz + 1.);

      // u
      mxp[i] = ppx;
      // v
      myp[i] = ppy;
      // w
      mzp[i] = ppz;

      double g = sqrt(1.0 + ppx*ppx + ppy*ppy + ppz*ppz);
      u[i] = mxp[i] / g;
      v[i] = myp[i] / g;
      w[i] = mzp[i] / g;

      // Charge (equal to the background)
      q[i] = qbg;

//cout << u[i] << " " << v[i] << " " << w[i] << endl;
    }
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

void Particles3D::get_Bl(const double weights[2][2][2], int ix, int iy, int iz, double& Bxl, double& Byl, double& Bzl, double*** Bx, double*** By, double*** Bz){

  Bxl = 0.0;
  Byl = 0.0;
  Bzl = 0.0;

  int l = 0;
  for (int i=0; i<=1; i++)
    for (int j=0; j<=1; j++)
      for (int k=0; k<=1; k++) {
        Bxl += weights[i][j][k] * Bx[ix-i][iy-j][iz-k];
        Byl += weights[i][j][k] * By[ix-i][iy-j][iz-k];
        Bzl += weights[i][j][k] * Bz[ix-i][iy-j][iz-k];
        l = l + 1;
      }
}

/** Relativistic ES mover: solve momentum equation explicitly */
int Particles3D::mover_relativistic_mom_ES(Grid * grid, VirtualTopology3D * vct, double*** Ex, double*** Ey, double*** Ez) {
  double start_mover_rel = MPI_Wtime();
  double weights[2][2][2];

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
  }                   // END OF ALL THE PARTICLES

}

/** Relativistic EM mover: solve momentum equation with a Newton iteration particle by particle */
int Particles3D::mover_relativistic_mom_EM(Grid * grid, VirtualTopology3D * vct, double*** Ex, double*** Ey, double*** Ez, double*** Bx, double*** By, double*** Bz) {
  double start_mover_rel = MPI_Wtime();
  double weights[2][2][2];

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
    get_El(weights, ix, iy, iz, Exl, Eyl, Ezl, Ex, Ey, Ez);
    get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz);

    // Precondition the newton variables (can be changed)
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
    // Jacobian is inverted manually for speed
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

    // Update the momentum
    mxp[rest] = pxk;
    myp[rest] = pyk;
    mzp[rest] = pzk;
  }                  // END OF ALL THE PARTICLES

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
