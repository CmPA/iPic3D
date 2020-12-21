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
#include "Field.h"
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
  delete[]q;
}

/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::initMaxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2 + ns);

  double harvest;
  double prob, theta, sign;
  long long counter = 0;
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
              q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOL());
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

/** Kappa random velocity (anisotropic) and uniform spatial distribution */
void Particles3D::initAnisotropicKappa(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);
  double harvest;
  double reverser;
  double prob, theta, sign;
  // the following variables should be coming from input file, eventually
  double k_kappa = 2.0;//has to be larger than 1.5
  double apar = vth*vth*2.0*(k_kappa-1.5)/k_kappa;
  double aperp = uth*uth*2.0*(k_kappa-1.5)/k_kappa;
  long long counter = 0;
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
              q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOL());
              // u

              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              harvest = rand() / (double) RAND_MAX;

              double vpar=sqrt(k_kappa*apar*(pow(harvest,-1.0/(k_kappa-0.5))-1.0))*cos(theta);

              harvest = rand() / (double) RAND_MAX;

              double vperp=sqrt(k_kappa*aperp*(1.0+vpar*vpar/(k_kappa*apar+1e-10))*(pow(1.0-harvest,-1.0/k_kappa)-1.0));

              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
 //             if(abs(vpar)>sqrt(apar)*10) cout << "vpar=" << vpar << endl;
 //             if(abs(vperp)>sqrt(aperp)*10) cout << "vperp=" << vperp << endl;

              // Assuming now perp to be x and z and parallel y
              u[counter] = u0 + vperp * cos(theta);

              w[counter] = w0 + vperp * sin(theta);

              v[counter] = v0 + vpar;

              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];


              counter++;
            }
}

/** Maxwell-Juttner velocity distribution and uniform spatial distribution */
void Particles3D::initMaxwellJuttner(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor and for different species */
  srand(vct->getCartesian_rank() + 2 + ns);

  // uth here is taken as the temperature, u0/v0/w0 as the shift gamma
  // note: the shift direction and magnitude is determined by preference: first x, second y, third z
  // and by direction: first +, second -
  double Tcrit = 0.1;
  double T = uth;
  double shift_gamma;
  int    shift_dir;
  if      (u0 > 1.) {shift_gamma=u0;  shift_dir=1;}
  else if (u0 < -1.) {shift_gamma=-u0; shift_dir=-1;}
  else if (v0 > 1.) {shift_gamma=v0;  shift_dir=2;}
  else if (v0 < -1.) {shift_gamma=-v0; shift_dir=-2;}
  else if (w0 > 1.) {shift_gamma=w0;  shift_dir=3;}
  else if (w0 < -1.) {shift_gamma=-w0; shift_dir=-3;}
  else {shift_gamma=1.; shift_dir=0;}
  if (shift_gamma < 0.)
    if (vct->getCartesian_rank()==0) {
      cout << "Shift gamma is <0! Aborting..." << endl;
      abort();
    }

  // If T<Tcrit, the Sobol method will likely take too long
  // Revert to tabulated Maxwellian instead
  int npoints = 2000;
  double u_t[npoints];
  double DF_t[npoints];
  if (T<Tcrit) {
    double u_max = 5.*sqrt(2.*T);
    for (int iT=0; iT<npoints; iT++) {
      double u_1 = u_max * (double)(iT - 1) / ((double) npoints);
      double u_2 = u_max * (double) iT / ((double) npoints);
      double df = (u_2 - u_1) * 0.5 * (exp(-u_1*u_1 * 0.5 / T) * u_1*u_1 + exp(-u_2*u_2 * 0.5 / T) * u_2*u_2);
      u_t[iT] = u_2;
      if (iT==0) DF_t[iT] = df;
      else DF_t[iT] = DF_t[iT-1] + df;
    }
    for (int iT=0; iT<npoints; iT++) DF_t[iT] /= DF_t[npoints-1];
  }

  double harvest;
  double prob, fi, sign, gammap, mu, ep, denom, numerator;
  long long counter = 0;
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
              q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOL());

              double U;
              if (T<Tcrit) { // Use tabulated Maxwellian
                double X3 = rand() / (double) RAND_MAX;
                for (int iT=0; iT<npoints; iT++) {
                  if (DF_t[iT] >= X3) {
                    double dx1, dx2;
                    if (iT > 0) {
                      dx1 = (DF_t[iT]-X3)/(DF_t[iT]-DF_t[iT-1]);
                      dx2 = (X3-DF_t[iT-1])/(DF_t[iT]-DF_t[iT-1]);
                      U = u_t[iT]*dx2 + u_t[iT-1]*dx1;
                    }
                    else {
                      dx2 = X3/DF_t[iT];
                      U = u_t[iT]*dx2;
                    }
//                    if (fabs(U)<1.)
                        U /= sqrt(1.-U*U);
//if (U != U) cout << counter << " " << U << endl;
                    break;
                  }
                }
              }
              else { // Sobol method to generate momentum
                int  flag = 0;
                double ETA, X4, X5, X6, X7;
                while (!(flag)) {
                  X4 = rand() / (double) RAND_MAX;
                  X5 = rand() / (double) RAND_MAX;
                  X6 = rand() / (double) RAND_MAX;
                  X7 = rand() / (double) RAND_MAX;
                  if (X4 * X5 * X6 * X7 == 0.) continue;
  
                  U = -T * log(X4 * X5 * X6);
                  ETA = -T * log(X4 * X5 * X6 * X7);
                  if (ETA*ETA - U*U > 1.) flag = 1;
                }
              }

              // Generate projections
              double X1 = rand() / (double) RAND_MAX;
              double X2 = rand() / (double) RAND_MAX;
              double ux = U * (2.0 * X1 - 1.0);
              double uy = 2.0 * U * sqrt(X1 * (1.0 - X1)) * cos(2.0 * M_PI * X2);
              double uz = 2.0 * U * sqrt(X1 * (1.0 - X1)) * sin(2.0 * M_PI * X2);

              // Shift Maxwellian
              if (shift_dir != 0) {
                double X8 = rand() / (double) RAND_MAX;
                double gamma = sqrt(1.0 + U*U);
                double BETA = sqrt(1.0 - 1.0 / (shift_gamma*shift_gamma));
                if (shift_dir == 1) { // +x
                  if (-BETA * ux / gamma > X8) ux = -ux;
                  ux = shift_gamma * (ux + BETA * sqrt(1 + U*U));
                }
                else if (shift_dir == -1) { // -x
                  BETA = -BETA;
                  if (-BETA * ux / gamma > X8) ux = -ux;
                  ux = shift_gamma * (ux + BETA * sqrt(1 + U*U));
                }
                else if (shift_dir == 2) { // +y
                  if (-BETA * uy / gamma > X8) uy = -uy;
                  uy = shift_gamma * (uy + BETA * sqrt(1 + U*U));
                }
                else if (shift_dir == -2) { // -y
                  BETA = -BETA;
                  if (-BETA * uy / gamma > X8) uy = -uy;
                  uy = shift_gamma * (uy + BETA * sqrt(1 + U*U));
                }
                else if (shift_dir == 3) { // +z
                  if (-BETA * uz / gamma > X8) uz = -uz;
                  uz = shift_gamma * (uz + BETA * sqrt(1 + U*U));
                }
                else if (shift_dir == -3) { // -z
                  BETA = -BETA;
                  if (-BETA * uz / gamma > X8) uz = -uz;
                  uz = shift_gamma * (uz + BETA * sqrt(1 + U*U));
                }
              }

              double g = sqrt(1. + ux*ux+uy*uy+uz*uz);
              u[counter] = ux/g;
              v[counter] = uy/g;
              w[counter] = uz/g;
//if (ux != ux || uy != uy || uz != uz) 
//cout << "init particle " << counter << " u,v,w " <<u[counter]<<" " << u[counter] << " " << u[counter] << endl;

              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

              counter++;
            }
}

/** For KAW Turbulence: Maxwellian random velocity, uniform spatial distribution BUT variable number os ppc */
/** Also a perturbation in velocity on top */
void Particles3D::KAWTurbulencePert(Grid * grid, Field * EMf, VirtualTopology3D * vct, double B0x, double mime, double TiTe, bool symmetric) {

  long long seed = (vct->getCartesian_rank() + 1)*20 + ns;
  srand(seed);
  srand48(seed);

  // Profile parameters
  double h = 0.2;
  double r = 10.;
  // Magnetic field parameters
  double B0 = B0x;
  double Bm = 2.*B0x;
  double alpha = (Bm-B0)*r/(2.*pow(2*h,r)*pow(1.+pow(2*h,-r),2));
  // Density parameters
  double betam = 0.5;
  double rhom = rhoINIT/4./M_PI;
  double vthi = sqrt(betam*Bm*Bm/2./4./M_PI/rhom);
  double vthe = vthi*sqrt(mime/TiTe);
  double Ptot = Bm*Bm/2./4./M_PI + rhom*(vthi*vthi + vthe*vthe/mime);
  // Perturbation parameters
  double a = Bm/10.;

  double harvest;
  double prob, theta, sign;
  long long counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++) {
	
        // Density in this cell
        double rho = fabs(EMf->getRHOcs(i, j, k, ns));
        // Calculate nppc in this cell: it will be scaled according to rho/rhom
        double fs=1.;
        if (symmetric) fs=2.;
        int npcelhere = (int) (double(npcel) * rho/rhom / fs);
        int npcelxhere = (int) sqrt(double(npcelhere));
        int npcelyhere = npcelhere/npcelxhere;
        int npcelzhere = 1;
        for (int ii = 0; ii < npcelxhere; ii++)
          for (int jj = 0; jj < npcelyhere; jj++)
            for (int kk = 0; kk < npcelzhere; kk++) {
              x[counter] = (ii + .5) * (dx / npcelxhere) + grid->getXN(i, j, k);
              y[counter] = (jj + .5) * (dy / npcelyhere) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelzhere) + grid->getZN(i, j, k);
              
              // q = charge
              q[counter] = (qom / fabs(qom)) * (rho / (double(npcelhere)*fs)) * (1.0 / invVOL);

              // Drfit velocities for this particle
              double curlB_z = (Bm-B0)*r*pow(y[counter]-Ly/2.,r-1.)/pow(Ly*h,r)/pow(1.+pow((y[counter]-Ly/2.)/Ly/h,r),2)
                               - 2.*alpha/(Ly/2.)*(y[counter]/(Ly/2.)-1);
              double wdrift = curlB_z/4./M_PI/rho*(qom/fabs(qom)) / (1.+1./mime)*fabs(qom)/mime;
              double curlB_y = a*2.*M_PI/Lx*sin(2.*M_PI/Lx*x[counter]);
              double vdrift = curlB_y/4./M_PI/rho*(qom/fabs(qom)) / (1.+1./mime)*fabs(qom)/mime;
              // Add perturbation
              wdrift -= a*cos(2.*M_PI/Lx*x[counter])/sqrt(4.*M_PI*(1.+1./mime)*rho);
              
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              u[counter] = uth * prob * cos(theta);
              // v
              v[counter] = vdrift + vth * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              w[counter] = wdrift + wth * prob * cos(theta);
              
              counter++;

              // Symmetric case: place another particle here
              // Same position, opposite thermal velocity
              if (symmetric) {
                x[counter] = x[counter-1];
                y[counter] = y[counter-1];
                z[counter] = z[counter-1];
                q[counter] = q[counter-1];
                
                u[counter] = -u[counter-1];
                v[counter] = 2.*vdrift - v[counter-1];
                w[counter] = 2.*wdrift - w[counter-1];
               
                counter++;
              }
            }
      }
  nop = counter;
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
        weights[ii][jj][kk] = xi[ii] * eta[jj] * zeta[kk] * invVOL;
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

/** relativistic mover: LM push */
int Particles3D::mover_relativistic(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
	  if (vct->getCartesian_rank() == 0) {
	    cout << "*** MOVER species " << ns << " ***"<< " with " << nop << " particles ***"  << NiterMover << " ITERATIONS   ****" << endl;
	  }
	  double start_mover_PC = MPI_Wtime();
	  double weights[2][2][2];
	  double ***Ex = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getExth());
	  double ***Ey = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEyth());
	  double ***Ez = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEzth());
	  double ***Bx = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx());
	  double ***By = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy());
	  double ***Bz = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz());

	  double ***Bx_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx_ext());
	  double ***By_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy_ext());
	  double ***Bz_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz_ext());

	  double Fext = EMf->getFext();

	  // const double dto2 = .5 * dt, qomdt2 = qom * dto2 / c;
	  // don't bother trying to push any particles simultaneously;
	  // MIC already does vectorization automatically, and trying
	  // to do it by hand only hurts performance.
	//#pragma omp parallel for
	//#pragma simd                    // this just slows things down (why?)

	  for (long long rest = 0; rest < nop; rest++) {
	    // copy the particle
	    double xp = x[rest];
	    double yp = y[rest];
	    double zp = z[rest];
	    double ux0 = u[rest];
	    double uy0 = v[rest];
	    double uz0 = w[rest];
	    double gamma0, gamma, gamma_new;
	    // Assuming up, vp, wp to be VELOCITY and

	    gamma0 = 1.0/sqrt(1.0 - ux0*ux0 - uy0*uy0 - uz0*uz0);
	    ux0 *= gamma0;
	    uy0 *= gamma0;
	    uz0 *= gamma0;
	    const double xp0 = x[rest];
	    const double yp0 = y[rest];
	    const double zp0 = z[rest];
	    double uxnew;
	    double uynew;
	    double uznew;

	    double Exl = 0.0;
	    double Eyl = 0.0;
	    double Ezl = 0.0;
	    double Bxl = 0.0;
	    double Byl = 0.0;
	    double Bzl = 0.0;
	    int ix;
	    int iy;
	    int iz;

	    // No subcyling in relativisitc case

//	    get_weights(grid, xp, yp, zp, ix, iy, iz, weights);
//	    get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);

	   // const double B_mag      = sqrt(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
	    //double       dt_sub     = M_PI*c/(4*abs(qom)*B_mag);
	    //const int    sub_cycles = int(dt/dt_sub) + 1;

	    //dt_sub = dt/double(sub_cycles);
	    double dt_sub = dt;
	    const int sub_cycles = 1;

	    const double dto2 = .5 * dt_sub, qomdt2 = qom * dto2 / c;

	    // if (sub_cycles>1) cout << " >> sub_cycles = " << sub_cycles << endl;

	    for (int cyc_cnt = 0; cyc_cnt < sub_cycles; cyc_cnt++) {

	      // Picard iteration between position and velocity
	      int nit = NiterMover;
	      if (sub_cycles > 2*NiterMover) nit = 1;

              // Initial guess: old velocity
              uxnew = ux0;
              uynew = uy0;
              uznew = uz0;
              gamma_new = gamma0;
	      for (int innter = 0; innter < nit; innter++) {

                // update position (mid of the time step)
	        xp = xp0 + (uxnew + ux0) /(gamma_new + gamma0) * dto2;
	        yp = yp0 + (uynew + uy0) /(gamma_new + gamma0) * dto2;
	        zp = zp0 + (uznew + uz0) /(gamma_new + gamma0) * dto2;

	        // interpolation G-->P
	        get_weights(grid, xp, yp, zp, ix, iy, iz, weights);
	        get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);
	        get_El(weights, ix, iy, iz, Exl, Eyl, Ezl, Ex, Ey, Ez);

	        // solve the momentum equation: 4th order polynomial
                long double epsx = qomdt2*Exl;
                long double epsy = qomdt2*Eyl;
                long double epsz = qomdt2*Ezl;
                long double betax = qomdt2*Bxl;
                long double betay = qomdt2*Byl;
                long double betaz = qomdt2*Bzl;
                long double beta2 = betax*betax+betay*betay+betaz*betaz;
                long double upx = ux0 + epsx;
                long double upy = uy0 + epsy;
                long double upz = uz0 + epsz;
//cout << "subcycle " << cyc_cnt <<" particle " << rest << " epsx, betax, upx " << epsx << " "<< betax << " " << upx << endl;
                
                // Polynomial coefficients
                long double updote = upx*epsx+upy*epsy+upz*epsz;
//                if (fabs(updote)<1.e-14) updote = 0.;
                long double bdote = betax*epsx+betay*epsy+betaz*epsz;
//                if (fabs(bdote)<1.e-14) bdote = 0.;
                long double updotb = upx*betax+upy*betay+upz*betaz;
//                if (fabs(updotb)<1.e-14) updotb = 0.;
                long double upcrossb_x = (upy*betaz-upz*betay);
                long double upcrossb_y = (-upx*betaz+upz*betax);
                long double upcrossb_z = (upx*betay-upy*betax);
//                if (fabs(updotb-sqrt(upx*upx+upy*upy+upz*upz))<1.e-14) {
//                  upcrossb_x = 0.;
//                  upcrossb_y = 0.;
//                  upcrossb_z = 0.;
//                }                
                long double aa = updote - beta2;
                long double bb = upcrossb_x*epsx+upcrossb_y*epsy+upcrossb_z*epsz+ gamma0*beta2;
                long double cc = updotb*bdote;
//cout << "subcycle " << cyc_cnt <<" particle " << rest << " coeffs0 " << aa << " "<< bb << " " << cc << endl;
                /* METHOD 1: DIRECT SOLVE */
                
                // Solution coefficients
                double AA = 2.*aa/3.+gamma0*gamma0/4.;
                double BB = 4.*aa*gamma0+8.*bb+gamma0*gamma0*gamma0;
                double DD = aa*aa-3.*bb*gamma0-12.*cc;
                double FF = -2.*aa*aa*aa+9.*aa*bb*gamma0-72.*aa*cc+27.*bb*bb-27.*cc*gamma0*gamma0;
                std::complex<double> GG = FF*FF-4.*DD*DD*DD;
//                if (fabs(GG)<1.e-14) GG = 0.;
//                std::complex<double> EE = cbrt((FF+sqrt(GG))/2.); // pow((FF+sqrt(GG))/2.,1./3.);
                std::complex<double> EE;
                if (std::real((FF+sqrt(GG))/2.)<0.) EE = -pow(-(FF+sqrt(GG))/2.,1./3.);
                else EE = pow((FF+sqrt(GG))/2.,1./3.);
                std::complex<double> CC = DD/(EE+1.e-20)/3.+EE/3.;
                // Solution
                std::complex<double> gbarc = gamma0/4.+sqrt(2.*AA+BB/4./sqrt(AA+CC+1.e-20)-CC)/2.+sqrt(AA+CC)/2.;
                double gbar = (double) std::real(gbarc);
//cout << "subcycle " << innter <<" particle " << rest << " coeffs " << AA << " "<< BB << " " << CC << " " << DD << " " << EE << " " << FF << " " << GG << endl;
                
                /* METHOD 2: PolyRoots */
                /*
                long double rrr[4];
                long double iii[4];
                long double ccc[5];
                ccc[0] = -1.; ccc[1] = gamma0; ccc[2] = aa; ccc[3] = bb; ccc[4] = cc;
if (ccc[4]<0.) ccc[4] = -ccc[4];
if (fabs(ccc[4])<1.e-13) ccc[4] = 0.;
                QuadCubicRoots(ccc, 4, rrr, iii);
                double gbar = 1.0;
                for (int im=0; im<4; im++)
                  if (iii[im] == 0. && rrr[im]>gbar) gbar = rrr[im];
cout << "co " << ccc[1] <<" " << ccc[2] << " " << ccc[3] << " " << ccc[4] << endl;
cout << "im " << iii[0] <<" " << iii[1] << " " << iii[2] << " " << iii[3] << endl;
cout << "re " << rrr[0] <<" " << rrr[1] << " " << rrr[2] << " " << rrr[3] << endl;
cout << "!!!!!!!!!!!!!!!!!!!! subcycle " << innter <<" particle " << rest << " gbar " << gbar << endl;
                */

//abort();
                double uxbar = (upx+(upx*betax+upy*betay+upz*betaz)*betax/(gbar*gbar)+(upy*betaz-upz*betay)/gbar)/(1.+beta2/gbar/gbar);
                double uybar = (upy+(upx*betax+upy*betay+upz*betaz)*betay/(gbar*gbar)+(-upx*betaz+upz*betax)/gbar)/(1.+beta2/gbar/gbar);
                double uzbar = (upz+(upx*betax+upy*betay+upz*betaz)*betaz/(gbar*gbar)+(upx*betay-upy*betax)/gbar)/(1.+beta2/gbar/gbar);
                
                uxnew = 2.*uxbar - ux0;
                uynew = 2.*uybar - uy0;
                uznew = 2.*uzbar - uz0;
                gamma_new = 2.*gbar - gamma0;

	      } // end of velocity iteration
	      // update the final position and velocity

              gamma_new = sqrt(1.+uxnew*uxnew+uynew*uynew+uznew*uznew); 
	      u[rest] = uxnew/gamma_new;
	      v[rest] = uynew/gamma_new;
	      w[rest] = uznew/gamma_new;

	      x[rest] = xp0 + (uxnew + ux0)/(gamma_new + gamma0) * dt;
	      y[rest] = yp0 + (uynew + uy0)/(gamma_new + gamma0) * dt;
	      z[rest] = zp0 + (uznew + uz0)/(gamma_new + gamma0) * dt;

	    } // END  OF SUBCYCLING LOOP
//cout << "final: particle " << rest << " ux " << u[rest] << endl;
//if (fabs(u[rest])>1. || fabs(v[rest])>1. || fabs(w[rest])>1.)
// cout << "moved particle " << rest << " u,v,w " <<u[rest]<<" " << v[rest] << " " << w[rest] << endl;
	  }                             // END OF ALL THE PARTICLES
//abort();
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
	  return (0);                   // exit succcesfully (hopefully)
}

int Particles3D::particle_repopulator(Grid* grid,VirtualTopology3D* vct, Field* EMf, int is){

  /* -- NOTE: Hardcoded option -- */
  enum {LINEAR,INITIAL,FFIELD};
  int rtype = FFIELD;
  /* -- END NOTE -- */

  if (vct->getCartesian_rank()==0){
    cout << "*** Repopulator species " << ns << " ***" << endl;
  }
  double  FourPI =16*atan(1.0);
  int avail;
  long long store_nop=nop;

  ////////////////////////
  // INJECTION FROM XLEFT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcPfaceXleft > 1){ // use Field topology in this case
    long long particles_index=0;
    long long nplast = nop-1;

    while (particles_index < nplast+1) {
      if (x[particles_index] < 3.0*dx ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }

    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< 4;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);

                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];


                particles_index++ ;
              }
    nop = particles_index;
  }

  store_nop = nop;

  ////////////////////////
  // INJECTION FROM YLEFT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getYleft_neighbor() == MPI_PROC_NULL  && bcPfaceYleft > 1)
  {
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (y[particles_index] < 3.0*dy ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=1; j< 4;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM ZLEFT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getZleft_neighbor() == MPI_PROC_NULL  && bcPfaceZleft >1)
  {
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (z[particles_index] < 3.0*dz ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< 4;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM XRIGHT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getXright_neighbor() == MPI_PROC_NULL  && bcPfaceXright > 1){
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (x[particles_index] > (Lx-3.0*dx) ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=(grid->getNXC()-4); i< grid->getNXC()-1;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM YRIGHT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getYright_neighbor() == MPI_PROC_NULL  && bcPfaceYright > 1)
  {
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (y[particles_index] > (Ly-3.0*dy) ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=(grid->getNYC()-4); j< grid->getNYC()-1;j++)
        for (int k=1; k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  ////////////////////////
  // INJECTION FROM ZRIGHT
  ////////////////////////
  srand (vct->getCartesian_rank()+1+ns+(int(MPI_Wtime()))%10000);
  if (vct->getZright_neighbor() == MPI_PROC_NULL  && bcPfaceZright > 1)
  {
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (z[particles_index] > (Lz-3.0*dz) ) {
        del_pack(particles_index,&nplast);
      } else {
        particles_index++;
      }
    }
    nop = nplast+1;
    particles_index = nop;
    double harvest;
    double prob, theta, sign;
    //particles_index;
    for (int i=1; i< grid->getNXC()-1;i++)
      for (int j=1; j< grid->getNYC()-1;j++)
        for (int k=(grid->getNZC()-4); k< grid->getNZC()-1;k++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++)
              for (int kk=0; kk < npcelz; kk++){
                harvest =   rand()/(double)RAND_MAX ;
                x[particles_index] = (ii + harvest)*(dx/npcelx) + grid->getXN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                y[particles_index] = (jj + harvest)*(dy/npcely) + grid->getYN(i,j,k);
                harvest =   rand()/(double)RAND_MAX ;
                z[particles_index] = (kk + harvest)*(dz/npcelz) + grid->getZN(i,j,k);
                // q = charge
                /* ATTENTION: OVther methods can be use, i.e. using the values close to the boundary: */
                   double rho = 1.0/FourPI;
                   if (rtype==FFIELD)  rho = EMf->getRHOcs(i,j,k,is);
                   if (rtype==LINEAR)  rho = (0.1 + 0.9*(grid->getXC(i, j, k)/Lx)) / FourPI;
                   if (rtype==INITIAL) rho = rhoINJECT/FourPI;
                   q[particles_index] = (qom / fabs(qom))*(fabs(rho)/npcel)*(1.0/grid->getInvVOL());
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[particles_index] = u0 + uth*prob*cos(theta);
                // v
                v[particles_index] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[particles_index] = w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[particles_index]= particles_index*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                particles_index++ ;
              }
    nop = particles_index;
  }

  if (vct->getCartesian_rank()==0){
    cout << "*** number of particles " << nop << " ***" << endl;
  }

  //********************//
  // COMMUNICATION
  // *******************//
  avail = communicate(vct);
  if (avail < 0) return(-1);

  MPI_Barrier(MPI_COMM_WORLD);

  // communicate again if particles are not in the correct domain
  while(isMessagingDone(vct) >0){
    // COMMUNICATION
    avail = communicate(vct);
    if (avail < 0)
      return(-1);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return(0); // exit succcesfully (hopefully)
}


int Particles3D::particle_reflector(Grid* grid,VirtualTopology3D* vct, Field* EMf, int is){

  /* -- NOTE: Hardcoded option -- */
  enum {LINEAR,INITIAL,FFIELD};
  int rtype = FFIELD;
  /* -- END NOTE -- */

  if (vct->getCartesian_rank()==0){
    cout << "*** Repopulator species " << ns << " ***" << endl;
  }
  double weights[2][2][2];
  double B_mag;
  double bdotn;
  double bxin ;
  double byin ;
  double bzin ;
  double vdotb;
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

  double Exl = 0.0;
  double Eyl = 0.0;
  double Ezl = 0.0;
  double Bxl = 0.0;
  double Byl = 0.0;
  double Bzl = 0.0;
  int ix;
  int iy;
  int iz;
  double  FourPI =16*atan(1.0);
  int avail;
  long long store_nop=nop;
  double vtestin, vtestout,vabsin,vabsout;

  ////////////////////////
  // INJECTION FROM XLEFT
  ////////////////////////

  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcPfaceXleft == 3){ // use Field topology in this case
    long long particles_index=0;
    long long nplast = nop-1;

    while (particles_index < nplast+1) {
      if (x[particles_index] > 3.0*dx && x[particles_index] < 4.0*dx) {
        // reflect particles
          //vtestin = u[particles_index];
          //vabsin = u[particles_index] * u[particles_index] + v[particles_index] * v[particles_index] + w[particles_index] * w[particles_index];
    	      get_weights(grid, x[particles_index], y[particles_index], z[particles_index], ix, iy, iz, weights);
    	      get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);

    	       // Compute the unit vector  b along B and pointing INSIDE the domain.
     	       // This is achieved by computing the dot product of B with the outgoing normal
     	       // And then imposing the vector b to point opposite to the outgoing normal

    	       B_mag      = sqrt(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
    	       bdotn = -Bxl /abs(Bxl+1e-10);
    	       bxin = -Bxl/(B_mag+1e-10)*bdotn;
    	       byin = -Byl/(B_mag+1e-10)*bdotn;
    	       bzin = -Bzl/(B_mag+1e-10)*bdotn;
    	       vdotb = u[particles_index]*bxin + v[particles_index]*byin + w[particles_index]*bzin;
    	       if(abs(vdotb)>uth){
    	      u[particles_index] += -bxin * ( vdotb - abs(vdotb));
    	      v[particles_index] += -byin * ( vdotb - abs(vdotb));
    	      w[particles_index] += -bzin * ( vdotb - abs(vdotb));
    	       }
      	     // vdotb = u[particles_index]*bxin + v[particles_index]*byin + w[particles_index]*bzin;
             // vtestout = u[particles_index];
             // vabsout = u[particles_index] * u[particles_index] + v[particles_index] * v[particles_index] + w[particles_index] * w[particles_index];
             // if(vtestin < 0.0 )
             // cout << "test reflector  " << "vtestin=" << vtestin << "  vtestout=" << vtestout<< "   bxin="<< bxin<< "   vdotb="<< vdotb<<endl;
              particles_index++;
      } else {
        particles_index++;
      }
    }

  }

  ////////////////////////
  // INJECTION FROM XRIGHT
  ////////////////////////

  if (vct->getXright_neighbor() == MPI_PROC_NULL  && bcPfaceXright == 3){
    long long particles_index=0;
    long long nplast = nop-1;
    while (particles_index < nplast+1) {
      if (x[particles_index] < (Lx-3.0*dx) && x[particles_index] > (Lx-4.0*dx)) {
          // reflect particles

      	      get_weights(grid, x[particles_index], y[particles_index], z[particles_index], ix, iy, iz, weights);
      	      get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);


      	       B_mag      = sqrt(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
     	       bdotn = Bxl /abs(Bxl+1e-10);
      	       bxin = -Bxl/(B_mag+1e-10)*bdotn;
      	       byin = -Byl/(B_mag+1e-10)*bdotn;
      	       bzin = -Bzl/(B_mag+1e-10)*bdotn;
      	       vdotb = u[particles_index]*bxin + v[particles_index]*byin + w[particles_index]*bzin;
      	      u[particles_index] += -bxin * ( vdotb - abs(vdotb));
      	      v[particles_index] += -byin * ( vdotb - abs(vdotb));
      	      w[particles_index] += -bzin * ( vdotb - abs(vdotb));

      	      particles_index++;
      } else {
        particles_index++;
      }
    }

  }
  if (vct->getCartesian_rank()==0){
    cout << "*** Reflector: number of particles " << nop << " ***" << endl;
  }

  //********************//
  // COMMUNICATION
  // *******************//
  avail = communicate(vct);
  if (avail < 0) return(-1);

  MPI_Barrier(MPI_COMM_WORLD);

  // communicate again if particles are not in the correct domain
  while(isMessagingDone(vct) >0){
    // COMMUNICATION
    avail = communicate(vct);
    if (avail < 0)
      return(-1);
    MPI_Barrier(MPI_COMM_WORLD);
  }

  return(0); // exit succcesfully (hopefully)
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
  double ep;
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

    ep = 0.5 / qom *( u[i]*u[i] + v[i]*v[i] + w[i]*w[i]);

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
    // EFx
    eqValue(0.0, temp, 2, 2, 2);
    addscale(u[i]*ep, temp, weight, 2, 2, 2);
    EMf->addEFx(temp, ix, iy, iz, ns);
    // EFy
    eqValue(0.0, temp, 2, 2, 2);
    addscale(v[i]*ep, temp, weight, 2, 2, 2);
    EMf->addEFy(temp, ix, iy, iz, ns);
    // EFz
    eqValue(0.0, temp, 2, 2, 2);
    addscale(w[i]*ep, temp, weight, 2, 2, 2);
    EMf->addEFz(temp, ix, iy, iz, ns);

  }
  // communicate contribution from ghost cells 
  EMf->communicateGhostP2G(ns, 0, 0, 0, 0, vct);
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

/*! Delete the particles inside the sphere with radius R and center x_center y_center and return the total charge removed */
double Particles3D::deleteParticlesInsideSphere(double R, double x_center, double y_center, double z_center){

  long long np_current = 0;
  long long nplast     = nop-1;

  while (np_current < nplast+1){

    double xd = x[np_current] - x_center;
    double yd = y[np_current] - y_center;
    double zd = z[np_current] - z_center;

    if ( (xd*xd+yd*yd+zd*zd) < R*R ){
      Q_removed += q[np_current];
      del_pack(np_current,&nplast);

    } else {
      np_current++;
    }
  }
  nop = nplast +1;
  return(Q_removed);
}


/** Delete the particles outer multx*dx, multy*dy, multz*dz celles */
double Particles3D::deleteParticlesOuterFrame(double multx, double multy, double multz){
	// calculate accumulated charge on the spacecraft
	long long np_current = 0, nplast = nop-1;
	while (np_current < nplast+1) {
		if (x[np_current] < multx*dx || x[np_current] > (Lx-multx*dx) ||
		    y[np_current] < multy*dy || y[np_current] > (Ly-multy*dx) ||
		    z[np_current] < multz*dz || z[np_current] > (Lz-multz*dx)) {
			Q_removed += q[np_current];
			del_pack(np_current,&nplast);
		} else {
			np_current++;
		}
	}

	nop = nplast +1;
	return(Q_removed);
}

/** Reflect particles in the outer shell towards the center */
double Particles3D::ReturnToCenterOuterFrame(double multx, double multy, double multz){
	// calculate accumulated charge on the spacecraft
	long long np_current = 0, nplast = nop-1;
	double r;
	double udotr;
	while (np_current < nplast+1) {
		if (x[np_current] < multx*dx || x[np_current] > (Lx-multx*dx) ||
		    y[np_current] < multy*dy || y[np_current] > (Ly-multy*dx) ||
		    z[np_current] < multz*dz || z[np_current] > (Lz-multz*dx)) {
			r = 1e-10+sqrt((x[np_current]-Lx/2.0)*(x[np_current]-Lx/2.0) +
					(y[np_current]-Ly/2.0)*(y[np_current]-Ly/2.0) +
					(z[np_current]-Lz/2.0)*(z[np_current]-Lz/2.0));
			udotr= (u[np_current] * (x[np_current]-Lx/2.0)+
					v[np_current] * (y[np_current]-Ly/2.0)+
					w[np_current] * (z[np_current]-Lz/2.0))/r;
		//	udotr += abs(udotr);
if(udotr>0){
			u[np_current] = u[np_current] -  2.0 *udotr * (x[np_current]-Lx/2.0)/r;
			v[np_current] = v[np_current] -  2.0* udotr * (y[np_current]-Ly/2.0)/r;
			w[np_current] = w[np_current] -  2.0 *udotr * (z[np_current]-Lz/2.0)/r; }
			np_current++;
		} else {
			np_current++;
		}
	}
	nop = nplast +1;
	return(0.0);
}
double Particles3D::ReturnToCenterCircle(){
	// calculate accumulated charge on the spacecraft
	long long np_current = 0, nplast = nop-1;
	double r;
	double udotr;
	double ukick = 3.0 *uth;
	while (np_current < nplast+1) {
		r = 1e-10+sqrt((x[np_current]-Lx/2.0)*(x[np_current]-Lx/2.0) +
							(y[np_current]-Ly/2.0)*(y[np_current]-Ly/2.0) +
							(z[np_current]-Lz/2.0)*(z[np_current]-Lz/2.0));
		if (r> L_outer) {
/*			udotr= (u[np_current] * (x[np_current]-Lx/2.0)+
					v[np_current] * (y[np_current]-Ly/2.0)+
					w[np_current] * (z[np_current]-Lz/2.0))/r;
		//	udotr += abs(udotr);
if(udotr>0){
			u[np_current] = u[np_current] -  2.0 *udotr * (x[np_current]-Lx/2.0)/r;
			v[np_current] = v[np_current] -  2.0* udotr * (y[np_current]-Ly/2.0)/r;
			w[np_current] = w[np_current] -  2.0 *udotr * (z[np_current]-Lz/2.0)/r; }
			*/
			u[np_current] = 0.0*u[np_current] -  ukick * (x[np_current]-Lx/2.0)/r;
			v[np_current] = 0.0*v[np_current] -  ukick * (y[np_current]-Ly/2.0)/r;
			w[np_current] = 0.0*w[np_current] -  ukick * (z[np_current]-Lz/2.0)/r;

			np_current++;
		} else {
			np_current++;
		}
	}
	nop = nplast +1;
	return(0.0);
}
double Particles3D::ReturnRegeneratedToCenterOuterFrame(double multx, double multy, double multz){
	// calculate accumulated charge on the spacecraft
	long long np_current = 0, nplast = nop-1;
	  double harvest;
	  double prob, theta, sign;
	double r;
	double udotr;
	while (np_current < nplast+1) {
		if (x[np_current] < multx*dx || x[np_current] > (Lx-multx*dx) ||
		    y[np_current] < multy*dy || y[np_current] > (Ly-multy*dx) ||
		    z[np_current] < multz*dz || z[np_current] > (Lz-multz*dx)) {
			r = 1e-10+sqrt((x[np_current]-Lx/2.0)*(x[np_current]-Lx/2.0) +
					(y[np_current]-Ly/2.0)*(y[np_current]-Ly/2.0) +
					(z[np_current]-Lz/2.0)*(z[np_current]-Lz/2.0));
			// regenerate particle speed from intial maxwellian
            // u
            harvest = rand() / (double) RAND_MAX;
            prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
            harvest = rand() / (double) RAND_MAX;
            theta = 2.0 * M_PI * harvest;
            u[np_current] = u0 + uth * prob * cos(theta);
            // v
            v[np_current] = v0 + vth * prob * sin(theta);
            // w
            harvest = rand() / (double) RAND_MAX;
            prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
            harvest = rand() / (double) RAND_MAX;
            theta = 2.0 * M_PI * harvest;
            w[np_current] = w0 + wth * prob * cos(theta);

            // turn it towards the center if it is not going there already
			udotr= (u[np_current] * (x[np_current]-Lx/2.0)+
					v[np_current] * (y[np_current]-Ly/2.0)+
					w[np_current] * (z[np_current]-Lz/2.0))/r;
		//	udotr += abs(udotr);
if(udotr>0){
			u[np_current] = u[np_current] -  2.0 *udotr * (x[np_current]-Lx/2.0)/r;
			v[np_current] = v[np_current] -  2.0* udotr * (y[np_current]-Ly/2.0)/r;
			w[np_current] = w[np_current] -  2.0 *udotr * (z[np_current]-Lz/2.0)/r; }
			np_current++;
		} else {
			np_current++;
		}
	}
	nop = nplast +1;
	return(0.0);
}

/** Delete the particles inside the sphere with radius R and center x_center y_center */
double Particles3D::deleteParticlesOutsideSphere(double R, double x_center, double y_center, double z_center){
	// calculate accumulated charge on the spacecraft
	long long np_current = 0, nplast = nop-1;
	while (np_current < nplast+1){
		if (sqrt( pow(x[np_current] - x_center,2) + pow(y[np_current] - y_center,2) + pow(z[np_current] - z_center,2) ) > R){
			// delete the particle and pack the particle array, the value of nplast changes
			Q_removed += q[np_current];
			del_pack(np_current,&nplast);

		} else {
			// particle is still in the domain, procede with the next particle
			np_current++;
		}
	}
	nop = nplast +1;
	return(Q_removed);
}

//** Delete the particles outside the cube with dimension L */
double Particles3D::deleteParticlesOutsideBox(double L){
	// calculate accumulated charge on the spacecraft
	long long np_current = 0, nplast = nop-1;
	while (np_current < nplast+1){
		if (x[np_current] > L || x[np_current] < 0 ||
			y[np_current] > L || y[np_current] < 0 ||
			z[np_current] > L || z[np_current] < 0){
			// delete the particle and pack the particle array, the value of nplast changes
			Q_removed += q[np_current];
			del_pack(np_current,&nplast);

		} else {
			// particle is still in the domain, procede with the next particle
			np_current++;
		}
	}
	nop = nplast +1;
	return(Q_removed);
}
/* Computes overalp of two intervals */
double Particles3D::interval_overlap(double xd0, double xd1, double xp0, double xp1, double& xstart_interval)
{
	double dx_inject;
	xstart_interval = max(xp0, xd0);
	double xend_interval = min(xp1,xd1);
	if(xstart_interval < xend_interval )
		dx_inject = xend_interval-xstart_interval;
	else
		dx_inject = 0.0;
	return(dx_inject);
}

/* Inject from a box, checking overal with the local processor */
int Particles3D::injector_rand_box(Grid* grid,VirtualTopology3D* vct, Field* EMf, double x_center_inject, double y_center_inject, double z_center_inject, double L_inject)
{
	double harvest;
	double prob, theta, sign;
	double r;
	int avail;
	long long store_nop=nop;
	long long counter=nop;

	//Vinj from the input file has hte number of particles to be injected
	//Ninj is the RHOinject from the input file
	double Iinj=Vinj;
	long long npinject=1;
	if (Iinj > 1)
		npinject = (int)floor(Iinj);

	double xstart_inject;
	double dx_inject = interval_overlap(x_center_inject-L_inject,x_center_inject+L_inject, xstart, xend, xstart_inject);
	double ystart_inject;
	double dy_inject = interval_overlap(y_center_inject-L_inject,y_center_inject+L_inject, ystart, yend, ystart_inject);
	double zstart_inject;
	double dz_inject = interval_overlap(z_center_inject-L_inject,z_center_inject+L_inject, zstart, zend, zstart_inject);

//	cout << xstart_inject  << "  " << ystart_inject << "   " << zstart_inject << endl;
	if (dx_inject * dy_inject * dz_inject > 0.0)
	{

//cout << dx_inject * dy_inject * dz_inject  << "  " << npinject << endl;

		for (int inject = 0; inject < npinject; inject++)
		{
			harvest =   rand()/(double)RAND_MAX ;
			x[nop] = xstart_inject+(harvest) * dx_inject;
			harvest =   rand()/(double)RAND_MAX ;
			y[nop] = ystart_inject+(harvest) * dy_inject;
			harvest =   rand()/(double)RAND_MAX ;
			z[nop] = zstart_inject+(harvest) * dz_inject;
			r = 1e-10+sqrt((x[nop]-Lx/2.0)*(x[nop]-Lx/2.0) +
					(y[nop]-Ly/2.0)*(y[nop]-Ly/2.0) +
					(z[nop]-Lz/2.0)*(z[nop]-Lz/2.0));

			q[nop] =  (qom/fabs(qom)) * rhoINJECT * dx_inject * dy_inject * dz_inject / npinject;
			// u
			harvest =   rand()/(double)RAND_MAX;
			prob  = sqrt(-2.0*log(1.0-.999999*harvest));
			harvest =   rand()/(double)RAND_MAX;
			theta = 2.0*M_PI*harvest;

			u[nop] = - v0 * (x[nop]-Lx/2.0)/r + uth*prob*cos(theta);
			v[nop] = - v0 * (y[nop]-Ly/2.0)/r + vth*prob*sin(theta);
			harvest =   rand()/(double)RAND_MAX;
			prob  = sqrt(-2.0*log(1.0-.999999*harvest));
			harvest =   rand()/(double)RAND_MAX;
			theta = 2.0*M_PI*harvest;
			w[nop] = - v0 * (z[nop]-Lz/2.0)/r + wth*prob*cos(theta);;


			if (TrackParticleID)
				ParticleID[nop]= nop*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

			//printf("*** DEBUG: Added a particle at (x,y,z) = (%f,%f,%f) with (vx,vy,vz) = (%f,%f,%f)\n", x[nop], y[nop], z[nop], u[nop], v[nop], w[nop]);

			nop++;
//			cout << nop << endl;
		}

	}
		// move particles inside the sim box
		for (int i=store_nop; i < nop; i++){
			x[i] += u[i]*dt;
			y[i] += v[i]*dt;
			z[i] += w[i]*dt;
		}


		// ******************** //
		// COMMUNICATION
		// ******************* //
	    avail = communicate(vct);
	    if (avail < 0){
		cout <<"Cannot communicate! Failing!!!!\n";
		return(-1);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
//	    cout << "Past barrier...\n";
	    // communicate again if particles are not in the correct domain
	    while(isMessagingDone(vct) >0){
			// COMMUNICATION
			avail = communicate(vct);

			if (avail < 0)
		        return(-1);
			MPI_Barrier(MPI_COMM_WORLD);
		}


//	cout << "Communicated!!\n";

	return(0); // exit succcesfully (hopefully)

}



int Particles3D::injector_rand_box_mono(Grid* grid,VirtualTopology3D* vct, Field* EMf)
{
	double harvest;
	double prob, theta, sign, fi;
	int avail;
	long long store_nop=nop;
	long long counter=nop;

	//double Iinj=0.002;
	double Iinj=Vinj;
	long long npinject=1;
	if (Iinj > 1)
		npinject = (int)floor(Iinj);

	if (vct->getCartesian_rank() == 0){
			cout << "*** Injector Rand Box species " << ns << " ***" << NiterMover <<" ITERATIONS : INJECTING " << npinject << " particles, each with a charge of " << (qom/fabs(qom))*(EMf->getRHOcs(1,1,1,ns)/npcel)*(1.0/grid->getInvVOL()) << " C  ****" << endl;


		for (long long inject = 0; inject < npinject; inject++)
		{
			harvest =   rand()/(double)RAND_MAX ;
			x[nop] = x_center+(harvest-0.5)*L_square;
			harvest =   rand()/(double)RAND_MAX ;
			y[nop] = y_center+(harvest-0.5)*L_square;
			harvest =   rand()/(double)RAND_MAX ;
			z[nop] = z_center; //z_center+(harvest-0.5)*L_square;
			//cout << "x=" <<x[nop] << "y=" << y[nop]<< "z=" << z[nop] <<endl;
			//this assigns charge the same as does maxwell_box
			q[nop] =  (qom/fabs(qom))*(EMf->getRHOcs(1,1,1,ns)/npcel)*pow(L_square,3);

			// random numbers for Maxwellian
			harvest =   rand()/(double)RAND_MAX;
			fi  = M_PI*harvest;
			harvest =   rand()/(double)RAND_MAX;
			theta = 2.0*M_PI*harvest;

			//u
			u[nop] = u0 + uth*cos(fi)*cos(theta);

			// v
			v[nop] = v0 + vth*cos(fi)*sin(theta);

			// w
			w[nop] = w0 + wth*sin(fi);

			if (TrackParticleID)
				ParticleID[nop]= nop*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

			//printf("*** DEBUG: Added a particle at (x,y,z) = (%f,%f,%f) with (vx,vy,vz) = (%f,%f,%f)\n", x[nop], y[nop], z[nop], u[nop], v[nop], w[nop]);

			nop++;
		}


		// move particles inside the sim box
		for (long long i=store_nop; i < nop; i++){
			x[i] += u[i]*dt;
			y[i] += v[i]*dt;
			z[i] += w[i]*dt;
		}
	}

		// ******************** //
		// COMMUNICATION
		// ******************* //
	    avail = communicate(vct);
	    if (avail < 0){
		cout <<"Cannot communicate! Failing!!!!\n";
		return(-1);
	    }
	    MPI_Barrier(MPI_COMM_WORLD);
//	    cout << "Past barrier...\n";
	    // communicate again if particles are not in the correct domain
	    while(isMessagingDone(vct) >0){
			// COMMUNICATION
			avail = communicate(vct);

			if (avail < 0)
		        return(-1);
			MPI_Barrier(MPI_COMM_WORLD);
		}


//	cout << "Communicated!!\n";

	return(0); // exit succcesfully (hopefully)

}




