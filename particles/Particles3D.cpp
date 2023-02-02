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
/** Empty particle distribution */
void Particles3D::empty(Grid* grid,Field* EMf,VirtualTopology3D* vct){
	nop = 0;
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
void Particles3D::MaxwellianFromFluid(Grid* grid,Field* EMf,VirtualTopology3D* vct,Collective *col, int is){
#ifdef BATSRUS
  /*
   * Constuctiong the distrebution function from a Fluid model
   */

  // loop over grid cells and set position, velociy and charge of all particles indexed by counter
  // there are multiple (27 or so) particles per grid cell.
  int i,j,k,counter=0;
  for (i=1; i< grid->getNXC()-1;i++)
    for (j=1; j< grid->getNYC()-1;j++)
      for (k=1; k< grid->getNZC()-1;k++)
        MaxwellianFromFluidCell(grid,col,is, i,j,k,counter,x,y,z,q,u,v,w,ParticleID);
#endif
}

void Particles3D::MaxwellianFromFluidCell(Grid* grid, Collective *col, int is, int i, int j, int k, int &ip, double *x, double *y, double *z, double *q, double *vx, double *vy, double *vz, unsigned long* ParticleID)
{
#ifdef BATSRUS
  /*
   * grid           : local grid object (in)
   * col            : collective (global) object (in)
   * is             : species index (in)
   * i,j,k          : grid cell index on proc (in)
   * ip             : particle number counter (inout)
   * x,y,z          : particle position (out)
   * q              : particle charge (out)
   * vx,vy,vz       : particle velocity (out)
   * ParticleID     : particle tracking ID (out)
   */

  double harvest;
  double prob, theta;

  // loop over particles inside grid cell i,j,k
  for (int ii=0; ii < npcelx; ii++)
    for (int jj=0; jj < npcely; jj++)
      for (int kk=0; kk < npcelz; kk++){
        // Assign particle positions: uniformly spaced. x_cellnode + dx_particle*(0.5+index_particle)
        x[ip] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);
        y[ip] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
        z[ip] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
        // q = charge
        q[ip] =  (qom/fabs(qom))*(col->getFluidRhoCenter(i,j,k,is)/npcel)*(1.0/grid->getInvVOL());
        // u = X velocity
        harvest =   rand()/(double)RAND_MAX;
        prob  = sqrt(-2.0*log(1.0-.999999*harvest));
        harvest =   rand()/(double)RAND_MAX;
        theta = 2.0*M_PI*harvest;
        u[ip] = col->getFluidUx(i,j,k,is) + col->getFluidUthx(i,j,k,is)*prob*cos(theta);
        // v = Y velocity
        v[ip] = col->getFluidUy(i,j,k,is) + col->getFluidUthy(i,j,k,is)*prob*sin(theta);
        // w = Z velocity
        harvest =   rand()/(double)RAND_MAX;
        prob  = sqrt(-2.0*log(1.0-.999999*harvest));
        harvest =   rand()/(double)RAND_MAX;
        theta = 2.0*M_PI*harvest;
        w[ip] = col->getFluidUz(i,j,k,is) + col->getFluidUthz(i,j,k,is)*prob*cos(theta);
        if (TrackParticleID)
          ParticleID[ip]= ip*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
        ip++ ;
      }
#endif
}

/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::MaxwellianFromFields(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);

  double ebc[3];
  double vec[3];

  double harvest;
  double prob, theta, sign;
  long long counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++) {

        // Rebuild the E field from Eb = -v x B
        double Bx = EMf->getBxc(i,j,k);
        double By = EMf->getByc(i,j,k);
        double Bz = EMf->getBzc(i,j,k);

        cross_product(u0,v0,w0,Bx,By,Bz,ebc);
        scale(ebc,-1.0,3);

        // Rebuild velocity using v = Eb x B / B**2
        cross_product(ebc[0],ebc[1],ebc[2],Bx,By,Bz,vec);
        scale(vec,1.0/(Bx*Bx+By*By+Bz*Bz),3);

        double rho = fabs(EMf->getRHOcs(i, j, k, ns));

        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {
              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);

              // q = charge
              q[counter] = (qom / fabs(qom)) * (rho / npcel) * (1.0 / grid->getInvVOL());
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              u[counter] = vec[0] + uth * prob * cos(theta);
              // v
              v[counter] = vec[1] + vth * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              w[counter] = vec[2] + wth * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];


              counter++;
            }
      }


}

/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);

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
/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::maxwellian_reversed(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);
  double harvest;
  double reverser;
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

              reverser = 1.0;
              if(y[counter] < Ly/2.0 ) reverser=-1.0;

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
              w[counter] = reverser * w0 + wth * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];


              counter++;
            }


}

/** Kappa random velocity and uniform spatial distribution */
void Particles3D::kappa(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);
  double harvest;
  double reverser;
  double prob, theta, sign;
  // the following variables should be coming form input file, eventually
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


/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::maxwellian_whistler(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);
  double harvest;
  double reverser;
  double prob, theta, sign;
  double uthx, uthy, uthz;
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

              reverser = 1.0;
              if(y[counter] < Ly/2.0 ) reverser=-1.0;

              // q = charge
              q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOL());
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              double x1 = Lx/2 - Lx/5;
              double x2 = Lx/2 + Lx/5;
              double shaper_x = (tanh(x[counter]- x1)+tanh(x2-x[counter]))/2.0;
              double y1 = Ly/2 - Lx/5;
              double y2 = Ly/2 + Lx/5;
              double shaper_y = (tanh(y[counter]- y1)+tanh(y2-y[counter]))/2.0;

        	  uthx = uth;
        	  //uthy = vth * shaper_x * shaper_y + uth * (1.0 - shaper_x * shaper_y);
        	  //Used for the runs after whistler7
        	  //uthy = vth * shaper_x * shaper_y + uth / 1.1 * (1.0 - shaper_x * shaper_y);
        	  uthy = vth;
// Used for whistler 4
        	  //     	  uthy = vth * shaper_x * shaper_y + (.8*uth+.2*vth) * (1.0 - shaper_x * shaper_y);
        	  uthz = wth;

        	  /**
              if((fabs(x[counter]-Lx/2.0)<Lx/5.0)&&(fabs(y[counter]-Ly/2.0)<Lx/5.0)){
            	  uthx = uth;
            	  uthy = vth;
            	  uthz = wth;
              }
              */

              u[counter] = u0 + uthx * prob * cos(theta);
              // v
              v[counter] = v0 + uthy * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              w[counter] = reverser * w0 + uthz * prob * cos(theta);

              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];


              counter++;
            }


}

/** Maxellian velocity from currents and prescribed spatial distribution */
void Particles3D::drift_maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct)
	{
	  double harvest, prob, theta;
	  long long counter = 0;


		/* initialize random generator with different seed on different processor */
		srand(vct->getCartesian_rank()+2);

		const double q_sgn = (qom / fabs(qom));
		const double q_factor =  q_sgn / invVOL / npcel;

		for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
		for (int k=1; k< grid->getNZC()-1;k++){

			double shaper_th = fabs(EMf->getRHOcs(i, j, k, ns));

			const double qpart = q_factor * EMf->getRHOcs(i, j, k, ns);

			// determine the drift velocity from current X
			u0 = EMf->getJxs(i,j,k,ns)/EMf->getRHOns(i,j,k,ns);
			if (u0 > c){
				cout << "DRIFT VELOCITY x > c : B init field too high!" << endl;
				MPI_Abort(MPI_COMM_WORLD,2);
			}
			// determine the drift velocity from current Y
			v0 = EMf->getJys(i,j,k,ns)/EMf->getRHOns(i,j,k,ns);
			if (v0 > c){
				cout << "DRIFT VELOCITY y > c : B init field too high!" << endl;
				MPI_Abort(MPI_COMM_WORLD,2);
			}
			// determine the drift velocity from current Z
			w0 = EMf->getJzs(i,j,k,ns)/EMf->getRHOns(i,j,k,ns);
			if (w0 > c){
				cout << "DRIFT VELOCITY z > c : B init field too high!" << endl;
				MPI_Abort(MPI_COMM_WORLD,2);
			}
			for (int ii=0; ii < npcelx; ii++)
			for (int jj=0; jj < npcely; jj++)
			for (int kk=0; kk < npcelz; kk++){

	            x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
	            y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
	            z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
	            // q = charge
	            q[counter] = qpart;
				harvest = rand() / (double) RAND_MAX;
				prob = shaper_th * sqrt(-2.0 * log(1.0 - .999999 * harvest));
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
				if (TrackParticleID)
				    ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

				counter++;
			}
		}
	}

/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::relativistic_maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);

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

              //Generate speed in rest frame using Boltzmann energy distribution

              // Generate first the particle energy
              // We assume a Boltzmann energy distribution
              // Only one thermal speed is considered: kT wich
              // reinterprets vth from the input file
              // This assumes uth is sqrt(kT / m_0 c^2)
              //

              harvest = rand() / (double) RAND_MAX;
              prob = log(1.0 - .999999 * harvest);
              gammap = 1.0 - uth * uth /abs(qom) * prob;

              // Generate solid angle
              harvest = rand() / (double) RAND_MAX;
              mu = 2.0 * harvest - 1.0; //this is cos theta
              harvest = rand() / (double) RAND_MAX;
              fi = 2.0 * M_PI * harvest;


              u[counter] = sqrt((1.0-1.0/gammap/gammap)) * sqrt(1-mu*mu) * cos(fi);
              // v
              v[counter] = sqrt((1.0-1.0/gammap/gammap)) * sqrt(1-mu*mu) * sin(fi);
              // w
              w[counter] = sqrt((1.0-1.0/gammap/gammap)) * mu;
              // w


              // Boost it from rest frame to lab frame
              // This procedure is not correct in the opinion of Melzani.
              // Further analysis is needed
              // x direction
              denom = (1.0 + u0 * u[counter]/c/c);
              numerator = sqrt(1.0 - u0 * u0/c/c);
              u[counter] = (u[counter] + u0)/ denom;
              v[counter] *=  numerator/ denom;
              w[counter] *=  numerator/ denom;

              // y direction
              denom = (1.0 + v0 * v[counter]/c/c);
              numerator = sqrt(1.0 - v0 * v0/c/c);
              v[counter] = (v[counter] + v0)/ denom;
              u[counter] *=  numerator/ denom;
              w[counter] *=  numerator/ denom;

              // z direction
              denom = (1.0 + w0 * w[counter]/c/c);
              numerator = sqrt(1.0 - w0 * w0/c/c);
              w[counter] = (w[counter] + w0)/ denom;
              u[counter] *=  numerator/ denom;
              v[counter] *=  numerator/ denom;


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

/** Maxwellian random velocity and uniform spatial distribution */
int Particles3D::maxwell_box(Grid* grid,Field* EMf,VirtualTopology3D* vct, double L_square, double x_center, double y_center, double z_center, double multiple){

/* initialize random generator with different seed on different processor */
 srand(vct->getCartesian_rank()+2);

  double harvest;
  double prob, theta, sign;
  long long counter=0;
  for (int i=1; i< grid->getNXC()-1;i++)
    for (int j=1; j< grid->getNYC()-1;j++)
      for (int k=1; k< grid->getNZC()-1;k++)
        for (int ii=0; ii < npcelx; ii++)
          for (int jj=0; jj < npcely; jj++)
            for (int kk=0; kk < npcelz; kk++){
              x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);   // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
              y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
              z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
              if (fabs(x[counter] - x_center) < L_square/2 && fabs(y[counter] - y_center) < L_square/2 && fabs(z[counter] - z_center) < L_square/2) {
                // q = charg
                q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,k,ns)/npcel)*(1.0/grid->getInvVOL());
                if (ii == 0 && jj == 0 && kk == 0 && i == 1 && j == 1 && k == 1)
                  cout << "Q at ("<< i<<","<<j<<","<<k<<") = "<<q[counter]<<endl;
                // u
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[counter] = multiple*u0 + uth*prob*cos(theta);
                // v
                v[counter] = multiple*v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[counter] = multiple*w0 + wth*prob*cos(theta);
                if (TrackParticleID)
                  ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
                counter++ ;
              }
            }
  nop = counter - 1;
  return nop;
}

/** Maxellian random velocity and uniform spatial distribution */
int Particles3D::maxwell_box_thin(Grid* grid,Field* EMf,VirtualTopology3D* vct, double L_square, double x_center, double y_center, double z_center, double multiple){

/* initialize random generator with different seed on different processor */
 srand(vct->getCartesian_rank()+2);

	double harvest;
	double prob, theta, sign;
	long long counter=0;
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int k=1; k< grid->getNZC()-1;k++)
				for (int ii=0; ii < npcelx; ii++)
					for (int jj=0; jj < npcely; jj++)
						for (int kk=0; kk < npcelz; kk++){
							x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);   // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
							y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
							z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
							if (fabs(x[counter] - x_center) < L_square/4 && fabs(y[counter] - y_center) < L_square/2 && fabs(z[counter] - z_center) < L_square/2){
							// q = charg
							q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,k,ns)/npcel)*(1.0/grid->getInvVOL());
							if (ii == 0 && jj == 0 && kk == 0 && i == 1 && j == 1 && k == 1)
								cout << "Q at ("<< i<<","<<j<<","<<k<<") = "<<q[counter]<<endl;
							// u
							harvest =   rand()/(double)RAND_MAX;
							prob  = sqrt(-2.0*log(1.0-.999999*harvest));
							harvest =   rand()/(double)RAND_MAX;
							theta = 2.0*M_PI*harvest;
							u[counter] = multiple*u0 + uth*prob*cos(theta);
							// v
							v[counter] = multiple*v0 + vth*prob*sin(theta);
							// w
							harvest =   rand()/(double)RAND_MAX;
							prob  = sqrt(-2.0*log(1.0-.999999*harvest));
							harvest =   rand()/(double)RAND_MAX;
							theta = 2.0*M_PI*harvest;
							w[counter] = multiple*w0 + wth*prob*cos(theta);
							if (TrackParticleID)
								ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];


							counter++ ;
							}
						}

	nop = counter - 1;
	return nop;
}

/** Maxwellian random velocity and uniform spatial distribution */
void Particles3D::dual_spark_plug(Grid* grid,Field* EMf,VirtualTopology3D* vct, double L_square, double x_center, double y_center, double z_center){
	long long storeNop = 0;
	nop = 0;

	maxwell_box_thin(grid, EMf, vct, L_square, L_square*(3.0/4.0), y_center, z_center);
	storeNop = nop;
	storeNop += maxwell_box_thin(grid, EMf, vct, L_square, Lx-(L_square*(3.0/4.0)), y_center, z_center, -1.0);
	nop = storeNop;
}

/** Monoenergetic random velocity and uniform spatial distribution */
int Particles3D::monoenergetic_box(Grid* grid,Field* EMf,VirtualTopology3D* vct, double L_square, double x_center, double y_center, double z_center, double multiple){

/* initialize random generator with different seed on different processor */
 srand(vct->getCartesian_rank()+2);

	double harvest;
	double prob, theta, sign;
	long long counter=0;
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int k=1; k< grid->getNZC()-1;k++)
				for (int ii=0; ii < npcelx; ii++)
					for (int jj=0; jj < npcely; jj++)
						for (int kk=0; kk < npcelz; kk++){
							x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);   // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
							y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
							z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
							if (fabs(x[counter] - x_center) < L_square/2 && fabs(y[counter] - y_center) < L_square/2 && fabs(z[counter] - z_center) < L_square/2){
							// q = charg
							q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,k,ns)/npcel)*(1.0/grid->getInvVOL());
							if (ii == 0 && jj == 0 && kk == 0 && i == 1 && j == 1 && k == 1)
								cout << "Q at ("<< i<<","<<j<<","<<k<<") = "<<q[counter]<<endl;
							// u
							harvest =   rand()/(double)RAND_MAX;
							theta = 2.0*M_PI*harvest;
							u[counter] = multiple*u0 + uth*cos(theta);
							// v
							v[counter] = multiple*v0 + vth*sin(theta);
							// w
							harvest =   rand()/(double)RAND_MAX;
							theta = 2.0*M_PI*harvest;
							w[counter] = multiple*w0 + wth*cos(theta);
							if (TrackParticleID)
								ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];


							counter++ ;
							}
						}

	nop = counter - 1;
	return nop;
}


/*
int Particles3D::getGlobalFlux(int i){
    int glob = 0;

    //we should check that i is not greater than nFluxLoops
    if (i > 1)
        return -1;

    MPI_Allreduce(&fluxCounter[i], &glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return(glob);
}

double Particles3D::getGlobalFluxEnergy(int i){
    double glob = 0;

    //we should check that i is not greater than nFluxLoops
    if (i > 1)
        return -1;

    MPI_Allreduce(&fluxEnergy[i], &glob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return(glob);
}


double getGlobalFluxEnergy(int i);
int getGlobalFlux(int i);
 double fluxCounter[1];
 double fluxEnergy[1];
void Particles3D::recordFlux(double oldX, double oldY, double oldZ, double newX, double newY, double newZ, int ptcl)
{
	double interpX, interpY, interpZ, centerX, centerY, centerZ;
	centerX = Lx/2.0;
	centerY = Ly/2.0;
	centerZ = Lz/2.0;

	int nFluxLoops = 1;
	int fluxPlane[1] = {1};
	double fluxCenterX[1], fluxCenterY[1], fluxCenterZ[1];
	fluxCenterY[0] = fluxCenterZ[0] = centerY;
	fluxCenterX[0] = centerX+30;
	double fluxRadius[1] = {12};

	for (int i = 0; i < nFluxLoops; i++)
	{

			switch(fluxPlane[i]){
				case 1:				//the case where the fluxLoop is normal to the x-plane

					//interpolate the coordinates up to the position of the flux-loop
					interpY = oldY + (newY-oldY)*(fluxCenterX[i]-oldX)/(newX-oldX);
					interpZ = oldZ + (newZ-oldZ)*(fluxCenterX[i]-oldX)/(newX-oldX);

					//check if the particle has gone from the center of the machine, though the loop, to the outside of the machine
					if (oldX-centerX < fluxCenterX[i]-centerX && newX-centerX > fluxCenterX[i]-centerX &&
						pow(interpY-fluxCenterY[i], 2.0) + pow(interpZ-fluxCenterZ[i], 2.0) < pow(fluxRadius[i], 2.0) )
					{
						fluxCounter[i]++;
						fluxEnergy[i] += 0.5*getQ(ptcl)*(pow(getU(ptcl), 2.0)+pow(getV(ptcl), 2.0)+pow(getW(ptcl), 2.0))/qom;
					}

					//check if the particle has gone from outside the center of the machine, though the loop, to the inside of the machine
					else if (oldX-centerX > fluxCenterX[i]-centerX && newX-centerX < fluxCenterX[i]-centerX &&
						pow(interpY-fluxCenterY[i], 2.0) + pow(interpZ-fluxCenterZ[i], 2.0) < pow(fluxRadius[i], 2.0) )
					{
						fluxCounter[i]--;
						fluxEnergy[i] -= 0.5*getQ(ptcl)*(pow(getU(ptcl), 2.0)+pow(getV(ptcl), 2.0)+pow(getW(ptcl), 2.0))/qom;
					}
				break;

				case 2:				//the case where the fluxLoop is normal to the y-plane

					//interpolate the coordinates up to the position of the flux-loop
					interpX = oldX + (newX-oldX)*(fluxCenterY[i]-oldY)/(newY-oldY);
					interpZ = oldZ + (newZ-oldZ)*(fluxCenterY[i]-oldY)/(newY-oldY);

					//check if the particle has gone from the center of the machine, though the loop, to the outside of the machine
					if (oldY-centerY < fluxCenterY[i]-centerY && newY-centerY > fluxCenterY[i]-centerY &&
						pow(interpX-fluxCenterX[i], 2.0) + pow(interpZ-fluxCenterZ[i], 2.0) < pow(fluxRadius[i], 2.0) )
					{
						fluxCounter[i]++;
						fluxEnergy[i] += 0.5*getQ(ptcl)*(pow(getU(ptcl), 2.0)+pow(getV(ptcl), 2.0)+pow(getW(ptcl), 2.0))/qom;
					}

					//check if the particle has gone from outside the center of the machine, though the loop, to the inside of the machine
					else if (oldY-centerY > fluxCenterY[i]-centerY && newY-centerY < fluxCenterY[i]-centerY &&
						pow(interpX-fluxCenterX[i], 2.0) + pow(interpZ-fluxCenterZ[i], 2.0) < pow(fluxRadius[i], 2.0) )
					{
						fluxCounter[i]--;
						fluxEnergy[i] -= 0.5*getQ(ptcl)*(pow(getU(ptcl), 2.0)+pow(getV(ptcl), 2.0)+pow(getW(ptcl), 2.0))/qom;
					}
				break;

				case 3:				//the case where the fluxLoop is normal to the z-plane

					//interpolate the coordinates up to the position of the flux-loop
					interpX = oldX + (newX-oldX)*(fluxCenterY[i]-oldY)/(newY-oldY);
					interpY = oldZ + (newZ-oldZ)*(fluxCenterY[i]-oldY)/(newY-oldY);

					//check if the particle has gone from the center of the machine, though the loop, to the outside of the machine
					if (oldZ-centerZ < fluxCenterZ[i]-centerZ && newZ-centerZ > fluxCenterZ[i]-centerZ &&
						pow(interpX-fluxCenterX[i], 2.0) + pow(interpY-fluxCenterY[i], 2.0) < pow(fluxRadius[i], 2.0) )
					{
						fluxCounter[i]++;
						fluxEnergy[i] += 0.5*getQ(ptcl)*(pow(getU(ptcl), 2.0)+pow(getV(ptcl), 2.0)+pow(getW(ptcl), 2.0))/qom;
					}

					//check if the particle has gone from outside the center of the machine, though the loop, to the inside of the machine
					else if (oldZ-centerZ > fluxCenterZ[i]-centerZ && newZ-centerZ < fluxCenterZ[i]-centerZ &&
						pow(interpX-fluxCenterX[i], 2.0) + pow(interpY-fluxCenterY[i], 2.0) < pow(fluxRadius[i], 2.0) )
					{
						fluxCounter[i]--;
						fluxEnergy[i] -= 0.5*getQ(ptcl)*(pow(getU(ptcl), 2.0)+pow(getV(ptcl), 2.0)+pow(getW(ptcl), 2.0))/qom;
					}
				break;

				default:
				break;
		}
	}

}

*/

/**Add a periodic perturbation in J exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
void Particles3D::AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0, Grid * grid) {

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

/** mover with a Predictor-Corrector scheme */
int Particles3D::mover_PC_old(Grid* grid,VirtualTopology3D* vct, Field* EMf){

  int avail;
  double dto2 = .5*dt, qomdt2 = qom*dto2/c;
  double omdtsq[P_SAME_TIME], denom[P_SAME_TIME], ut[P_SAME_TIME], vt[P_SAME_TIME], wt[P_SAME_TIME], udotb[P_SAME_TIME];
  double Exl[P_SAME_TIME], Eyl[P_SAME_TIME], Ezl[P_SAME_TIME], Bxl[P_SAME_TIME], Byl[P_SAME_TIME], Bzl[P_SAME_TIME];
  double Exlp[P_SAME_TIME], Eylp[P_SAME_TIME], Ezlp[P_SAME_TIME], Bxlp[P_SAME_TIME], Bylp[P_SAME_TIME], Bzlp[P_SAME_TIME];
  double xptilde[P_SAME_TIME], yptilde[P_SAME_TIME], zptilde[P_SAME_TIME], uptilde[P_SAME_TIME], vptilde[P_SAME_TIME], wptilde[P_SAME_TIME];
  double xp[P_SAME_TIME], yp[P_SAME_TIME], zp[P_SAME_TIME], up[P_SAME_TIME], vp[P_SAME_TIME], wp[P_SAME_TIME];
  double ixd[P_SAME_TIME], iyd[P_SAME_TIME], izd[P_SAME_TIME];
  int ix[P_SAME_TIME], iy[P_SAME_TIME], iz[P_SAME_TIME];
  double weight[P_SAME_TIME][2][2][2];
  double xi[2]; double eta[2]; double zeta[2];
  double inv_dx = 1.0/dx, inv_dy = 1.0/dy, inv_dz = 1.0/dz;
  double Fext = EMf->getFext();
  
  if(Gravity) Fext /= qom;
  
  if (vct->getCartesian_rank()==0) {
    cout << "*** MOVER species " << ns << " ***" << NiterMover <<" ITERATIONS   ****" << "Fext=" << Fext << endl;
  }
  
  // move each particle with new fields: MOVE P_SAME_TIME PARTICLES AT THE SAME TIME TO ALLOW AUTOVECTORIZATION
  int i;
  for (i=0; i <  (nop-(P_SAME_TIME-1)); i+=P_SAME_TIME) {
    // copy x, y, z
    for (int p = 0; p < P_SAME_TIME; p++) {
      xp[p] = x[i+p];   yp[p] = y[i+p];   zp[p] = z[i+p];
      up[p] = u[i+p];   vp[p] = v[i+p];   wp[p] = w[i+p];
    }
    for (int p = 0; p < P_SAME_TIME; p++) { // VECTORIZED
      xptilde[p] = xp[p]; yptilde[p] = yp[p]; zptilde[p] = zp[p];
    }
    // calculate the average velocity iteratively
    for (int innter=0; innter < NiterMover; innter++) {
      // interpolation G-->P
      for(int p=0; p < P_SAME_TIME; p++) ixd[p] = floor((xp[p] - xstart)*inv_dx);  // VECTORIZED
      for(int p=0; p < P_SAME_TIME; p++) iyd[p] = floor((yp[p] - ystart)*inv_dy);  // VECTORIZED
      for(int p=0; p < P_SAME_TIME; p++) izd[p] = floor((zp[p] - zstart)*inv_dz);  // VECTORIZED
      for(int p=0; p < P_SAME_TIME; p++){ ix[p] = 2 + int(ixd[p]); iy[p] = 2 + int(iyd[p]); iz[p] = 2 + int(izd[p]);}
      // check if they are out of the boundary
      for(int p=0; p < P_SAME_TIME; p++){ if(ix[p] < 1) ix[p] = 1;}
      for(int p=0; p < P_SAME_TIME; p++){ if(iy[p] < 1) iy[p] = 1;}
      for(int p=0; p < P_SAME_TIME; p++){ if(iz[p] < 1) iz[p] = 1;}
      for(int p=0; p < P_SAME_TIME; p++){ if(ix[p] > nxn-1) ix[p] = nxn-1;}
      for(int p=0; p < P_SAME_TIME; p++){ if(iy[p] > nyn-1) iy[p] = nyn-1;}
      for(int p=0; p < P_SAME_TIME; p++){ if(iz[p] > nzn-1) iz[p] = nzn-1;}
      
      // CALCULATE WEIGHTS
      for(int p=0; p < P_SAME_TIME; p++){
        xi[0]   = xp[p] - grid->getXN(ix[p]-1,iy[p],iz[p]); eta[0]  = yp[p] - grid->getYN(ix[p],iy[p]-1,iz[p]); zeta[0] = zp[p] - grid->getZN(ix[p],iy[p],iz[p]-1);
        xi[1]   = grid->getXN(ix[p],iy[p],iz[p]) - xp[p];   eta[1]  = grid->getYN(ix[p],iy[p],iz[p]) - yp[p];   zeta[1] = grid->getZN(ix[p],iy[p],iz[p]) - zp[p];
        for (int ii=0; ii < 2; ii++)
          for (int jj=0; jj < 2; jj++)
            for(int kk=0; kk < 2; kk++)
        	    weight[p][ii][jj][kk] = xi[ii]*eta[jj]*zeta[kk]*invVOL;
      }
      // clear the electric and the magnetic field field acting on the particles
      for (int p = 0; p < P_SAME_TIME; p++){Exl[p]=0.0; Eyl[p] = 0.0; Ezl[p] = 0.0; Bxl[p] = 0.0; Byl[p] = 0.0; Bzl[p] = 0.0;}
      // calculate fields acting on the particles
      for (int ii=0; ii < 2; ii++)
        for (int jj=0; jj < 2; jj++)
      	for(int kk=0; kk < 2; kk++) {
            for(int p = 0; p < P_SAME_TIME; p++) {
              Exlp[p] = EMf->getEx(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              Eylp[p] = EMf->getEy(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              Ezlp[p] = EMf->getEz(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              Bxlp[p] = EMf->getBx(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              Bylp[p] = EMf->getBy(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              Bzlp[p] = EMf->getBz(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              // Add external fields
              Exlp[p] += Fext * EMf->getEx_ext(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              Eylp[p] += Fext * EMf->getEy_ext(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              Ezlp[p] += Fext * EMf->getEz_ext(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              Bxlp[p] += Fext * EMf->getBx_ext(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              Bylp[p] += Fext * EMf->getBy_ext(ix[p] - ii,iy[p] -jj,iz[p] - kk);
              Bzlp[p] += Fext * EMf->getBz_ext(ix[p] - ii,iy[p] -jj,iz[p] - kk);
            }
            for(int p = 0; p < P_SAME_TIME; p++) { // VECTORIZED
              Exlp[p] = weight[p][ii][jj][kk]*Exlp[p];
              Eylp[p] = weight[p][ii][jj][kk]*Eylp[p];
              Ezlp[p] = weight[p][ii][jj][kk]*Ezlp[p];
              Bxlp[p] = weight[p][ii][jj][kk]*Bxlp[p];
              Bylp[p] = weight[p][ii][jj][kk]*Bylp[p];
              Bzlp[p] = weight[p][ii][jj][kk]*Bzlp[p];
            }
            // finished with the two particles: add the contributions
            for(int p = 0; p < P_SAME_TIME; p++){
              Exl[p] += Exlp[p]; Eyl[p] += Eylp[p]; Ezl[p] += Ezlp[p];
              Bxl[p] += Bxlp[p]; Byl[p] += Bylp[p]; Bzl[p] += Bzlp[p];
            }
          }
      // end interpolation
      for(int p = 0; p < P_SAME_TIME; p++) { // PARTIALLY VECTORIZED
        omdtsq[p] = qomdt2*qomdt2*(Bxl[p]*Bxl[p]+Byl[p]*Byl[p]+Bzl[p]*Bzl[p]);
        denom[p] = 1.0/(1.0 + omdtsq[p]);
        // solve the position equation
        ut[p]= up[p] + qomdt2*Exl[p];
        vt[p]= vp[p] + qomdt2*Eyl[p];
        wt[p]= wp[p] + qomdt2*Ezl[p];
        udotb[p] = ut[p]*Bxl[p] + vt[p]*Byl[p] + wt[p]*Bzl[p];
        // solve the velocity equation
        uptilde[p] = (ut[p]+qomdt2*(vt[p]*Bzl[p] -wt[p]*Byl[p] + qomdt2*udotb[p]*Bxl[p]))*denom[p];
        vptilde[p] = (vt[p]+qomdt2*(wt[p]*Bxl[p] -ut[p]*Bzl[p] + qomdt2*udotb[p]*Byl[p]))*denom[p];
        wptilde[p] = (wt[p]+qomdt2*(ut[p]*Byl[p] -vt[p]*Bxl[p] + qomdt2*udotb[p]*Bzl[p]))*denom[p];
        // update position
        xp[p] = xptilde[p] + uptilde[p]*dto2;
        yp[p] = yptilde[p] + vptilde[p]*dto2;
        zp[p] = zptilde[p] + wptilde[p]*dto2;
      }
    } // end of iteration

    // update the final position and velocity
    for(int p=0; p < P_SAME_TIME; p++) { // VECTORIZED
       up[p]= 2.0*uptilde[p] - up[p];
       vp[p]= 2.0*vptilde[p] - vp[p];
       wp[p]= 2.0*wptilde[p] - wp[p];
       xp[p] = xptilde[p] + uptilde[p]*dt;
       yp[p] = yptilde[p] + vptilde[p]*dt;
       zp[p] = zptilde[p] + wptilde[p]*dt;
    }
    // copy back the particles in the array
    for(int p=0; p < P_SAME_TIME; p++){
      x[i+p]   = xp[p]; y[i+p]   = yp[p]; z[i+p]   = zp[p];
      u[i+p]   = up[p]; v[i+p]   = vp[p]; w[i+p]   = wp[p];
    }
  }

  // FINISH WITH PARTICLE LEFT, IF ANY (EVEN NUMBER OF PARTICLES)
  // move each particle with new fields
  for (int rest = (i+1); rest <  nop; rest++){
      // copy the particle
  	xp[0] = x[rest];   yp[0] = y[rest];   zp[0] = z[rest];   up[0] = u[rest];   vp[0] = v[rest];   wp[0] = w[rest];
  	xptilde[0] = x[rest];
  	yptilde[0] = y[rest];
  	zptilde[0] = z[rest];
  	// calculate the average velocity iteratively
  	for(int innter=0; innter < 1; innter++){
  		// interpolation G-->P
  		ixd[0] = floor((xp[0]-xstart)*inv_dx);
  		iyd[0] = floor((yp[0]-ystart)*inv_dy);
  		izd[0] = floor((zp[0]-zstart)*inv_dz);
  		ix[0] = 2 +  int(ixd[0]);
  		iy[0] = 2 +  int(iyd[0]);
  		iz[0] = 2 +  int(izd[0]);
  		if(ix[0] < 1) ix[0] = 1;
  		if(iy[0] < 1) iy[0] = 1;
  		if(iz[0] < 1) iz[0] = 1;
  		if(ix[0] > nxn-1) ix[0] = nxn-1;
  		if(iy[0] > nyn-1) iy[0] = nyn-1;
  		if(iz[0] > nzn-1) iz[0] = nzn-1;
  
  		xi[0]   = xp[0] - grid->getXN(ix[0]-1,iy[0],iz[0]); eta[0]  = yp[0] - grid->getYN(ix[0],iy[0]-1,iz[0]); zeta[0] = zp[0] - grid->getZN(ix[0],iy[0],iz[0]-1);
      xi[1]   = grid->getXN(ix[0],iy[0],iz[0]) - xp[0];   eta[1]  = grid->getYN(ix[0],iy[0],iz[0]) - yp[0];   zeta[1] = grid->getZN(ix[0],iy[0],iz[0]) - zp[0];
      for (int ii=0; ii < 2; ii++)
  	      for (int jj=0; jj < 2; jj++)
  		     for(int kk=0; kk < 2; kk++)
  			   weight[0][ii][jj][kk] = xi[ii]*eta[jj]*zeta[kk]*invVOL;
  
  	Exl[0]=0.0, Eyl[0] = 0.0, Ezl[0] = 0.0, Bxl[0] = 0.0, Byl[0] = 0.0, Bzl[0] = 0.0;
      for (int ii=0; ii < 2; ii++)
  			for (int jj=0; jj < 2; jj++)
  				for(int kk=0; kk < 2; kk++){
  					Exlp[0] = weight[0][ii][jj][kk]*EMf->getEx(ix[0] - ii,iy[0] -jj,iz[0]- kk );
  					Eylp[0] = weight[0][ii][jj][kk]*EMf->getEy(ix[0] - ii,iy[0] -jj,iz[0]- kk );
  					Ezlp[0] = weight[0][ii][jj][kk]*EMf->getEz(ix[0] - ii,iy[0] -jj,iz[0] -kk );
  					Bxlp[0] = weight[0][ii][jj][kk]*EMf->getBx(ix[0] - ii,iy[0] -jj,iz[0] -kk );
  					Bylp[0] = weight[0][ii][jj][kk]*EMf->getBy(ix[0] - ii,iy[0] -jj,iz[0] -kk );
  					Bzlp[0] = weight[0][ii][jj][kk]*EMf->getBz(ix[0] - ii,iy[0] -jj,iz[0] -kk );
  
  					Exlp[0] += Fext * weight[0][ii][jj][kk]*EMf->getEx_ext(ix[0] - ii,iy[0] -jj,iz[0]- kk );
  				    Eylp[0] += Fext * weight[0][ii][jj][kk]*EMf->getEy_ext(ix[0] - ii,iy[0] -jj,iz[0]- kk );
  					Ezlp[0] += Fext * weight[0][ii][jj][kk]*EMf->getEz_ext(ix[0] - ii,iy[0] -jj,iz[0] -kk );
  					Bxlp[0] += Fext * weight[0][ii][jj][kk]*EMf->getBx_ext(ix[0] - ii,iy[0] -jj,iz[0] -kk );
  					Bylp[0] += Fext * weight[0][ii][jj][kk]*EMf->getBy_ext(ix[0] - ii,iy[0] -jj,iz[0] -kk );
  					Bzlp[0] += Fext * weight[0][ii][jj][kk]*EMf->getBz_ext(ix[0] - ii,iy[0] -jj,iz[0] -kk );
  
  					Exl[0] += Exlp[0];
  					Eyl[0] += Eylp[0];
  					Ezl[0] += Ezlp[0];
  					Bxl[0] += Bxlp[0];
  					Byl[0] += Bylp[0];
  					Bzl[0] += Bzlp[0];
  				}
  		// end interpolation
      	omdtsq[0] = qomdt2*qomdt2*(Bxl[0]*Bxl[0]+Byl[0]*Byl[0]+Bzl[0]*Bzl[0]);
  		denom[0] = 1.0/(1.0 + omdtsq[0]);
  		// solve the position equation
  		ut[0]= up[0] + qomdt2*Exl[0];
  		vt[0]= vp[0] + qomdt2*Eyl[0];
  		wt[0]= wp[0] + qomdt2*Ezl[0];
  		udotb[0] = ut[0]*Bxl[0] + vt[0]*Byl[0] + wt[0]*Bzl[0];
  		// solve the velocity equation
  		uptilde[0] = (ut[0]+qomdt2*(vt[0]*Bzl[0] -wt[0]*Byl[0] + qomdt2*udotb[0]*Bxl[0]))*denom[0];
  		vptilde[0] = (vt[0]+qomdt2*(wt[0]*Bxl[0] -ut[0]*Bzl[0] + qomdt2*udotb[0]*Byl[0]))*denom[0];
  		wptilde[0] = (wt[0]+qomdt2*(ut[0]*Byl[0] -vt[0]*Bxl[0] + qomdt2*udotb[0]*Bzl[0]))*denom[0];
  		// update position
  		xp[0] = xptilde[0] + uptilde[0]*dto2;
  		yp[0] = yptilde[0] + vptilde[0]*dto2;
  		zp[0] = zptilde[0] + wptilde[0]*dto2;
  	} // end of iteration
  	// update the final position and velocity
  	up[0]= 2.0*uptilde[0] - u[rest];
  	vp[0]= 2.0*vptilde[0] - v[rest];
  	wp[0]= 2.0*wptilde[0] - w[rest];
  	xp[0] = xptilde[0] + uptilde[0]*dt;
  	yp[0] = yptilde[0] + vptilde[0]*dt;
  	zp[0] = zptilde[0] + wptilde[0]*dt;
  	x[rest]   = xp[0]; y[rest]   = yp[0]; z[rest]   = zp[0]; u[rest]   = up[0]; v[rest]   = vp[0]; w[rest]   = wp[0];
  } // END OF ALL THE PARTICLES
  
  //********************//
  // COMMUNICATION
  // *******************//
    avail = communicate(vct);
    if (avail < 0)
		return(-1);
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

/** mover with a Predictor-Corrector scheme */
int Particles3D::mover_PC(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER species " << ns << " ***" << " with " << nop << " particles ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  double start_mover_PC = MPI_Wtime();
  double ***Ex = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEx());
  double ***Ey = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEy());
  double ***Ez = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEz());
  double ***Ex_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEx_ext());
  double ***Ey_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEy_ext());
  double ***Ez_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEz_ext());
  double ***Bx = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx());
  double ***By = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy());
  double ***Bz = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz());
  double ***Bx_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx_ext());
  double ***By_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy_ext());
  double ***Bz_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz_ext());
  double Fext = EMf->getFext();
  
  if(Gravity) Fext /= qom;
  
  const double dto2 = .5 * dt, qomdt2 = qom * dto2 / c;
  const double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy, inv_dz = 1.0 / dz;
  // don't bother trying to push any particles simultaneously;
  // MIC already does vectorization automatically, and trying
  // to do it by hand only hurts performance.
#pragma omp parallel for
#pragma simd                    // this just slows things down (why?)
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
      
      // MIC refuses to vectorize this ...
      //
      // double weight[2][2][2];
      // for (int ii = 0; ii < 2; ii++)
      // for (int jj = 0; jj < 2; jj++)
      // for (int kk = 0; kk < 2; kk++)
      // weight[ii][jj][kk] = xi[ii] * eta[jj] * zeta[kk] * invVOL;
      // for (int ii = 0; ii < 2; ii++)
      // for (int jj = 0; jj < 2; jj++)
      // for (int kk = 0; kk < 2; kk++) {
      // const double Exlp = weight[ii][jj][kk] * Ex.get(ix - ii, iy - jj, iz - kk);
      // const double Eylp = weight[ii][jj][kk] * Ey.get(ix - ii, iy - jj, iz - kk);
      // const double Ezlp = weight[ii][jj][kk] * Ez.get(ix - ii, iy - jj, iz - kk);
      // const double Bxlp = weight[ii][jj][kk] * Bx.get(ix - ii, iy - jj, iz - kk);
      // const double Bylp = weight[ii][jj][kk] * By.get(ix - ii, iy - jj, iz - kk);
      // const double Bzlp = weight[ii][jj][kk] * Bz.get(ix - ii, iy - jj, iz - kk);
      // Exl += Exlp;
      // Eyl += Eylp;
      // Ezl += Ezlp;
      // Bxl += Bxlp;
      // Byl += Bylp;
      // Bzl += Bzlp;
      // }
      
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
      //
      Exl += weight000 * (Ex[ix][iy][iz] 			 + Fext * Ex_ext[ix][iy][iz]);
      Exl += weight001 * (Ex[ix][iy][iz - 1] 		 + Fext * Ex_ext[ix][iy][iz - 1]);
      Exl += weight010 * (Ex[ix][iy - 1][iz] 		 + Fext * Ex_ext[ix][iy - 1][iz]);
      Exl += weight011 * (Ex[ix][iy - 1][iz - 1] 	 + Fext * Ex_ext[ix][iy - 1][iz - 1]);
      Exl += weight100 * (Ex[ix - 1][iy][iz] 		 + Fext * Ex_ext[ix - 1][iy][iz]);
      Exl += weight101 * (Ex[ix - 1][iy][iz - 1] 	 + Fext * Ex_ext[ix - 1][iy][iz - 1]);
      Exl += weight110 * (Ex[ix - 1][iy - 1][iz] 	 + Fext * Ex_ext[ix - 1][iy - 1][iz]);
      Exl += weight111 * (Ex[ix - 1][iy - 1][iz - 1] + Fext * Ex_ext[ix - 1][iy - 1][iz - 1]);
      //
      Eyl += weight000 * (Ey[ix][iy][iz] 			 + Fext * Ey_ext[ix][iy][iz]);
      Eyl += weight001 * (Ey[ix][iy][iz - 1] 		 + Fext * Ey_ext[ix][iy][iz - 1]);
      Eyl += weight010 * (Ey[ix][iy - 1][iz] 		 + Fext * Ey_ext[ix][iy - 1][iz]);
      Eyl += weight011 * (Ey[ix][iy - 1][iz - 1] 	 + Fext * Ey_ext[ix][iy - 1][iz - 1]);
      Eyl += weight100 * (Ey[ix - 1][iy][iz] 		 + Fext * Ey_ext[ix - 1][iy][iz]);
      Eyl += weight101 * (Ey[ix - 1][iy][iz - 1] 	 + Fext * Ey_ext[ix - 1][iy][iz - 1]);
      Eyl += weight110 * (Ey[ix - 1][iy - 1][iz] 	 + Fext * Ey_ext[ix - 1][iy - 1][iz]);
      Eyl += weight111 * (Ey[ix - 1][iy - 1][iz - 1] + Fext * Ey_ext[ix - 1][iy - 1][iz - 1]);
      
      Ezl += weight000 * (Ez[ix][iy][iz] 			 + Fext * Ez_ext[ix][iy][iz]);
      Ezl += weight001 * (Ez[ix][iy][iz - 1] 		 + Fext * Ez_ext[ix][iy][iz - 1]);
      Ezl += weight010 * (Ez[ix][iy - 1][iz] 		 + Fext * Ez_ext[ix][iy - 1][iz]);
      Ezl += weight011 * (Ez[ix][iy - 1][iz - 1] 	 + Fext * Ez_ext[ix][iy - 1][iz - 1]);
      Ezl += weight100 * (Ez[ix - 1][iy][iz] 		 + Fext * Ez_ext[ix - 1][iy][iz]);
      Ezl += weight101 * (Ez[ix - 1][iy][iz - 1] 	 + Fext * Ez_ext[ix - 1][iy][iz - 1]);
      Ezl += weight110 * (Ez[ix - 1][iy - 1][iz] 	 + Fext * Ez_ext[ix - 1][iy - 1][iz]);
      Ezl += weight111 * (Ez[ix - 1][iy - 1][iz - 1] + Fext * Ez_ext[ix - 1][iy - 1][iz - 1]);
      
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
  // timeTasks.start_communicate();
  const int avail = communicate(vct);
  if (avail < 0) return (-1);
  MPI_Barrier(MPI_COMM_WORLD);
  // communicate again if particles are not in the correct domain
  while (isMessagingDone(vct) > 0) {
    // COMMUNICATION
    const int avail = communicate(vct);
    if (avail < 0) return (-1);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // timeTasks.addto_communicate();
  return (0);                   // exit succcesfully (hopefully)
}

/** relativistic mover with a Predictor-Corrector scheme */
int Particles3D::mover_PC_rel(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER species " << ns << " ***" << " with " << nop << " particles ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  double start_mover_PC = MPI_Wtime();
  double ***Ex = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEx());
  double ***Ey = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEy());
  double ***Ez = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEz());
  double ***Ex_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEx_ext());
  double ***Ey_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEy_ext());
  double ***Ez_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getEz_ext());
  double ***Bx = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx());
  double ***By = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy());
  double ***Bz = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz());
  double ***Bx_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBx_ext());
  double ***By_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBy_ext());
  double ***Bz_ext = asgArr3(double, grid->getNXN(), grid->getNYN(), grid->getNZN(), EMf->getBz_ext());
  double Fext = EMf->getFext();
  
  if(Gravity) Fext /= qom;
  
  const double dto2 = .5 * dt, qomdt2 = qom * dto2 / c;
  const double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy, inv_dz = 1.0 / dz;
  // don't bother trying to push any particles simultaneously;
  // MIC already does vectorization automatically, and trying
  // to do it by hand only hurts performance.
#pragma omp parallel for
#pragma simd                    // this just slows things down (why?)
  for (long long rest = 0; rest < nop; rest++) {
    // copy the particle
    double xp = x[rest];
    double yp = y[rest];
    double zp = z[rest];
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];
    double gp = sqrt(1.+up*up+vp*vp+wp*wp);
    const double xpold = xp;
    const double ypold = yp;
    const double zpold = zp;
    const double upold = up;
    const double vpold = vp;
    const double wpold = wp;
    const double gpold = gp;
    double uxbar, uybar, uzbar, gbar;

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
      //
      Exl += weight000 * (Ex[ix][iy][iz] 			 + Fext * Ex_ext[ix][iy][iz]);
      Exl += weight001 * (Ex[ix][iy][iz - 1] 		 + Fext * Ex_ext[ix][iy][iz - 1]);
      Exl += weight010 * (Ex[ix][iy - 1][iz] 		 + Fext * Ex_ext[ix][iy - 1][iz]);
      Exl += weight011 * (Ex[ix][iy - 1][iz - 1] 	 + Fext * Ex_ext[ix][iy - 1][iz - 1]);
      Exl += weight100 * (Ex[ix - 1][iy][iz] 		 + Fext * Ex_ext[ix - 1][iy][iz]);
      Exl += weight101 * (Ex[ix - 1][iy][iz - 1] 	 + Fext * Ex_ext[ix - 1][iy][iz - 1]);
      Exl += weight110 * (Ex[ix - 1][iy - 1][iz] 	 + Fext * Ex_ext[ix - 1][iy - 1][iz]);
      Exl += weight111 * (Ex[ix - 1][iy - 1][iz - 1] + Fext * Ex_ext[ix - 1][iy - 1][iz - 1]);
      //
      Eyl += weight000 * (Ey[ix][iy][iz] 			 + Fext * Ey_ext[ix][iy][iz]);
      Eyl += weight001 * (Ey[ix][iy][iz - 1] 		 + Fext * Ey_ext[ix][iy][iz - 1]);
      Eyl += weight010 * (Ey[ix][iy - 1][iz] 		 + Fext * Ey_ext[ix][iy - 1][iz]);
      Eyl += weight011 * (Ey[ix][iy - 1][iz - 1] 	 + Fext * Ey_ext[ix][iy - 1][iz - 1]);
      Eyl += weight100 * (Ey[ix - 1][iy][iz] 		 + Fext * Ey_ext[ix - 1][iy][iz]);
      Eyl += weight101 * (Ey[ix - 1][iy][iz - 1] 	 + Fext * Ey_ext[ix - 1][iy][iz - 1]);
      Eyl += weight110 * (Ey[ix - 1][iy - 1][iz] 	 + Fext * Ey_ext[ix - 1][iy - 1][iz]);
      Eyl += weight111 * (Ey[ix - 1][iy - 1][iz - 1] + Fext * Ey_ext[ix - 1][iy - 1][iz - 1]);
      
      Ezl += weight000 * (Ez[ix][iy][iz] 			 + Fext * Ez_ext[ix][iy][iz]);
      Ezl += weight001 * (Ez[ix][iy][iz - 1] 		 + Fext * Ez_ext[ix][iy][iz - 1]);
      Ezl += weight010 * (Ez[ix][iy - 1][iz] 		 + Fext * Ez_ext[ix][iy - 1][iz]);
      Ezl += weight011 * (Ez[ix][iy - 1][iz - 1] 	 + Fext * Ez_ext[ix][iy - 1][iz - 1]);
      Ezl += weight100 * (Ez[ix - 1][iy][iz] 		 + Fext * Ez_ext[ix - 1][iy][iz]);
      Ezl += weight101 * (Ez[ix - 1][iy][iz - 1] 	 + Fext * Ez_ext[ix - 1][iy][iz - 1]);
      Ezl += weight110 * (Ez[ix - 1][iy - 1][iz] 	 + Fext * Ez_ext[ix - 1][iy - 1][iz]);
      Ezl += weight111 * (Ez[ix - 1][iy - 1][iz - 1] + Fext * Ez_ext[ix - 1][iy - 1][iz - 1]);
      
      // end interpolation
      double epsx = qomdt2*Exl;
      double epsy = qomdt2*Eyl;
      double epsz = qomdt2*Ezl;
      double betax = qomdt2*Bxl;
      double betay = qomdt2*Byl;
      double betaz = qomdt2*Bzl;
      double beta2 = betax*betax+betay*betay+betaz*betaz;
      double upx = upold + epsx;
      double upy = vpold + epsy;
      double upz = wpold + epsz;

      // Polynomial coefficients
      double updote = upx*epsx+upy*epsy+upz*epsz;
      double bdote = betax*epsx+betay*epsy+betaz*epsz;
      double updotb = upx*betax+upy*betay+upz*betaz;
      double upcrossb_x = (upy*betaz-upz*betay);
      double upcrossb_y = (-upx*betaz+upz*betax);
      double upcrossb_z = (upx*betay-upy*betax);
      double aa = updote - beta2;
      double bb = upcrossb_x*epsx+upcrossb_y*epsy+upcrossb_z*epsz+ gpold*beta2;
      double cc = updotb*bdote;
      
      // Solution coefficients
      double AA = 2.*aa/3.+gpold*gpold/4.;
      double BB = 4.*aa*gpold+8.*bb+gpold*gpold*gpold;
      double DD = aa*aa-3.*bb*gpold-12.*cc;
      double FF = -2.*aa*aa*aa+9.*aa*bb*gpold-72.*aa*cc+27.*bb*bb-27.*cc*gpold*gpold;
      std::complex<double> GG = FF*FF-4.*DD*DD*DD;
      std::complex<double> EE;
      if (std::real((FF+sqrt(GG))/2.)<0.) EE = -pow(-(FF+sqrt(GG))/2.,1./3.);
      else EE = pow((FF+sqrt(GG))/2.,1./3.);
      std::complex<double> CC = DD/(EE+1.e-20)/3.+EE/3.;
      // Solution
      std::complex<double> gbarc = gpold/4.+sqrt(2.*AA+BB/4./sqrt(AA+CC+1.e-20)-CC)/2.+sqrt(AA+CC)/2.;
      gbar = (double) std::real(gbarc);
      
      uxbar = (upx+(upx*betax+upy*betay+upz*betaz)*betax/(gbar*gbar)+(upy*betaz-upz*betay)/gbar)/(1.+beta2/gbar/gbar);
      uybar = (upy+(upx*betax+upy*betay+upz*betaz)*betay/(gbar*gbar)+(-upx*betaz+upz*betax)/gbar)/(1.+beta2/gbar/gbar);
      uzbar = (upz+(upx*betax+upy*betay+upz*betaz)*betaz/(gbar*gbar)+(upx*betay-upy*betax)/gbar)/(1.+beta2/gbar/gbar);
      
      // update position
      xp = xpold + uxbar/gbar * dto2;
      yp = ypold + uybar/gbar * dto2;
      zp = zpold + uzbar/gbar * dto2;
    } // end of PC iteration

    // update the final velocity and position
    up = 2.*uxbar - upold;
    vp = 2.*uybar - vpold;
    wp = 2.*uzbar - wpold;
    xp = xpold + uxbar/gbar * dt;
    yp = ypold + uybar/gbar * dt;
    zp = zpold + uzbar/gbar * dt;
    x[rest] = xp;
    y[rest] = yp;
    z[rest] = zp;
    u[rest] = up;
    v[rest] = vp;
    w[rest] = wp;
  } // END OF ALL THE PARTICLES
  
  // ********************//
  // COMMUNICATION
  // *******************//
  // timeTasks.start_communicate();
  const int avail = communicate(vct);
  if (avail < 0) return (-1);
  MPI_Barrier(MPI_COMM_WORLD);
  // communicate again if particles are not in the correct domain
  while (isMessagingDone(vct) > 0) {
    // COMMUNICATION
    const int avail = communicate(vct);
    if (avail < 0) return (-1);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // timeTasks.addto_communicate();
  return (0);                   // exit succcesfully (hopefully)
}

/** mover with a Predictor-Corrector scheme */
int Particles3D::mover_PC_sub(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER species " << ns << " ***"<< " with " << nop << " particles ***"  << NiterMover << " ITERATIONS   ****" << endl;
  }
  double start_mover_PC = MPI_Wtime();
  double weights[2][2][2];
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
  if(Gravity)
		Fext /= qom;

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
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];
    const double xptilde = x[rest];
    const double yptilde = y[rest];
    const double zptilde = z[rest];
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
      xp = xptilde + uptilde * dt_sub;
      yp = yptilde + vptilde * dt_sub;
      zp = zptilde + wptilde * dt_sub;
      x[rest] = xp;
      y[rest] = yp;
      z[rest] = zp;
      u[rest] = up;
      v[rest] = vp;
      w[rest] = wp;
    } // END  OF SUBCYCLING LOOP
  }                             // END OF ALL THE PARTICLES

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

/** mover with a Predictor-Corrector scheme for 2D cylindrical symmetric systems */
int Particles3D::mover_PC_sub_cyl(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER CYLINDRICAL species " << ns << " ***"<< " with " << nop << " particles ***"  << NiterMover << " ITERATIONS   ****" << endl;
  }
  double start_mover_PC = MPI_Wtime();
  double weights[2][2][2];
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
  if(Gravity)
		Fext /= qom;

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
    double up = u[rest];
    double vp = v[rest];
    double wp = w[rest];
    const double xptilde = x[rest];
    const double yptilde = y[rest];
    const double zptilde = z[rest];
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
      xp = xptilde + uptilde * dt_sub;
      yp = yptilde + vptilde * dt_sub;
      zp = zptilde + wptilde * dt_sub;

      // Apply cylindrical correction
      // This assumes x= r; y= z; z= theta
      double gamma = atan(zp/xp);
      xp=xp-Lx/2.0;
      xp=xp*cos(gamma)+zp*sin(gamma);
      xp=xp+Lx/2.0;
      zp=0.0;
      double uprot=up*cos(gamma)+wp*sin(gamma);
      wp=-up*sin(gamma)+wp*cos(gamma);
      up=uprot;

      x[rest] = xp;
      y[rest] = yp;
      z[rest] = zp;
      u[rest] = up;
      v[rest] = vp;
      w[rest] = wp;
    } // END  OF SUBCYCLING LOOP
  }                             // END OF ALL THE PARTICLES

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


/** 1D Two-Stream instability and uniform spatial distribution */
void Particles3D::twostream1D(Grid * grid, VirtualTopology3D * vct, int mode) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2 + 10*ns);

  double mxp, myp, mzp;
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
              mxp = u0 / sqrt(1.0 - u0*u0) + uth / sqrt(1.0 - uth*uth) * prob * cos(theta);
              // v
              myp = v0 / sqrt(1.0 - v0*v0) + vth / sqrt(1.0 - vth*vth) * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              mzp = w0 / sqrt(1.0 - w0*w0) + wth / sqrt(1.0 - wth*wth) * prob * cos(theta);

              mxp += uth * sin(x[counter] * 2.0 * M_PI / Lx * mode);
              double g = sqrt(1.0 + mxp*mxp + myp*myp + mzp*mzp);
              u[counter] = mxp / g;
              v[counter] = myp / g;
              w[counter] = mzp / g;
              if ( counter % 2 ==0)
            	  u[counter] = - u[counter];
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

              counter++;
            }


}


/** relativistic with a Predictor Corrector based directly
 * on Boris-like rotation  */
int Particles3D::mover_relativistic(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
	  if (vct->getCartesian_rank() == 0) {
	    cout << "*** MOVER species " << ns << " ***"<< " with " << nop << " particles ***"  << NiterMover << " ITERATIONS   ****" << endl;
	  }
	  double start_mover_PC = MPI_Wtime();
	  double weights[2][2][2];
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
	  if(Gravity)
			Fext /= qom;

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
	    double uxp = u[rest];
	    double uyp = v[rest];
	    double uzp = w[rest];
	    double gamma0, gamma, gamma_new;
	    // Assuming up, vp, wp to be VELOCITY and

	    gamma0 = 1.0/sqrt(1.0 - uxp*uxp - uyp*uyp - uzp*uzp);
	    uxp *= gamma0;
	    uyp *= gamma0;
	    uzp *= gamma0;
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

	    get_weights(grid, xp, yp, zp, ix, iy, iz, weights);
	    get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);

	   // const double B_mag      = sqrt(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
	    //double       dt_sub     = M_PI*c/(4*abs(qom)*B_mag);
	    //const int    sub_cycles = int(dt/dt_sub) + 1;

	    //dt_sub = dt/double(sub_cycles);
	    double dt_sub = dt;
	    const int sub_cycles = 1;

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

	        //relativistic mover

	        // solve the position equation
	        const double wx = uxp + qomdt2 * Exl;
	        const double wy = uyp + qomdt2 * Eyl;
	        const double wz = uzp + qomdt2 * Ezl;

	        gamma = sqrt(1.0 + wx *wx + wy *wy + wz *wz);

			  Bxl *=qomdt2 / gamma;
			  Byl *=qomdt2 / gamma;
			  Bzl *=qomdt2 / gamma;

	        const double h2 = Bxl * Bxl + Byl * Byl + Bzl * Bzl;

	        // solve the velocity equation (Relativistic Boris method)
	       /* uxnew = -(pow(Byl,2)*wx) - pow(Bzl,2)*wx + Bxl*Byl*wy +
	             2*Bzl*wy - 2*Byl*wz + Bxl*Bzl*wz,
		    uynew = Bxl*Byl*wx - 2*Bzl*wx - pow(Bxl,2)*wy -
	             pow(Bzl,2)*wy + 2*Bxl*wz + Byl*Bzl*wz;
            uznew = 2*Byl*wx + Bxl*Bzl*wx - 2*Bxl*wy + Byl*Bzl*wy -
	             pow(Bxl,2)*wz - pow(Byl,2)*wz;
            */
            uxnew = -(pow(Byl,2)*wx) - pow(Bzl,2)*wx + Bxl*Byl*wy +
            	             Bzl*wy - Byl*wz + Bxl*Bzl*wz,
            uynew = Bxl*Byl*wx - Bzl*wx - pow(Bxl,2)*wy -
            	             pow(Bzl,2)*wy + Bxl*wz + Byl*Bzl*wz;
            uznew = Byl*wx + Bxl*Bzl*wx - Bxl*wy + Byl*Bzl*wy -
            	             pow(Bxl,2)*wz - pow(Byl,2)*wz;

            uxnew = wx + uxnew *2.0/(1+h2) + qomdt2 * Exl;
            uynew = wy + uynew *2.0/(1+h2) + qomdt2 * Eyl;
            uznew = wz + uznew *2.0/(1+h2) + qomdt2 * Ezl;

            gamma_new = sqrt(1.0 + uxnew *uxnew + uynew *uynew + uznew *uznew);

            // update position (mid of the time step)
	        xp = xp0 + (uxnew + uxp) /(gamma_new + gamma0) * dto2;
	        yp = yp0 + (uynew + uyp) /(gamma_new + gamma0) * dto2;
	        zp = zp0 + (uznew + uzp) /(gamma_new + gamma0) * dto2;
	      }                           // end of iteration
	      // update the final position and velocity


	      u[rest] = uxnew/ gamma_new;
	      v[rest] = uynew/ gamma_new;
	      w[rest] = uznew/ gamma_new;


	      x[rest] = xp0 +  (uxnew + uxp) /(gamma_new + gamma0) * dt;
	      y[rest] = yp0 +  (uynew + uyp) /(gamma_new + gamma0) * dt;
	      z[rest] = zp0 +  (uznew + uzp) /(gamma_new + gamma0) * dt;

	    } // END  OF SUBCYCLING LOOP
	  }                             // END OF ALL THE PARTICLES

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

/** relativistic mover with a scheme from the old but glorious Celeste3D */
int Particles3D::mover_relativistic_celeste(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
	  if (vct->getCartesian_rank() == 0) {
	    cout << "*** MOVER species " << ns << " ***"<< " with " << nop << " particles ***"  << NiterMover << " ITERATIONS   ****" << endl;
	  }
	  double start_mover_PC = MPI_Wtime();
	  double weights[2][2][2];
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
	  if(Gravity)
			Fext /= qom;
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
	    double pxp = u[rest];
	    double pyp = v[rest];
	    double pzp = w[rest];
	    // Assuming up, vp, wp to be VELOCITY and

	    double gamma0 = 1.0/sqrt(1.0 - (pxp*pxp - pyp*pyp - pzp*pzp)/c/c );
	    pxp *= gamma0;
	    pyp *= gamma0;
	    pzp *= gamma0;
	    const double xptilde = x[rest];
	    const double yptilde = y[rest];
	    const double zptilde = z[rest];
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

	    // No subcyling in relativisitc case

	    get_weights(grid, xp, yp, zp, ix, iy, iz, weights);
	    get_Bl(weights, ix, iy, iz, Bxl, Byl, Bzl, Bx, By, Bz, Bx_ext, By_ext, Bz_ext, Fext);

	   // const double B_mag      = sqrt(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
	    //double       dt_sub     = M_PI*c/(4*abs(qom)*B_mag);
	    //const int    sub_cycles = int(dt/dt_sub) + 1;

	    //dt_sub = dt/double(sub_cycles);
	    double dt_sub = dt;
	    const int sub_cycles = 1;

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
	        double denom = 1.0 / (1.0 + omdtsq);

	        //relativistic correction

			  Bxl /=gamma0;
			  Byl /=gamma0;
			  Bzl /=gamma0;
			  denom /=gamma0;

	        // solve the position equation
	        const double ut = pxp + qomdt2 * Exl;
	        const double vt = pyp + qomdt2 * Eyl;
	        const double wt = pzp + qomdt2 * Ezl;
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



	      // This needs to be revisited because it assumes c=1
	      	double u02=pxp*pxp+pyp*pyp+pzp*pzp;
	      	double v2=uptilde*uptilde+vptilde*vptilde+wptilde*wptilde;
	      	double vdu=pxp*uptilde+pyp*vptilde+pzp*wptilde;

	      	double cfa=1.0-v2;
	      	double cfb=-2.0*(-vdu+gamma0*v2);
	      	double cfc=-1.0-gamma0*gamma0*v2+2.0*gamma0*vdu-u02;

	      	double delta=cfb*cfb-4.0*cfa*cfc;

	      	if(delta<0.0) {
	      		// this in pronciple should never happen. Yeah, right!
	      		cout << "Relativity violated: gamma0=" << gamma0 << ",  vavg_sq=" << v2 << endl;
	      		//u[rest] = (gamma0*2.)*uptilde - u[rest]*gamma0;
	      		//v[rest] = (gamma0*2.)*vptilde - v[rest]*gamma0;
	      		//w[rest] = (gamma0*2.)*wptilde - w[rest]*gamma0;
	      	}
	      	else{
	      		double gamma=(-cfb+sqrt(delta))/2.0/cfa;
	      		double unew = ( (gamma+gamma0)*uptilde - pxp )/gamma;
	      		v[rest] = ( (gamma+gamma0)*vptilde - pyp )/gamma;
	      		w[rest] = ( (gamma+gamma0)*wptilde - pzp )/gamma;
	      		if(abs(unew)>1.0){
	      				cout << "Relativity violated: u>c=" << u[rest] << endl;
	      		}
	      		else
	      			{
	      			u[rest] = unew;
	      			}
	      	}

	      x[rest] = xptilde +  dt * uptilde;
	      y[rest] = yptilde +  dt * vptilde;
	      z[rest] = zptilde +  dt * wptilde;

	    } // END  OF SUBCYCLING LOOP
	  }                             // END OF ALL THE PARTICLES

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
    cout << "*** Reflector species " << ns << " ***" << endl;
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




