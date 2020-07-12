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
	      /*      #ifdef EB
	      if (nyn==4){ // meaning it's 1D; now particles are populated around y=0
		y[counter]-=Ly/2;
	      }
	      if (nzn==4){ // meaning it's 1D
		z[counter]-=Lz/2;
	      }
	      #endif*/
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

/**** maxwellian expansion deformed *****/
void Particles3D::maxwellian_EBDeformationWithoutInteraction(Grid * grid, Field * EMf, VirtualTopology3D * vct) {
  
  
  /* the reason for this: 
     assuming absence of interaction, VDF gets deformed during expansion, apart from the expected 1/R
     in the perpendicular direction;
     deformation is expected larger for sub-sonic species (in our case, electron), see Landi et al, 2012, Fig 2
     so here I want to start from the deformed VDF expected at a certain, hardcoded R, and see if a field is
     generated such that the deformation is dissipated 

     how to do it:
     particle per particle, first I generate a "regular maxwellian"
     then I impose (copied from the matlab script)
     ---vvy= v0./R_R0;
     ---vvz= w0./R_R0;

     conservation of energy, keeping into account the UEB
     it's important here to select the right solution when solving the equation
     the correct solution is one that does not artificially scatter the particle parallely

     ---vperp_0= sqrt(v0.^2 + w0.^2);
     ---vperp= sqrt(vvy.^2 +vvz.^2 );

     ---%%% just substitute v before with sqrt(v^2 + w^2)

     ---a= 1;
     ---b= 2*U0;
     ---c= vperp.^2 + U0^2- ( (U0 + u0).^2 + vperp_0.^2);

     ---ff= sqrt(b.^2 -4.*a.*c);

     ---vvx1=( -b + ff)./ (2*a);
     ---vvx2=( -b - ff)./ (2*a);

     ---q= ones(size(vvy));
     ---%%% select the right solution, avoiding paralleling scattering

     ---%%% NB: if i do not add the +U0, I am selecting the wrong solution for
     ---%%% part of the population
     ---cond1= find(u0+ U0>=0);
     ---cond2= find(u0+ U0<0);


     ---vx1=vvx1(cond1);
     ---vx2=vvx2(cond2);

     ---vy1=vvy(cond1);
     ---vy2=vvy(cond2);

     ---vz1=vvz(cond1);
     ---vz2=vvz(cond2);

     ---q1=q(cond1);
     ---q2=q(cond2); 


     ---vx= [vx1; vx2];
     ---vy= [vy1 ;vy2];
     ---vz= [vz1 ;vz2];
        
     ---q= [q1; q2];


     Warning: this thing is not correct if there are transverse B components, because it is
     assumed that the yz directions (transverse directions) are also perpendicular
     also, I am assuming there are no drift velocities, I did not think of what happens 
     with drifts

     /*
  
  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);

     double R= 196;
     double R_R0= R/REB;
     cout << "In maxwellian_EBDeformationWithoutInteraction, I am using an hardcoded value of R, R=" << 196 <<endl;
     cout << "R/R_0= " << R_R0 << ", U0_EB" << UEB << endl;


  double harvest;
  double prob, theta, sign;
  long long counter = 0;
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++)
        for (int ii = 0; ii < npcelx; ii++)
          for (int jj = 0; jj < npcely; jj++)
            for (int kk = 0; kk < npcelz; kk++) {

	      double u_tmp;
	      double v_tmp;
	      double w_tmp;

              x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);  // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
              y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
              z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
	      
              q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOL());
              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              //u[counter] = u0 + uth * prob * cos(theta);
	      u_tmp= u0 + uth * prob * cos(theta);
              // v
              //v[counter] = v0 + vth * prob * sin(theta);
	      v_tmp= v0 + vth * prob * sin(theta);  
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              //w[counter] = w0 + wth * prob * cos(theta);
	      w_tmp = w0 + wth * prob * cos(theta); 


	      double vperp0= sqrt(v_tmp*v_tmp + w_tmp*w_tmp);
	      v[counter]= v_tmp/R_R0;
	      w[counter]= w_tmp/R_R0;

	      double vperp= sqrt(v[counter]* v[counter]+ w[counter]*w[counter]);

	      double a=1;
	      double b=2*UEB;
	      double c= vperp*vperp + UEB*UEB - ((UEB + u_tmp)* (UEB + u_tmp) + vperp0*vperp0  );
	      double ff= sqrt(b*b - 4 *a*c);

	      double sol1= (-b + ff)/ (2*a);
	      double sol2= (-b - ff)/ (2*a);

	      if (u_tmp+ UEB > 0 ){
		u[counter]= sol1;
		  }else{
		u[counter]= sol2;
	      }

              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];

	      // now, fake expansion

              counter++;
            }


}

/**** end maxwellian expansion-deformed ****/

/** maxwellian with perpendicular current oscillations for electrons supporting whistler **/
void Particles3D::maxwellian_WhistlerCurrent(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* to be used with initBy for fields;
     deltaB to be read from inputfile */

  if (qom<0){
    cout << "Particle init: imposing perp current oscillations, deltaB= " << deltaB <<endl;
  }

  // specific of whistler init
  double k0= 2* M_PI/ Lx;
  double V0= - k0* deltaB/ fabs(rhoINIT)/ ElectronSpNumber;
  cout << "species " << ns << " rhoINIT " << rhoINIT << " V0 " << V0 << "=(- k0 delta B/ |n0|)/" <<ElectronSpNumber <<", with " << ElectronSpNumber << "the # of current-carrying electron species" <<endl;
  double COS, SIN, Vy, Vz;

  // I have to distribute the current among the current-carrying species, i.e. electrons
  // (I may have more than one species if i do tracking)
  

  // end specific of whistler init 

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);

  double harvest;
  double prob, theta, sign;
  long long counter = 0;
  
  for (int i = 1; i < grid->getNXC() - 1; i++){
    COS= cos(k0* grid->getXN(i, 0, 0));
    SIN= sin(k0* grid->getXN(i, 0, 0));
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
	      /* the shape of the perp electron current is given through this */
	      if (qom< 0){
		Vy= v0+ V0*SIN;
		Vz= w0+ V0*COS;
	      }else{
		Vy= v0;
		Vz= w0;
	      }

              // u
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              u[counter] = u0 + uth * prob * cos(theta);
              // v
              v[counter] = Vy + vth * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              w[counter] = Vz + wth * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];


              counter++;
	    }
  }
  

}

/* end maxwellian with perpendicular current oscillations for electrons supporting whistler */

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
int Particles3D::mover_PC(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  double start_mover_PC = MPI_Wtime();
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


/*** mover, EB **/
/** mover with a Predictor-Corrector scheme, with expanding box */
int Particles3D::mover_PC_EB(Grid * grid, VirtualTopology3D * vct, Field * EMf, EBox * ebox){
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER !!!EXPANDING BOX!!! species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
  }
  /* this is the intermediate position of the grid, 
     still add intermediate position of the particle the particle */
  double R_nth= ebox->getR_nth();
  if (vct->getCartesian_rank() == 0 and ns==0) {
    cout << "R_nth in mover: " << R_nth << endl;
  }
  // R0
  double R0_EB= ebox->getREB_0();
  // U0
  double U0_EB= ebox->getUEB_0();

  double start_mover_PC = MPI_Wtime();
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

      // end interpolation
      const double omdtsq = qomdt2 * qomdt2 * (Bxl * Bxl + Byl * Byl + Bzl * Bzl);
      // solve the position equation
      const double ut = up + qomdt2 * Exl;
      const double vt = vp + qomdt2 * Eyl;
      const double wt = wp + qomdt2 * Ezl;
      const double udotb = ut * Bxl + vt * Byl + wt * Bzl;
      // solve the velocity equation 
      /* No EB
      const double denom = 1.0 / (1.0 + omdtsq);   
      uptilde = (ut + qomdt2 * (vt * Bzl - wt * Byl + qomdt2 * udotb * Bxl)) * denom;
      vptilde = (vt + qomdt2 * (wt * Bxl - ut * Bzl + qomdt2 * udotb * Byl)) * denom;
      wptilde = (wt + qomdt2 * (ut * Byl - vt * Bxl + qomdt2 * udotb * Bzl)) * denom;*/
      /* EB: extra term at the denom of perpendicular terms, U0/R= U0/ (R0 +x) */
      const double denom_par = 1.0 / (1.0 + omdtsq);
      // Rnth= intermediate position grid (R_nth) + intermediate position particle (xp)
      //const double denom_perp = 1.0 / (1.0 + omdtsq + dto2* U0_EB/ (R_nth+xp));
      // for the moment, neglect local displacement within the grid
      const double denom_perp = 1.0 / (1.0 + omdtsq + dto2* U0_EB/ (R_nth)); 
      uptilde = (ut + qomdt2 * (vt * Bzl - wt * Byl + qomdt2 * udotb * Bxl)) * denom_par;  
      vptilde = (vt + qomdt2 * (wt * Bxl - ut * Bzl + qomdt2 * udotb * Byl)) * denom_perp; 
      wptilde = (wt + qomdt2 * (ut * Byl - vt * Bxl + qomdt2 * udotb * Bzl)) * denom_perp;
      // update position
      xp = xptilde + uptilde * dto2;                                                  
      // Rn+th = intermediate position of grid (R_nth) + intermediate position of particle  
      /*yp = yptilde + vptilde * dto2* (R0_EB/ (R_nth + xptilde));                         
	zp = zptilde + wptilde * dto2* (R0_EB/ (R_nth + xptilde));*/
      // for the moment, neglect local displacement within the grid
      yp = yptilde + vptilde * dto2* (R0_EB/ (R_nth ));                                                              
      zp = zptilde + wptilde * dto2* (R0_EB/ (R_nth ));
      // update position, expanding box (R_0/R) in perp direction
    }                           // end of iteration
    // update the final position and velocity
    up = 2.0 * uptilde - u[rest];
    vp = 2.0 * vptilde - v[rest];
    wp = 2.0 * wptilde - w[rest];
    /*xp = xptilde + uptilde * dt;
    yp = yptilde + vptilde * dt;
    zp = zptilde + wptilde * dt;*/
    double x_int= xp; // intermediate x position, needed for Px
    xp = xptilde + uptilde * dt;
    // Rn+th = intermediate position of grid (R_nth) + intermediate position of particle
    /*yp = yptilde + vptilde * dt* (R0_EB/ (R_nth + x_int)); 
      zp = zptilde + wptilde * dt* (R0_EB/ (R_nth + x_int));*/
    // for the moment, neglect local displacement within the grod
    yp = yptilde + vptilde * dt* (R0_EB/ (R_nth ));                                                                    
    zp = zptilde + wptilde * dt* (R0_EB/ (R_nth ));
    /* EB: add the R0_R term in the perp direction */
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

/** end mover, EB **/



/** mover with a Predictor-Corrector scheme */
int Particles3D::mover_PC_sub(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  if (vct->getCartesian_rank() == 0) {
    cout << "*** MOVER species " << ns << " ***" << NiterMover << " ITERATIONS   ****" << endl;
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

/** relativistic mover with a Predictor-Corrector scheme */
int Particles3D::mover_relativistic(Grid * grid, VirtualTopology3D * vct, Field * EMf) {
  return (0);
}

int Particles3D::particle_repopulator(Grid* grid,VirtualTopology3D* vct, Field* EMf, int is){

  /* -- NOTE: Hardcoded option -- */
  enum {LINEAR,INITIAL,FFIELD};
  int rtype = INITIAL;
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
  if (vct->getXleft_neighbor() == MPI_PROC_NULL && bcPfaceXleft == 2){ // use Field topology in this case
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
  if (vct->getYleft_neighbor() == MPI_PROC_NULL  && bcPfaceYleft == 2)
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
  if (vct->getZleft_neighbor() == MPI_PROC_NULL  && bcPfaceZleft == 2)
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
  if (vct->getXright_neighbor() == MPI_PROC_NULL  && bcPfaceXright == 2){
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
  if (vct->getYright_neighbor() == MPI_PROC_NULL  && bcPfaceYright == 2)
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
  if (vct->getZright_neighbor() == MPI_PROC_NULL  && bcPfaceZright == 2)
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

/** Maxellian random velocity and uniform spatial distribution 
 double periodic, current sign switched in the two half planes*/
void Particles3D::maxwellian_HarrisDoublePeriodic(Grid * grid, Field * EMf, VirtualTopology3D * vct) {

  /* initialize random generator with different seed on different processor */
  srand(vct->getCartesian_rank() + 2);

  double harvest;
  double prob, theta, sign;
  double flvx, flvy, flvz;
  double shaperx, shapery, shaperz;
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

	      shaperz= 1.0;
	      shapery= 1.0 ;
              shaperz= -tanh((grid->getYN(i, j, k) - Ly/2)/delta) ;

	      flvx =u0*shaperx;
	      flvy =v0*shapery;
	      flvz =w0*shaperz;

              u[counter] = flvx + uth * prob * cos(theta);
              // v
              v[counter] = flvy + vth * prob * sin(theta);
              // w
              harvest = rand() / (double) RAND_MAX;
              prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
              harvest = rand() / (double) RAND_MAX;
              theta = 2.0 * M_PI * harvest;
              w[counter] = flvz + wth * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];


              counter++;
            }


}
// Papini's turbulence initialization, in the xy plane
void Particles3D::alfredoturbulence(Grid * grid, Field * EMf, VirtualTopology3D * vct, Collective* col) {

  const double coarsedx = grid->getDX();
  const double coarsedy = grid->getDY();
  double globalxC,globalyC;

  double prob, theta, sign;
  long long counter = 0;

  double uu0;
  double vv0;

  const double pi = 3.141592653589793;
  const double Lx = col->getLx();
  const double Ly = col->getLy();
  const double amp0 = .06 ;
  const double slope = 0.;
  const double rkx0 = 2*pi/Lx;
  const double rky0 = 2*pi/Ly;
  const int    mkx0 = -3;
  const int    mkx1 = 3;
  const int    mky0 = -3;
  const int    mky1 = 3;
  int    idummy=0;
  long int myseed=vct->getCartesian_rank()+47 + ns*vct->getNprocs();
  double B0z;
  double rphv[(mkx1-mkx0+1)*(mky1-mky0+1)];
  //double* rphv = new double[mkx1-mkx0+1,mky1-mky0+1];

  double fact;
  int ikx,iky;
  double totalrms=0.0;

  //B0z = EMf->getB0z();
  B0z= col->getB0z(); 
  //totalrms=EMf->getBrms();

  //cout << "rank, seed, ns " <<vct->getCartesian_rank() << ", "<< myseed << ", " <<ns << endl;


  srand48(1);
  for ( ikx = mkx0; ikx <= mkx1; ikx++)
    for ( iky = mky0; iky <= mky1; iky++)
      {
	rphv[idummy] = drand48()*2*pi;
	idummy++;
	//rphv[ikx-mkx0,iky-mky0] = drand48()*2*pi;
      }


  srand48(myseed);
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++){
	//need variables in the center of the cell to calculate fluctuations
	globalxC= grid->getXN(i, j, k) + coarsedx/2;
	globalyC= grid->getYN(i, j, k) + coarsedy/2;
	uu0=0.;
	vv0=0.;
	// initialize alfven like fluctuations in uthp and vthp       
	idummy=0;
	for ( ikx = mkx0; ikx <= mkx1; ikx++)
	  for ( iky = mky0; iky <= mky1; iky++)
	    {
	      //skip the k=0 point
	      if ((iky != 0) || (ikx != 0)) {
		fact = B0z*amp0*pow(sqrt(pow(rkx0*ikx,2) + pow(rky0*iky,2) ),slope-1.);

		uu0 += + fact*rky0*iky*cos( rkx0*ikx*globalxC + rky0*iky*globalyC +rphv[idummy]);
		vv0 += - fact*rkx0*ikx*cos( rkx0*ikx*globalxC + rky0*iky*globalyC +rphv[idummy]);
	      }
	      idummy++;
	    }

	for (int ii = 0; ii < npcelx; ii++)
	  for (int jj = 0; jj < npcely; jj++)
	    for (int kk = 0; kk < npcelz; kk++) {



	      //x[counter] = (ii + .5) * (dx / npcelx) + grid->getXN(i, j, k);
	      //y[counter] = (jj + .5) * (dy / npcely) + grid->getYN(i, j, k);
	      //z[counter] = (kk + .5) * (dz / npcelz) + grid->getZN(i, j, k);
	      x[counter] = (1.0 - .999999 * drand48())* dx  + grid->getXN(i, j, k);
	      y[counter] = (1.0 - .999999 * drand48())* dy  + grid->getYN(i, j, k);
	      z[counter] = (1.0 - .999999 * drand48())* dz  + grid->getZN(i, j, k);
	      // q = charge
	      q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOL());
	      // u
	      prob = sqrt(-2.0 * log(1.0 - .999999 * drand48()));
	      //theta = 2.0 * M_PI * drand48();
	      theta = 2.0 * pi * drand48();
	      u[counter] = u0 + uu0 + uth * prob * cos(theta);
	      // v
	      v[counter] = v0 + vv0 + vth * prob * sin(theta);
	      // w
	      prob = sqrt(-2.0 * log(1.0 - .999999 * drand48()));
	      //theta = 2.0 * M_PI * drand48();
	      theta = 2.0 * pi * drand48();
	      w[counter] = w0 + wth * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];
	      counter++;
	    }
      }


  nop = counter;
}

// from here
// Papini's turbulence initialization, in the yz plane
void Particles3D::alfredoturbulence_yz(Grid * grid, Field * EMf, VirtualTopology3D * vct, Collective* col) {

  const double coarsedx = grid->getDX();
  const double coarsedy = grid->getDY();
  const double coarsedz = grid->getDZ();
  double globalxC,globalyC,globalzC ;

  double prob, theta, sign;
  long long counter = 0;

  double uu0;
  double vv0;

  const double pi = 3.141592653589793;
  const double Lx = col->getLx();
  const double Ly = col->getLy();
  const double amp0 = .06 ;
  const double slope = 0.;
  /*const double rkx0 = 2*pi/Lx;
    const double rky0 = 2*pi/Ly;*/
  // inverted here
  const double rkx0 = 2*pi/Ly;
  const double rky0 = 2*pi/Lz;
  // inverted here
  const int    mkx0 = -3;
  const int    mkx1 = 3;
  const int    mky0 = -3;
  const int    mky1 = 3;
  int    idummy=0;
  long int myseed=vct->getCartesian_rank()+47 + ns*vct->getNprocs();
  //double B0z;
  double B0x;
  double rphv[(mkx1-mkx0+1)*(mky1-mky0+1)];
  //double* rphv = new double[mkx1-mkx0+1,mky1-mky0+1];

  double fact;
  int ikx,iky;
  double totalrms=0.0;

  //B0z= col->getB0z();
  B0x= col->getB0x(); 

  srand48(1);
  for ( ikx = mkx0; ikx <= mkx1; ikx++)
    for ( iky = mky0; iky <= mky1; iky++)
      {
	rphv[idummy] = drand48()*2*pi;
	idummy++;
	//rphv[ikx-mkx0,iky-mky0] = drand48()*2*pi;
      }


  srand48(myseed);
  for (int i = 1; i < grid->getNXC() - 1; i++)
    for (int j = 1; j < grid->getNYC() - 1; j++)
      for (int k = 1; k < grid->getNZC() - 1; k++){
	//need variables in the center of the cell to calculate fluctuations
	/*globalxC= grid->getXN(i, j, k) + coarsedx/2;
	  globalyC= grid->getYN(i, j, k) + coarsedy/2;*/
	// inverted here
	globalxC= grid->getYN(i, j, k) + coarsedy/2;
	globalyC= grid->getZN(i, j, k) + coarsedz/2;
	uu0=0.;
	vv0=0.;
	// initialize alfven like fluctuations in uthp and vthp       
	idummy=0;
	for ( ikx = mkx0; ikx <= mkx1; ikx++)
	  for ( iky = mky0; iky <= mky1; iky++)
	    {
	      //skip the k=0 point
	      if ((iky != 0) || (ikx != 0)) {
		//fact = B0z*amp0*pow(sqrt(pow(rkx0*ikx,2) + pow(rky0*iky,2) ),slope-1.);
		fact = B0x*amp0*pow(sqrt(pow(rkx0*ikx,2) + pow(rky0*iky,2) ),slope-1.);
		// uu0, vv0 are actually in y and z
		uu0 += + fact*rky0*iky*cos( rkx0*ikx*globalxC + rky0*iky*globalyC +rphv[idummy]);
		vv0 += - fact*rkx0*ikx*cos( rkx0*ikx*globalxC + rky0*iky*globalyC +rphv[idummy]);
	      }
	      idummy++;
	    }

	for (int ii = 0; ii < npcelx; ii++)
	  for (int jj = 0; jj < npcely; jj++)
	    for (int kk = 0; kk < npcelz; kk++) {

	      x[counter] = (1.0 - .999999 * drand48())* dx  + grid->getXN(i, j, k);
	      y[counter] = (1.0 - .999999 * drand48())* dy  + grid->getYN(i, j, k);
	      z[counter] = (1.0 - .999999 * drand48())* dz  + grid->getZN(i, j, k);
	      // q = charge
	      q[counter] = (qom / fabs(qom)) * (fabs(EMf->getRHOcs(i, j, k, ns)) / npcel) * (1.0 / grid->getInvVOL());
	      // u
	      prob = sqrt(-2.0 * log(1.0 - .999999 * drand48()));
	      //theta = 2.0 * M_PI * drand48();
	      theta = 2.0 * pi * drand48();
	      //u[counter] = u0 + uu0 + uth * prob * cos(theta);
	      u[counter] = u0 + uth * prob * cos(theta);
	      // v
	      //v[counter] = v0 + vv0 + vth * prob * sin(theta);
	      v[counter] = v0 + uu0 + vth * prob * sin(theta); 
	      // w
	      prob = sqrt(-2.0 * log(1.0 - .999999 * drand48()));
	      //theta = 2.0 * M_PI * drand48();
	      theta = 2.0 * pi * drand48();
	      //w[counter] = w0 + wth * prob * cos(theta);
	      w[counter] = w0 + vv0 + wth * prob * cos(theta);
              if (TrackParticleID)
                ParticleID[counter] = counter * (unsigned long) pow(10.0, BirthRank[1]) + BirthRank[0];
	      counter++;
	    }
      }


  nop = counter;
}

// to here

