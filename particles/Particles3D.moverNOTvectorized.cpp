/*******************************************************************************************
Particles3D.cpp  -  Class for particles of the same species, in a 3D space and 3component velocity
-------------------
developers: Stefano Markidis, Giovanni Lapenta
	********************************************************************************************/

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


#include "hdf5.h"
#include <complex>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define MIN_VAL   1E-16
/**
* 
 * Class for particles of the same species, in a 2D space and 3component velocity
 * @date Fri Jun 4 2009
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */

/** constructor */
Particles3D::Particles3D(){
	// see allocate(int species, CollectiveIO* col, VirtualTopology3D* vct, Grid* grid)
	
}
/** deallocate particles */
Particles3D::~Particles3D(){
    delete[] x;
    delete[] y;
	delete[] z;
    delete[] u;
    delete[] v;
    delete[] w;
    delete[] q;
}

/** particles are uniformly distributed with zero velocity   */
void Particles3D::uniform_background(Grid* grid,Field* EMf){
	int counter=0;
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int k=1; k< grid->getNZC()-1;k++)
				for (int ii=0; ii < npcelx; ii++)
					for (int jj=0; jj < npcely; jj++)
						for (int kk=0; kk < npcelz; kk++){
							x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);
							y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
							z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
							u[counter] = 0.0;
							v[counter] = 0.0;
							w[counter] = 0.0;
							q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,k,ns)/npcel)*(1.0/grid->getInvVOL());
							if (TrackParticleID)	      
								ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
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
void Particles3D::constantVelocity(double vel, int dim,Grid* grid,Field* EMf){
	switch(dim)
	{
		case 0:
			for (int i=0; i < nop; i++)
				u[i] = vel, v[i] = 0.0, w[i] = 0.0;
			break;
		case 1:
			for (register int i=0; i < nop; i++)
				u[i] = 0.0,v[i] = vel,w[i] = 0.0;
			break;
		case 2:
			for (register int i=0; i < nop; i++)
				u[i] = 0.0, v[i] = 0.0, w[i] = vel;
			break;
			
	}
	
}

/** alternative routine maxellian random velocity and uniform spatial distribution */
void Particles3D::alt_maxwellian(Grid* grid,Field* EMf,VirtualTopology3D* vct){
	
	
}

/** Maxellian random velocity and uniform spatial distribution */
void Particles3D::maxwellian(Grid* grid,Field* EMf,VirtualTopology3D* vct){
	
	
	double harvest;
	double prob, theta, sign;
	int counter=0;
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int k=1; k< grid->getNZC()-1;k++)
				for (int ii=0; ii < npcelx; ii++)
					for (int jj=0; jj < npcely; jj++)
						for (int kk=0; kk < npcelz; kk++){
							x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);   // x[i] = xstart + (xend-xstart)/2.0 + harvest1*((xend-xstart)/4.0)*cos(harvest2*2.0*M_PI);
							y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
							z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
							// q = charge
							q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,k,ns)/npcel)*(1.0/grid->getInvVOL());
							// u
							harvest =   rand()/(double)RAND_MAX;
							prob  = sqrt(-2.0*log(1.0-.999999*harvest));
							harvest =   rand()/(double)RAND_MAX;
							theta = 2.0*M_PI*harvest;
							u[counter] = u0 + uth*prob*cos(theta);
							// v
							v[counter] = v0 + vth*prob*sin(theta);
							// w
							harvest =   rand()/(double)RAND_MAX;
							prob  = sqrt(-2.0*log(1.0-.999999*harvest));
							harvest =   rand()/(double)RAND_MAX;
							theta = 2.0*M_PI*harvest;
							w[counter] = w0 + wth*prob*cos(theta);
							if (TrackParticleID)	      
								ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
							
							
							counter++ ;
						}
							
							
}

/** Force Free initialization (JxB=0) for particles */
void Particles3D::force_free(Grid* grid,Field* EMf,VirtualTopology3D* vct){
	
	
	double harvest, prob, theta;
	int counter=0;
	double shaperx, shapery, shaperz;
	double flvx=1.0,flvy=1.0,flvz=1.0;
	
	
	/* initialize random generator */
	srand (vct->getCartesian_rank()+1+ns);
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int k=1; k< grid->getNZC()-1;k++)
				for (int ii=0; ii < npcelx; ii++)
					for (int jj=0; jj < npcely; jj++)
						for (int kk=0; kk < npcelz; kk++){
							flvx = 1.0;
							flvy = 1.0;
							flvz = 1.0;
							x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,k);   
							y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,k);
							z[counter] = (kk + .5)*(dz/npcelz) + grid->getZN(i,j,k);
							// q = charge
							q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,k,ns)/npcel)*(1.0/invVOL);
							shaperx = tanh((y[counter] - Ly/2)/delta)/cosh((y[counter] - Ly/2)/delta)/delta;
							shaperz = 1.0/(cosh((y[counter] - Ly/2)/delta)*cosh((y[counter] - Ly/2)/delta))/delta;
							shapery = shapery;
							// new drift velocity to satisfy JxB=0
							flvx =u0*flvx*shaperx;
							flvz =w0*flvz*shaperz;
							flvy =v0*flvy*shapery;
							u[counter] = c; v[counter] = c; w[counter] = c;
							while ((fabs(u[counter])>=c) | (fabs(v[counter])>=c) | (fabs(w[counter])>=c)){
								harvest =   rand()/(double)RAND_MAX;
								prob  = sqrt(-2.0*log(1.0-.999999*harvest));
								harvest =   rand()/(double)RAND_MAX;
								theta = 2.0*M_PI*harvest;
								u[counter] = flvx + uth*prob*cos(theta);
								// v
								v[counter] = flvy + vth*prob*sin(theta);
								// w
								harvest =   rand()/(double)RAND_MAX;
								prob  = sqrt(-2.0*log(1.0-.999999*harvest));
								harvest =   rand()/(double)RAND_MAX;
								theta = 2.0*M_PI*harvest;
								w[counter] = flvz + wth*prob*cos(theta);}
							if (TrackParticleID)	      
								ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
							
							counter++ ;
						}
							
}

/**Add a periodic perturbation in J exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
inline void Particles3D::AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0, Grid* grid){
	
	// rescaling of amplitudes according to deltaBoB //
	double alpha;
	alpha=deltaBoB*B0/sqrt(Bx_mod*Bx_mod+By_mod*By_mod+Bz_mod*Bz_mod);
	jx_mod *=alpha;
	jy_mod *=alpha;
	jz_mod *=alpha;
	for (register int i=0; i<nop; i++){
		u[i] += jx_mod/q[i]/npcel/invVOL* cos(kx*x[i] + ky*y[i] + jx_phase);
		v[i] += jy_mod/q[i]/npcel/invVOL* cos(kx*x[i] + ky*y[i] + jy_phase);
		w[i] += jz_mod/q[i]/npcel/invVOL* cos(kx*x[i] + ky*y[i] + jz_phase);}
}


/** explicit mover */
void Particles3D::mover_explicit(Grid* grid,VirtualTopology3D* vct, Field* EMf){
 	// to be implemented
	
}
/** mover with a Predictor-Corrector scheme */
int Particles3D::mover_PC(Grid* grid,VirtualTopology3D* vct, Field* EMf){
	if (vct->getCartesian_rank()==0){
		cout << "*** MOVER species " << ns << " ***" << endl;
		cout << "*** 3 ITERATIONS   ****"  <<  endl;
	}
	int avail;
	double dto2 = .5*dt, qomdt2 = qom*dto2/c, omdtsq, denom, ut, vt, wt, udotb;
	double Exl=0.0, Eyl=0.0, Ezl=0.0, Bxl=0.0, Byl=0.0, Bzl=0.0;
	double inv_dx = 1.0/dx, inv_dy = 1.0/dy, inv_dz = 1.0/dz;
	int ix,iy,iz;
	double xptilde, yptilde, zptilde, uptilde, vptilde, wptilde;
	double weight[2][2][2]; //double*** weight = newArr3(double,2,2,2);
	// move each particle with new fields
        for (int i=0; i <  nop; i++){
		xptilde = x[i];
		yptilde = y[i];
		zptilde = z[i];
		// calculate the average velocity iteratively
		for(int innter=0; innter < 1; innter++){	  
			// interpolation G-->P
			ix = 2 +  int(floor((x[i]-xstart)*inv_dx));
			iy = 2 +  int(floor((y[i]-ystart)*inv_dy));
			iz = 2 +  int(floor((z[i]-zstart)*inv_dz));
			calculateWeights(weight,x[i],y[i],z[i],ix,iy,iz,grid);
			Exl=0.0, Eyl = 0.0, Ezl = 0.0, Bxl = 0.0, Byl = 0.0, Bzl = 0.0;
            for (int ii=0; ii < 2; ii++)
				for (int jj=0; jj < 2; jj++)
					for(int kk=0; kk < 2; kk++){
						Exl += weight[ii][jj][kk]*(EMf->getEx(ix - ii,iy -jj,iz- kk ));
						Eyl += weight[ii][jj][kk]*(EMf->getEy(ix - ii,iy -jj,iz- kk ));
						Ezl += weight[ii][jj][kk]*(EMf->getEz(ix - ii,iy -jj,iz -kk ));
						Bxl += weight[ii][jj][kk]*(EMf->getBx(ix - ii,iy -jj,iz -kk ));
						Byl += weight[ii][jj][kk]*(EMf->getBy(ix - ii,iy -jj,iz -kk ));
						Bzl += weight[ii][jj][kk]*(EMf->getBz(ix - ii,iy -jj,iz -kk ));
					}
			// end interpolation
			omdtsq = qomdt2*qomdt2*(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
			denom = 1.0/(1.0 + omdtsq);
			// solve the position equation
			ut= u[i] + qomdt2*Exl;
			vt= v[i] + qomdt2*Eyl;
			wt= w[i] + qomdt2*Ezl;
			udotb = ut*Bxl + vt*Byl + wt*Bzl;
			// solve the velocity equation 
			uptilde = (ut+qomdt2*(vt*Bzl -wt*Byl + qomdt2*udotb*Bxl))*denom; 
			vptilde = (vt+qomdt2*(wt*Bxl -ut*Bzl + qomdt2*udotb*Byl))*denom;
			wptilde = (wt+qomdt2*(ut*Byl -vt*Bxl + qomdt2*udotb*Bzl))*denom;
			// update position
			x[i] = xptilde + uptilde*dto2;
			y[i] = yptilde + vptilde*dto2;
			z[i] = zptilde + wptilde*dto2;
		} // end of iteration
		// update the final position and velocity
		u[i]= 2.0*uptilde - u[i];
		v[i]= 2.0*vptilde - v[i];
		w[i]= 2.0*wptilde - w[i];
		x[i] = xptilde + uptilde*dt;
		y[i] = yptilde + vptilde*dt;
		z[i] = zptilde + wptilde*dt;
		// apply the BC
		if (x[i] < 0 && vct->getXleft_neighbor() == MPI_PROC_NULL)
			BCpart(&x[i],&u[i],&v[i],&w[i],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft);
		else if (x[i] > Lx && vct->getXright_neighbor() == MPI_PROC_NULL)
			BCpart(&x[i],&u[i],&v[i],&w[i],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft); 
		if (y[i] < 0 && vct->getYleft_neighbor() == MPI_PROC_NULL)  
			BCpart(&y[i],&v[i],&u[i],&w[i],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft);
		else if (y[i] > Ly && vct->getYright_neighbor() == MPI_PROC_NULL) 
			BCpart(&y[i],&v[i],&u[i],&w[i],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft); 
		if (z[i] < 0 && vct->getZleft_neighbor() == MPI_PROC_NULL)  
			BCpart(&z[i],&v[i],&u[i],&w[i],Lz,vth,uth,wth,bcPfaceZright,bcPfaceZleft);
		else if (z[i] > Lz && vct->getZright_neighbor() == MPI_PROC_NULL)
		    BCpart(&z[i],&v[i],&u[i],&w[i],Lz,vth,uth,wth,bcPfaceZright,bcPfaceZleft);
	} // end of particles
	// COMMUNICATION
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
	// deallocate
	//delArr3(weight,2,2);
	return(0); // exit succcesfully	
}



/** relativistic mover with a Predictor-Corrector scheme */
int Particles3D::mover_relativistic(Grid* grid,VirtualTopology3D* vct, Field* EMf){
	return(0);
}


/** interpolation Particle->Grid only for pressure tensor */    
void Particles3D::interpP2G_onlyP(Field* EMf, Grid *grid, VirtualTopology3D* vct){
	double weight[2][2][2];
	double temp[2][2][2];
	int ix, iy, iz, temp1, temp2, temp3;
	for (register int i=0; i < nop; i++){
		ix = 2 +  int(floor((x[i]-grid->getXstart())/grid->getDX()))  ;
		iy = 2 +  int(floor((y[i]-grid->getYstart())/grid->getDY()))  ;
		iz = 2 +  int(floor((z[i]-grid->getZstart())/grid->getDZ()))  ;
		calculateWeights(weight,x[i],y[i],z[i],ix,iy,iz,grid);
		scale(weight,q[i],2,2,2);
		//Pxx
		eqValue(0.0,temp,2,2,2);
		addscale(u[i]*u[i],temp,weight,2,2,2);
		EMf->addPxx(temp,ix,iy,iz,ns);
		// Pxy
		eqValue(0.0,temp,2,2,2);
		addscale(u[i]*v[i],temp,weight,2,2,2);
		EMf->addPxy(temp,ix,iy,iz,ns);
		// Pxz
		eqValue(0.0,temp,2,2,2);
		addscale(u[i]*w[i],temp,weight,2,2,2);
		EMf->addPxz(temp,ix,iy,iz,ns);
		// Pyy
		eqValue(0.0,temp,2,2,2);
		addscale(v[i]*v[i],temp,weight,2,2,2);
		EMf->addPyy(temp,ix,iy,iz,ns);
		// Pyz
		eqValue(0.0,temp,2,2,2);
		addscale(v[i]*w[i],temp,weight,2,2,2);
		EMf->addPyz(temp,ix,iy,iz,ns);
		// Pzz
		eqValue(0.0,temp,2,2,2);
		addscale(w[i]*w[i],temp,weight,2,2,2);
		EMf->addPzz(temp,ix,iy,iz,ns);
	}
	//delArr3(weight,2,2);
	//delArr3(temp,2,2);
}
/** interpolation Particle->Grid only charge density, current */    
void Particles3D::interpP2G_notP(Field* EMf, Grid *grid, VirtualTopology3D* vct){
	double weight[2][2][2];
	double temp[2][2][2];
	int ix,iy, iz, temp2,temp1,temp3;
	for (register int i=0; i < nop; i++){
		ix = 2 +  int(floor((x[i]-grid->getXstart())/grid->getDX()))  ;
		iy = 2 +  int(floor((y[i]-grid->getYstart())/grid->getDY()))  ;
		iz = 2 +  int(floor((z[i]-grid->getZstart())/grid->getDZ()))  ;
		temp1 = (int) min(ix, nxn-2);
		temp2 = (int) min(iy, nyn-2);
		temp3 = (int) min(iz, nzn-2);
		ix = (int) max(temp1,2);
		iy = (int) max(temp2,2);
		iz = (int) max(temp3,2);
		calculateWeights(weight,x[i],y[i],z[i],ix,iy,iz,grid);
		scale(weight,q[i],2,2,2);
		// rho
		EMf->addRho(weight,ix,iy,iz,ns);
		// Jx
		eqValue(0.0,temp,2,2,2);
		addscale(u[i],temp,weight,2,2,2);
		EMf->addJx(temp,ix,iy,iz,ns);
		// Jy
		eqValue(0.0,temp,2,2,2);
		addscale(v[i],temp,weight,2,2,2);
		EMf->addJy(temp,ix,iy,iz,ns);
		// Jz
		eqValue(0.0,temp,2,2,2);
		addscale(w[i],temp,weight,2,2,2);
		EMf->addJz(temp,ix,iy,iz,ns);
		
	}
	// communicate contribution from ghost cells     
	EMf->communicateGhostP2G(ns,0,0,0,0,vct);
	//delArr3(weight,2,2);
	//delArr3(temp,2,2);
	
}
/** apply a linear perturbation to particle distribution */
void Particles3D::linear_perturbation(double deltaBoB, double kx, double ky, double angle, double omega_r,     double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Grid* grid,Field* EMf,VirtualTopology3D* vct){
	
	double value1=0.0,value2=0.0,max_value=0.0,min_value=0.0,phi,n;
	int counter=0, total_generated=0;
	bool rejected;
	double harvest, prob, theta;
	// rescaling of amplitudes according to deltaBoB //
	double alpha;
	double integral=0.0;
	
	alpha=deltaBoB*sqrt(EMf->getBx(1,1,0)*EMf->getBx(1,1,0)+EMf->getBy(1,1,0)*EMf->getBy(1,1,0)+EMf->getBz(1,1,0)*EMf->getBz(1,1,0))/sqrt(Bx_mod*Bx_mod+By_mod*By_mod+Bz_mod*Bz_mod);
	
	Ex_mod *= alpha;
	Ey_mod *= alpha;
	Ez_mod *= alpha;
	Bx_mod *= alpha;
	By_mod *= alpha;
	Bz_mod *= alpha;
	
	
	
	// find the maximum value of f=1+delta_f/f0
	for (register double vpar=-2*uth; vpar<=2*uth; vpar += 0.0005)
		for (register double vperp=1e-10; vperp<=2*vth; vperp += 0.0005)
			for (register double X=xstart; X<=xend; X += 2*grid->getDX())
				for (register double Y=ystart; Y<=yend; Y += 2*grid->getDY()){
					value1=1+delta_f(vpar,vperp,0.0, X, Y, kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle,  EMf)/f0(vpar,vperp);
					
					if (value1>max_value)
						max_value=value1;
					
					
				}
					
					
					
					max_value *=3.2;phi=1.48409;n=2.948687; // security factor...
					if (ns==1){
						max_value *=3.0;phi=-1.65858;n=2.917946;} // security factor...
					cout<<"max-value="<<max_value<<" min-value="<<min_value<<endl;
					
					/* initialize random generator */
					srand (vct->getCartesian_rank()+2);
					
					for (int i=1; i< grid->getNXC()-1;i++)
						for (int j=1; j< grid->getNYC()-1;j++)
							for (int ii=0; ii < npcelx+(int)(2*n*(cos(2*M_PI*0.4125*grid->getXN(i,j,0)+phi))); ii++)
								for (int jj=0; jj < npcely; jj++){
									x[counter] = (ii + .5)*(dx/(npcelx+(int)(2*n*(cos(2*M_PI*0.4125*grid->getXN(i,j,0)+phi))))) + grid->getXN(i,j,0);   
									y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,0);
									q[counter] =  (qom/fabs(qom))*((0.19635)/npcel)*(1.0/invVOL);
									
									// apply rejection method in velocity space
									rejected=true;
									while (rejected){
										total_generated ++;
										harvest =   rand()/(double)RAND_MAX;
										prob  = sqrt(-2.0*log(1.0-.999999*harvest));
										harvest =   rand()/(double)RAND_MAX;
										theta = 2.0*M_PI*harvest;
										// u
										u[counter] = u0 + uth*prob*cos(theta);
										// v
										v[counter] = v0 + vth*prob*sin(theta);
										// w
										harvest =   rand()/(double)RAND_MAX;
										prob  = sqrt(-2.0*log(1.0-.999999*harvest));
										harvest =   rand()/(double)RAND_MAX;
										theta = 2.0*M_PI*harvest;
										w[counter] = w0 + wth*prob*cos(theta);
										
										// test: if rand < (1+delta_f/f0)/max_value --> accepted
										if ( rand()/(double)RAND_MAX <= (1+ delta_f(u[counter],v[counter],w[counter],x[counter], y[counter], kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase,Ez_mod, Ez_phase, angle,  EMf)/f0(u[counter],sqrt(v[counter]*v[counter]+w[counter]*w[counter])))/max_value)
											rejected=false;
										
									}
									if (TrackParticleID)	      
										ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
									counter++ ;
								}
									nop=counter+1;
					//		     if (vct->getCartesian_rank()==0)
					cout<<"Rejection method: "<<(counter+1)/ double(total_generated) * 100<< " % of particles are accepted for species "<<ns<<" counter="<<counter<<endl; 
}

/** Linear delta f for bi-maxwellian plasma */
double Particles3D::delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase,double Ez_mod, double Ez_phase, double theta, Field* EMf){
	const complex<double> I(0.0,1.0);
	const double vperp = sqrt(v*v+w*w);
	const double vpar = u;
	const double kpar= kx;
	double kperp;
	if (ky==0.0) // because this formula is not valid for exactly parallel
		kperp = 1e-9;
	else kperp = ky;
	const double om_c = qom/c*sqrt(EMf->getBx(1,1,0)*EMf->getBx(1,1,0)+EMf->getBy(1,1,0)*EMf->getBy(1,1,0))/2/M_PI;
	const double phi = atan2(w,v);
	const double lambda = kperp*vperp/om_c;
	const complex<double> omega (omega_re,omega_i);
	
	const int lmax=5; // sum from -lmax to lmax
	
	double bessel_Jn_array[ lmax + 2 ];
	double bessel_Jn_prime_array[ lmax +1 ];
	complex <double> a1[2*lmax+1], a2[2*lmax+1],a3[2*lmax+1] ;
	complex <double> factor, deltaf; 
	
	// rotation of x,y
	double temp;
	temp=x;
	x=x*cos(theta)-y*sin(theta);
	y=temp*sin(theta)+y*cos(theta);
	
	
	/** for compilation issues comment this part: PUT in the math stuff */
	//calc_bessel_Jn_seq(lambda, lmax, bessel_Jn_array, bessel_Jn_prime_array);
	factor = (kpar*vperp/omega*df0_dvpar(vpar,vperp) + (1.0-(kpar*vpar/omega))*df0_dvperp(vpar,vperp) );
	for (register int l=-lmax; l<0; l++){ // negative index
		a1[l+lmax] = factor/lambda*pow(-1.0,-l)*bessel_Jn_array[-l];
		a1[l+lmax] *= (double) l;
		a2[l+lmax] = factor* I* 0.5*pow(-1.0,-l)*(bessel_Jn_array[-l-1]-bessel_Jn_array[-l+1]);
		a3[l+lmax] = kperp/omega*(vpar*df0_dvperp(vpar,vperp) - vperp*df0_dvpar(vpar,vperp))/lambda*pow(-1.0,-l)*bessel_Jn_array[-l];
		a3[l+lmax] *= (double) l;
		a3[l+lmax] += df0_dvpar(vpar,vperp)*pow(-1.0,-l)*bessel_Jn_array[-l];
	}
	
	for (register int l=0; l<lmax+1; l++){ //positive index
		a1[l+lmax] = factor/lambda*bessel_Jn_array[l];
		a1[l+lmax] *= (double) l;
		a2[l+lmax] = factor* I* bessel_Jn_prime_array[l];
		a3[l+lmax] = kperp/omega*(vpar*df0_dvperp(vpar,vperp) - vperp*df0_dvpar(vpar,vperp))/lambda*bessel_Jn_array[l];
		a3[l+lmax] *= (double) l;
		a3[l+lmax] += df0_dvpar(vpar,vperp)*bessel_Jn_array[l];
	}
	
	deltaf=(0.0,0.0);
	for (register int l=-lmax; l<lmax+1; l++){
		deltaf += (a3[l+lmax]*Ex_mod*exp(I*Ex_phase) + a1[l+lmax]*Ey_mod*exp(I*Ey_phase) + a2[l+lmax]*Ez_mod*exp(I*Ez_phase))/(kpar*vpar+l*om_c-omega)*exp(-I*phi*(double)l);}
	deltaf *= I*qom*exp(I*lambda*sin(phi))*exp(I*(2*M_PI*kx*x+2*M_PI*ky*y));
	
	return(real(deltaf));
}

double Particles3D::df0_dvpar(double vpar,double vperp){
	double result;
	result = -2*(vpar-u0)/uth/uth*exp(-(vperp*vperp/vth/vth + (vpar-u0) * (vpar-u0) /uth/uth));
	result *= 3.92e6/pow(M_PI,3/2)/vth/vth/uth;
	return(result);
}

double Particles3D::df0_dvperp(double vpar,double vperp){
	double result;
	result = -2*(vperp)/vth/vth*exp(-(vperp*vperp/vth/vth + (vpar-u0) * (vpar-u0) /uth/uth));
	result *= 3.92e6/pow(M_PI,3/2)/vth/vth/uth;
	return(result);
}

double Particles3D::f0(double vpar,double vperp){
	double result;
	result = exp(-(vperp*vperp/vth/vth + (vpar-u0) * (vpar-u0) /uth/uth));
	result *= 3.92e6/pow(M_PI,3/2)/vth/vth/uth;
	return(result);
}

void Particles3D::RotatePlaneXY(double theta){
	double temp,temp2;
	for (register int s=0; s < nop; s++){
		temp=u[s];temp2=v[s];
		u[s]=temp*cos(theta)+v[s]*sin(theta);
		v[s]=-temp*sin(theta)+temp2*cos(theta);
	}
}




