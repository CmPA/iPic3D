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
        cout << "*** MOVER species " << ns << " ***" << NiterMover <<" ITERATIONS   ****" << endl;
    }
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
    // move each particle with new fields: MOVE P_SAME_TIME PARTICLES AT THE SAME TIME TO ALLOW AUTOVECTORIZATION
    int i;
    for (i=0; i <  (nop-(P_SAME_TIME-1)); i+=P_SAME_TIME){ 
        // copy x, y, z
        for (int p = 0; p < P_SAME_TIME; p++){
            xp[p] = x[i+p];   yp[p] = y[i+p];   zp[p] = z[i+p];
            up[p] = u[i+p];   vp[p] = v[i+p];   wp[p] = w[i+p];
        }
        for (int p = 0; p < P_SAME_TIME; p++){ // VECTORIZED
            xptilde[p] = xp[p]; yptilde[p] = yp[p]; zptilde[p] = zp[p];
        }
        // calculate the average velocity iteratively
        for(int innter=0; innter < NiterMover; innter++){	  
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
                    for(int kk=0; kk < 2; kk++){
                        for(int p = 0; p < P_SAME_TIME; p++){
                            Exlp[p] = EMf->getEx(ix[p] - ii,iy[p] -jj,iz[p] - kk);  
                            Eylp[p] = EMf->getEy(ix[p] - ii,iy[p] -jj,iz[p] - kk);  
                            Ezlp[p] = EMf->getEz(ix[p] - ii,iy[p] -jj,iz[p] - kk);  
                            Bxlp[p] = EMf->getBx(ix[p] - ii,iy[p] -jj,iz[p] - kk);  
                            Bylp[p] = EMf->getBy(ix[p] - ii,iy[p] -jj,iz[p] - kk);  
                            Bzlp[p] = EMf->getBz(ix[p] - ii,iy[p] -jj,iz[p] - kk); 
                        }
                        for(int p = 0; p < P_SAME_TIME; p++){ // VECTORIZED
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
            for(int p = 0; p < P_SAME_TIME; p++){ // PARTIALLY VECTORIZED
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
        for(int p=0; p < P_SAME_TIME; p++){ // VECTORIZED
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




