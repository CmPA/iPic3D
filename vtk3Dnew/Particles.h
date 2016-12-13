/*
 *  Particles.h
 *  
 *
 *  Created by Giovanni Lapenta on 7/29/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "Alloc.h"
#include "math.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

void load(double* xp, double* vp);
void mover(double* xp, double* vp);


int np;
double Dt;
double qoms = 1;
double *xp;
double *yp;
double *zp;

double *up;
double *vp;
double *wp;

void load(double* xp, double* vp){
	double u0 = 0.0;
	double v0 = 0.0;
	double w0 = 0.0;
	double uth = 1.0;
	double vth = 1.0;
	double wth = 1.0;	
	double harvest;
	double prob, theta, sign;
	for (int ip=0; ip<np; ip++){
		harvest = rand() / (double) RAND_MAX;
		xp[ip] = harvest;
		harvest = rand() / (double) RAND_MAX;
		yp[ip] = harvest;
		harvest = rand() / (double) RAND_MAX;
		zp[ip] = harvest;
		
		harvest = rand() / (double) RAND_MAX;
		prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
		harvest = rand() / (double) RAND_MAX;
		theta = 2.0 * M_PI * harvest;
		up[ip] = u0 + uth * prob * cos(theta);
		// v
		vp[ip] = v0 + vth * prob * sin(theta);
		// w
		harvest = rand() / (double) RAND_MAX;
		prob = sqrt(-2.0 * log(1.0 - .999999 * harvest));
		harvest = rand() / (double) RAND_MAX;
		theta = 2.0 * M_PI * harvest;
		wp[ip] = w0 + wth * prob * cos(theta);
  		
     }

}

void mover(double* xp, double* vp){
	
	for (int ip=0; ip<np; ip++){
		xp[ip] += up[ip] * Dt;
		yp[ip] += vp[ip] * Dt;
		zp[ip] += wp[ip] * Dt;

		up[ip] += qoms * Dt * 0.0;
		vp[ip] += qoms * Dt * 0.0;
		wp[ip] += qoms * Dt * 0.0;
		
		cout <<"xp= " << xp[ip] << endl; 
	}
	
}