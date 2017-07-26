/*
 *  Manipulator.h
 *  
 *
 *  Created by Giovanni Lapenta on 7/29/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "hdf5.h"
#include "Alloc.h"
#include "math.h"
#include "Basic.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

extern "C" void dgeev_( char* jobvl, char* jobvr, int* n, double* a,
                int* lda, double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info );

using namespace std;
using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;

// Manipulations

void sort(int n, double* vet);

void mat2vet(int n, double** mat, double* vet);

double ppar(int n, double** p, double* b);

double norm(int n,  double** p);

double norm(int n,  double* b);

void vet2mat(int n, double** mat, double* vet, double* b, double* dot);

void eq(double ***vect1, double ***vect2, int nx, int ny, int nz);

void smooth(int Nvolte, double*** A, double*** B, int nx, int ny, int nz);

void smooth(int Nvolte, double*** B, int nx, int ny, int nz);

void average(double** EXB, double*** EX);

void averageSTD(double** EXB, double** EXSTD, double*** EX);

void average(double** EXB, double** EYB, double** EZB, double*** EX, double*** EY, double*** EZ);

void averagedot(double** EdotJX, double** EdotJY, double** EdotJZ,
		   double** dEdotJX, double** dEdotJY, double** dEdotJZ,
		   double** EXB, double** EYB, double** EZB,
		   double** VXB, double** VYB, double** VZB,
		   double*** EX, double*** EY, double*** EZ, 
		   double*** VX, double*** VY, double*** VZ);

void averagedot(double** EdotJX, double** EdotJY, double** EdotJZ,
                double** dEdotJX, double** dEdotJY, double** dEdotJZ,
                double** sEdotJX, double** sEdotJY, double** sEdotJZ,
                double** EXB, double** EYB, double** EZB,
                double** VXB, double** VYB, double** VZB,
                double*** EX, double*** EY, double*** EZ,
                double*** VX, double*** VY, double*** VZ);

void averagecross(double** CX, double** CY, double** CZ,
				double** dCX, double** dCY, double** dCZ,
				double** EXB, double** EYB, double** EZB,
				double** BXB, double** BYB, double** BZB,
				double*** EX, double*** EY, double*** EZ, 
				double*** BX, double*** BY, double*** BZ);

void averagecross(double** CX, double** CY, double** CZ,
                  double** dCX, double** dCY, double** dCZ,
                  double** sCX, double** sCY, double** sCZ,
                  double** EXB, double** EYB, double** EZB,
                  double** BXB, double** BYB, double** BZB,
                  double*** EX, double*** EY, double*** EZ,
                  double*** BX, double*** BY, double*** BZ);

 void average_energy(double** AVGENTH, double** AVGENIN, double** AVGJDOTE,
		double** STDENTH, double** STDENIN, double** STDJDOTE,
		double*** VDIVP, double*** DIVJV, double*** JDOTE);

 void averageY(double** EXB, double*** EX);

 void averageYSTD(double** EXB, double** EXSTD, double*** EX);

 void averageY(double** EXB, double** EYB, double** EZB, double*** EX, double*** EY, double*** EZ);

 void averageYdot(double** EdotJX, double** EdotJY, double** EdotJZ,
 		   double** dEdotJX, double** dEdotJY, double** dEdotJZ,
 		   double** EXB, double** EYB, double** EZB,
 		   double** VXB, double** VYB, double** VZB,
 		   double*** EX, double*** EY, double*** EZ,
 		   double*** VX, double*** VY, double*** VZ);

 void averageYdot(double** EdotJX, double** EdotJY, double** EdotJZ,
                 double** dEdotJX, double** dEdotJY, double** dEdotJZ,
                 double** sEdotJX, double** sEdotJY, double** sEdotJZ,
                 double** EXB, double** EYB, double** EZB,
                 double** VXB, double** VYB, double** VZB,
                 double*** EX, double*** EY, double*** EZ,
                 double*** VX, double*** VY, double*** VZ);

 void averageYcross(double** CX, double** CY, double** CZ,
 				double** dCX, double** dCY, double** dCZ,
 				double** EXB, double** EYB, double** EZB,
 				double** BXB, double** BYB, double** BZB,
 				double*** EX, double*** EY, double*** EZ,
 				double*** BX, double*** BY, double*** BZ);

 void averageYcross(double** CX, double** CY, double** CZ,
                   double** dCX, double** dCY, double** dCZ,
                   double** sCX, double** sCY, double** sCZ,
                   double** EXB, double** EYB, double** EZB,
                   double** BXB, double** BYB, double** BZB,
                   double*** EX, double*** EY, double*** EZ,
                   double*** BX, double*** BY, double*** BZ);

  void averageY_energy(double** AVGENTH, double** AVGENIN, double** AVGJDOTE,
 		double** STDENTH, double** STDENIN, double** STDJDOTE,
 		double*** VDIVP, double*** DIVJV, double*** JDOTE);


void extract_pressure(double qom, double*** BX, double*** BY, double*** BZ,
					  double*** VX, double*** VY, double*** VZ,
					  double*** N, double*** pXX, double*** pXY,
					  double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ,
					  double*** pPAR, double*** pPER1, double*** pPER2, double*** EPS);
/* Computes Bulk Energy Flux */
void bulk_energy_flux(double*** Ubulk,
					  double*** QXbulk, double*** QYbulk, double*** QZbulk,
					  double*** VX, double*** VY, double*** VZ);
/* Computes Intenal Energy Flux */
void internal_energy_flux(double*** Uth,
		              double*** QX, double*** QY, double*** QZ,
					  double*** VX, double*** VY, double*** VZ,
					  double*** pXX, double*** pXY,
					  double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ);
// Computes P:gradU
void pgradU(double*** PgradU, double*** VX, double*** VY, double*** VZ,
		double*** WX, double*** WY, double*** WZ,
		double*** TXX, double***  TXY, double*** TXZ,
		double*** TYY, double*** TYZ, double*** TZZ);

void compute_energy(double qom, double*** Ubulk, double*** Uth,
		  double*** JX, double*** JY, double*** JZ,
		  double*** N, double*** pXX, double*** pYY, double*** pZZ);

void gradient(double*** EX, double*** EY, double*** EZ, double*** phi);

void gradient_par(double*** BX, double*** BY, double*** BZ,
				  double*** phi, double*** out);

void gradient_per1(double*** BX, double*** BY, double*** BZ,
				  double*** phi, double*** out);

void gradient_per2(double*** BX, double*** BY, double*** BZ,
				  double*** phi, double*** out);

void divergence(double*** div, double*** EX, double*** EY, double*** EZ);

// Computes div (J/N) where J is a vector JX, JY, JZ
void divergenceN(double*** div, double*** JX, double*** JY, double*** JZ, double*** N);

// Computes div (J*P/N) where J is a vector JX, JY, JZ
void divergenceNP(double*** div, double*** JX, double*** JY, double*** JZ, double*** N, double*** P);

// Computes  U div (P) where P is hte pressure tensor
void udivP(double*** O, double*** OX, double*** OY, double*** OZ,
		double*** JX, double*** JY, double*** JZ,
		double*** pXX, double*** pXY,
				  double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ);

// Computes  Agyrotropy according to Scudder
void agyro(double*** agyro_scudder, double*** agyro_aunai, double*** nongyro_swisdak, double*** align,
		double*** BX, double*** BY, double*** BZ,
		double*** pXX, double*** pXY,
		double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ);

void curl(double*** curlx, double*** curly, double*** curlz, double*** EX, double*** EY, double*** EZ);

void vdotgrad(double*** vdotgrad, double*** PHI, double*** RHO,
                         double*** JX, double*** JY, double*** JZ);

void vdiv(double*** vdiv, double*** PXX, double*** PYY, double*** PZZ,
                         double*** PXY, double*** PXZ, double*** PYZ,
                         double*** RHO, 
                         double*** JX, double*** JY, double*** JZ);

//Computes P_hi *d(phi)/dx_i in cartesian coordinates, output is OX, OY, OZ
void p_dot_gradient(double*** OX, double*** OY, double*** OZ,
		  double*** pXX, double*** pXY,
		  double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ,
		  double*** phi);

//Computes P_hi *dU_k/dx_i in cartesian coordinates
void pdv(double*** pdvX, double*** pdvY, double*** pdvZ,
						 double*** PXX, double*** PYY, double*** PZZ,
                         double*** PXY, double*** PXZ, double*** PYZ,
						 double*** BX, double*** BY, double*** BZ,
                         double*** RHO,
                         double*** JX, double*** JY, double*** JZ,
						 double*** V,
						 double*** VX, double*** VY, double*** VZ);

//Computes P_hi *dU_k/dx_i in magnetic coordinates (par, perp1, perp2)
void pdv_par_per(double*** pdvpar, double*** pdvper1, double*** pdvper2,
						 double*** PXX, double*** PYY, double*** PZZ,
                         double*** PXY, double*** PXZ, double*** PYZ,
						 double*** pPAR, double*** pPER1, double*** pPER2,
						 double*** BX, double*** BY, double*** BZ,
                         double*** RHO,
                         double*** JX, double*** JY, double*** JZ,
						 double*** V,
						 double*** VX, double*** VY, double*** VZ);

// Compute divergence of enthalpy component flux
void divPv(double*** OX, double*** OY, double*** OZ,
		  double*** PXX, double*** PYY, double*** PZZ,
		  double*** PXY, double*** PXZ, double*** PYZ,
          double*** JX, double*** JY, double*** JZ, double*** N,
		  double*** WX, double*** WY, double*** WZ);

void divj(double qom, double*** divj,  
          double*** EX, double*** EY, double*** EZ, double VSQ);
          
void dot(double*** dot, 
          double*** AX, double*** AY, double*** AZ, 
          double*** BX, double*** BY, double*** BZ, double*** RHO);
          
void dot(double*** dot, 
          double*** AX, double*** AY, double*** AZ, 
          double*** BX, double*** BY, double*** BZ);

void adddot(double*** dot,
          double*** AX, double*** AY, double*** AZ,
          double*** BX, double*** BY, double*** BZ);

// Computation of ExBfor 3D arrays
void cross(double*** EX, double*** EY, double*** EZ,
		   double*** BX, double*** BY, double*** BZ,
		   double*** VXBX, double*** VXBY, double*** VXBZ);


// Computation of ExB for 1D array A and 2D array B
void cross(int m, double* A, double** B, double** AXB);

// Computation of ExB for 1D array A and 2D array B transposed
void cross_trans(int m, double* A, double** B, double** AXB);


// Computation of ExB/B^2 for 3D arrays
void crossoverb(double*** EX, double*** EY, double*** EZ,
		   double*** BX, double*** BY, double*** BZ,
		   double*** VXBX, double*** VXBY, double*** VXBZ);

// Computation of ExB for 2D arrays
void cross(double** EX, double** EY, double** EZ,
		   double** BX, double** BY, double** BZ,
		   double** VXBX, double** VXBY, double** VXBZ);

//Computes S=ExB and its divergence
void div_cross(double*** EX, double*** EY, double*** EZ,
		   double*** BX, double*** BY, double*** BZ,
		   double*** SX, double*** SY, double*** SZ,
		   double*** divS);

void vnonfrozen(double*** N, double*** BX, double*** BY, double*** BZ, 
				double*** EX, double*** EY, double*** EZ, 
				double*** VX, double*** VY, double*** VZ);

void ohm(double*** N, double*** BX, double*** BY, double*** BZ, 
		 double*** EX, double*** EY, double*** EZ, 
		 double*** VX, double*** VY, double*** VZ,
		 double*** OHMX, double*** OHMY, double*** OHMZ);

void zenitani(double*** NE, double*** NI, 
 		 double*** BX, double*** BY, double*** BZ, 
		 double*** EX, double*** EY, double*** EZ, 
		 double*** JEX, double*** JEY, double*** JEZ,
		 double*** JIX, double*** JIY, double*** JIZ,
		 double*** DE, double*** DI);

// B x curl (Eparallel)
void biskamp(double*** BX, double*** BY, double*** BZ,
		 double*** EX, double*** EY, double*** EZ,
		 double*** WX, double*** WY, double*** WZ,
		 double*** BiskampX, double*** BiskampY, double*** BiskampZ);

// Lagrangian Energy Balance TErm: S.gradient(log(B))
void LEBT(double*** BX, double*** BY, double*** BZ,
		 double*** SX, double*** SY, double*** SZ,
		 double*** WX, double*** WY, double*** WZ,
		 double*** psi);

void mult(double a, double*** E);

// C++ implementation of the circshift MATLAB function
void circshift(double ***vect, double ***tmp, int xshift, int yshift, int zshift);


//Compute A = A *B /C
void divmult(double*** A, double*** B, double*** C);

double vecpot(double*** AZ, double*** BX, double*** BY, double*** BZ);

int nxn, nyn, nzn;
int nodes_z;
int nxc, nyc, nzc;
int XLEN, YLEN, ZLEN;
int nproc;
int ns;
double* qom;

double Lx, Ly, Lz;
double Dx, Dy, Dz;

void average(double** EXB, double** EYB, double** EZB, double*** EX, double*** EY, double*** EZ){
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][jj] = 0.0;
			EYB[ii][jj] = 0.0;
			EZB[ii][jj] = 0.0;
		}			
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				EXB[ii][jj] += EX[ii][jj][kk];				
				EYB[ii][jj] += EY[ii][jj][kk];
				EZB[ii][jj] += EZ[ii][jj][kk];			
			}
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][jj] = EXB[ii][jj] / nzn/ZLEN;
			EYB[ii][jj] = EYB[ii][jj] / nzn/ZLEN;
			EZB[ii][jj] = EZB[ii][jj] / nzn/ZLEN;
		}	
}

void average(double** EXB, double*** EX){
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][jj] = 0.0;

		}			
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				EXB[ii][jj] += EX[ii][jj][kk];				
			
			}
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][jj] = EXB[ii][jj] / nzn/ZLEN;

		}	
}


void averageSTD(double** EXB, double** EXSTD, double*** EX){
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][jj] = 0.0;
            
		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				EXB[ii][jj] += EX[ii][jj][kk];
                
			}
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][jj] = EXB[ii][jj] / nzn/ZLEN;
            
		}	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXSTD[ii][jj] = 0.0;
            
		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				EXSTD[ii][jj] += pow(EX[ii][jj][kk]-EXB[ii][jj],2);
            }
    for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXSTD[ii][jj] = sqrt(EXSTD[ii][jj] / nzn/ZLEN);
            
		}
                
}

void averagedot(double** EdotJX, double** EdotJY, double** EdotJZ,
				double** dEdotJX, double** dEdotJY, double** dEdotJZ,
				double** EXB, double** EYB, double** EZB,
				double** VXB, double** VYB, double** VZB,
				double*** EX, double*** EY, double*** EZ, 
				double*** VX, double*** VY, double*** VZ){
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EdotJX[ii][jj] = EXB[ii][jj] * VXB[ii][jj];
			EdotJY[ii][jj] = EYB[ii][jj] * VYB[ii][jj];
			EdotJZ[ii][jj] = EZB[ii][jj] * VZB[ii][jj];
			dEdotJX[ii][jj] = 0.0;
			dEdotJY[ii][jj] = 0.0;
			dEdotJZ[ii][jj] = 0.0;		
		}			
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				dEdotJX[ii][jj] += (EX[ii][jj][kk]-EXB[ii][jj]) * (VX[ii][jj][kk]-VXB[ii][jj]);				
				dEdotJY[ii][jj] += (EY[ii][jj][kk]-EYB[ii][jj]) * (VY[ii][jj][kk]-VYB[ii][jj]);		
				dEdotJZ[ii][jj] += (EZ[ii][jj][kk]-EZB[ii][jj]) * (VZ[ii][jj][kk]-VZB[ii][jj]);							
			}
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dEdotJX[ii][jj] = dEdotJX[ii][jj] / nzn/ZLEN;
			dEdotJY[ii][jj] = dEdotJY[ii][jj] / nzn/ZLEN;
			dEdotJZ[ii][jj] = dEdotJZ[ii][jj] / nzn/ZLEN;
		}	
	
}

void averagedot(double** EdotJX, double** EdotJY, double** EdotJZ,
				double** dEdotJX, double** dEdotJY, double** dEdotJZ,
				double** sEdotJX, double** sEdotJY, double** sEdotJZ,
				double** EXB, double** EYB, double** EZB,
				double** VXB, double** VYB, double** VZB,
				double*** EX, double*** EY, double*** EZ,
				double*** VX, double*** VY, double*** VZ){
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EdotJX[ii][jj] = EXB[ii][jj] * VXB[ii][jj];
			EdotJY[ii][jj] = EYB[ii][jj] * VYB[ii][jj];
			EdotJZ[ii][jj] = EZB[ii][jj] * VZB[ii][jj];
			dEdotJX[ii][jj] = 0.0;
			dEdotJY[ii][jj] = 0.0;
			dEdotJZ[ii][jj] = 0.0;
			sEdotJX[ii][jj] = 0.0;
			sEdotJY[ii][jj] = 0.0;
			sEdotJZ[ii][jj] = 0.0;
		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				dEdotJX[ii][jj] += (EX[ii][jj][kk]-EXB[ii][jj]) * (VX[ii][jj][kk]-VXB[ii][jj]);
				dEdotJY[ii][jj] += (EY[ii][jj][kk]-EYB[ii][jj]) * (VY[ii][jj][kk]-VYB[ii][jj]);
				dEdotJZ[ii][jj] += (EZ[ii][jj][kk]-EZB[ii][jj]) * (VZ[ii][jj][kk]-VZB[ii][jj]);
			}
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dEdotJX[ii][jj] = dEdotJX[ii][jj] / nzn/ZLEN;
			dEdotJY[ii][jj] = dEdotJY[ii][jj] / nzn/ZLEN;
			dEdotJZ[ii][jj] = dEdotJZ[ii][jj] / nzn/ZLEN;
		}
    for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				sEdotJX[ii][jj] += pow(EX[ii][jj][kk] * VX[ii][jj][kk]
                                    - EdotJX[ii][jj] - dEdotJX[ii][jj], 2);
				sEdotJY[ii][jj] += pow(EY[ii][jj][kk] * VY[ii][jj][kk]
                                       - EdotJY[ii][jj] - dEdotJY[ii][jj], 2);
				sEdotJZ[ii][jj] += pow(EZ[ii][jj][kk] * VZ[ii][jj][kk]
                                       - EdotJZ[ii][jj] - dEdotJZ[ii][jj], 2);
            }
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			sEdotJX[ii][jj] = sqrt(sEdotJX[ii][jj] / nzn/ZLEN);
			sEdotJY[ii][jj] = sqrt(sEdotJY[ii][jj] / nzn/ZLEN);
			sEdotJZ[ii][jj] = sqrt(sEdotJZ[ii][jj] / nzn/ZLEN);
		}

}


void averagecross(double** CX, double** CY, double** CZ,
				double** dCX, double** dCY, double** dCZ,
				double** EXB, double** EYB, double** EZB,
				double** BXB, double** BYB, double** BZB,
				double*** EX, double*** EY, double*** EZ, 
				double*** BX, double*** BY, double*** BZ){
	
	cross(EXB,EYB,EZB,BXB,BYB,BZB,CX,CY,CZ);

	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dCX[ii][jj] = 0.0;
			dCY[ii][jj] = 0.0;
			dCZ[ii][jj] = 0.0;		
		}			
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){

				dCX[ii][jj] += (EY[ii][jj][kk] - EYB[ii][jj]) *
                                               (BZ[ii][jj][kk] - BZB[ii][jj]) - 
                                               (EZ[ii][jj][kk] - EZB[ii][jj]) *
					       (BY[ii][jj][kk] - BYB[ii][jj]);
				
				dCY[ii][jj] += (EZ[ii][jj][kk] - EZB[ii][jj]) *
                                               (BX[ii][jj][kk] - BXB[ii][jj]) - 
                                               (EX[ii][jj][kk] - EXB[ii][jj]) *
					       (BZ[ii][jj][kk] - BZB[ii][jj]);

				dCZ[ii][jj] += (EX[ii][jj][kk] - EXB[ii][jj]) *
                                               (BY[ii][jj][kk] - BYB[ii][jj]) - 
                                               (EY[ii][jj][kk] - EYB[ii][jj]) *
					       (BX[ii][jj][kk] - BXB[ii][jj]);
			}
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dCX[ii][jj] = dCX[ii][jj] / nzn/ZLEN;
			dCY[ii][jj] = dCY[ii][jj] / nzn/ZLEN;
			dCZ[ii][jj] = dCZ[ii][jj] / nzn/ZLEN;
		}	
	
}	


void averagecross(double** CX, double** CY, double** CZ,
                  double** dCX, double** dCY, double** dCZ,
                  double** sCX, double** sCY, double** sCZ,
                  double** EXB, double** EYB, double** EZB,
                  double** BXB, double** BYB, double** BZB,
                  double*** EX, double*** EY, double*** EZ,
                  double*** BX, double*** BY, double*** BZ){
	
	cross(EXB,EYB,EZB,BXB,BYB,BZB,CX,CY,CZ);
    
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dCX[ii][jj] = 0.0;
			dCY[ii][jj] = 0.0;
			dCZ[ii][jj] = 0.0;
            sCX[ii][jj] = 0.0;
			sCY[ii][jj] = 0.0;
			sCZ[ii][jj] = 0.0;
		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
                
				dCX[ii][jj] += (EY[ii][jj][kk] - EYB[ii][jj]) *
                (BZ[ii][jj][kk] - BZB[ii][jj]) -
                (EZ[ii][jj][kk] - EZB[ii][jj]) *
                (BY[ii][jj][kk] - BYB[ii][jj]);
				
				dCY[ii][jj] += (EZ[ii][jj][kk] - EZB[ii][jj]) *
                (BX[ii][jj][kk] - BXB[ii][jj]) -
                (EX[ii][jj][kk] - EXB[ii][jj]) *
                (BZ[ii][jj][kk] - BZB[ii][jj]);
                
				dCZ[ii][jj] += (EX[ii][jj][kk] - EXB[ii][jj]) *
                (BY[ii][jj][kk] - BYB[ii][jj]) -
                (EY[ii][jj][kk] - EYB[ii][jj]) *
                (BX[ii][jj][kk] - BXB[ii][jj]);
			}
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dCX[ii][jj] = dCX[ii][jj] / nzn/ZLEN;
			dCY[ii][jj] = dCY[ii][jj] / nzn/ZLEN;
			dCZ[ii][jj] = dCZ[ii][jj] / nzn/ZLEN;
		}
    for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
                
				sCX[ii][jj] += pow(EY[ii][jj][kk] * BZ[ii][jj][kk]  -
                                   EZ[ii][jj][kk] * BY[ii][jj][kk] - CX[ii][jj]
                                   -dCX[ii][jj],2);

                sCY[ii][jj] += pow(EZ[ii][jj][kk] * BX[ii][jj][kk]  -
                                   EX[ii][jj][kk] * BZ[ii][jj][kk] - CY[ii][jj]
                                   -dCY[ii][jj],2);

                sCZ[ii][jj] += pow(EX[ii][jj][kk] * BY[ii][jj][kk]  -
                                   EY[ii][jj][kk] * BX[ii][jj][kk] - CZ[ii][jj]
                                   -dCZ[ii][jj],2);
                
			}
    for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			sCX[ii][jj] = sqrt(sCX[ii][jj] / nzn/ZLEN);
			sCY[ii][jj] = sqrt(sCY[ii][jj] / nzn/ZLEN);
			sCZ[ii][jj] = sqrt(sCZ[ii][jj] / nzn/ZLEN);
		}

	
}


void average_energy(double** AVGENTH, double** AVGENIN, double** AVGJDOTE,
		double** STDENTH, double** STDENIN, double** STDJDOTE,
		double*** VDIVP, double*** DIVJV, double*** JDOTE){

	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			AVGJDOTE[ii][jj] = 0.0;
			STDJDOTE[ii][jj] = 0.0;
			AVGENTH[ii][jj] = 0.0;
			AVGENIN[ii][jj] = 0.0;
			STDENTH[ii][jj] = 0.0;
			STDENIN[ii][jj] = 0.0;
		}

    for (int kk=1; kk < nzn*ZLEN-1;kk++)
		for (int jj=1; jj < nyn*YLEN-1;jj++)
			for (int ii=1; ii < nxn*XLEN-1;ii++){

				AVGENTH[ii][jj] += VDIVP[ii][jj][kk]  ;
				AVGENIN[ii][jj] += (JDOTE[ii][jj][kk] - VDIVP[ii][jj][kk])  ;
				AVGJDOTE[ii][jj] += JDOTE[ii][jj][kk]   ;
			}

    for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			AVGENTH[ii][jj] = AVGENTH[ii][jj] / (nzn*ZLEN-2);
			AVGENIN[ii][jj] = AVGENIN[ii][jj] / (nzn*ZLEN-2);
			AVGJDOTE[ii][jj] = AVGJDOTE[ii][jj] / (nzn*ZLEN-2);
		}


    for (int kk=1; kk < nzn*ZLEN-1;kk++)
		for (int jj=1; jj < nyn*YLEN-1;jj++)
			for (int ii=1; ii < nxn*XLEN-1;ii++){

				STDENTH[ii][jj] += pow(VDIVP[ii][jj][kk]  - AVGENTH[ii][jj],2);
				STDENIN[ii][jj] += pow(JDOTE[ii][jj][kk] - VDIVP[ii][jj][kk]  - AVGENIN[ii][jj],2);
				STDJDOTE[ii][jj] += pow(JDOTE[ii][jj][kk]  - AVGJDOTE[ii][jj],2);

			}
    for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			STDENTH[ii][jj] = sqrt(STDENTH[ii][jj] / (nzn*ZLEN-2));
			STDENIN[ii][jj] = sqrt(STDENIN[ii][jj] / (nzn*ZLEN-2));
			STDJDOTE[ii][jj] = sqrt(STDJDOTE[ii][jj] / (nzn*ZLEN-2));
		}

}

// Averaging in Y

void averageY(double** EXB, double** EYB, double** EZB, double*** EX, double*** EY, double*** EZ){

	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][kk] = 0.0;
			EYB[ii][kk] = 0.0;
			EZB[ii][kk] = 0.0;
		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				EXB[ii][kk] += EX[ii][jj][kk];
				EYB[ii][kk] += EY[ii][jj][kk];
				EZB[ii][kk] += EZ[ii][jj][kk];
			}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][kk] = EXB[ii][kk] / nyn/YLEN;
			EYB[ii][kk] = EYB[ii][kk] / nyn/YLEN;
			EZB[ii][kk] = EZB[ii][kk] / nyn/YLEN;
		}
}

void averageY(double** EXB, double*** EX){

	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][kk] = 0.0;

		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				EXB[ii][kk] += EX[ii][jj][kk];

			}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][kk] = EXB[ii][kk] / nyn/YLEN;

		}
}


void averageYSTD(double** EXB, double** EXSTD, double*** EX){

	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][kk] = 0.0;

		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				EXB[ii][kk] += EX[ii][jj][kk];

			}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXB[ii][kk] = EXB[ii][kk] / nyn/YLEN;

		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXSTD[ii][kk] = 0.0;

		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				EXSTD[ii][kk] += pow(EX[ii][jj][kk]-EXB[ii][kk],2);
            }
    for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EXSTD[ii][kk] = sqrt(EXSTD[ii][kk] / nyn/YLEN);

		}

}

void averageYdot(double** EdotJX, double** EdotJY, double** EdotJZ,
				double** dEdotJX, double** dEdotJY, double** dEdotJZ,
				double** EXB, double** EYB, double** EZB,
				double** VXB, double** VYB, double** VZB,
				double*** EX, double*** EY, double*** EZ,
				double*** VX, double*** VY, double*** VZ){

	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EdotJX[ii][kk] = EXB[ii][kk] * VXB[ii][kk];
			EdotJY[ii][kk] = EYB[ii][kk] * VYB[ii][kk];
			EdotJZ[ii][kk] = EZB[ii][kk] * VZB[ii][kk];
			dEdotJX[ii][kk] = 0.0;
			dEdotJY[ii][kk] = 0.0;
			dEdotJZ[ii][kk] = 0.0;
		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				dEdotJX[ii][kk] += (EX[ii][jj][kk]-EXB[ii][kk]) * (VX[ii][jj][kk]-VXB[ii][kk]);
				dEdotJY[ii][kk] += (EY[ii][jj][kk]-EYB[ii][kk]) * (VY[ii][jj][kk]-VYB[ii][kk]);
				dEdotJZ[ii][kk] += (EZ[ii][jj][kk]-EZB[ii][kk]) * (VZ[ii][jj][kk]-VZB[ii][kk]);
			}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dEdotJX[ii][kk] = dEdotJX[ii][kk] / nyn/YLEN;
			dEdotJY[ii][kk] = dEdotJY[ii][kk] / nyn/YLEN;
			dEdotJZ[ii][kk] = dEdotJZ[ii][kk] / nyn/YLEN;
		}

}

void averageYdot(double** EdotJX, double** EdotJY, double** EdotJZ,
				double** dEdotJX, double** dEdotJY, double** dEdotJZ,
				double** sEdotJX, double** sEdotJY, double** sEdotJZ,
				double** EXB, double** EYB, double** EZB,
				double** VXB, double** VYB, double** VZB,
				double*** EX, double*** EY, double*** EZ,
				double*** VX, double*** VY, double*** VZ){

	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			EdotJX[ii][kk] = EXB[ii][kk] * VXB[ii][kk];
			EdotJY[ii][kk] = EYB[ii][kk] * VYB[ii][kk];
			EdotJZ[ii][kk] = EZB[ii][kk] * VZB[ii][kk];
			dEdotJX[ii][kk] = 0.0;
			dEdotJY[ii][kk] = 0.0;
			dEdotJZ[ii][kk] = 0.0;
			sEdotJX[ii][kk] = 0.0;
			sEdotJY[ii][kk] = 0.0;
			sEdotJZ[ii][kk] = 0.0;
		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				dEdotJX[ii][kk] += (EX[ii][jj][kk]-EXB[ii][kk]) * (VX[ii][jj][kk]-VXB[ii][kk]);
				dEdotJY[ii][kk] += (EY[ii][jj][kk]-EYB[ii][kk]) * (VY[ii][jj][kk]-VYB[ii][kk]);
				dEdotJZ[ii][kk] += (EZ[ii][jj][kk]-EZB[ii][kk]) * (VZ[ii][jj][kk]-VZB[ii][kk]);
			}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dEdotJX[ii][kk] = dEdotJX[ii][kk] / nyn/YLEN;
			dEdotJY[ii][kk] = dEdotJY[ii][kk] / nyn/YLEN;
			dEdotJZ[ii][kk] = dEdotJZ[ii][kk] / nyn/YLEN;
		}
    for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				sEdotJX[ii][kk] += pow(EX[ii][jj][kk] * VX[ii][jj][kk]
                                    - EdotJX[ii][kk] - dEdotJX[ii][kk], 2);
				sEdotJY[ii][kk] += pow(EY[ii][jj][kk] * VY[ii][jj][kk]
                                       - EdotJY[ii][kk] - dEdotJY[ii][kk], 2);
				sEdotJZ[ii][kk] += pow(EZ[ii][jj][kk] * VZ[ii][jj][kk]
                                       - EdotJZ[ii][kk] - dEdotJZ[ii][kk], 2);
            }
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			sEdotJX[ii][kk] = sqrt(sEdotJX[ii][kk] / nyn/YLEN);
			sEdotJY[ii][kk] = sqrt(sEdotJY[ii][kk] / nyn/YLEN);
			sEdotJZ[ii][kk] = sqrt(sEdotJZ[ii][kk] / nyn/YLEN);
		}

}


void averageYcross(double** CX, double** CY, double** CZ,
				double** dCX, double** dCY, double** dCZ,
				double** EXB, double** EYB, double** EZB,
				double** BXB, double** BYB, double** BZB,
				double*** EX, double*** EY, double*** EZ,
				double*** BX, double*** BY, double*** BZ){

	cross(EXB,EYB,EZB,BXB,BYB,BZB,CX,CY,CZ);

	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dCX[ii][kk] = 0.0;
			dCY[ii][kk] = 0.0;
			dCZ[ii][kk] = 0.0;
		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){

				dCX[ii][kk] += (EY[ii][jj][kk] - EYB[ii][kk]) *
                                               (BZ[ii][jj][kk] - BZB[ii][kk]) -
                                               (EZ[ii][jj][kk] - EZB[ii][kk]) *
					       (BY[ii][jj][kk] - BYB[ii][kk]);

				dCY[ii][kk] += (EZ[ii][jj][kk] - EZB[ii][kk]) *
                                               (BX[ii][jj][kk] - BXB[ii][kk]) -
                                               (EX[ii][jj][kk] - EXB[ii][kk]) *
					       (BZ[ii][jj][kk] - BZB[ii][kk]);

				dCZ[ii][kk] += (EX[ii][jj][kk] - EXB[ii][kk]) *
                                               (BY[ii][jj][kk] - BYB[ii][kk]) -
                                               (EY[ii][jj][kk] - EYB[ii][kk]) *
					       (BX[ii][jj][kk] - BXB[ii][kk]);
			}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dCX[ii][kk] = dCX[ii][kk] / nyn/YLEN;
			dCY[ii][kk] = dCY[ii][kk] / nyn/YLEN;
			dCZ[ii][kk] = dCZ[ii][kk] / nyn/YLEN;
		}

}


void averageYcross(double** CX, double** CY, double** CZ,
                  double** dCX, double** dCY, double** dCZ,
                  double** sCX, double** sCY, double** sCZ,
                  double** EXB, double** EYB, double** EZB,
                  double** BXB, double** BYB, double** BZB,
                  double*** EX, double*** EY, double*** EZ,
                  double*** BX, double*** BY, double*** BZ){

	cross(EXB,EYB,EZB,BXB,BYB,BZB,CX,CY,CZ);

	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dCX[ii][kk] = 0.0;
			dCY[ii][kk] = 0.0;
			dCZ[ii][kk] = 0.0;
            sCX[ii][kk] = 0.0;
			sCY[ii][kk] = 0.0;
			sCZ[ii][kk] = 0.0;
		}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){

				dCX[ii][kk] += (EY[ii][jj][kk] - EYB[ii][kk]) *
                (BZ[ii][jj][kk] - BZB[ii][kk]) -
                (EZ[ii][jj][kk] - EZB[ii][kk]) *
                (BY[ii][jj][kk] - BYB[ii][kk]);

				dCY[ii][kk] += (EZ[ii][jj][kk] - EZB[ii][kk]) *
                (BX[ii][jj][kk] - BXB[ii][kk]) -
                (EX[ii][jj][kk] - EXB[ii][kk]) *
                (BZ[ii][jj][kk] - BZB[ii][kk]);

				dCZ[ii][kk] += (EX[ii][jj][kk] - EXB[ii][kk]) *
                (BY[ii][jj][kk] - BYB[ii][kk]) -
                (EY[ii][jj][kk] - EYB[ii][kk]) *
                (BX[ii][jj][kk] - BXB[ii][kk]);
			}
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			dCX[ii][kk] = dCX[ii][kk] / nyn/YLEN;
			dCY[ii][kk] = dCY[ii][kk] / nyn/YLEN;
			dCZ[ii][kk] = dCZ[ii][kk] / nyn/YLEN;
		}
    for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){

				sCX[ii][kk] += pow(EY[ii][jj][kk] * BZ[ii][jj][kk]  -
                                   EZ[ii][jj][kk] * BY[ii][jj][kk] - CX[ii][kk]
                                   -dCX[ii][kk],2);

                sCY[ii][jj] += pow(EZ[ii][jj][kk] * BX[ii][jj][kk]  -
                                   EX[ii][jj][kk] * BZ[ii][jj][kk] - CY[ii][kk]
                                   -dCY[ii][kk],2);

                sCZ[ii][jj] += pow(EX[ii][jj][kk] * BY[ii][jj][kk]  -
                                   EY[ii][jj][kk] * BX[ii][jj][kk] - CZ[ii][kk]
                                   -dCZ[ii][kk],2);

			}
    for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			sCX[ii][kk] = sqrt(sCX[ii][kk] / nyn/YLEN);
			sCY[ii][kk] = sqrt(sCY[ii][kk] / nyn/YLEN);
			sCZ[ii][kk] = sqrt(sCZ[ii][kk] / nyn/YLEN);
		}


}


void averageY_energy(double** AVGENTH, double** AVGENIN, double** AVGJDOTE,
		double** STDENTH, double** STDENIN, double** STDJDOTE,
		double*** VDIVP, double*** DIVJV, double*** JDOTE){

	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			AVGJDOTE[ii][kk] = 0.0;
			STDJDOTE[ii][kk] = 0.0;
			AVGENTH[ii][kk] = 0.0;
			AVGENIN[ii][kk] = 0.0;
			STDENTH[ii][kk] = 0.0;
			STDENIN[ii][kk] = 0.0;
		}

    for (int kk=1; kk < nzn*ZLEN-1;kk++)
		for (int jj=1; jj < nyn*YLEN-1;jj++)
			for (int ii=1; ii < nxn*XLEN-1;ii++){

				AVGENTH[ii][kk] += VDIVP[ii][jj][kk]  ;
				AVGENIN[ii][kk] += (JDOTE[ii][jj][kk] - VDIVP[ii][jj][kk])  ;
				AVGJDOTE[ii][kk] += JDOTE[ii][jj][kk]   ;
			}

    for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			AVGENTH[ii][kk] = AVGENTH[ii][kk] / (nyn*YLEN-2);
			AVGENIN[ii][kk] = AVGENIN[ii][kk] / (nyn*YLEN-2);
			AVGJDOTE[ii][kk] = AVGJDOTE[ii][kk] / (nyn*YLEN-2);
		}


    for (int kk=1; kk < nzn*ZLEN-1;kk++)
		for (int jj=1; jj < nyn*YLEN-1;jj++)
			for (int ii=1; ii < nxn*XLEN-1;ii++){

				STDENTH[ii][kk] += pow(VDIVP[ii][jj][kk]  - AVGENTH[ii][jj],2);
				STDENIN[ii][kk] += pow(JDOTE[ii][jj][kk] - VDIVP[ii][jj][kk]  - AVGENIN[ii][kk],2);
				STDJDOTE[ii][kk] += pow(JDOTE[ii][jj][kk]  - AVGJDOTE[ii][kk],2);

			}
    for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			STDENTH[ii][kk] = sqrt(STDENTH[ii][kk] / (nyn*YLEN-2));
			STDENIN[ii][kk] = sqrt(STDENIN[ii][kk] / (nyn*YLEN-2));
			STDJDOTE[ii][kk] = sqrt(STDJDOTE[ii][kk] / (nyn*YLEN-2));
		}

}

void cross(double*** EX, double*** EY, double*** EZ,
		   double*** BX, double*** BY, double*** BZ,
		   double*** VXBX, double*** VXBY, double*** VXBZ){
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){

				VXBX[ii][jj][kk] = (EY[ii][jj][kk]*BZ[ii][jj][kk] - EZ[ii][jj][kk]*BY[ii][jj][kk]);
				VXBY[ii][jj][kk] = (EZ[ii][jj][kk]*BX[ii][jj][kk] - EX[ii][jj][kk]*BZ[ii][jj][kk]);
				VXBZ[ii][jj][kk] = (EX[ii][jj][kk]*BY[ii][jj][kk] - EY[ii][jj][kk]*BX[ii][jj][kk]);

			}
}

void crossoverb(double*** EX, double*** EY, double*** EZ,
		   double*** BX, double*** BY, double*** BZ,
		   double*** VXBX, double*** VXBY, double*** VXBZ){
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				double b = 1e-10 + BX[ii][jj][kk]*BX[ii][jj][kk] + 
				BY[ii][jj][kk]*BY[ii][jj][kk] +
				BZ[ii][jj][kk]*BZ[ii][jj][kk];
				VXBX[ii][jj][kk] = (EY[ii][jj][kk]*BZ[ii][jj][kk] - EZ[ii][jj][kk]*BY[ii][jj][kk])/b;
				VXBY[ii][jj][kk] = (EZ[ii][jj][kk]*BX[ii][jj][kk] - EX[ii][jj][kk]*BZ[ii][jj][kk])/b;
				VXBZ[ii][jj][kk] = (EX[ii][jj][kk]*BY[ii][jj][kk] - EY[ii][jj][kk]*BX[ii][jj][kk])/b;
				
			}
}
void div_cross(double*** EX, double*** EY, double*** EZ,
		   double*** BX, double*** BY, double*** BZ,
		   double*** SX, double*** SY, double*** SZ,
		   double*** divS){
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				SX[ii][jj][kk] = (EY[ii][jj][kk]*BZ[ii][jj][kk] - EZ[ii][jj][kk]*BY[ii][jj][kk]);
				SY[ii][jj][kk] = (EZ[ii][jj][kk]*BX[ii][jj][kk] - EX[ii][jj][kk]*BZ[ii][jj][kk]);
				SZ[ii][jj][kk] = (EX[ii][jj][kk]*BY[ii][jj][kk] - EY[ii][jj][kk]*BX[ii][jj][kk]);

			}
	divergence(divS, SX, SY, SZ);
}

void cross(double** EX, double** EY, double** EZ,
		   double** BX, double** BY, double** BZ,
		   double** VXBX, double** VXBY, double** VXBZ){
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				VXBX[ii][jj] = (EY[ii][jj]*BZ[ii][jj] - EZ[ii][jj]*BY[ii][jj]);
				VXBY[ii][jj] = (EZ[ii][jj]*BX[ii][jj] - EX[ii][jj]*BZ[ii][jj]);
				VXBZ[ii][jj] = (EX[ii][jj]*BY[ii][jj] - EY[ii][jj]*BX[ii][jj]);
				
			}
}

void cross(int n, double* EX, double* EY, double* EZ,
		   double* BX, double* BY, double* BZ,
		   double* VXBX, double* VXBY, double* VXBZ){
		for (int ii=0; ii < n;ii++){
				VXBX[ii] = (EY[ii]*BZ[ii] - EZ[ii]*BY[ii]);
				VXBY[ii] = (EZ[ii]*BX[ii] - EX[ii]*BZ[ii]);
				VXBZ[ii] = (EX[ii]*BY[ii] - EY[ii]*BX[ii]);

			}
}

void cross(int m, double* A, double** B, double** AXB){
	for (int jj=0; jj< m;jj++ ){
				AXB[jj][0] = (A[1]*B[jj][2] - A[2]*B[jj][1]);
				AXB[jj][1] = (A[2]*B[jj][0] - A[0]*B[jj][2]);
				AXB[jj][2] = (A[0]*B[jj][1] - A[1]*B[jj][0]);

			}
}

void cross_trans(int m, double* A, double** B, double** AXB){
	for (int jj=0; jj< m;jj++ ){
				AXB[0][jj] = (A[1]*B[2][jj] - A[2]*B[1][jj]);
				AXB[1][jj] = (A[2]*B[0][jj] - A[0]*B[2][jj]);
				AXB[2][jj] = (A[0]*B[1][jj] - A[1]*B[0][jj]);

			}
}

void vnonfrozen(double*** N, double*** BX, double*** BY, double*** BZ, 
				double*** EX, double*** EY, double*** EZ, 
				double*** VX, double*** VY, double*** VZ){
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				double b = 1e-10 + BX[ii][jj][kk] * BX[ii][jj][kk] + 
				BY[ii][jj][kk] * BY[ii][jj][kk] +
				BZ[ii][jj][kk] * BZ[ii][jj][kk];
				
				double vdotB = (VX[ii][jj][kk] * BX[ii][jj][kk] +
						VY[ii][jj][kk] * BY[ii][jj][kk] +
						VZ[ii][jj][kk] * BZ[ii][jj][kk]) / N[ii][jj][kk];

				VX[ii][jj][kk] = VX[ii][jj][kk] / N[ii][jj][kk] - 
				(EY[ii][jj][kk]*BZ[ii][jj][kk] - EZ[ii][jj][kk]*BY[ii][jj][kk])/b;
				
				VY[ii][jj][kk] = VY[ii][jj][kk] / N[ii][jj][kk] -
				(EZ[ii][jj][kk]*BX[ii][jj][kk] - EX[ii][jj][kk]*BZ[ii][jj][kk])/b;
				
				VZ[ii][jj][kk] = VZ[ii][jj][kk] / N[ii][jj][kk] -
				(EX[ii][jj][kk]*BY[ii][jj][kk] - EY[ii][jj][kk]*BX[ii][jj][kk])/b;
// Removal of Vpar Added on Dec 30, 2014. Prior it was not subtracted.
				VX[ii][jj][kk] = VX[ii][jj][kk] - vdotB * BX[ii][jj][kk] / b;
				VY[ii][jj][kk] = VY[ii][jj][kk] - vdotB * BY[ii][jj][kk] / b;
				VZ[ii][jj][kk] = VZ[ii][jj][kk] - vdotB * BZ[ii][jj][kk] / b;
				
			}
}

void ohm(double*** N, double*** BX, double*** BY, double*** BZ, 
		 double*** EX, double*** EY, double*** EZ, 
		 double*** VX, double*** VY, double*** VZ,
		 double*** OHMX, double*** OHMY, double*** OHMZ){
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				
				OHMX[ii][jj][kk] = EX[ii][jj][kk]  + 
				(VY[ii][jj][kk]*BZ[ii][jj][kk] - VZ[ii][jj][kk]*BY[ii][jj][kk]) / N[ii][jj][kk];
				
				OHMY[ii][jj][kk] = EY[ii][jj][kk]  +
				(VZ[ii][jj][kk]*BX[ii][jj][kk] - VX[ii][jj][kk]*BZ[ii][jj][kk]) / N[ii][jj][kk];
				
				OHMZ[ii][jj][kk] = EZ[ii][jj][kk]  +
				(VX[ii][jj][kk]*BY[ii][jj][kk] - VY[ii][jj][kk]*BX[ii][jj][kk]) / N[ii][jj][kk];
				
			}
}

void zenitani(double*** NE, double*** NI, 
 		 double*** BX, double*** BY, double*** BZ, 
		 double*** EX, double*** EY, double*** EZ, 
		 double*** JEX, double*** JEY, double*** JEZ,
		 double*** JIX, double*** JIY, double*** JIZ,
		 double*** DE, double*** DI){
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				
// Compute E + ve x B
				double OHMEX = EX[ii][jj][kk]  + 
				(JEY[ii][jj][kk]*BZ[ii][jj][kk] - JEZ[ii][jj][kk]*BY[ii][jj][kk]) / NE[ii][jj][kk];
				
				double OHMEY = EY[ii][jj][kk]  +
				(JEZ[ii][jj][kk]*BX[ii][jj][kk] - JEX[ii][jj][kk]*BZ[ii][jj][kk]) / NE[ii][jj][kk];
				
				double OHMEZ = EZ[ii][jj][kk]  +
				(JEX[ii][jj][kk]*BY[ii][jj][kk] - JEY[ii][jj][kk]*BX[ii][jj][kk]) / NE[ii][jj][kk];

// compute E + vi x B
				double OHMIX = EX[ii][jj][kk]  + 
				(JIY[ii][jj][kk]*BZ[ii][jj][kk] - JIZ[ii][jj][kk]*BY[ii][jj][kk]) / NI[ii][jj][kk];
				
				double OHMIY = EY[ii][jj][kk]  +
				(JIZ[ii][jj][kk]*BX[ii][jj][kk] - JIX[ii][jj][kk]*BZ[ii][jj][kk]) / NI[ii][jj][kk];
				
				double OHMIZ = EZ[ii][jj][kk]  +
				(JIX[ii][jj][kk]*BY[ii][jj][kk] - JIY[ii][jj][kk]*BX[ii][jj][kk]) / NI[ii][jj][kk];

// Compute rho_net
				double rhonet= NE[ii][jj][kk] + NI[ii][jj][kk];

// Compute De
				DE[ii][jj][kk] = (JEX[ii][jj][kk] + JIX[ii][jj][kk]) * OHMEX +
						 (JEY[ii][jj][kk] + JIY[ii][jj][kk]) * OHMEY +
						 (JEZ[ii][jj][kk] + JIZ[ii][jj][kk]) * OHMEZ -
						 rhonet * ( JEX[ii][jj][kk] * EX[ii][jj][kk] + JEY[ii][jj][kk] * EY[ii][jj][kk] + JEZ[ii][jj][kk] * EZ[ii][jj][kk]) / NE[ii][jj][kk]; 
				
// Compute Di
				DI[ii][jj][kk] = (JEX[ii][jj][kk] + JIX[ii][jj][kk]) * OHMIX +
						 (JEY[ii][jj][kk] + JIY[ii][jj][kk]) * OHMIY +
						 (JEZ[ii][jj][kk] + JIZ[ii][jj][kk]) * OHMIZ -
						 rhonet * ( JIX[ii][jj][kk] * EX[ii][jj][kk] + JIY[ii][jj][kk] * EY[ii][jj][kk] + JIZ[ii][jj][kk] * EZ[ii][jj][kk]) / NI[ii][jj][kk]; 

			}
}




void extract_pressure(double qom, double*** BX, double*** BY, double*** BZ,
					  double*** VX, double*** VY, double*** VZ,
					  double*** N, double*** pXX, double*** pXY,
					  double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ,
					  double*** pPAR, double*** pPER1, double*** pPER2, double*** EPS) {
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				if(N[ii][jj][kk]!=0.0)
					
					pXX[ii][jj][kk] = pXX[ii][jj][kk]-VX[ii][jj][kk]*VX[ii][jj][kk]/N[ii][jj][kk];
				pXX[ii][jj][kk] = pXX[ii][jj][kk] / qom;
				if(N[ii][jj][kk]!=0.0)
					pXY[ii][jj][kk] = pXY[ii][jj][kk]-VX[ii][jj][kk]*VY[ii][jj][kk]/N[ii][jj][kk];
				pXY[ii][jj][kk] = pXY[ii][jj][kk] / qom;
				if(N[ii][jj][kk]!=0.0)
					pXZ[ii][jj][kk] = pXZ[ii][jj][kk]-VX[ii][jj][kk]*VZ[ii][jj][kk]/N[ii][jj][kk];
				pXZ[ii][jj][kk] = pXZ[ii][jj][kk] / qom;
				if(N[ii][jj][kk]!=0.0)
					pYY[ii][jj][kk] = pYY[ii][jj][kk]-VY[ii][jj][kk]*VY[ii][jj][kk]/N[ii][jj][kk];
				pYY[ii][jj][kk] = pYY[ii][jj][kk] / qom;
				if(N[ii][jj][kk]!=0.0)
					pYZ[ii][jj][kk] = pYZ[ii][jj][kk]-VY[ii][jj][kk]*VZ[ii][jj][kk]/N[ii][jj][kk];
				pYZ[ii][jj][kk] = pYZ[ii][jj][kk] / qom;
				if(N[ii][jj][kk]!=0.0)
					pZZ[ii][jj][kk] = pZZ[ii][jj][kk]-VZ[ii][jj][kk]*VZ[ii][jj][kk]/N[ii][jj][kk];
				pZZ[ii][jj][kk] = pZZ[ii][jj][kk] / qom;
				
				double b2D = 1e-10 + BX[ii][jj][kk]*BX[ii][jj][kk] + BY[ii][jj][kk]*BY[ii][jj][kk];
				double b = b2D + BZ[ii][jj][kk]*BZ[ii][jj][kk];
				double perp2x = BZ[ii][jj][kk]*BX[ii][jj][kk] /sqrt(b*b2D);
				double perp2y = BZ[ii][jj][kk]*BY[ii][jj][kk] /sqrt(b*b2D);
				double perp2z = -sqrt(b2D/b);
				
				pPAR[ii][jj][kk] = BX[ii][jj][kk]*pXX[ii][jj][kk]*BX[ii][jj][kk] + 
				2*BX[ii][jj][kk]*pXY[ii][jj][kk]*BY[ii][jj][kk] + 
				2*BX[ii][jj][kk]*pXZ[ii][jj][kk]*BZ[ii][jj][kk];
				pPAR[ii][jj][kk]+= BY[ii][jj][kk]*pYY[ii][jj][kk]*BY[ii][jj][kk] + 
				2*BY[ii][jj][kk]*pYZ[ii][jj][kk]*BZ[ii][jj][kk];
				pPAR[ii][jj][kk]+= BZ[ii][jj][kk]*pZZ[ii][jj][kk]*BZ[ii][jj][kk];
				
				pPAR[ii][jj][kk] = pPAR[ii][jj][kk]/b;
				
				pPER1[ii][jj][kk] = BY[ii][jj][kk]*pXX[ii][jj][kk]*BY[ii][jj][kk] - 
				2*BY[ii][jj][kk]*pXY[ii][jj][kk]*BX[ii][jj][kk];
				pPER1[ii][jj][kk]+= BX[ii][jj][kk]*pYY[ii][jj][kk]*BX[ii][jj][kk];
				
				pPER1[ii][jj][kk] = pPER1[ii][jj][kk]/b2D;
				
				pPER2[ii][jj][kk] = perp2x*pXX[ii][jj][kk]*perp2x + 2*perp2x*pXY[ii][jj][kk]*perp2y + 2*perp2x*pXZ[ii][jj][kk]*perp2z;
				pPER2[ii][jj][kk]+= perp2y*pYY[ii][jj][kk]*perp2y + 2*perp2y*pYZ[ii][jj][kk]*perp2z;
				pPER2[ii][jj][kk]+= perp2z*pZZ[ii][jj][kk]*perp2z;  
				
				EPS[ii][jj][kk] = 1.0 - 4.0* M_PI * ( pPAR[ii][jj][kk] - sqrt(pPER1[ii][jj][kk] * pPER2[ii][jj][kk] ) ) / b;
				
			}
}

void bulk_energy_flux(double*** Ubulk,
		  	  	  	  double*** QXbulk, double*** QYbulk, double*** QZbulk,
					  double*** VX, double*** VY, double*** VZ){

	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){


				QXbulk[ii][jj][kk] = VX[ii][jj][kk]  * Ubulk[ii][jj][kk];

				QYbulk[ii][jj][kk] = VY[ii][jj][kk]  * Ubulk[ii][jj][kk];

				QZbulk[ii][jj][kk] = VZ[ii][jj][kk]  * Ubulk[ii][jj][kk];

			}
}

void internal_energy_flux(double*** Uth,
					  double*** QX, double*** QY, double*** QZ,
					  double*** VX, double*** VY, double*** VZ,
					  double*** pXX, double*** pXY,
					  double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ){

	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){

				QX[ii][jj][kk] = VX[ii][jj][kk] * Uth[ii][jj][kk];
				QX[ii][jj][kk] += VX[ii][jj][kk] * pXX[ii][jj][kk] + VY[ii][jj][kk] *pXY[ii][jj][kk] + VZ[ii][jj][kk] *pXZ[ii][jj][kk];

				QY[ii][jj][kk] = VY[ii][jj][kk] * Uth[ii][jj][kk] ;
				QY[ii][jj][kk] += VX[ii][jj][kk] * pXY[ii][jj][kk] + VY[ii][jj][kk] *pYY[ii][jj][kk] + VZ[ii][jj][kk] *pYZ[ii][jj][kk];

				QZ[ii][jj][kk] = VZ[ii][jj][kk] * Uth[ii][jj][kk] ;
				QZ[ii][jj][kk] += VX[ii][jj][kk] * pXZ[ii][jj][kk] + VY[ii][jj][kk] *pYZ[ii][jj][kk] + VZ[ii][jj][kk] *pZZ[ii][jj][kk];

			}
}

void compute_energy(double qom, double*** Ubulk, double*** Uth,
		  double*** VX, double*** VY, double*** VZ,
		  double*** N, double*** pXX, double*** pYY, double*** pZZ){
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){

				Ubulk[ii][jj][kk] = 0.5 *(
						VX[ii][jj][kk] * VX[ii][jj][kk] +
						VY[ii][jj][kk] * VY[ii][jj][kk] +
						VZ[ii][jj][kk] * VZ[ii][jj][kk] ) * N[ii][jj][kk] /qom;

				Uth[ii][jj][kk] = 0.5 *(
						pXX[ii][jj][kk] + pYY[ii][jj][kk] + pZZ[ii][jj][kk]);

			}
}

void divergence(double*** div, double*** EX, double*** EY, double*** EZ){
  
 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			    
		div[ii][jj][kk] = 0.0;
			}
  
 if(nzn*ZLEN==1) {
     int kk=0;
    	for (int jj=1; jj < nyn*YLEN-1;jj++)
    		for (int ii=1; ii < nxn*XLEN-1;ii++){

    		div[ii][jj][kk] = (EX[ii+1][jj][kk] - EX[ii-1][jj][kk])/2.0/Dx +
    			              (EY[ii][jj+1][kk] - EY[ii][jj-1][kk])/2.0/Dy ;
    			}
 }
 else
 {
 for (int kk=1; kk < nzn*ZLEN-1;kk++)
	for (int jj=1; jj < nyn*YLEN-1;jj++)
		for (int ii=1; ii < nxn*XLEN-1;ii++){
			    
		div[ii][jj][kk] = (EX[ii+1][jj][kk] - EX[ii-1][jj][kk])/2.0/Dx +
			              (EY[ii][jj+1][kk] - EY[ii][jj-1][kk])/2.0/Dy +
				          (EZ[ii][jj][kk+1] - EZ[ii][jj][kk-1])/2.0/Dz;
			}
 }
}  

void divergenceN(double*** div, double*** EX, double*** EY, double*** EZ, double*** N){

 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){

		div[ii][jj][kk] = 0.0;
			}

 if(nzn*ZLEN==1) {
     int kk=0;
    	for (int jj=1; jj < nyn*YLEN-1;jj++)
    		for (int ii=1; ii < nxn*XLEN-1;ii++){

    		div[ii][jj][kk] = (EX[ii+1][jj][kk] /N[ii+1][jj][kk] - EX[ii-1][jj][kk]/N[ii-1][jj][kk])/2.0/Dx +
    			              (EY[ii][jj+1][kk]/ N[ii][jj+1][kk] - EY[ii][jj-1][kk]/ N[ii][jj-1][kk])/2.0/Dy ;
    			}
 }
 else
 {
 for (int kk=1; kk < nzn*ZLEN-1;kk++)
	for (int jj=1; jj < nyn*YLEN-1;jj++)
		for (int ii=1; ii < nxn*XLEN-1;ii++){

		div[ii][jj][kk] = (EX[ii+1][jj][kk]/N[ii+1][jj][kk] - EX[ii-1][jj][kk]/N[ii-1][jj][kk])/2.0/Dx +
			              (EY[ii][jj+1][kk]/N[ii][jj+1][kk] - EY[ii][jj-1][kk]/N[ii][jj-1][kk])/2.0/Dy +
				          (EZ[ii][jj][kk+1]/N[ii][jj][kk+1] - EZ[ii][jj][kk-1]/N[ii][jj][kk-1])/2.0/Dz;
			}
 }
}


void divergenceNP(double*** div, double*** EX, double*** EY, double*** EZ, double*** N, double*** P){

 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){

		div[ii][jj][kk] = 0.0;
			}

 if(nzn*ZLEN==1) {
     int kk=0;
    	for (int jj=1; jj < nyn*YLEN-1;jj++)
    		for (int ii=1; ii < nxn*XLEN-1;ii++){

    		div[ii][jj][kk] = (EX[ii+1][jj][kk] * P[ii+1][jj][kk] / N[ii+1][jj][kk] -
    				           EX[ii-1][jj][kk] * P[ii-1][jj][kk] / N[ii-1][jj][kk])/2.0/Dx +
    			              (EY[ii][jj+1][kk] * P[ii][jj+1][kk] / N[ii][jj+1][kk] -
    			               EY[ii][jj-1][kk] * P[ii][jj-1][kk] / N[ii][jj-1][kk])/2.0/Dy ;
    			}
 }
 else
 {
 for (int kk=1; kk < nzn*ZLEN-1;kk++)
	for (int jj=1; jj < nyn*YLEN-1;jj++)
		for (int ii=1; ii < nxn*XLEN-1;ii++){

		div[ii][jj][kk] = (EX[ii+1][jj][kk]*P[ii+1][jj][kk]/N[ii+1][jj][kk] -
								  EX[ii-1][jj][kk]*P[ii-1][jj][kk]/N[ii-1][jj][kk])/2.0/Dx +
			              (EY[ii][jj+1][kk]*P[ii][jj+1][kk]/N[ii][jj+1][kk] -
			            		  EY[ii][jj-1][kk]*P[ii][jj-1][kk]/N[ii][jj-1][kk])/2.0/Dy +
				          (EZ[ii][jj][kk+1]*P[ii][jj][kk+1]/N[ii][jj][kk+1] -
				        		  EZ[ii][jj][kk-1]*P[ii][jj][kk-1]/N[ii][jj][kk-1])/2.0/Dz;
			}
 }
}

/* Computing U dot div P */
void udivP(double*** O, double*** OX, double*** OY, double*** OZ,
		double*** VX, double*** VY, double*** VZ,
		double*** pXX, double*** pXY,
				  double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ){

		divergence(OX, pXX, pXY, pXZ);
		prod(O, 1.0, OX, VX, nxn*XLEN, nyn*YLEN, nzn*ZLEN);

		divergence(OY, pXY, pYY, pYZ);
		addprod(O, 1.0, OY, VY, nxn*XLEN, nyn*YLEN, nzn*ZLEN);

		divergence(OZ, pXZ, pYZ, pZZ);
		addprod(O, 1.0, OZ, VZ, nxn*XLEN, nyn*YLEN, nzn*ZLEN);

}

/* Compute P:gradU */
void pgradU(double*** PgradU, double*** VX, double*** VY, double*** VZ,
		double*** WX, double*** WY, double*** WZ,
		double*** TXX, double***  TXY, double*** TXZ,
		double*** TYY, double*** TYZ, double*** TZZ){

	gradient(WX, WY, WZ, VX);

	dot(PgradU, WX, WY, WZ, TXX, TXY, TXZ);


	gradient(WX, WY, WZ, VY);

	adddot(PgradU, WX, WY, WZ, TXY, TYY, TYZ);


	gradient(WX, WY, WZ, VZ);

	adddot(PgradU, WX, WY, WZ, TXZ, TYZ, TZZ);

}


void LEBT(double*** BX, double*** BY, double*** BZ,
		 double*** SX, double*** SY, double*** SZ,
		 double*** WX, double*** WY, double*** WZ,
		 double*** psi){
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){

			psi[ii][jj][kk] = 0.0;

				}

	// Compute log Bnorm
	 for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){

		    double bnorm2 = BX[ii][jj][kk] * BX[ii][jj][kk] +
		    			   BY[ii][jj][kk] * BY[ii][jj][kk] +
		    			   BZ[ii][jj][kk] * BZ[ii][jj][kk];



			psi[ii][jj][kk] = log (bnorm2);

			}
	// Take the gradient and dot it
	   	   gradient(WX, WY, WZ, psi);
	   	   dot(psi, WX, WY, WZ, SX, SY, SZ);
}

void biskamp(double*** BX, double*** BY, double*** BZ,
		 double*** EX, double*** EY, double*** EZ,
		 double*** WX, double*** WY, double*** WZ,
		 double*** BiskampX, double*** BiskampY, double*** BiskampZ){

 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){

		BiskampX[ii][jj][kk] = 0.0;
		BiskampY[ii][jj][kk] = 0.0;
		BiskampZ[ii][jj][kk] = 0.0;
			}

// Compute the vector E parallel
 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){

	    double bnorm2 = BX[ii][jj][kk] * BX[ii][jj][kk] +
	    			   BY[ii][jj][kk] * BY[ii][jj][kk] +
	    			   BZ[ii][jj][kk] * BZ[ii][jj][kk];

	    double EdotB = BX[ii][jj][kk] * EX[ii][jj][kk] +
 			   BY[ii][jj][kk] * EY[ii][jj][kk] +
 			   BZ[ii][jj][kk] * EZ[ii][jj][kk];

		BiskampX[ii][jj][kk] = BX[ii][jj][kk] * EdotB / bnorm2;
		BiskampY[ii][jj][kk] = BY[ii][jj][kk] * EdotB / bnorm2;
		BiskampZ[ii][jj][kk] = BZ[ii][jj][kk] * EdotB / bnorm2;

		}
// Take the curl of Eparallel
   	   curl(WX, WY, WZ, BiskampX, BiskampY, BiskampZ);
   	   cross(BX, BY, BZ, WX, WY, WZ, BiskampX, BiskampY, BiskampZ);
}


void gradient(double*** EX, double*** EY, double*** EZ, double*** phi){

 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){

		EX[ii][jj][kk] = 0.0;
		EY[ii][jj][kk] = 0.0;
		EZ[ii][jj][kk] = 0.0;
			}

 if(nzn*ZLEN==1) {
     int kk=0;
    	for (int jj=1; jj < nyn*YLEN-1;jj++)
    		for (int ii=1; ii < nxn*XLEN-1;ii++){

    			EX[ii][jj][kk] = (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx;

    			EY[ii][jj][kk] = (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy;

    			}
 }
 else
 {
 for (int kk=1; kk < nzn*ZLEN-1;kk++)
	for (int jj=1; jj < nyn*YLEN-1;jj++)
		for (int ii=1; ii < nxn*XLEN-1;ii++){

		EX[ii][jj][kk] = (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx;

		EY[ii][jj][kk] = (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy;

		EZ[ii][jj][kk] = (phi[ii][jj][kk+1] - phi[ii][jj][kk-1])/2.0/Dz;

			}
 }
}

void p_dot_gradient(double*** EX, double*** EY, double*** EZ,
		  double*** pXX, double*** pXY,
		  double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ,
		  double*** phi){

 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){

		EX[ii][jj][kk] = 0.0;
		EY[ii][jj][kk] = 0.0;
		EZ[ii][jj][kk] = 0.0;
			}

 if(nzn*ZLEN==1) {
     int kk=0;
    	for (int jj=1; jj < nyn*YLEN-1;jj++)
    		for (int ii=1; ii < nxn*XLEN-1;ii++){

    			EX[ii][jj][kk] = pXX[ii][jj][kk] * (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx +
    					         pXY[ii][jj][kk] * (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy ;

    			EY[ii][jj][kk] = pXY[ii][jj][kk] * (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx +
    					         pYY[ii][jj][kk] * (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy ;

    			EZ[ii][jj][kk] = pXZ[ii][jj][kk] * (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx +
    					         pYZ[ii][jj][kk] * (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy ;
    			}
 }
 else
 {
 for (int kk=1; kk < nzn*ZLEN-1;kk++)
	for (int jj=1; jj < nyn*YLEN-1;jj++)
		for (int ii=1; ii < nxn*XLEN-1;ii++){

		EX[ii][jj][kk] = pXX[ii][jj][kk] * (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx +
				         pXY[ii][jj][kk] * (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy +
						 pXZ[ii][jj][kk] * (phi[ii][jj][kk+1] - phi[ii][jj][kk-1])/2.0/Dz;

		EY[ii][jj][kk] = pXY[ii][jj][kk] * (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx +
				         pYY[ii][jj][kk] * (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy +
						 pYZ[ii][jj][kk] * (phi[ii][jj][kk+1] - phi[ii][jj][kk-1])/2.0/Dz;

		EZ[ii][jj][kk] = pXZ[ii][jj][kk] * (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx +
				         pYZ[ii][jj][kk] * (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy +
						 pZZ[ii][jj][kk] * (phi[ii][jj][kk+1] - phi[ii][jj][kk-1])/2.0/Dz;
			}
 }
}


void gradient_par(double*** BX, double*** BY, double*** BZ,
		  double*** phi, double*** out){


 if(nzn*ZLEN==1) {
     int kk=0;
    	for (int jj=1; jj < nyn*YLEN-1;jj++)
    		for (int ii=1; ii < nxn*XLEN-1;ii++){

    			double EX = (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx;

    			double EY = (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy;

    			out[ii][jj][kk] = EX * BX[ii][jj][kk] + EY * BY[ii][jj][kk] ;
    			out[ii][jj][kk] = out[ii][jj][kk] /  sqrt(BX[ii][jj][kk] * BX[ii][jj][kk] +
    					          BY[ii][jj][kk] * BY[ii][jj][kk]);
    			}
 }
 else
 {
 for (int kk=1; kk < nzn*ZLEN-1;kk++)
	for (int jj=1; jj < nyn*YLEN-1;jj++)
		for (int ii=1; ii < nxn*XLEN-1;ii++){

		double EX = (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx;

		double EY = (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy;

		double EZ = (phi[ii][jj][kk+1] - phi[ii][jj][kk-1])/2.0/Dz;

		out[ii][jj][kk] = EX * BX[ii][jj][kk] + EY * BY[ii][jj][kk] + EZ * BZ[ii][jj][kk];
		out[ii][jj][kk] = out[ii][jj][kk] /  sqrt(BX[ii][jj][kk] * BX[ii][jj][kk] +
				          BY[ii][jj][kk] * BY[ii][jj][kk] +
						  BZ[ii][jj][kk] * BZ[ii][jj][kk]);

			}
 }
}



void gradient_per1(double*** BX, double*** BY, double*** BZ,
		  double*** phi, double*** out){


 if(nzn*ZLEN==1) {
     int kk=0;
    	for (int jj=1; jj < nyn*YLEN-1;jj++)
    		for (int ii=1; ii < nxn*XLEN-1;ii++){

    			double EX = (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx;

    			double EY = (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy;

    			out[ii][jj][kk] = EX * BY[ii][jj][kk] - EY * BX[ii][jj][kk] ;
    			out[ii][jj][kk] = out[ii][jj][kk] /  sqrt(BX[ii][jj][kk] * BX[ii][jj][kk] +
    					          BY[ii][jj][kk] * BY[ii][jj][kk]);
    			}
 }
 else
 {
 for (int kk=1; kk < nzn*ZLEN-1;kk++)
	for (int jj=1; jj < nyn*YLEN-1;jj++)
		for (int ii=1; ii < nxn*XLEN-1;ii++){

		double EX = (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx;

		double EY = (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy;

		double EZ = (phi[ii][jj][kk+1] - phi[ii][jj][kk-1])/2.0/Dz;

		out[ii][jj][kk] = EX * BY[ii][jj][kk] - EY * BX[ii][jj][kk];
		out[ii][jj][kk] = out[ii][jj][kk] /  sqrt(BX[ii][jj][kk] * BX[ii][jj][kk] +
				          BY[ii][jj][kk] * BY[ii][jj][kk]);

			}
 }
}

void gradient_per2(double*** BX, double*** BY, double*** BZ,
		  double*** phi, double*** out){


 if(nzn*ZLEN==1) {
     int kk=0;
    	for (int jj=1; jj < nyn*YLEN-1;jj++)
    		for (int ii=1; ii < nxn*XLEN-1;ii++){

    			double EX = (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx;

    			double EY = (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy;

    			out[ii][jj][kk] = EX * BX[ii][jj][kk] + EY * BY[ii][jj][kk] ;
    			out[ii][jj][kk] = out[ii][jj][kk] /  sqrt(BX[ii][jj][kk] * BX[ii][jj][kk] +
    					          BY[ii][jj][kk] * BY[ii][jj][kk]);
    			}
 }
 else
 {
 for (int kk=1; kk < nzn*ZLEN-1;kk++)
	for (int jj=1; jj < nyn*YLEN-1;jj++)
		for (int ii=1; ii < nxn*XLEN-1;ii++){

		double EX = (phi[ii+1][jj][kk] - phi[ii-1][jj][kk])/2.0/Dx;

		double EY = (phi[ii][jj+1][kk] - phi[ii][jj-1][kk])/2.0/Dy;

		double EZ = (phi[ii][jj][kk+1] - phi[ii][jj][kk-1])/2.0/Dz;

		out[ii][jj][kk] = EX * BX[ii][jj][kk] + EY * BY[ii][jj][kk] + EZ * BZ[ii][jj][kk];
		out[ii][jj][kk] = out[ii][jj][kk] /  sqrt(BX[ii][jj][kk] * BX[ii][jj][kk] +
				          BY[ii][jj][kk] * BY[ii][jj][kk] +
						  BZ[ii][jj][kk] * BZ[ii][jj][kk]);

			}
 }
}
void curl(double*** curlx, double*** curly, double*** curlz, double*** EX, double*** EY, double*** EZ){

 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){

		curlx[ii][jj][kk] = 0.0;
		curly[ii][jj][kk] = 0.0;
		curlz[ii][jj][kk] = 0.0;
			}

 if(nzn*ZLEN==1) {
     int kk=0;
    	for (int jj=1; jj < nyn*YLEN-1;jj++)
    		for (int ii=1; ii < nxn*XLEN-1;ii++){

    			curlx[ii][jj][kk] = (EZ[ii][jj+1][kk] - EZ[ii][jj-1][kk])/2.0/Dy;
    			curly[ii][jj][kk] = -(EZ[ii+1][jj][kk] - EZ[ii-1][jj][kk])/2.0/Dx ;
    			curlz[ii][jj][kk] = (EY[ii+1][jj][kk] - EY[ii-1][jj][kk])/2.0/Dx -
    								(EX[ii][jj+1][kk] - EX[ii][jj-1][kk])/2.0/Dy;    			}
 }
 else
 {
 for (int kk=1; kk < nzn*ZLEN-1;kk++)
	for (int jj=1; jj < nyn*YLEN-1;jj++)
		for (int ii=1; ii < nxn*XLEN-1;ii++){

		curlx[ii][jj][kk] = (EZ[ii][jj+1][kk] - EZ[ii][jj-1][kk])/2.0/Dy -
							(EY[ii][jj][kk+1] - EY[ii][jj][kk-1])/2.0/Dz;
		curly[ii][jj][kk] = (EX[ii][jj][kk+1] - EX[ii][jj][kk-1])/2.0/Dz -
							(EZ[ii+1][jj][kk] - EZ[ii-1][jj][kk])/2.0/Dx ;
		curlz[ii][jj][kk] = (EY[ii+1][jj][kk] - EY[ii-1][jj][kk])/2.0/Dx -
							(EX[ii][jj+1][kk] - EX[ii][jj-1][kk])/2.0/Dy;

			}
 }
}


void divPv(double*** OX, double*** OY, double*** OZ,
		  double*** PXX, double*** PYY, double*** PZZ,
		  double*** PXY, double*** PXZ, double*** PYZ,
          double*** JX, double*** JY, double*** JZ, double*** N,
		  double*** WX, double*** WY, double*** WZ){

	for (int kk=1; kk < nzn*ZLEN-1;kk++)
		for (int jj=1; jj < nyn*YLEN-1;jj++)
			for (int ii=1; ii < nxn*XLEN-1;ii++){
				WX[ii][jj][kk] = PXX[ii][jj][kk] * JX[ii][jj][kk] / N[ii][jj][kk];
				WY[ii][jj][kk] = PXY[ii][jj][kk] * JX[ii][jj][kk] / N[ii][jj][kk];
				WZ[ii][jj][kk] = PXZ[ii][jj][kk] * JX[ii][jj][kk] / N[ii][jj][kk];
			}
	divergence(OX,WX, WY, WZ);

	for (int kk=1; kk < nzn*ZLEN-1;kk++)
		for (int jj=1; jj < nyn*YLEN-1;jj++)
			for (int ii=1; ii < nxn*XLEN-1;ii++){
				WX[ii][jj][kk] = PXY[ii][jj][kk] * JY[ii][jj][kk] / N[ii][jj][kk];
				WY[ii][jj][kk] = PYY[ii][jj][kk] * JY[ii][jj][kk] / N[ii][jj][kk];
				WZ[ii][jj][kk] = PYZ[ii][jj][kk] * JY[ii][jj][kk] / N[ii][jj][kk];
			}
	divergence(OY, WX, WY, WZ);

	for (int kk=1; kk < nzn*ZLEN-1;kk++)
		for (int jj=1; jj < nyn*YLEN-1;jj++)
			for (int ii=1; ii < nxn*XLEN-1;ii++){
				WX[ii][jj][kk] = PXZ[ii][jj][kk] * JZ[ii][jj][kk] / N[ii][jj][kk];
				WY[ii][jj][kk] = PYZ[ii][jj][kk] * JZ[ii][jj][kk] / N[ii][jj][kk];
				WZ[ii][jj][kk] = PZZ[ii][jj][kk] * JZ[ii][jj][kk] / N[ii][jj][kk];
			}
	divergence(OZ, WX, WY, WZ);

}

// Computes div( E VSQ /2)
void divj(double qom, double*** divj,
          double*** EX, double*** EY, double*** EZ, double*** VSQ){
  
 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			    
		divj[ii][jj][kk] = 0.0;
			}
  
 if(nzn*ZLEN==1) {
     int kk=0;
    	for (int jj=1; jj < nyn*YLEN-1;jj++)
    		for (int ii=1; ii < nxn*XLEN-1;ii++){


    		divj[ii][jj][kk] = (EX[ii+1][jj][kk] * VSQ[ii+1][jj][kk] - EX[ii-1][jj][kk] * VSQ[ii-1][jj][kk])/ 4.0/ Dx/ qom +
    			               (EY[ii][jj+1][kk] * VSQ[ii][jj+1][kk] - EY[ii][jj-1][kk] * VSQ[ii][jj-1][kk])/ 4.0/ Dy/ qom ;
    			}
     }
     else
     {
 for (int kk=1; kk < nzn*ZLEN-1;kk++)
	for (int jj=1; jj < nyn*YLEN-1;jj++)
		for (int ii=1; ii < nxn*XLEN-1;ii++){
			    
		
		divj[ii][jj][kk] = (EX[ii+1][jj][kk] * VSQ[ii+1][jj][kk] - EX[ii-1][jj][kk] * VSQ[ii-1][jj][kk])/ 4.0/ Dx/ qom +
			               (EY[ii][jj+1][kk] * VSQ[ii][jj+1][kk] - EY[ii][jj-1][kk] * VSQ[ii][jj-1][kk])/ 4.0/ Dy/ qom +
				           (EZ[ii][jj][kk+1] * VSQ[ii][jj][kk+1] - EZ[ii][jj][kk-1] * VSQ[ii][jj][kk-1])/ 4.0/ Dz/ qom;
			}
     }
}  

void adddot(double*** dot,
          double*** AX, double*** AY, double*** AZ,
          double*** BX, double*** BY, double*** BZ){


 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){


		dot[ii][jj][kk] += AX[ii][jj][kk] * BX[ii][jj][kk] +
						  AY[ii][jj][kk] * BY[ii][jj][kk] +
						  AZ[ii][jj][kk] * BZ[ii][jj][kk] ;
			}

}
void dot(double*** dot,  
          double*** AX, double*** AY, double*** AZ, 
          double*** BX, double*** BY, double*** BZ){
  
 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			    
		dot[ii][jj][kk] = 0.0;
			}
  
 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			    
		
		dot[ii][jj][kk] = AX[ii][jj][kk] * BX[ii][jj][kk] + 
						  AY[ii][jj][kk] * BY[ii][jj][kk] + 
						  AZ[ii][jj][kk] * BZ[ii][jj][kk] ;
			}
  
}  

void dot(double*** dot,  
          double*** AX, double*** AY, double*** AZ, 
          double*** BX, double*** BY, double*** BZ, double*** RHO){
  
 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			    
		dot[ii][jj][kk] = 0.0;
			}
  
 for (int kk=0; kk < nzn*ZLEN;kk++)
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			    
		
		dot[ii][jj][kk] = AX[ii][jj][kk] * BX[ii][jj][kk] /pow(RHO[ii][jj][kk],2)+ 
						  AY[ii][jj][kk] * BY[ii][jj][kk] /pow(RHO[ii][jj][kk],2)+ 
						  AZ[ii][jj][kk] * BZ[ii][jj][kk] /pow(RHO[ii][jj][kk],2);
			}
  
}  

void vdotgrad(double*** vdotgrad, double*** PHI, double*** RHO,
                         double*** JX, double*** JY, double*** JZ){
    for (int kk=0; kk < nzn*ZLEN;kk++)
        for (int jj=0; jj < nyn*YLEN;jj++)
            for (int ii=0; ii < nxn*XLEN;ii++){

                vdotgrad[ii][jj][kk] = 0.0;
			}
    if(nzn*ZLEN==1) {
        int kk=0;
            for (int jj=1; jj < nyn*YLEN-1;jj++)
                for (int ii=1; ii < nxn*XLEN-1;ii++){

                    vdotgrad[ii][jj][kk] = JX[ii][jj][kk] / (RHO[ii][jj][kk]+1e-10) *
                    		(PHI[ii+1][jj][kk] - PHI[ii-1][jj][kk])/2.0/Dx;

                    vdotgrad[ii][jj][kk] += JY[ii][jj][kk] / (RHO[ii][jj][kk]+1e-10) *
                                        		(PHI[ii][jj+1][kk] - PHI[ii][jj-1][kk])/2.0/Dy;
    			}

    }
    else
    {
    for (int kk=1; kk < nzn*ZLEN-1;kk++)
        for (int jj=1; jj < nyn*YLEN-1;jj++)
            for (int ii=1; ii < nxn*XLEN-1;ii++){

                vdotgrad[ii][jj][kk] = JX[ii][jj][kk] / (RHO[ii][jj][kk]+1e-10) *
                		(PHI[ii+1][jj][kk] - PHI[ii-1][jj][kk])/2.0/Dx;

                vdotgrad[ii][jj][kk] += JY[ii][jj][kk] / (RHO[ii][jj][kk]+1e-10) *
                                    		(PHI[ii][jj+1][kk] - PHI[ii][jj-1][kk])/2.0/Dy;


                vdotgrad[ii][jj][kk] += JZ[ii][jj][kk] / (RHO[ii][jj][kk]+1e-10) *
                                    		(PHI[ii][jj][kk+1] - PHI[ii][jj][kk-1])/2.0/Dz;
            }
    }
}

void vdiv(double*** vdiv, double*** PXX, double*** PYY, double*** PZZ,
                         double*** PXY, double*** PXZ, double*** PYZ,
                         double*** RHO, 
                         double*** JX, double*** JY, double*** JZ){
    
    for (int kk=0; kk < nzn*ZLEN;kk++)
        for (int jj=0; jj < nyn*YLEN;jj++)
            for (int ii=0; ii < nxn*XLEN;ii++){
			    
                vdiv[ii][jj][kk] = 0.0;
			}
    if(nzn*ZLEN==1) {
        int kk=0;
            for (int jj=1; jj < nyn*YLEN-1;jj++)
                for (int ii=1; ii < nxn*XLEN-1;ii++){

                    vdiv[ii][jj][kk] = ((PXX[ii+1][jj][kk] - PXX[ii-1][jj][kk])/2.0/Dx +
                                        (PXY[ii][jj+1][kk] - PXY[ii][jj-1][kk])/2.0/Dy ) *
                                        JX[ii][jj][kk] / RHO[ii][jj][kk];


                    vdiv[ii][jj][kk] += ((PXY[ii+1][jj][kk] - PXY[ii-1][jj][kk])/2.0/Dx +
                                        (PYY[ii][jj+1][kk] - PYY[ii][jj-1][kk])/2.0/Dy ) *
                                        JY[ii][jj][kk] / RHO[ii][jj][kk];


                    vdiv[ii][jj][kk] += ((PXZ[ii+1][jj][kk] - PXZ[ii-1][jj][kk])/2.0/Dx +
                                        (PYZ[ii][jj+1][kk] - PYZ[ii][jj-1][kk])/2.0/Dy ) *
                                        JZ[ii][jj][kk] / RHO[ii][jj][kk];
    			}

    }
    else
    {
    for (int kk=1; kk < nzn*ZLEN-1;kk++)
        for (int jj=1; jj < nyn*YLEN-1;jj++)
            for (int ii=1; ii < nxn*XLEN-1;ii++){
			    
                vdiv[ii][jj][kk] = ((PXX[ii+1][jj][kk] - PXX[ii-1][jj][kk])/2.0/Dx +
                                    (PXY[ii][jj+1][kk] - PXY[ii][jj-1][kk])/2.0/Dy +
                                    (PXZ[ii][jj][kk+1] - PXZ[ii][jj][kk-1])/2.0/Dz) *
                                    JX[ii][jj][kk] / RHO[ii][jj][kk];
                                    
                                    
                vdiv[ii][jj][kk] += ((PXY[ii+1][jj][kk] - PXY[ii-1][jj][kk])/2.0/Dx +
                                    (PYY[ii][jj+1][kk] - PYY[ii][jj-1][kk])/2.0/Dy +
                                    (PYZ[ii][jj][kk+1] - PYZ[ii][jj][kk-1])/2.0/Dz) *
                                    JY[ii][jj][kk] / RHO[ii][jj][kk];
                                    
                                    
                vdiv[ii][jj][kk] += ((PXZ[ii+1][jj][kk] - PXZ[ii-1][jj][kk])/2.0/Dx +
                                    (PYZ[ii][jj+1][kk] - PYZ[ii][jj-1][kk])/2.0/Dy +
                                    (PZZ[ii][jj][kk+1] - PZZ[ii][jj][kk-1])/2.0/Dz) *
                                    JZ[ii][jj][kk] / RHO[ii][jj][kk];
			}
    }
}  

void pdv_par_per(double*** pdvpar, double*** pdvper1, double*** pdvper2,
						 double*** PXX, double*** PYY, double*** PZZ,
                         double*** PXY, double*** PXZ, double*** PYZ,
						 double*** pPAR, double*** pPER1, double*** pPER2,
						 double*** BX, double*** BY, double*** BZ,
                         double*** RHO,
                         double*** JX, double*** JY, double*** JZ,
						 double*** V,
						 double*** VX, double*** VY, double*** VZ){


    for (int kk=0; kk < nzn*ZLEN;kk++)
        for (int jj=0; jj < nyn*YLEN;jj++)
            for (int ii=0; ii < nxn*XLEN;ii++){

                pdvpar[ii][jj][kk] = 0.0;
                pdvper1[ii][jj][kk] = 0.0;
                pdvper2[ii][jj][kk] = 0.0;
			}
    if(nzn*ZLEN==1) {
        int kk=0;
            for (int jj=1; jj < nyn*YLEN-1;jj++)
                for (int ii=1; ii < nxn*XLEN-1;ii++){


    			}

    }
    else
    {
    for (int kk=1; kk < nzn*ZLEN-1;kk++)
        for (int jj=1; jj < nyn*YLEN-1;jj++)
            for (int ii=1; ii < nxn*XLEN-1;ii++){

            	V[ii][jj][kk] = (JX[ii][jj][kk] / RHO[ii][jj][kk] *BX[ii][jj][kk] +
            	                		JY[ii][jj][kk] / RHO[ii][jj][kk] *BY[ii][jj][kk] +
            							JZ[ii][jj][kk] / RHO[ii][jj][kk] *BZ[ii][jj][kk]);
            	V[ii][jj][kk] = V[ii][jj][kk] / sqrt(BX[ii][jj][kk] * BX[ii][jj][kk] +
				          BY[ii][jj][kk] * BY[ii][jj][kk] +
						  BZ[ii][jj][kk] * BZ[ii][jj][kk]);
            }
            	p_dot_gradient(VX, VY, VZ, PXX, PXY, PXZ, PYY, PYZ, PZZ, V );

                for (int kk=1; kk < nzn*ZLEN-1;kk++)
                    for (int jj=1; jj < nyn*YLEN-1;jj++)
                        for (int ii=1; ii < nxn*XLEN-1;ii++){

                        	pdvpar[ii][jj][kk] = (VX[ii][jj][kk] *BX[ii][jj][kk] +
        	                		VY[ii][jj][kk] *BY[ii][jj][kk] +
        							VZ[ii][jj][kk] *BZ[ii][jj][kk]) / sqrt(BX[ii][jj][kk] * BX[ii][jj][kk] +
        							          BY[ii][jj][kk] * BY[ii][jj][kk] +
        									  BZ[ii][jj][kk] * BZ[ii][jj][kk]);
                        }

            	//Vperp1
                for (int kk=1; kk < nzn*ZLEN-1;kk++)
                    for (int jj=1; jj < nyn*YLEN-1;jj++)
                        for (int ii=1; ii < nxn*XLEN-1;ii++){
                        	double b2D = 1e-10 + BX[ii][jj][kk]*BX[ii][jj][kk] + BY[ii][jj][kk]*BY[ii][jj][kk];
                        	V[ii][jj][kk] = (JX[ii][jj][kk] / RHO[ii][jj][kk] *BY[ii][jj][kk] -
            	           	     JY[ii][jj][kk] / RHO[ii][jj][kk] *BX[ii][jj][kk])/ sqrt(b2D);
                        }

            	p_dot_gradient(VX, VY, VZ, PXX, PXY, PXZ, PYY, PYZ, PZZ, V );

                for (int kk=1; kk < nzn*ZLEN-1;kk++)
                    for (int jj=1; jj < nyn*YLEN-1;jj++)
                        for (int ii=1; ii < nxn*XLEN-1;ii++){

                        	double b2D = 1e-10 + BX[ii][jj][kk]*BX[ii][jj][kk] + BY[ii][jj][kk]*BY[ii][jj][kk];
                        	pdvper1[ii][jj][kk] = (VX[ii][jj][kk] *BY[ii][jj][kk] -
                        	            	           	     VY[ii][jj][kk] *BX[ii][jj][kk])/ sqrt(b2D);
                        }

            	//Vperp2
      for (int kk=1; kk < nzn*ZLEN-1;kk++)
        for (int jj=1; jj < nyn*YLEN-1;jj++)
            for (int ii=1; ii < nxn*XLEN-1;ii++){


            	double b2D = 1e-10 + BX[ii][jj][kk]*BX[ii][jj][kk] + BY[ii][jj][kk]*BY[ii][jj][kk];
            	double b = b2D + BZ[ii][jj][kk]*BZ[ii][jj][kk];
            	double perp2x = BZ[ii][jj][kk]*BX[ii][jj][kk] /sqrt(b*b2D);
            	double perp2y = BZ[ii][jj][kk]*BY[ii][jj][kk] /sqrt(b*b2D);
            	double perp2z = -sqrt(b2D/b);

            	V[ii][jj][kk] = (JX[ii][jj][kk] / RHO[ii][jj][kk] * perp2x +
            	            	        JY[ii][jj][kk] / RHO[ii][jj][kk] * perp2y +
										JZ[ii][jj][kk] / RHO[ii][jj][kk] * perp2z);
            }

      	  	  p_dot_gradient(VX, VY, VZ, PXX, PXY, PXZ, PYY, PYZ, PZZ, V );

                for (int kk=1; kk < nzn*ZLEN-1;kk++)
                    for (int jj=1; jj < nyn*YLEN-1;jj++)
                        for (int ii=1; ii < nxn*XLEN-1;ii++){
                        	double b2D = 1e-10 + BX[ii][jj][kk]*BX[ii][jj][kk] + BY[ii][jj][kk]*BY[ii][jj][kk];
                        	double b = b2D + BZ[ii][jj][kk]*BZ[ii][jj][kk];
                        	double perp2x = BZ[ii][jj][kk]*BX[ii][jj][kk] /sqrt(b*b2D);
                        	double perp2y = BZ[ii][jj][kk]*BY[ii][jj][kk] /sqrt(b*b2D);
                        	double perp2z = -sqrt(b2D/b);

                        	pdvper2[ii][jj][kk] = (VX[ii][jj][kk]  * perp2x +
        	            	        VY[ii][jj][kk]  * perp2y +
									VZ[ii][jj][kk]  * perp2z);
                        }

    }

		divergenceN(V,JX, JY, JZ, RHO);

		for (int kk=1; kk < nzn*ZLEN-1;kk++)
		                    for (int jj=1; jj < nyn*YLEN-1;jj++)
		                        for (int ii=1; ii < nxn*XLEN-1;ii++){

		                        	pdvpar[ii][jj][kk] = 2.0*pdvpar[ii][jj][kk] + V[ii][jj][kk] * pPAR[ii][jj][kk];
		                        	pdvper1[ii][jj][kk] = 2.0*pdvper1[ii][jj][kk] + V[ii][jj][kk] * pPER1[ii][jj][kk];
		                        	pdvper2[ii][jj][kk] = 2.0*pdvper2[ii][jj][kk] + V[ii][jj][kk] * pPER2[ii][jj][kk];
		                        }


}
void pdv(double*** pdvX, double*** pdvY,double*** pdvZ,
						 double*** PXX, double*** PYY, double*** PZZ,
                         double*** PXY, double*** PXZ, double*** PYZ,
						 double*** BX, double*** BY, double*** BZ,
                         double*** RHO,
                         double*** JX, double*** JY, double*** JZ,
						 double*** V,
						 double*** VX, double*** VY, double*** VZ){


    for (int kk=0; kk < nzn*ZLEN;kk++)
        for (int jj=0; jj < nyn*YLEN;jj++)
            for (int ii=0; ii < nxn*XLEN;ii++){

                pdvX[ii][jj][kk] = 0.0;
                pdvY[ii][jj][kk] = 0.0;
                pdvZ[ii][jj][kk] = 0.0;
			}
    if(nzn*ZLEN==1) {
        int kk=0;
            for (int jj=1; jj < nyn*YLEN-1;jj++)
                for (int ii=1; ii < nxn*XLEN-1;ii++){


    			}

    }
    else
    {
    for (int kk=1; kk < nzn*ZLEN-1;kk++)
        for (int jj=1; jj < nyn*YLEN-1;jj++)
            for (int ii=1; ii < nxn*XLEN-1;ii++){

            	V[ii][jj][kk] = JX[ii][jj][kk] / RHO[ii][jj][kk] ;

            }
            	p_dot_gradient(pdvX, VY, VZ, PXX, PXY, PXZ, PYY, PYZ, PZZ, V);



    for (int kk=1; kk < nzn*ZLEN-1;kk++)
         for (int jj=1; jj < nyn*YLEN-1;jj++)
             for (int ii=1; ii < nxn*XLEN-1;ii++){

             	V[ii][jj][kk] = JY[ii][jj][kk] / RHO[ii][jj][kk] ;

             }
             	p_dot_gradient(VX, pdvY, VZ, PXX, PXY, PXZ, PYY, PYZ, PZZ, V);



for (int kk=1; kk < nzn*ZLEN-1;kk++)
     for (int jj=1; jj < nyn*YLEN-1;jj++)
         for (int ii=1; ii < nxn*XLEN-1;ii++){

         	V[ii][jj][kk] = JZ[ii][jj][kk] / RHO[ii][jj][kk] ;

         }
         	p_dot_gradient(VX, VY, pdvZ, PXX, PXY, PXZ, PYY, PYZ, PZZ, V);


 }
    /*
		divergenceN(V,JX, JY, JZ, RHO);

		for (int kk=1; kk < nzn*ZLEN-1;kk++)
		                    for (int jj=1; jj < nyn*YLEN-1;jj++)
		                        for (int ii=1; ii < nxn*XLEN-1;ii++){

		                        	pdvX[ii][jj][kk] = 2.0*pdvX[ii][jj][kk] + V[ii][jj][kk] * PXX[ii][jj][kk];
		                        	pdvY[ii][jj][kk] = 2.0*pdvY[ii][jj][kk] + V[ii][jj][kk] * PYY[ii][jj][kk];
		                        	pdvZ[ii][jj][kk] = 2.0*pdvZ[ii][jj][kk] + V[ii][jj][kk] * PZZ[ii][jj][kk];
		                        }
*/

}

double vecpot(double*** AZ, double*** BX, double*** BY, double*** BZ){
    for (int kk=0; kk < nzn*ZLEN;kk++){
    	AZ[0][0][kk] = 0.0;
        for (int jj=0; jj < nyn*YLEN;jj++)
        {
            for (int ii=1; ii < nxn*XLEN;ii++){

                AZ[ii][jj][kk] = AZ[ii-1][jj][kk] - (BY[ii][jj][kk] + BY[ii][jj][kk])*Dx/2.0;
			}
       if(jj<nyn*YLEN-1)  AZ[0][jj+1][kk] = AZ[0][jj][kk] + (BX[0][jj][kk] + BX[0][jj+1][kk])*Dy/2.0;
    }
}
    double maxmax = -1.0e10;
    double minmax = 1.0e10;
    for (int kk=0; kk < nzn*ZLEN;kk++)
            for (int ii=1; ii < nxn*XLEN;ii++){
            	double maxAz =0.0;
            	for (int jj=0; jj < nyn*YLEN;jj++){
            		if (AZ[ii][jj][kk] > maxAz) maxAz = AZ[ii][jj][kk];
			}
         if(maxAz < minmax) minmax = maxAz;
         if(maxAz > maxmax) maxmax = maxAz;
     }
    double recflux = maxmax-minmax;
    cout << "RecFlux = " << recflux << endl;
    //cout << "RecFlux = " << maxmax << "    " << minmax << endl;
    return(recflux);
}

void mult(double a, double*** E){
    for (int kk=1; kk < nzn*ZLEN-1;kk++)
        for (int jj=1; jj < nyn*YLEN-1;jj++)
            for (int ii=1; ii < nxn*XLEN-1;ii++){
            E[ii][jj][kk] *=a;
            }
}

void divmult(double*** A, double*** B, double*** C){
    for (int kk=1; kk < nzn*ZLEN-1;kk++)
        for (int jj=1; jj < nyn*YLEN-1;jj++)
            for (int ii=1; ii < nxn*XLEN-1;ii++){
            A[ii][jj][kk] *= B[ii][jj][kk] / C[ii][jj][kk];
            }
}

void smooth(int Nvolte, double*** B, int nx, int ny, int nz){

	double*** C = newArr3(double,nx, ny, nz);

	for (int ivolte=0; ivolte<Nvolte; ivolte++){

    if(nx>2){
    for (int kk=0; kk < nz;kk++)
        for (int jj=0; jj < ny;jj++){
        	C[0][jj][kk] = .25 * (2*B[0][jj][kk] + B[1][jj][kk] + B[2][jj][kk]);
         	C[nx-1][jj][kk]= .25*(2*B[nx-1][jj][kk]+B[nx-2][jj][kk]+B[nx-3][jj][kk]);
            for (int ii=1; ii < nx-1;ii++){
            C[ii][jj][kk] = .25 *(B[ii+1][jj][kk] + B[ii-1][jj][kk] +2.0 * B[ii][jj][kk]) ;
            }}
    }
    eq(B, C, nx, ny, nz);
    if(ny>2){
    for (int kk=0; kk < nz;kk++)
    	for (int ii=0; ii < nx;ii++){
        	C[ii][0][kk] = .25*(2*B[ii][0][kk]+B[ii][1][kk]+B[ii][2][kk]);
         	C[ii][ny-1][kk] = .25*(2*B[ii][ny-1][kk]+B[ii][ny-2][kk]+B[ii][ny-3][kk]);
         	for (int jj=1; jj < ny-1;jj++){
            C[ii][jj][kk] = .25 *(B[ii][jj+1][kk] + B[ii][jj-1][kk] +2.0 * B[ii][jj][kk]) ;
            }}
    }
    eq(B, C, nx, ny, nz);
    if(nz>2){
        for (int jj=0; jj < ny;jj++)
            for (int ii=0; ii < nx;ii++) {
            	C[ii][jj][0] = .25*(2*B[ii][jj][0]+B[ii][jj][1]+B[ii][jj][2]);
             	C[ii][jj][nz-1] = .25*(2*B[ii][jj][nz-1]+B[ii][jj][nz-2]+B[ii][jj][nz-3]);
            	for (int kk=1; kk < nz-1;kk++){
            C[ii][jj][kk] = .25 *(B[ii][jj][kk+1] + B[ii][jj][kk-1] +2.0 * B[ii][jj][kk]) ;
            }}
    }
    eq(B, C, nx, ny, nz);

	}
    delArr3(C,nx,ny);
}

void smooth(int Nvolte, double*** A, double*** B, int nx, int ny, int nz){

	double*** C = newArr3(double,nx, ny, nz);

	for (int ivolte=0; ivolte<Nvolte; ivolte++){


    eq(B, A, nx, ny, nz);
    if(nx>2){
    for (int kk=0; kk < nz;kk++)
        for (int jj=0; jj < ny;jj++){
        	C[0][jj][kk] = .25 * (2*B[0][jj][kk] + B[1][jj][kk] + B[2][jj][kk]);
         	C[nx-1][jj][kk]= .25*(2*B[nx-1][jj][kk]+B[nx-2][jj][kk]+B[nx-3][jj][kk]);
            for (int ii=1; ii < nx-1;ii++){
            C[ii][jj][kk] = .25 *(B[ii+1][jj][kk] + B[ii-1][jj][kk] +2.0 * B[ii][jj][kk]) ;
            }}
    }
    eq(B, C, nx, ny, nz);
    if(ny>2){
    for (int kk=0; kk < nz;kk++)
    	for (int ii=0; ii < nx;ii++){
        	C[ii][0][kk] = .25*(2*B[ii][0][kk]+B[ii][1][kk]+B[ii][2][kk]);
         	C[ii][ny-1][kk] = .25*(2*B[ii][ny-1][kk]+B[ii][ny-2][kk]+B[ii][ny-3][kk]);
         	for (int jj=1; jj < ny-1;jj++){
            C[ii][jj][kk] = .25 *(B[ii][jj+1][kk] + B[ii][jj-1][kk] +2.0 * B[ii][jj][kk]) ;
            }}
    }
    eq(B, C, nx, ny, nz);
    if(nz>2){
        for (int jj=0; jj < ny;jj++)
            for (int ii=0; ii < nx;ii++) {
            	C[ii][jj][0] = .25*(2*B[ii][jj][0]+B[ii][jj][1]+B[ii][jj][2]);
             	C[ii][jj][nz-1] = .25*(2*B[ii][jj][nz-1]+B[ii][jj][nz-2]+B[ii][jj][nz-3]);
            	for (int kk=1; kk < nz-1;kk++){
            C[ii][jj][kk] = .25 *(B[ii][jj][kk+1] + B[ii][jj][kk-1] +2.0 * B[ii][jj][kk]) ;
            }}
    }
    eq(B, C, nx, ny, nz);

	}
    delArr3(C,nx,ny);
}


void agyro(double*** agyro_scudder, double*** agyro_aunai, double*** nongyro_swisdak, double*** align,
		double*** BX, double*** BY, double*** BZ,
		double*** pXX, double*** pXY,
		double*** pXZ, double*** pYY, double*** pYZ, double*** pZZ){

    int n = 3, lda = 3, ldvl = 3, ldvr = 3, info, lwork;
    double wkopt;
    double* work;
    double wr[n], wi[n], vl[ldvl*n], vr[ldvr*n];

	double** p = newArr2(double,n,n);
	double** p1 = newArr2(double,n,n);
	double** p2 = newArr2(double,n,n);
	double* b = new double[n];
	double* dot = new double[n];
	double* a = new double[n*n];


	 for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				p[0][0] = pXX[ii][jj][kk]+1e-10;
				p[0][1] = pXY[ii][jj][kk];
				p[0][2] = pXZ[ii][jj][kk];
				p[1][1] = pYY[ii][jj][kk]+1e-10;
				p[1][2] = pYZ[ii][jj][kk];
				p[2][2] = pZZ[ii][jj][kk]+1e-10;
				p[1][0] = p[0][1];
				p[2][0] = p[0][2];
				p[2][1] = p[1][2];
				b[0] = BX[ii][jj][kk];
				b[1] = BY[ii][jj][kk];
				b[2] = BZ[ii][jj][kk];


	 double bnorm=norm(n, b);
			 b[0]=b[0]/(bnorm+1e-10);
			 b[1]=b[1]/(bnorm+1e-10);
			 b[2]=b[2]/(bnorm+1e-10);

	// Scudder Agyro
    cross(n,b,p,p1);
    cross_trans(n,b,p1,p2);
    mat2vet(n,p2,a);

    lwork = -1;
    dgeev_( (char*)"N", (char*)"N", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
     &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = new double[lwork];
    /* Solve eigenproblem */
    dgeev_( (char*)"N", (char*)"N", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
     work, &lwork, &info );

    sort(n,wr);
    //cout << wr[0]<< "  "<< wr[1]<< "  "<< wr[2]<< endl;
    agyro_scudder[ii][jj][kk] = 2.0* fabs(wr[2]-wr[1])/(wr[2]+wr[1]);

    //Non-alignment
    mat2vet(n,p,a);
	delete[] work;

    lwork = -1;
    dgeev_( (char*)"N", (char*)"V", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
     &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = new double[lwork];
    /* Solve eigenproblem */
    dgeev_( (char*)"N", (char*)"V", &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr,
     work, &lwork, &info );

    vet2mat(n, p2, vr, b, dot);
    sort(n, dot); //put the minimum first, maximum last

    align[ii][jj][kk] = -log(1.0-.9999*fabs(dot[2]));


    // Aunai Agyro
    double p_par = ppar(n, p, b);
    double Tr = p[0][0] + p[1][1] + p[2][2];
    double p_per = (Tr-p_par)/2.0;

    for (int i=0; i<n; i++)
    	for (int j=0; j<n; j++){
    		p2[i][j] = (p_par-p_per) * b[i]*b[j];
    		if(i==j) p2[i][j] += p_per;
    		p1[i][j] = p[i][j] - p2[i][j];
    	}

    agyro_aunai[ii][jj][kk] = norm(n, p1)/Tr;

    // Swisdak Nongyro

    double I2=p[0][0]*p[1][1]+p[0][0]*p[2][2]+p[1][1]*p[2][2];
    I2=I2-(p[0][1]*p[0][1] + p[0][2]*p[0][2] + p[1][2] * p[1][2]);
    double Q=1.0-4.0*I2/((Tr-p_par)*(Tr+3.0*p_par));
    nongyro_swisdak[ii][jj][kk] = sqrt(fabs(Q));
}
	delArr2(p,n);
	delArr2(p1,n);
	delArr2(p2,n);
	delete[] a;
	delete[] b;
	delete[] dot;
	delete[] work;
}
void mat2vet(int n, double** mat, double* vet){
	int counter =0;

	for (int j=0; j<n; j++)
		for (int i=0; i<n; i++){
		vet[counter] = mat[i][j];
		counter++;
	}
}

void vet2mat(int n, double** mat, double* vet, double* b, double* dot){
	int counter =0;




	for (int j=0; j<n; j++)
		for (int i=0; i<n; i++){
		mat[i][j] =  vet[counter];
		counter++;
	}
// compute the dot product of the eigenvectors with b
	for (int j=0; j<n; j++){
		dot[j] = 0.0;
	    double tmp = 0.0;
		for (int i=0; i<n; i++){
				dot[j] += mat[i][j] *b[i];
				tmp += mat[i][j] * mat[i][j];
			}
        //now it normalizes the dot product of the eigenvectors
        //with b using the norm of teh eigenvectors
        //(in case LAPACK does not give unit egenvectors)
		dot[j] = dot[j]/sqrt(tmp);
	}

}
void sort(int n, double* vet){
	for (int i=0; i<n-1;i++)
	   for(int j=i+1;j<n;j++){
		if(fabs(vet[j])<fabs(vet[i])){
			double temp=vet[i];
			vet[i] = vet[j];
			vet[j] = temp;
		}

	}
}

double ppar(int n, double** p, double* b){
	double ppar = 0.0;
	for (int i=0; i<n; i++)
		for(int j=0; j<n; j++){
			ppar += b[i]*p[i][j]*b[j];
		}
	return(ppar);
}


double norm(int n,  double** p){
	double norm = 0.0;
	for (int i=0; i<n;i++)
		for(int j=0 ;j<n;j++){
			norm += p[i][j]*p[i][j];
		}
	norm= sqrt(norm);
	return( norm ) ;
}


double norm(int n,  double* b){
	double norm = 0.0;
	for (int i=0; i<n;i++){
			norm += b[i]*b[i];
		}
	norm= sqrt(norm);
	return( norm ) ;
}

void circshift(double ***vect, double ***tmp, int xshift, int yshift, int zshift)
{
	int xdim = nxn*XLEN;
	int ydim = nyn*YLEN;
	int zdim = nzn*ZLEN;
   for (int i=0; i < xdim ;i++){
   int ii = (i + xshift) % xdim;
   if (ii<0) ii = xdim + ii;
   for (int j=0; j < ydim;j++){
     int jj = (j + yshift) % ydim;
     if (jj<0) jj = ydim + jj;
     for (int k=0; k < zdim;k++){
    	 int kk = (k + zshift) % zdim;
    	 if (kk<0) kk = zdim +kk;
             tmp[ii][jj][kk] = vect[i][j][k];
   }}}
   for (int i=0; i < xdim ;i++)
   for (int j=0; j < ydim;j++)
   for (int k=0; k < zdim;k++){
    vect[i][j][k] = tmp[i][j][k];
   }
}


