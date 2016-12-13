
/*
 *  ASCIIStuff.h
 *  
 *
 *  Created by Giovanni Lapenta on 7/29/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "hdf5.h"
#include "Alloc.h"
#include "math.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;


void readASCIIvect(string filename, double ***BX, double ***BY,double ***BZ,
			double ***VX, double ***VY,double ***VZ, 
			double ***JX, double ***JY,double ***JZ,
			double ***RHO, double ***P);


// Various useful variables about the dataset




void readASCIIvect(string filename, double*** BX, double*** BY,double*** BZ, 
			double*** VX, double*** VY,double*** VZ, 
			double*** JX, double*** JY,double*** JZ,
			double*** RHO, double*** P) {
	string line;
	string s;
	ifstream myfile(filename.c_str());
	
	double x, y, z;
	double bx, by, bz;
	double vx, vy, vz;
	double jx, jy, jz;
	double density;
	double pressure;

    double xmin, xmax, dx, kx;
    
    // x-direction
	myfile >> s;
	myfile >> s;
	myfile >> s;
    myfile >> s;
    myfile >> xmin;
    myfile >> xmax;
    myfile >> dx;
    myfile >> kx;
    // y-direction
    myfile >> s;
    myfile >> s;
    myfile >> s;
    myfile >> s;
    myfile >> xmin;
    myfile >> xmax;
    myfile >> dx;
    myfile >> kx;
    // z-direction
    myfile >> s;
    myfile >> s;
    myfile >> s;
    myfile >> s;
    myfile >> xmin;
    myfile >> xmax;
    myfile >> dx;
    myfile >> kx;


// Data structure:	
// x(RE) y(RE) z(RE) bx(nT) by(nT) bz(nT) vx(km/s) vy(km/s) vz(km/s) 
// density(cm-3) pressure(pPa) jx(nA/m2) jy(nA/m2) jz(nA/m2)	
     	
	for(int i = 0; i < nxc; i++) 
		for(int j = 0; j < nyc; j++)
			for(int k = 0; k < nzc; k++){    
    // read first line
	myfile >> x; 
	myfile >> y; 
	myfile >> z; 
	
				myfile >> BX[i][j][k];
				myfile >> BY[i][j][k];
				myfile >> BZ[i][j][k];
	
				myfile >> VX[i][j][k];
				myfile >> VY[i][j][k];
				myfile >> VZ[i][j][k];

				myfile >> RHO[i][j][k];
				myfile >> P[i][j][k];
	
				myfile >> JX[i][j][k];
				myfile >> JY[i][j][k];
				myfile >> JZ[i][j][k];
				
		cout << "writing   " << x << "   BX" << BX[i][j][k] <<endl;
}
    myfile.close();	
}	
