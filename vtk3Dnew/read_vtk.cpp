#include "hdf5.h"
#include "Alloc.h"
#include "math.h"
#include "Manipulator.h"
#include "VtkStuff.h"
#include "Particles.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


int main (int argc, char **argv) {
    

	readVTKpreamble("B1.vtk");
	
	double*** BX = newArr3(double,nxc,nyc,nzc);
	double*** BY = newArr3(double,nxc,nyc,nzc);
	double*** BZ = newArr3(double,nxc,nyc,nzc);
	
	double*** EX = newArr3(double,nxc,nyc,nzc);
	double*** EY = newArr3(double,nxc,nyc,nzc);
	double*** EZ = newArr3(double,nxc,nyc,nzc);

		np=100;
	xp = new double[np];
	yp = new double[np];
	zp = new double[np];

	up = new double[np];
	vp = new double[np];
	wp = new double[np];

	
	readVTKvect("B25.vtk", BX, BY, BZ);
	readVTKvect("E25.vtk", BX, BY, BZ);
	
/*	for(int i = 0; i < nxc; i++) 
		for(int j = 0; j < nyc; j++)
				for(int k = 0; k < nzc; k++){
	cout <<"bx= " << BX[i][j][k] << " by=" << BY[i][j][k] << " bz=" << BZ[i][j][k] << endl; 
		}
*/	

	load(xp,vp);
	Dt = .1;
	mover(xp,vp);
	

	delArr3(BX,nxc,nyc);
	delArr3(BY,nxc,nyc);
	delArr3(BZ,nxc,nyc);
	delArr3(EX,nxc,nyc);
	delArr3(EY,nxc,nyc);
	delArr3(EZ,nxc,nyc);	
}
