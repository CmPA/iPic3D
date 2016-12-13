/***************************************************************************
read_ASCII.cpp  -  Convert program to open Parsek Output
 -------------------
 begin                : Jun 2008
 copyright            : (C) 2004 by Stefano Markidis, Giovanni Lapenta
 ************************************************************************** */

#include "hdf5.h"
#include "Alloc.h"
#include "math.h"
#include "Manipulator.h"
#include "VtkStuff.h"
#include "Particles.h"
#include "ASCIIStuff.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


int main (int argc, char **argv) {
    
	// tail
    nxc = 151;
    nyc = 61;
    nzc = nyc;
	//pause
	nxc = 201;
	nyc = 401;
	nzc = 9;
    
    XLEN=1;
    YLEN=1;
    ZLEN=1;
    nxn=nxc;
    nyn=nyc;
    nzn=nzc;
    // tail
    Lx=.2*nxc;
    Ly=.2*nyc;
    Lz=.2*nzc;
    Dx=0.02;
    Dy=0.02;
    Dz=0.02;
    //pause
//    Lx=.1*nxc;
//    Ly=.1*nyc;
 //   Lz=.1*nzc;
 //   Dx=0.1;
//    Dy=0.1;
 //   Dz=0.1;
    
    
	double*** BX = newArr3(double,nxc,nyc,nzc);
	double*** BY = newArr3(double,nxc,nyc,nzc);
	double*** BZ = newArr3(double,nxc,nyc,nzc);
	double*** VX = newArr3(double,nxc,nyc,nzc);
	double*** VY = newArr3(double,nxc,nyc,nzc);
	double*** VZ = newArr3(double,nxc,nyc,nzc);
	double*** JX = newArr3(double,nxc,nyc,nzc);
	double*** JY = newArr3(double,nxc,nyc,nzc);
	double*** JZ = newArr3(double,nxc,nyc,nzc);
	double*** RHO = newArr3(double,nxc,nyc,nzc);
	double*** P = newArr3(double,nxc,nyc,nzc);
    
// Tail
//	readASCIIvect("feb1508cgsec.box.0357.dat", BX, BY, BZ,
//	readASCIIvect("feb1508cgsec.box.0348.dat", BX, BY, BZ,
//	readASCIIvect("feb1508cgsec.box.035100.dat", BX, BY, BZ,
//readASCIIvect("feb1508gsm.zmin-0870zmax+0330.034800UT.dat"
	// New Tail in code coordinates
//		readASCIIvect("feb1508iPIC.035100UT.dat"
// MAgnetopause
	readASCIIvect("psy400.iPic_mp_box_ttcor.0015.dat"
				, BX, BY, BZ,
                  VX, VY, VZ,
                  JX, JY, JZ,
                  RHO, P);
    
    
    writeVTKvect(0,"Bmhd", "", BX, BY, BZ);
    writeVTKvect(0,"Vmhd", "", VX, VY, VZ);
    writeVTKvect(0,"Jmhd", "", JX, JY, JZ);
    writeVTKscalar(0,"RHOmhd", "", RHO);
    writeVTKscalar(0,"Pmhd", "", P);
    
    
	delArr3(BX,nxc,nyc);
	delArr3(BY,nxc,nyc);
	delArr3(BZ,nxc,nyc);
	delArr3(VX,nxc,nyc);
	delArr3(VY,nxc,nyc);
	delArr3(VZ,nxc,nyc);
	delArr3(JX,nxc,nyc);
	delArr3(JY,nxc,nyc);
	delArr3(JZ,nxc,nyc);	
	delArr3(RHO,nxc,nyc);
	delArr3(P,nxc,nyc);
    
    return(0);
}
