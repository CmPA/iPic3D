/***************************************************************************
 convHDF5.cpp  -  Convert program to open Parsek Output
 -------------------
 begin                : Jun 2008
 copyright            : (C) 2004 by Stefano Markidis, Giovanni Lapenta
 ************************************************************************** */

#include "hdf5.h"
#include "mpi.h"
#include "Alloc.h"
#include "math.h"
#include "Manipulator.h"
#include "VtkStuff.h"
#include "GdaStuff.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>


using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;


// binary files 
FILE *fp;

int main (int argc, char **argv) {
	// cycle we want to open
		sscanf(argv[1],"%d",&InitT);
	sscanf(argv[2],"%d",&MaxLevel);
	sscanf(argv[3],"%d",&DeltaT);


			NdimCode = 3;
			xshift = 0;
			yshift = 0;
			zshift = 0;
	if(argc>3) {
			sscanf(argv[4], "%d", &NdimCode);
			sscanf(argv[5],"%d",&xshift);
			sscanf(argv[6],"%d",&yshift);
			sscanf(argv[7],"%d",&zshift);
			grid_str = argv[8];
		}

	printf( "Shift in x is set to %d ", xshift);
	printf( "Shift in y is set to %d ", yshift);
	printf( "Shift in z is set to %d ", zshift);
	printf( "grid_Str %s\n ", grid_str);

     int rank, size;
    
    MPI_Init (&argc, &argv);	/* starts MPI */
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */
    printf( "Hello world from process %d of %d\n", rank, size );
    
	nlevels = MaxLevel/DeltaT;
	initlevel = InitT/DeltaT;

    int out;
	out = readsettings();
	if(out<0)
		return -1;

	nxn = nxc/XLEN;
	nyn = nyc/YLEN;
	nzn = nzc/ZLEN;




	temp_storageX = new double[(nxn+1)*(nyn+1)*(nzn+1)];
	temp_storageY = new double[(nxn+1)*(nyn+1)*(nzn+1)];
    temp_storageZ = new double[(nxn+1)*(nyn+1)*(nzn+1)];

	double*** EX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** EY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** EZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** BX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** BY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** BZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** VX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** VY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** VZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** NE = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** NI = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);

	double*** TXX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TXY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TXZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TYY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TYZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TZZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TPAR  = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TPER1 = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** TPER2 = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** EPS = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
/*
		double*** pdvX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** pdvY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** pdvZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		// The following are work space variables
		double*** WX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** WY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
		double*** WZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
*/
	//Electric field
	for (int it = initlevel+rank; it < nlevels+1; it=it+size){
    printf( "Process %d doing level %d\n", rank, it );
        
    readvect(it, "/fields/","E", EX, EY, EZ);
    writeGDAvect(it, "E", "", EX, EY, EZ);

	//Magnetic field
    readvect(it,"/fields/","B", BX, BY, BZ);
    writeGDAvect(it,"B", "", BX, BY, BZ);
/*
    //Divergence Poynting
    div_cross(EX, EY, EZ, BX, BY, BZ, TXX, TYY, TZZ, TXY);
    // Note: SX,SY,SZ have the Poynting flux, but not yet divided by 4pi,
    // same for divergence
    double invFourPi = 1.0/4.0/M_PI;
    mult(invFourPi, TXY);
    mult(invFourPi, TXX);
    mult(invFourPi, TYY);
    mult(invFourPi, TZZ);
    writeGDAscalar(it,"divS","", TXY);
    writeGDAvect(it,"S","", TXX, TYY, TZZ);
*/
/*
    // Lagrangian Energy Balance Term
    // To save space now we overwrite divS and V that will now be,
    // respectively a work variable and the Lagrangian Energy Balance Term
    LEBT(BX, BY, BZ, TXX, TYY, TZZ, VX, VY, VZ, TXY);
    writeGDAscalar(it,"LagrEBTerm","", TXY);

    // Reconnection measure by Biskmap
    // To save space now we overwrite S and V that will now be,
    // respectively a work variable and the Biskamp measure
    biskamp(BX, BY, BZ, EX, EY, EZ, VX, VY, VZ, TXX, TYY, TZZ);
    writeGDAvect(it,"Biskamp","", TXX, TYY, TZZ);

	//Compute ExB
    cross(BX, BY, BZ, EX, EY, EZ, VXBX, VXBY, VXBZ)
*/
	//Rho by species
    readscalar(it,"/moments/species_0/","rho",  NE);
    readscalar(it,"/moments/species_1/","rho",  NI);
    if (ns >2) addreadscalar(it,"/moments/species_2/","rho",  NE);
    if (ns >2) addreadscalar(it,"/moments/species_3/","rho",  NI);
    writeGDAscalar_species(it,"rho", NE, NI);
   // writeVTKscalar("rho", "e", VX);


	//Currents species0
    readvect(it,"/moments/species_0/","J",  VX, VY, VZ);
    if (ns >2) addreadvect(it,"/moments/species_2/","J",  VX, VY, VZ);
    writeGDAvect(it,"J", "e", VX, VY, VZ);
    
    //Pressure tensor species 0
     readscalar(it,"/moments/species_0/","pXX",  TXX);
     readscalar(it,"/moments/species_0/","pXY",  TXY);
     readscalar(it,"/moments/species_0/","pXZ",  TXZ);
     readscalar(it,"/moments/species_0/","pYY",  TYY);
     readscalar(it,"/moments/species_0/","pYZ",  TYZ);
     readscalar(it,"/moments/species_0/","pZZ",  TZZ);

     if (ns >2) addreadscalar(it,"/moments/species_2/","pXX",  TXX);
     if (ns >2) addreadscalar(it,"/moments/species_2/","pXY",  TXY);
     if (ns >2) addreadscalar(it,"/moments/species_2/","pXZ",  TXZ);
     if (ns >2) addreadscalar(it,"/moments/species_2/","pYY",  TYY);
     if (ns >2) addreadscalar(it,"/moments/species_2/","pYZ",  TYZ);
     if (ns >2) addreadscalar(it,"/moments/species_2/","pZZ",  TZZ);

     extract_pressure(qom[0], BX, BY, BZ, VX, VY, VZ, NE,
           TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2, EPS);
     writeGDAtensor(it, "P", "e", TXX, TXY, TXZ, TYY, TYZ, TZZ,
     	  TPAR, TPER1, TPER2, EPS);

    /*
    // To save space now we overwrite EPS that is now a service variable,
 	pdv(pdvX, pdvY, pdvZ, TXX, TXY, TXZ, TYY, TYZ, TZZ,
 			BX, BY, BZ, NE, VX, VY, VZ, EPS, WX, WY, WZ);
 	writeGDAvect(it,"PdV", "e", pdvX, pdvY, pdvZ);

  	// We reuse pvpar etc to save space
 	divergenceNP(pdvX, VX, VY, VZ, NE, TXX);
 	divergenceNP(pdvY, VX, VY, VZ, NE, TYY);
 	divergenceNP(pdvZ, VX, VY, VZ, NE, TZZ);
 	writeGDAvect(it,"divPu", "e", pdvX, pdvY, pdvZ);

 	cross(TXX, TXY, TXZ, BX, BY, BZ, pdvX, WY, WZ);
 	cross(TXY, TYY, TYZ, BX, BY, BZ, WX, pdvY, WZ);
 	cross(TXZ, TYZ, TZZ, BX, BY, BZ, WX, WY, pdvZ);
 	mult(qom[0],pdvX);
 	mult(qom[0],pdvY);
 	mult(qom[0],pdvY);
 	writeGDAvect(it,"PxB", "e", pdvX, pdvY, pdvZ);

 	 // div Enthalpy component flux (reuses pdv vectors)

 	 	divPv(pdvX, pdvY, pdvZ, TXX, TXY, TXZ, TYY, TYZ, TZZ,
 	 			VX, VY, VZ, NE,  WX, WY, WZ);
 	 	writeGDAvect(it,"Enth", "e", pdvX, pdvY, pdvZ);

 	 // Computes Work by pressure


 	 	udivP(pdvX, pdvY, pdvZ, VX, VY, VZ, NE,
 	 			TXX, TXY, TXZ, TYY, TYZ, TZZ);
 	 	writeGDAvect(it,"WorkP", "e", pdvX, pdvY, pdvZ);

     vnonfrozen(NE, BX, BY, BZ, EX, EY, EZ, VX, VY, VZ);
     writeGDAvect(it,"VNF", "e", VX, VY, VZ);
*/
    //Currents species1
    readvect(it,"/moments/species_1/","J",  VX, VY, VZ);
    if (ns >2) addreadvect(it,"/moments/species_3/","J",  VX, VY, VZ);
    writeGDAvect(it,"J", "i", VX, VY, VZ);

    //Pressure tensor species 1
     readscalar(it,"/moments/species_1/","pXX",  TXX);
     readscalar(it,"/moments/species_1/","pXY",  TXY);
     readscalar(it,"/moments/species_1/","pXZ",  TXZ);
     readscalar(it,"/moments/species_1/","pYY",  TYY);
     readscalar(it,"/moments/species_1/","pYZ",  TYZ);
     readscalar(it,"/moments/species_1/","pZZ",  TZZ);

     if (ns >2) addreadscalar(it,"/moments/species_3/","pXX",  TXX);
     if (ns >2) addreadscalar(it,"/moments/species_3/","pXY",  TXY);
     if (ns >2) addreadscalar(it,"/moments/species_3/","pXZ",  TXZ);
     if (ns >2) addreadscalar(it,"/moments/species_3/","pYY",  TYY);
     if (ns >2) addreadscalar(it,"/moments/species_3/","pYZ",  TYZ);
     if (ns >2) addreadscalar(it,"/moments/species_3/","pZZ",  TZZ);

     extract_pressure(qom[1], BX, BY, BZ, VX, VY, VZ, NI,
           TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2, EPS);
     writeGDAtensor(it, "P", "i", TXX, TXY, TXZ, TYY, TYZ, TZZ,
     	  TPAR, TPER1, TPER2, EPS);
/*
     // To save space now we overwrite EPS that is now a service variable,
  	pdv(pdvX, pdvY, pdvZ, TXX, TXY, TXZ, TYY, TYZ, TZZ,
  			BX, BY, BZ, NI, VX, VY, VZ, EPS, WX, WY, WZ);
  	writeGDAvect(it,"PdV", "i", pdvX, pdvY, pdvZ);

  	// We reuse pvpar etc to save space
 	divergenceNP(pdvX, VX, VY, VZ, NI, TXX);
 	divergenceNP(pdvY, VX, VY, VZ, NI, TYY);
 	divergenceNP(pdvZ, VX, VY, VZ, NI, TZZ);
 	writeGDAvect(it,"divPu", "i", pdvX, pdvY, pdvZ);

 	cross(TXX, TXY, TXZ, BX, BY, BZ, pdvX, WY, WZ);
 	cross(TXY, TYY, TYZ, BX, BY, BZ, WX, pdvY, WZ);
 	cross(TXZ, TYZ, TZZ, BX, BY, BZ, WX, WY, pdvZ);
 	mult(qom[1],pdvX);
 	mult(qom[1],pdvY);
 	mult(qom[1],pdvY);
 	writeGDAvect(it,"PxB", "i", pdvX, pdvY, pdvZ);


 // div Enthalpy component flux (reuses pdv vectors)

 	divPv(pdvX, pdvY, pdvZ, TXX, TXY, TXZ, TYY, TYZ, TZZ,
 			VX, VY, VZ, NI,  WX, WY, WZ);
 	writeGDAvect(it,"Enth", "i", pdvX, pdvY, pdvZ);

 // Computes Work by pressure


 	udivP(pdvX, pdvY, pdvZ, VX, VY, VZ, NI,
 			TXX, TXY, TXZ, TYY, TYZ, TZZ);
 	writeGDAvect(it,"WorkP", "i", pdvX, pdvY, pdvZ);

     vnonfrozen(NI, BX, BY, BZ, EX, EY, EZ, VX, VY, VZ);
     writeGDAvect(it,"VNF", "i", VX, VY, VZ);
*/
}
	delArr3(BX,nxn*XLEN,nyn*YLEN);
	delArr3(BY,nxn*XLEN,nyn*YLEN);
	delArr3(BZ,nxn*XLEN,nyn*YLEN);
	delArr3(EX,nxn*XLEN,nyn*YLEN);
	delArr3(EY,nxn*XLEN,nyn*YLEN);
	delArr3(EZ,nxn*XLEN,nyn*YLEN);
	delArr3(VX,nxn*XLEN,nyn*YLEN);
	delArr3(VY,nxn*XLEN,nyn*YLEN);
	delArr3(VZ,nxn*XLEN,nyn*YLEN);

	delArr3(TXX,nxn*XLEN,nyn*YLEN);
	delArr3(TXY,nxn*XLEN,nyn*YLEN);
	delArr3(TXZ,nxn*XLEN,nyn*YLEN);
	delArr3(TYY,nxn*XLEN,nyn*YLEN);
	delArr3(TYZ,nxn*XLEN,nyn*YLEN);
	delArr3(TZZ,nxn*XLEN,nyn*YLEN);
	delArr3(TPAR,nxn*XLEN,nyn*YLEN);
	delArr3(TPER1,nxn*XLEN,nyn*YLEN);
	delArr3(TPER2,nxn*XLEN,nyn*YLEN);
	delArr3(EPS,nxn*XLEN,nyn*YLEN);
/*
	delArr3(pdvX,nxn*XLEN,nyn*YLEN);
	delArr3(pdvY,nxn*XLEN,nyn*YLEN);
	delArr3(pdvZ,nxn*XLEN,nyn*YLEN);

	delArr3(WX,nxn*XLEN,nyn*YLEN);
	delArr3(WY,nxn*XLEN,nyn*YLEN);
	delArr3(WZ,nxn*XLEN,nyn*YLEN);
*/
	delete[] temp_storageX;
	delete[] temp_storageY;
	delete[] temp_storageZ;
	
    
    MPI_Finalize();
    return(0);
}


