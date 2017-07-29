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
#include "Basic.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

int main (int argc, char **argv) {
	// cycle we want to open
		sscanf(argv[1],"%d",&InitT);
	sscanf(argv[2],"%d",&MaxLevel);
	sscanf(argv[3],"%d",&DeltaT);

	if(argc>2) {
			sscanf(argv[4], "%d", &NdimCode);
		}
		else
		{
			NdimCode = 3;
		}

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

	double*** BX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** BY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** BZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** VX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** VY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** VZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
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
	double*** NE = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** NI = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** QX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** QY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** QZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** QXbulk = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** QYbulk = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** QZbulk = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** divQ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** divQbulk = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** divU = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** UdivP = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** Uth = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** Ubulk = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** PgradU = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);


	for (int it = initlevel+rank; it < nlevels+1; it=it+size){
	     printf( "Process %d doing level %d\n", rank, it );

    readvect(it, "/fields/","B", BX, BY, BZ);

		//Rho by species
		readscalar(it,"/moments/species_0/","rho",  NE);
		readscalar(it,"/moments/species_1/","rho",  NI);
		if (ns >2) addreadscalar(it,"/moments/species_2/","rho",  NE);
		if (ns >2) addreadscalar(it,"/moments/species_3/","rho",  NI);

		// ===========================================================
		// SPECIES 0 : electrons
		// ===========================================================
		
		//Currents species0
		readvect(it,"/moments/species_0/","J",  VX, VY, VZ);
		if (ns >2) addreadvect(it,"/moments/species_2/","J",  VX, VY, VZ);
		
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
		
		// Note: this call expects V to be the current for compatibility with other postprocessing tools
		extract_pressure(qom[0], BX, BY, BZ, VX, VY, VZ, NE,
						 TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2, EPS);

		//Compute velocity from Current
		div(VX, NE, nxn*XLEN, nyn*YLEN, nzn*ZLEN);
		div(VY, NE, nxn*XLEN, nyn*YLEN, nzn*ZLEN);
		div(VZ, NE, nxn*XLEN, nyn*YLEN, nzn*ZLEN);

		// Computes Energies
		compute_energy(qom[0], Ubulk, Uth, VX, VY, VZ, NE,
				TXX, TYY, TZZ);
		// Bulk Energy flux
		bulk_energy_flux(Ubulk,  QXbulk, QYbulk, QZbulk, VX, VY, VZ);
		writeVTKvect_binary(it, "Qbulk", "e", QXbulk, QYbulk, QZbulk);
		// Computes its divergence
		divergence(divQbulk, QXbulk, QYbulk, QZbulk);
	    writeVTKscalar_binary(it, "divQbulk", "e", divQbulk);

		// Intenal Energy flux
		internal_enthalpy_flux(Uth,  QX, QY, QZ, VX, VY, VZ,
					TXX, TXY, TXZ, TYY, TYZ, TZZ);
		writeVTKvect_binary(it, "Qenth", "e", QX, QY, QZ);
		// Computes its divergence
		divergence(divQ, QX, QY, QZ);
	    writeVTKscalar_binary(it, "divQenth", "e", divQ);

	    // Computes u dot div P
	    // ATTENTION it overwrites the internal energy flux using it as service variable
	    udivP(UdivP,  QX, QY, QZ,
	    		VX, VY, VZ,
				TXX, TXY, TXZ, TYY, TYZ, TZZ);
	    writeVTKscalar_binary(it, "UdivP", "e", UdivP);


		//Compute Lagrangian terms

	    // Compute div U
	    divergence(divU, VX, VY, VZ);
	    writeVTKscalar_binary(it, "divU", "e", divU);

	    //Computes  P:gradU
		pgradU(PgradU, VX, VY, VZ, QX, QY, QZ,
					TXX, TXY, TXZ, TYY, TYZ, TZZ);
	    writeVTKscalar_binary(it, "PgradU", "e", divU);

		// ===========================================================
		// SPECIES 1 : ions
		// ===========================================================

		//Currents species1
		readvect(it,"/moments/species_1/","J",  VX, VY, VZ);
		if (ns >2) addreadvect(it,"/moments/species_3/","J",  VX, VY, VZ);

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
		
		// Note: this call expects V to be the current for compatibility with other postprocessing tools
		extract_pressure(qom[0], BX, BY, BZ, VX, VY, VZ, NI,
						 TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2, EPS);

		//Compute velocity from Current
		div(VX, NI, nxn*XLEN, nyn*YLEN, nzn*ZLEN);
		div(VY, NI, nxn*XLEN, nyn*YLEN, nzn*ZLEN);
		div(VZ, NI, nxn*XLEN, nyn*YLEN, nzn*ZLEN);

		// Computes Energies
		compute_energy(qom[1], Ubulk, Uth, VX, VY, VZ, NE,
				TXX, TYY, TZZ);
		// Bulk Energy flux
		bulk_energy_flux(Ubulk,  QXbulk, QYbulk, QZbulk, VX, VY, VZ);
		writeVTKvect_binary(it, "Qbulk", "i", QXbulk, QYbulk, QZbulk);
		// Computes its divergence
		divergence(divQbulk, QXbulk, QYbulk, QZbulk);
	    writeVTKscalar_binary(it, "divQbulk", "i", divQbulk);

		// Intenal Energy flux
		internal_enthalpy_flux(Uth,  QX, QY, QZ, VX, VY, VZ,
					TXX, TXY, TXZ, TYY, TYZ, TZZ);
		writeVTKvect_binary(it, "Qenth", "i", QX, QY, QZ);
		// Computes its divergence
		divergence(divQ, QX, QY, QZ);
	    writeVTKscalar_binary(it, "divQenth", "i", divQ);

	    // Computes u dot div P
	    // ATTENTION it overwrites the internal energy flux using it as service variable
	    udivP(UdivP,  QX, QY, QZ,
	    		VX, VY, VZ,
				TXX, TXY, TXZ, TYY, TYZ, TZZ);
	    writeVTKscalar_binary(it, "UdivP", "i", UdivP);


		//Compute Lagrangian terms

	    // Compute div U
	    divergence(divU, VX, VY, VZ);
	    writeVTKscalar_binary(it, "divU", "i", divU);

	    //Computes  P:gradU
		pgradU(PgradU, VX, VY, VZ, QX, QY, QZ,
					TXX, TXY, TXZ, TYY, TYZ, TZZ);
	    writeVTKscalar_binary(it, "PgradU", "i", divU);

}


	delArr3(BX,nxn*XLEN,nyn*YLEN);
	delArr3(BY,nxn*XLEN,nyn*YLEN);
	delArr3(BZ,nxn*XLEN,nyn*YLEN);
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
	delArr3(NE,nxn*XLEN,nyn*YLEN);			
	delArr3(NI,nxn*XLEN,nyn*YLEN);
	delArr3(QX,nxn*XLEN,nyn*YLEN);
	delArr3(QY,nxn*XLEN,nyn*YLEN);
	delArr3(QZ,nxn*XLEN,nyn*YLEN);
	delArr3(QXbulk,nxn*XLEN,nyn*YLEN);
	delArr3(QYbulk,nxn*XLEN,nyn*YLEN);
	delArr3(QZbulk,nxn*XLEN,nyn*YLEN);
	delArr3(divQ,nxn*XLEN,nyn*YLEN);
	delArr3(divQbulk,nxn*XLEN,nyn*YLEN);
	delArr3(divU,nxn*XLEN,nyn*YLEN);
	delArr3(UdivP,nxn*XLEN,nyn*YLEN);
	delArr3(Uth,nxn*XLEN,nyn*YLEN);
	delArr3(Ubulk,nxn*XLEN,nyn*YLEN);
	delArr3(PgradU,nxn*XLEN,nyn*YLEN);

	delete[] temp_storageX;
	delete[] temp_storageY;
	delete[] temp_storageZ;
	
	
	return(0);
}

