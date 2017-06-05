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


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

int main (int argc, char **argv) {
	// cycle we want to open
	sscanf(argv[1],"%d",&InitT);
	sscanf(argv[2],"%d",&MaxLevel);
	sscanf(argv[3],"%d",&DeltaT);
	sscanf(argv[4],"%d",&Species);

	argc = argc -1; // the name of the program is added by the compiler as extra element in the array

	if(argc>4) sscanf(argv[5], "%d", &NdimCode);
	else NdimCode = 3;

	if(argc>5) sscanf(argv[6],"%d",&xshift);
		else xshift=0;

	if(argc>6) sscanf(argv[7],"%d",&yshift);
		else yshift=0;

	if(argc>7) sscanf(argv[8],"%d",&zshift);
		else zshift=0;

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

	double*** NS = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);


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
	double*** agyro_scudder = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** agyro_aunai = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** nongyro_swisdak = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** align = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);

	string species_id;

    stringstream species_name;
    species_name << Species;


	for (int it = initlevel+rank; it < nlevels+1; it=it+size){
	     printf( "Process %d doing level %d\n", rank, it );

	//Electric field

    readvect(it, "/fields/","E", EX, EY, EZ);
    writeVTKvect_binary(it, "E", "", EX, EY, EZ);

	//Magnetic field
    readvect(it,"/fields/","B", BX, BY, BZ);
    writeVTKvect_binary(it,"B", "", BX, BY, BZ);


	//Read Rho
    species_id = "/moments/species_" + species_name.str() +"/";
    cout << species_id <<endl;

    readscalar(it,species_id,"rho",  NS);
    writeVTKscalar_binary(it, "rho", species_name.str(), NS);

	//Currents species0
    readvect(it,species_id,"J",  VX, VY, VZ);
    writeVTKvect_binary(it,"J", species_name.str(), VX, VY, VZ);


    //Pressure tensor species 0
    readscalar(it,species_id,"pXX",  TXX);
    readscalar(it,species_id,"pXY",  TXY);
    readscalar(it,species_id,"pXZ",  TXZ);
    readscalar(it,species_id,"pYY",  TYY);
    readscalar(it,species_id,"pYZ",  TYZ);
    readscalar(it,species_id,"pZZ",  TZZ);

    
    extract_pressure(qom[Species], BX, BY, BZ, VX, VY, VZ, NS,
          TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2, EPS);
    writeVTKtensor_binary(it, "P", species_name.str(), TXX, TXY, TXZ, TYY, TYZ, TZZ,
    	  TPAR, TPER1, TPER2, EPS);


    int Nsmooth = 5;
    smooth(Nsmooth, TXX, nxn*XLEN, nyn*YLEN, nzn*ZLEN);
    smooth(Nsmooth, TXY, nxn*XLEN, nyn*YLEN, nzn*ZLEN);
    smooth(Nsmooth, TXZ, nxn*XLEN, nyn*YLEN, nzn*ZLEN);
    smooth(Nsmooth, TYY, nxn*XLEN, nyn*YLEN, nzn*ZLEN);
    smooth(Nsmooth, TYZ, nxn*XLEN, nyn*YLEN, nzn*ZLEN);
    smooth(Nsmooth, TZZ, nxn*XLEN, nyn*YLEN, nzn*ZLEN);
    agyro(agyro_scudder, agyro_aunai, nongyro_swisdak, align,
    	BX, BY, BZ, TXX, TXY, TXZ, TYY, TYZ, TZZ);


    writeVTKscalar_binary(it, "agyro_scudder", species_name.str(), agyro_scudder);
    writeVTKscalar_binary(it, "agyro_aunai", species_name.str(), agyro_aunai);
    writeVTKscalar_binary(it, "nongyro_swisdak", species_name.str(), nongyro_swisdak);
    writeVTKscalar_binary(it, "align", species_name.str(), align);


}


	delArr3(BX,nxn*XLEN,nyn*YLEN);
	delArr3(BY,nxn*XLEN,nyn*YLEN);
	delArr3(BZ,nxn*XLEN,nyn*YLEN);
	delArr3(EX,nxn*XLEN,nyn*YLEN);
	delArr3(EY,nxn*XLEN,nyn*YLEN);
	delArr3(EZ,nxn*XLEN,nyn*YLEN);
	delArr3(NS,nxn*XLEN,nyn*YLEN);
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
	delArr3(agyro_scudder,nxn*XLEN,nyn*YLEN);
	delArr3(agyro_aunai,nxn*XLEN,nyn*YLEN);
	delArr3(nongyro_swisdak,nxn*XLEN,nyn*YLEN);
	delArr3(align,nxn*XLEN,nyn*YLEN);

	delete[] temp_storageX;
	delete[] temp_storageY;
	delete[] temp_storageZ;
	
	
    MPI_Finalize();
    return(0);
}

