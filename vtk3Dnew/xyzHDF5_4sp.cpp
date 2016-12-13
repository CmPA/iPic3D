/***************************************************************************
 convHDF5.cpp  -  Convert program to open Parsek Output
 -------------------
 begin                : Jun 2008
 copyright            : (C) 2004 by Stefano Markidis, Giovanni Lapenta
 ************************************************************************** */

#include "hdf5.h"
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
	double*** OHMX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** OHMY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** OHMZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** NE = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** NI = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** AZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	

	string temp;

		temp = "recflux.txt";
		ofstream my_file(temp.c_str());

for (int it = initlevel; it < nlevels+1; it++){

	//Electric field

    readvect(it, "/fields/","E", EX, EY, EZ);
    writeVTKvect(it, "E", "", EX, EY, EZ);

	//Magnetic field
    readvect(it,"/fields/","B", BX, BY, BZ);
    writeVTKvect(it,"B", "", BX, BY, BZ);

    double recflux = 0.0;
    recflux = vecpot(AZ, BX, BY, BZ);
    writeVTKscalar(it, "Az", "", AZ);
    cout << recflux << endl;

    my_file << it * DeltaT << "   " << recflux << endl;

	//Compute ExB
    //cross(BX, BY, BZ, EX, EY, EZ, VXBX, VXBY, VXBZ)

	//Rho by species
    readscalar(it,"/moments/species_0/","rho",  NE);
    readscalar(it,"/moments/species_1/","rho",  NI);
    if (ns >2) addreadscalar(it,"/moments/species_2/","rho",  NE);
    if (ns >2) addreadscalar(it,"/moments/species_3/","rho",  NI);
    writeVTKscalar_species(it,"rho", NE, NI);
   // writeVTKscalar("rho", "e", VX);

	//Currents species0
    readvect(it,"/moments/species_0/","J",  VX, VY, VZ);
    if (ns >2) addreadvect(it,"/moments/species_2/","J",  VX, VY, VZ);
    writeVTKvect(it,"J", "e", VX, VY, VZ);
    
	//Evaluate and save E+VxB
    ohm(NE, BX, BY, BZ, EX, EY, EZ, VX, VY, VZ, OHMX, OHMY, OHMZ);
    writeVTKvect(it,"EMF", "e", VX, VY, VZ);
    vnonfrozen(NE, BX, BY, BZ, EX, EY, EZ, VX, VY, VZ);
    writeVTKvect(it,"VNF", "e", VX, VY, VZ);

    //Currents species1
    readvect(it,"/moments/species_1/","J",  VX, VY, VZ);
    if (ns >2) addreadvect(it,"/moments/species_3/","J",  VX, VY, VZ);
    writeVTKvect(it,"J", "i", VX, VY, VZ);

	//Evaluate and save E+VxB
    ohm(NI, BX, BY, BZ, EX, EY, EZ, VX, VY, VZ, OHMX, OHMY, OHMZ);
    writeVTKvect(it,"EMF", "i", VX, VY, VZ);
    vnonfrozen(NI, BX, BY, BZ, EX, EY, EZ, VX, VY, VZ);
    writeVTKvect(it,"VNF", "i", VX, VY, VZ);

}

	my_file.close();

	delArr3(BX,nxn*XLEN,nyn*YLEN);
	delArr3(BY,nxn*XLEN,nyn*YLEN);
	delArr3(BZ,nxn*XLEN,nyn*YLEN);
	delArr3(EX,nxn*XLEN,nyn*YLEN);
	delArr3(EY,nxn*XLEN,nyn*YLEN);
	delArr3(EZ,nxn*XLEN,nyn*YLEN);
	delArr3(VX,nxn*XLEN,nyn*YLEN);
	delArr3(VY,nxn*XLEN,nyn*YLEN);
	delArr3(VZ,nxn*XLEN,nyn*YLEN);
	delArr3(AZ,nxn*XLEN,nyn*YLEN);
	delete[] temp_storageX;
	delete[] temp_storageY;
	delete[] temp_storageZ;
	
	
	return(0);
}

