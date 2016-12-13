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
	double*** JEX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** JEY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** JEZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** JIX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** JIY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** JIZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** DE = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** DI = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** NE = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** NI = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
    double** DEB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** DIB = newArr2(double,nxn*XLEN,nyn*YLEN);
    double** DESTD = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** DISTD = newArr2(double,nxn*XLEN,nyn*YLEN);
	
	//Electric field
	for (int it = initlevel; it < nlevels+1; it++){
    readvect(it, "/fields/","E", EX, EY, EZ);

	//Magnetic field
    readvect(it,"/fields/","B", BX, BY, BZ);

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
    readvect(it,"/moments/species_0/","J",  JEX, JEY, JEZ);
    if (ns >2) addreadvect(it,"/moments/species_2/","J",  JEX, JEY, JEZ);

    //Currents species1
    readvect(it,"/moments/species_1/","J",  JIX, JIY, JIZ);
    if (ns >2) addreadvect(it,"/moments/species_3/","J",  JIX, JIY, JIZ);

    
	//Evaluate and save E+VxB
    zenitani(NE, NI, BX, BY, BZ, EX, EY, EZ, JEX, JEY, JEZ, JIX, JIY, JIZ, DE, DI);
    writeVTKscalar_species(it,"Zenitani", DE, DI);

    averageSTD(DEB,DESTD,DE);
    averageSTD(DIB,DISTD,DI);
    writeVTKscalar_species(it,"AVGZenitani", DEB, DIB);
    writeVTKscalar_species(it,"STDZenitani", DESTD, DISTD);


}
	delArr3(BX,nxn*XLEN,nyn*YLEN);
	delArr3(BY,nxn*XLEN,nyn*YLEN);
	delArr3(BZ,nxn*XLEN,nyn*YLEN);
	delArr3(EX,nxn*XLEN,nyn*YLEN);
	delArr3(EY,nxn*XLEN,nyn*YLEN);
	delArr3(EZ,nxn*XLEN,nyn*YLEN);
	delArr3(JEX,nxn*XLEN,nyn*YLEN);
	delArr3(JEY,nxn*XLEN,nyn*YLEN);
	delArr3(JEZ,nxn*XLEN,nyn*YLEN);
	delArr3(JIX,nxn*XLEN,nyn*YLEN);
	delArr3(JIY,nxn*XLEN,nyn*YLEN);
	delArr3(JIZ,nxn*XLEN,nyn*YLEN);
	delArr3(DE,nxn*XLEN,nyn*YLEN);
	delArr3(DI,nxn*XLEN,nyn*YLEN);
    delArr2(DEB,nxn*XLEN);
	delArr2(DIB,nxn*XLEN);
    delArr2(DESTD,nxn*XLEN);
	delArr2(DISTD,nxn*XLEN);
	delete[] temp_storageX;
	delete[] temp_storageY;
	delete[] temp_storageZ;
	
	
	return(0);
}

