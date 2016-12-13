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

	if(argc>3) {
			sscanf(argv[4], "%d", &NdimCode);
		}
		else
		{
			NdimCode = 3;
		}

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
	double*** TXX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** TXY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** TXZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** TYY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** TYZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** TZZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** TPAR  = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** TPER1 = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** TPER2 = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** VDIVP = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** DIVJV = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** JDOTE = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** VSQ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);	
	double*** EPS = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** NE = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** NI = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** pdvpar = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** pdvper1 = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** pdvper2 = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	// workspace
	double*** WX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** WY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** WZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);

	double** AVGENTH = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** AVGENIN = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** AVGJDOTE = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** STDENTH = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** STDENIN = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** STDJDOTE = newArr2(double,nxn*XLEN,nyn*YLEN);
	
	for (int it = initlevel; it < nlevels+1; it++){

	//Electric field
    readvect(it,"/fields/","E", EX, EY, EZ);
    //writeVTKvect(it,"E", "", EX, EY, EZ);
    
	//Magnetic field
    readvect(it,"/fields/","B", BX, BY, BZ);
    //writeVTKvect(it,"B", "", BX, BY, BZ);

	//Rho by species
    readscalar(it,"/moments/species_0/","rho",  NE);
    readscalar(it,"/moments/species_1/","rho",  NI);
    if (ns >2) addreadscalar(it,"/moments/species_2/","rho",  NE);
    if (ns >2) addreadscalar(it,"/moments/species_3/","rho",  NI);
    //writeVTKscalar_species(it,"rho", NE, NI);
   // writeVTKscalar("rho", "e", VX);

	//Currents species0
    readvect(it,"/moments/species_0/","J",  VX, VY, VZ);
    if (ns >2) addreadvect(it,"/moments/species_2/","J",  VX, VY, VZ);
    //writeVTKvect(it,"J", "e", VX, VY, VZ);
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
    writeVTKtensor(it, "P", "e", TXX, TXY, TXZ, TYY, TYZ, TZZ, 
    	  TPAR, TPER1, TPER2, EPS);
        
    dot(VSQ, VX, VY, VZ, VX, VY, VZ, NE);
    vdiv(VDIVP, TXX, TYY, TZZ, TXY, TXZ, TYZ, NE, VX, VY, VZ);
	divj(qom[0], DIVJV, VX, VY, VZ, VSQ);
	dot(JDOTE, VX, VY, VZ, EX, EY, EZ);
	writeVTKvect(it,"Energy", "e", VDIVP, DIVJV, JDOTE);
	
	pdv(pdvpar, pdvper1, pdvper2, TXX, TXY, TXZ, TYY, TYZ, TZZ,
			BX, BY, BZ, NE, VX, VY, VZ, VDIVP, WX, WY, WZ);
	writeVTKvect(it,"PdV", "e", pdvpar, pdvper1, pdvper2);

	average_energy(AVGENTH, AVGENIN, AVGJDOTE,
			STDENTH, STDENIN, STDJDOTE,
			VDIVP, DIVJV, JDOTE);
    writeVTKvect(it,"AVGEnergy","e", AVGENTH, AVGENIN, AVGJDOTE);
    writeVTKvect(it,"STDEnergy","e", STDENTH, STDENIN, STDJDOTE);

    //Currents species1
    readvect(it,"/moments/species_1/","J",  VX, VY, VZ);
    if (ns >2) addreadvect(it,"/moments/species_3/","J",  VX, VY, VZ);
    //writeVTKvect(it,"J", "i", VX, VY, VZ);
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
    writeVTKtensor(it, "P", "i", TXX, TXY, TXZ, TYY, TYZ, TZZ, 
    	  TPAR, TPER1, TPER2, EPS);
    	  
    dot(VSQ, VX, VY, VZ, VX, VY, VZ, NI);
    vdiv(VDIVP, TXX, TYY, TZZ, TXY, TXZ, TYZ, NI, VX, VY, VZ);
	divj(qom[1], DIVJV, VX, VY, VZ, VSQ);
	dot(JDOTE, VX, VY, VZ, EX, EY, EZ);
	writeVTKvect(it,"Energy", "i", VDIVP, DIVJV, JDOTE);

  	pdv(pdvpar, pdvper1, pdvper2, TXX, TXY, TXZ, TYY, TYZ, TZZ,
  			BX, BY, BZ, NI, VX, VY, VZ, VDIVP, WX, WY, WZ);
	writeVTKvect(it,"PdV", "i", pdvpar, pdvper1, pdvper2);

	average_energy(AVGENTH, AVGENIN, AVGJDOTE,
			STDENTH, STDENIN, STDJDOTE,
			VDIVP, DIVJV, JDOTE);
    writeVTKvect(it,"AVGEnergy","i", AVGENTH, AVGENIN, AVGJDOTE);
    writeVTKvect(it,"STDEnergy","i", STDENTH, STDENIN, STDJDOTE);

    //Potential
//    readpotential(it,"/potentials/","phi",  NI);
//    writeVTKscalar(it,"phi","", NI);
}
	delArr3(EX,nxn*XLEN,nyn*YLEN);
	delArr3(EY,nxn*XLEN,nyn*YLEN);
	delArr3(EZ,nxn*XLEN,nyn*YLEN);
	
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
	
	delArr3(VSQ,nxn*XLEN,nyn*YLEN);
	delArr3(VDIVP,nxn*XLEN,nyn*YLEN);
	delArr3(DIVJV,nxn*XLEN,nyn*YLEN);
	delArr3(JDOTE,nxn*XLEN,nyn*YLEN);
	
	delArr3(pdvpar,nxn*XLEN,nyn*YLEN);
	delArr3(pdvper1,nxn*XLEN,nyn*YLEN);
	delArr3(pdvper2,nxn*XLEN,nyn*YLEN);

	delArr3(WX,nxn*XLEN,nyn*YLEN);
	delArr3(WY,nxn*XLEN,nyn*YLEN);
	delArr3(WZ,nxn*XLEN,nyn*YLEN);

	delArr2(AVGENTH, nxn*XLEN);
	delArr2(AVGENIN, nxn*XLEN);
	delArr2(AVGJDOTE, nxn*XLEN);
	delArr2(STDENTH, nxn*XLEN);
	delArr2(STDENIN, nxn*XLEN);
	delArr2(STDJDOTE, nxn*XLEN);

	delete[] temp_storageX;
	delete[] temp_storageY;
	delete[] temp_storageZ;
	
	
	return(0);
}
