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

	temp_storageX = new double[(nxn+1)*(nyn+1)*(nzn+1)];
	temp_storageY = new double[(nxn+1)*(nyn+1)*(nzn+1)];
    temp_storageZ = new double[(nxn+1)*(nyn+1)*(nzn+1)];

// 3D Arrays
	double*** EX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** EY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** EZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** curlEX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** curlEY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** curlEZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** BX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** BY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** BZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	// V is a service array used for many things, not just speed
	double*** VX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** VY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** VZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** SX = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** SY = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** SZ = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	double*** divS = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
//	double*** dB2dt = newArr3(double,nxn*XLEN,nyn*YLEN,nzn*ZLEN);
	
	// 2D Arrays for avergaes along z indicated with last letter B meaning "bar".
	// Do not confuse it with magnetic field B
	double** EXB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** EYB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** EZB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** curlEXB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** curlEYB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** curlEZB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** BXB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** BYB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** BZB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** VXB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** VYB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** VZB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** divSB = newArr2(double,nxn*XLEN,nyn*YLEN);
//	double** dB2dtB = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** EdotJX = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** EdotJY = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** EdotJZ = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** dEdotJX = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** dEdotJY = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** dEdotJZ = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** sEdotJX = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** sEdotJY = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** sEdotJZ = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** PX = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** PY = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** PZ = newArr2(double,nxn*XLEN,nyn*YLEN);	
	double** dPX = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** dPY = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** dPZ = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** sPX = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** sPY = newArr2(double,nxn*XLEN,nyn*YLEN);
	double** sPZ = newArr2(double,nxn*XLEN,nyn*YLEN);
	
	//Electric field
	for (int it = initlevel; it < nlevels+1; it++){
    readvect(it, "/fields/","E", EX, EY, EZ);
		average(EXB,EYB,EZB,EX,EY,EZ);
		writeVTKvect(it,"E","", EXB, EYB, EZB);


	// Magnetic field
    readvect(it, "/fields/","B", BX, BY, BZ);
		average(BXB,BYB,BZB,BX,BY,BZ);
		writeVTKvect(it,"B","", BXB, BYB, BZB);

		// Magnetic energy change dB^2/dt=2 B.curlE
			curl(curlEX, curlEY, curlEZ, EX, EY, EZ);

	//		dot(dB2dt, BX, BY, BZ, curlEX, curlEY, curlEZ)
			//mult(2.0,curlEX);
			//mult(2.0,curlEY);
			//mult(2.0,curlEZ);
			average(curlEXB,curlEYB,curlEZB,curlEX,curlEY,curlEZ);
			averagedot(EdotJX,EdotJY,EdotJZ,dEdotJX,dEdotJY,dEdotJZ,
			                   sEdotJX,sEdotJY,sEdotJZ,
							   BXB,BYB,BZB,curlEXB,curlEYB,curlEZB,
							   BX,BY,BZ,curlEX,curlEY,curlEZ);
			writeVTKvect(it,"AVGdB2dt","", EdotJX, EdotJY, EdotJZ);
			writeVTKvect(it,"AVGdeldB2dt","", dEdotJX, dEdotJY, dEdotJZ);
			writeVTKvect(it,"STDdB2dt","", sEdotJX, sEdotJY, sEdotJZ);


	//Poynting Flux
    averagecross(PX, PY, PZ, dPX, dPY, dPZ, sPX, sPY, sPZ, EXB, EYB, EZB, BXB, BYB, BZB, EX, EY, EZ, BX, BY, BZ);
    writeVTKvect(it,"AVGPoy","", PX, PY, PZ);
    writeVTKvect(it,"AVGdelPoy","", dPX, dPY, dPZ);
    writeVTKvect(it,"STDPoy","", sPX, sPY, sPZ);

    //Divergence Poynting
    div_cross(EX, EY, EZ, BX, BY, BZ, SX, SY, SZ, divS);
    // Note: SX,SY,SZ have the Poynting flux, but not yet divided by 4pi,
    // same for divergence
    double invFourPi = 1.0/4.0/M_PI;
    mult(invFourPi, divS);
    mult(invFourPi, SX);
    mult(invFourPi, SY);
    mult(invFourPi, SZ);
    average(divSB,divS);
    writeVTKscalar(it,"divS","", divS);
    writeVTKvect(it,"S","", SX, SY, SZ);
    writeVTKscalar(it,"AVGdivS","", divSB);

    // Lagrangian Energy Balance Term
    // To save space now we overwrite divS and V that will now be,
    // respectively a work variable and the Lagrangian Energy Balance Term
    LEBT(BX, BY, BZ, SX, SY, SZ, VX, VY, VZ, divS);
    writeVTKscalar(it,"LagrEBTerm","", divS);

    // Reconnection measure by Biskmap
    // To save space now we overwrite S and V that will now be,
    // respectively a work variable and the Biskamp measure
    biskamp(BX, BY, BZ, EX, EY, EZ, VX, VY, VZ, SX, SY, SZ);
    writeVTKvect(it,"Biskamp","", SX, SY, SZ);

	//Currents species0
    readvect(it,"/moments/species_0/","J",  VX, VY, VZ);
    if (ns >2) addreadvect(it,"/moments/species_2/","J",  VX, VY, VZ);
		average(VXB,VYB,VZB,VX,VY,VZ);
		averagedot(EdotJX,EdotJY,EdotJZ,dEdotJX,dEdotJY,dEdotJZ,
                   sEdotJX,sEdotJY,sEdotJZ,
				   EXB,EYB,EZB,VXB,VYB,VZB,EX,EY,EZ,VX,VY,VZ);

		writeVTKvect(it,"J","e", VXB, VYB, VZB);
				writeVTKvect(it,"AVGJdotE","e", EdotJX, EdotJY, EdotJZ);
				writeVTKvect(it,"AVGdelJdotE","e", dEdotJX, dEdotJY, dEdotJZ);
				writeVTKvect(it,"STDJdotE","e", sEdotJX, sEdotJY, sEdotJZ);
		
		//Currents species1
		readvect(it,"/moments/species_1/","J",  VX, VY, VZ);
		if (ns >2) addreadvect(it,"/moments/species_3/","J",  VX, VY, VZ);
		average(VXB,VYB,VZB,VX,VY,VZ);
		averagedot(EdotJX,EdotJY,EdotJZ,dEdotJX,dEdotJY,dEdotJZ,
                   sEdotJX,sEdotJY,sEdotJZ,
				   EXB,EYB,EZB,VXB,VYB,VZB,EX,EY,EZ,VX,VY,VZ);
		
		writeVTKvect(it,"J","i", VXB, VYB, VZB);
		writeVTKvect(it,"AVGJdotE","i", EdotJX, EdotJY, EdotJZ);
		writeVTKvect(it,"AVGdelJdotE","i", dEdotJX, dEdotJY, dEdotJZ);
        writeVTKvect(it,"STDJdotE","i", sEdotJX, sEdotJY, sEdotJZ);
		
		
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
	delArr3(SX,nxn*XLEN,nyn*YLEN);
	delArr3(SY,nxn*XLEN,nyn*YLEN);
	delArr3(SZ,nxn*XLEN,nyn*YLEN);
	
	delArr2(EXB,nxn*XLEN);
	delArr2(EYB,nxn*XLEN);
	delArr2(EZB,nxn*XLEN);
	delArr2(VXB,nxn*XLEN);
	delArr2(VYB,nxn*XLEN);
	delArr2(VZB,nxn*XLEN);
	delArr2(EdotJX,nxn*XLEN);
	delArr2(EdotJY,nxn*XLEN);
	delArr2(EdotJZ,nxn*XLEN);
	delArr2(dEdotJX,nxn*XLEN);
	delArr2(dEdotJY,nxn*XLEN);
	delArr2(dEdotJZ,nxn*XLEN);
	delArr2(sEdotJX,nxn*XLEN);
	delArr2(sEdotJY,nxn*XLEN);
	delArr2(sEdotJZ,nxn*XLEN);
	delArr2(PX,nxn*XLEN);
	delArr2(PY,nxn*XLEN);
	delArr2(PZ,nxn*XLEN);
	delArr2(dPX,nxn*XLEN);
	delArr2(dPY,nxn*XLEN);
	delArr2(dPZ,nxn*XLEN);
	delArr2(sPX,nxn*XLEN);
	delArr2(sPY,nxn*XLEN);
	delArr2(sPZ,nxn*XLEN);
	delete[] temp_storageX;
	delete[] temp_storageY;
	delete[] temp_storageZ;
	
	
	return(0);
}

