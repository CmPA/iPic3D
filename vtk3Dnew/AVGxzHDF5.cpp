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
	double** EXB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** EYB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** EZB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** curlEXB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** curlEYB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** curlEZB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** BXB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** BYB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** BZB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** VXB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** VYB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** VZB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** divSB = newArr2(double,nxn*XLEN,nzn*ZLEN);
//	double** dB2dtB = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** EdotJX = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** EdotJY = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** EdotJZ = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** dEdotJX = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** dEdotJY = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** dEdotJZ = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** sEdotJX = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** sEdotJY = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** sEdotJZ = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** PX = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** PY = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** PZ = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** dPX = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** dPY = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** dPZ = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** sPX = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** sPY = newArr2(double,nxn*XLEN,nzn*ZLEN);
	double** sPZ = newArr2(double,nxn*XLEN,nzn*ZLEN);
	
	//Electric field
	for (int it = initlevel; it < nlevels+1; it++){
    readvect(it, "/fields/","E", EX, EY, EZ);
		averageY(EXB,EYB,EZB,EX,EY,EZ);
		writeVTKvect(it,"E","", EXB, EYB, EZB);


	// Magnetic field
    readvect(it, "/fields/","B", BX, BY, BZ);
		averageY(BXB,BYB,BZB,BX,BY,BZ);
		writeVTKvect(it,"B","", BXB, BYB, BZB);

		// Magnetic energy change dB^2/dt=2 B.curlE
			curl(curlEX, curlEY, curlEZ, EX, EY, EZ);

	//		dot(dB2dt, BX, BY, BZ, curlEX, curlEY, curlEZ)
			//mult(2.0,curlEX);
			//mult(2.0,curlEY);
			//mult(2.0,curlEZ);
			averageY(curlEXB,curlEYB,curlEZB,curlEX,curlEY,curlEZ);
			averageYdot(EdotJX,EdotJY,EdotJZ,dEdotJX,dEdotJY,dEdotJZ,
			                   sEdotJX,sEdotJY,sEdotJZ,
							   BXB,BYB,BZB,curlEXB,curlEYB,curlEZB,
							   BX,BY,BZ,curlEX,curlEY,curlEZ);
			writeVTKvect(it,"AVG_XZ_dB2dt","", EdotJX, EdotJY, EdotJZ);
			writeVTKvect(it,"AVG_XZ_deldB2dt","", dEdotJX, dEdotJY, dEdotJZ);
			writeVTKvect(it,"STDdB2dt","", sEdotJX, sEdotJY, sEdotJZ);


	//Poynting Flux
    averageYcross(PX, PY, PZ, dPX, dPY, dPZ, sPX, sPY, sPZ, EXB, EYB, EZB, BXB, BYB, BZB, EX, EY, EZ, BX, BY, BZ);
    writeVTKvect(it,"AVG_XZ_Poy","", PX, PY, PZ);
    writeVTKvect(it,"AVG_XZ_delPoy","", dPX, dPY, dPZ);
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
    averageY(divSB,divS);
    writeVTKscalar(it,"divS","", divS);
    writeVTKvect(it,"S","", SX, SY, SZ);
    writeVTKscalar(it,"AVG_XZ_divS","", divSB);

	//Currents species0
    readvect(it,"/moments/species_0/","J",  VX, VY, VZ);
    if (ns >2) addreadvect(it,"/moments/species_2/","J",  VX, VY, VZ);
		averageY(VXB,VYB,VZB,VX,VY,VZ);
		averageYdot(EdotJX,EdotJY,EdotJZ,dEdotJX,dEdotJY,dEdotJZ,
                   sEdotJX,sEdotJY,sEdotJZ,
				   EXB,EYB,EZB,VXB,VYB,VZB,EX,EY,EZ,VX,VY,VZ);

		writeVTKvect(it,"J","e", VXB, VYB, VZB);
				writeVTKvect(it,"AVG_XZ_JdotE","e", EdotJX, EdotJY, EdotJZ);
				writeVTKvect(it,"AVG_XZ_delJdotE","e", dEdotJX, dEdotJY, dEdotJZ);
				writeVTKvect(it,"STDJdotE","e", sEdotJX, sEdotJY, sEdotJZ);
		
		//Currents species1
		readvect(it,"/moments/species_1/","J",  VX, VY, VZ);
		if (ns >2) addreadvect(it,"/moments/species_3/","J",  VX, VY, VZ);
		averageY(VXB,VYB,VZB,VX,VY,VZ);
		averageYdot(EdotJX,EdotJY,EdotJZ,dEdotJX,dEdotJY,dEdotJZ,
                   sEdotJX,sEdotJY,sEdotJZ,
				   EXB,EYB,EZB,VXB,VYB,VZB,EX,EY,EZ,VX,VY,VZ);
		
		writeVTKvect(it,"J","i", VXB, VYB, VZB);
		writeVTKvect(it,"AVG_XZ_JdotE","i", EdotJX, EdotJY, EdotJZ);
		writeVTKvect(it,"AVG_XZ_delJdotE","i", dEdotJX, dEdotJY, dEdotJZ);
        writeVTKvect(it,"STDJdotE","i", sEdotJX, sEdotJY, sEdotJZ);
		
		
}

	delArr3(EX,nxn*XLEN,nzn*ZLEN);
	delArr3(EY,nxn*XLEN,nzn*ZLEN);
	delArr3(EZ,nxn*XLEN,nzn*ZLEN);
	delArr3(BX,nxn*XLEN,nzn*ZLEN);
	delArr3(BY,nxn*XLEN,nzn*ZLEN);
	delArr3(BZ,nxn*XLEN,nzn*ZLEN);
	delArr3(VX,nxn*XLEN,nzn*ZLEN);
	delArr3(VY,nxn*XLEN,nzn*ZLEN);
	delArr3(VZ,nxn*XLEN,nzn*ZLEN);
	delArr3(SX,nxn*XLEN,nzn*ZLEN);
	delArr3(SY,nxn*XLEN,nzn*ZLEN);
	delArr3(SZ,nxn*XLEN,nzn*ZLEN);
	
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

