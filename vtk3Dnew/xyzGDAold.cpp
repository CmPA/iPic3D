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


using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;


void writeGDAvect(int it, string vectname, string addname, double ***EX, double ***EY,double ***EZ);

void writeGDAscalar(int it, string scalarname, string addname, double ***f);

void writeGDAscalar_species(int it, string scalarname,  double ***f1, double ***f2);

void writeGDAscalar_species(int it, string scalarname,  double ***f1, double ***f2, double ***f3);

void writeGDAtensor(int it, string vectname, string addname,
		double*** TXX, double*** TXY, double*** TXZ,
		double*** TYY, double*** TYZ, double*** TZZ,
		double*** TPAR, double*** TPER1, double*** TPER2, double*** EPS);

// binary files 
FILE *fp;

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
	
	//Electric field
	for (int it = initlevel; it < nlevels+1; it++){
    readvect(it, "/fields/","E", EX, EY, EZ);
    writeGDAvect(it, "E", "", EX, EY, EZ);

	//Magnetic field
    readvect(it,"/fields/","B", BX, BY, BZ);
    writeGDAvect(it,"B", "", BX, BY, BZ);

	//Compute ExB
    //cross(BX, BY, BZ, EX, EY, EZ, VXBX, VXBY, VXBZ)

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

     vnonfrozen(NE, BX, BY, BZ, EX, EY, EZ, VX, VY, VZ);
     writeGDAvect(it,"VNF", "e", VX, VY, VZ);

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

     vnonfrozen(NI, BX, BY, BZ, EX, EY, EZ, VX, VY, VZ);
     writeGDAvect(it,"VNF", "i", VX, VY, VZ);

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

	delete[] temp_storageX;
	delete[] temp_storageY;
	delete[] temp_storageZ;
	
	
	return(0);
}



void writeGDAvect(int it, string vectname, string addname, double*** EX, double*** EY,double*** EZ) {

	

    int cycle = it * DeltaT;
	stringstream stringcycle;
	stringcycle << cycle;
	string temp;
	temp = vectname +addname +"_x_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo.write((char*)&EX[ii][jj][kk],sizeof(double));

			}
	foo.close();
	
	temp = vectname +addname +"_y_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo2(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo2.write((char*)&EY[ii][jj][kk],sizeof(double));

			}
	foo2.close();
	
	temp = vectname +addname +"_z_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo3(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo3.write((char*)&EZ[ii][jj][kk],sizeof(double));

			}
	foo3.close();
}


void writeGDAscalar(int it, string scalarname, string addname, double*** EX) {


    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +addname +"_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo.write((char*)&EX[ii][jj][kk],sizeof(double));
			}
	foo.close();
}


void writeGDAscalar_species(int it, string scalarname, double*** EX, double*** EY) {


    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_0_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo.write((char*)&EX[ii][jj][kk],sizeof(double));
			}
	foo.close();
	temp = scalarname +"_1_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo2(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo2.write((char*)&EY[ii][jj][kk],sizeof(double));
			}
	foo2.close();
	
	
}

void writeGDAscalar_species(int it, string scalarname, double*** EX, double*** EY, double*** EZ) {
    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_0_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo.write((char*)&EX[ii][jj][kk],sizeof(double));
			}
	foo.close();
	temp = scalarname +"_1_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo2(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo2.write((char*)&EY[ii][jj][kk],sizeof(double));
			}
	foo2.close();
	temp = scalarname +"_2_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo3(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo3.write((char*)&EZ[ii][jj][kk],sizeof(double));
			}
	foo3.close();
	
}


void writeGDAtensor(int it, string vectname, string addname,
		double*** TXX, double*** TXY,double*** TXZ,
		double*** TYY, double*** TYZ,double*** TZZ,
		double*** TPAR, double*** TPER1,double*** TPER2, double*** EPS) {



    int cycle = it * DeltaT;
	stringstream stringcycle;
	stringcycle << cycle;
	string temp;
	temp = vectname +addname +"_xx_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo.write((char*)&TXX[ii][jj][kk],sizeof(double));

			}
	foo.close();

	temp = vectname +addname +"_xy_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo2(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo2.write((char*)&TXY[ii][jj][kk],sizeof(double));

			}
	foo2.close();

	temp = vectname +addname +"_xz_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo3(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo3.write((char*)&TXZ[ii][jj][kk],sizeof(double));

			}
	foo3.close();

	temp = vectname +addname +"_yy_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo4(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo4.write((char*)&TYY[ii][jj][kk],sizeof(double));

			}
	foo4.close();

	temp = vectname +addname +"_yz_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo5(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo5.write((char*)&TYZ[ii][jj][kk],sizeof(double));

			}
	foo5.close();

	temp = vectname +addname +"_zz_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo6(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo6.write((char*)&TZZ[ii][jj][kk],sizeof(double));

			}
	foo6.close();

	temp = vectname +addname +"_par_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo7(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo7.write((char*)&TPAR[ii][jj][kk],sizeof(double));

			}
	foo7.close();

	temp = vectname +addname +"_per1_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo8(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo8.write((char*)&TPER1[ii][jj][kk],sizeof(double));

			}
	foo8.close();

	temp = vectname +addname +"_per2_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo9(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo9.write((char*)&TPER2[ii][jj][kk],sizeof(double));

			}
	foo9.close();

	temp = vectname +addname +"_eps_cycle" +stringcycle.str();
	temp += ".gda";
	ofstream foo10(temp.c_str(),ios::out | ios::binary);
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				foo10.write((char*)&EPS[ii][jj][kk],sizeof(double));

			}
	foo10.close();

}
