/*
 *  VtkStuff3D.h
 *
 *
 *  Created by Giovanni Lapenta on 7/29/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#include "hdf5.h"
#include "Alloc.h"
#include "math.h"



#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
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


void writeGDAvect(string vectname, string addname, double **EX, double **EY,double **EZ);

void writeGDAscalar(string scalarname, string addname, double **f);

void writeGDAscalar_species(string scalarname,  double **f1, double **f2);

void writeGDAtensor(string tensorname, string addname, double** pXX, double** pXY,
		double** pXZ, double** pYY, double** pYZ, double** pZZ,
		double** pPAR, double** pPER1, double** pPER2);




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



void writeGDAvect(string vectname, string addname, double** EX, double** EY,double** EZ) {



	string temp;
	temp = vectname +addname;
	temp += "_x.gda";
	ofstream foo(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo.write((char*)&EX[jj][it],sizeof(double));

			}
	foo.close();

	temp = vectname +addname;
	temp += "_y.gda";
	ofstream foo2(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo2.write((char*)&EY[jj][it],sizeof(double));

			}
	foo2.close();

	temp = vectname +addname;
	temp += "_z.gda";
	ofstream foo3(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo3.write((char*)&EZ[jj][it],sizeof(double));

			}
	foo3.close();
}


void writeGDAscalar(string scalarname, string addname, double** EX) {


	string temp;
	temp = scalarname +addname;
	temp += ".gda";
	ofstream foo(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo.write((char*)&EX[jj][it],sizeof(double));
			}
	foo.close();
}


void writeGDAscalar_species(int it, string scalarname, double** EX, double** EY) {


	string temp;
	temp = scalarname;
	temp += "_0.gda";
	ofstream foo(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo.write((char*)&EX[jj][it],sizeof(double));
			}
	foo.close();
	temp = scalarname ;
	temp += "_1.gda";
	ofstream foo2(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo2.write((char*)&EY[jj][it],sizeof(double));
			}
	foo2.close();


}

void writeGDAscalar_species(string scalarname, double** EX, double** EY, double** EZ) {
	string temp;
	temp = scalarname ;
	temp += "_0.gda";
	ofstream foo(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo.write((char*)&EX[jj][it],sizeof(double));
			}
	foo.close();
	temp = scalarname ;
	temp += "_1.gda";
	ofstream foo2(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo2.write((char*)&EY[jj][it],sizeof(double));
			}
	foo2.close();
	temp = scalarname ;
	temp += "_2.gda";
	ofstream foo3(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo3.write((char*)&EZ[jj][it],sizeof(double));
			}
	foo3.close();

}


void writeGDAtensor(string vectname, string addname,
		double** TXX, double** TXY,double** TXZ,
		double** TYY, double** TYZ,double** TZZ,
		double** TPAR, double** TPER1,double** TPER2, double** EPS) {

	string temp;
	temp = vectname +addname;
	temp += "_xx.gda";
	ofstream foo(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo.write((char*)&TXX[jj][it],sizeof(double));

			}
	foo.close();

	temp = vectname +addname;
	temp += "_xy.gda";
	ofstream foo2(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo2.write((char*)&TXY[jj][it],sizeof(double));

			}
	foo2.close();

	temp = vectname +addname;
	temp += "_xz.gda";
	ofstream foo3(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo3.write((char*)&TXZ[jj][it],sizeof(double));

			}
	foo3.close();

	temp = vectname +addname;
	temp += "_yy.gda";
	ofstream foo4(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo4.write((char*)&TYY[jj][it],sizeof(double));

			}
	foo4.close();

	temp = vectname +addname;
	temp += "_yz.gda";
	ofstream foo5(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo5.write((char*)&TYZ[jj][it],sizeof(double));

			}
	foo5.close();

	temp = vectname +addname;
	temp += "_zz.gda";
	ofstream foo6(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo6.write((char*)&TZZ[jj][it],sizeof(double));

			}
	foo6.close();

	temp = vectname +addname;
	temp += "_par.gda";
	ofstream foo7(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo7.write((char*)&TPAR[jj][it],sizeof(double));

			}
	foo7.close();

	temp = vectname +addname;
	temp += "_per1.gda";
	ofstream foo8(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo8.write((char*)&TPER1[jj][it],sizeof(double));

			}
	foo8.close();

	temp = vectname +addname;
	temp += "_per2.gda";
	ofstream foo9(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo9.write((char*)&TPER2[jj][it],sizeof(double));

			}
	foo9.close();

	temp = vectname +addname;
	temp += "_eps.gda";
	ofstream foo10(temp.c_str(),ios::out | ios::binary);
	for (int it=0; it < nlevels;it++)
		for (int jj=0; jj < nxn*XLEN;jj++){
				foo10.write((char*)&EPS[jj][it],sizeof(double));

			}
	foo10.close();

}


