/***************************************************************************
 convHDF5.cpp  -  Convert program to open Parsek Output
 -------------------
 begin                : Jun 2008
 copyright            : (C) 2004 by Stefano Markidis, Giovanni Lapenta
 ************************************************************************** */

#include "hdf5.h"
#include "Alloc.h"
#include "math.h"
#include "GdaStuff.h"
#include "timeplot.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using std::string;
using std::stringstream;
using std::ofstream;
using std::cout;
using std::endl;

int readsettings();

void readvect(string campo, string vectname, double **EX, double **EY,double **EZ);

void readscalar(string campo, string scalarname, double **EX);

void addreadvect(string campo, string vectname, double **EX, double **EY,double **EZ);

void addreadscalar(string campo, string scalarname, double **EX);

void readpotential(string campo, string scalarname, double **POT);

void extract_pressure(double qom, 
		double** BX, double** BY, double** BZ,
		double** VX, double** VY, double** VZ,
		double** N, double** pXX, double** pXY,
		double** pXZ, double** pYY, double** pYZ, double** pZZ,
		double** pPAR, double** pPER1, double** pPER2);


// x position in processor topology. it will read the first cell in x in that processor
int jproc;
// z position in processor topology. it will read the first cell in x in that processor
int kproc;



// binary files
FILE *fp;

int main (int argc, char **argv) {
	// cycle we want to open
	sscanf(argv[3],"%d",&MaxLevel);
	sscanf(argv[2],"%d",&DeltaT);
	sscanf(argv[1],"%d",&InitLevel);
	sscanf(argv[4],"%d",&jproc);
	sscanf(argv[5],"%d",&kproc);
	nlevels = (MaxLevel-InitLevel)/DeltaT+1;

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

	double** BX = newArr2(double,nxn*XLEN,nlevels);
	double** BY = newArr2(double,nxn*XLEN,nlevels);
	double** BZ = newArr2(double,nxn*XLEN,nlevels);
	double** VX = newArr2(double,nxn*XLEN,nlevels);
	double** VY = newArr2(double,nxn*XLEN,nlevels);
	double** VZ = newArr2(double,nxn*XLEN,nlevels);
	double** TXX= newArr2(double,nxn*XLEN,nlevels);
	double** TXY = newArr2(double,nxn*XLEN,nlevels);
	double** TXZ = newArr2(double,nxn*XLEN,nlevels);
	double** TYY = newArr2(double,nxn*XLEN,nlevels);
	double** TYZ = newArr2(double,nxn*XLEN,nlevels);
	double** TZZ = newArr2(double,nxn*XLEN,nlevels);
	double** TPAR= newArr2(double,nxn*XLEN,nlevels);
	double** TPER1 = newArr2(double,nxn*XLEN,nlevels);
	double** TPER2 = newArr2(double,nxn*XLEN,nlevels);
	double** NE = newArr2(double,nxn*XLEN,nlevels);
	double** NI = newArr2(double,nxn*XLEN,nlevels);
	
	//Electric field
    readvect("/fields/","E", VX, VY, VZ);
    writeGDAvect("E", "", VX, VY, VZ);

	//Magnetic field
    readvect("/fields/","B", BX, BY, BZ);
    writeGDAvect("B", "", BX, BY, BZ);

	//Rho by species
    readscalar("/moments/species_0/","rho",  NE);
    readscalar("/moments/species_1/","rho",  NI);
    if (ns >2) addreadscalar("/moments/species_2/","rho",  NE);
    if (ns >3) addreadscalar("/moments/species_3/","rho",  NI);
    writeGDAscalar_species("rho", NE, NI);
   // writeGDAscalar("rho", "e", VX);

	//Currents species0
    readvect("/moments/species_0/","J",  VX, VY, VZ);
    if (ns >2) addreadvect("/moments/species_2/","J",  VX, VY, VZ);
    writeGDAvect("J", "e", VX, VY, VZ);
    //Pressure tensor species 0
    readscalar("/moments/species_0/","pXX",  TXX);
    readscalar("/moments/species_0/","pXY",  TXY);
    readscalar("/moments/species_0/","pXZ",  TXZ);
    readscalar("/moments/species_0/","pYY",  TYY);
    readscalar("/moments/species_0/","pYZ",  TYZ);
    readscalar("/moments/species_0/","pZZ",  TZZ);
    
    if (ns >2) addreadscalar("/moments/species_2/","pXX",  TXX);
    if (ns >2) addreadscalar("/moments/species_2/","pXY",  TXY);
    if (ns >2) addreadscalar("/moments/species_2/","pXZ",  TXZ);
    if (ns >2) addreadscalar("/moments/species_2/","pYY",  TYY);
    if (ns >2) addreadscalar("/moments/species_2/","pYZ",  TYZ);
    if (ns >2) addreadscalar("/moments/species_2/","pZZ",  TZZ);
    
    extract_pressure(qom[0], BX, BY, BZ, VX, VY, VZ, NE, 
          TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2);
    writeGDAtensor("P", "e", TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2);

    //Currents species1
    readvect("/moments/species_1/","J",  VX, VY, VZ);
    if (ns >3) addreadvect("/moments/species_3/","J",  VX, VY, VZ);
    writeGDAvect("J", "i", VX, VY, VZ);
    //Pressure tensor species 1
    readscalar("/moments/species_1/","pXX",  TXX);
    readscalar("/moments/species_1/","pXY",  TXY);
    readscalar("/moments/species_1/","pXZ",  TXZ);
    readscalar("/moments/species_1/","pYY",  TYY);
    readscalar("/moments/species_1/","pYZ",  TYZ);
    readscalar("/moments/species_1/","pZZ",  TZZ);
    
    if (ns >3) addreadscalar("/moments/species_3/","pXX",  TXX);
    if (ns >3) addreadscalar("/moments/species_3/","pXY",  TXY);
    if (ns >3) addreadscalar("/moments/species_3/","pXZ",  TXZ);
    if (ns >3) addreadscalar("/moments/species_3/","pYY",  TYY);
    if (ns >3) addreadscalar("/moments/species_3/","pYZ",  TYZ);
    if (ns >3) addreadscalar("/moments/species_3/","pZZ",  TZZ);
    

    extract_pressure(qom[1], BX, BY, BZ, VX, VY, VZ, NI, 
          TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2);
    writeGDAtensor("P", "i", TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2);
    

    //Potential
/*    readpotential("/potentials/","phi",  NI);
    writeGDAscalar("phi","", NI);
*/

	delArr2(BX,nxn*XLEN);
	delArr2(BY,nxn*XLEN);
	delArr2(BZ,nxn*XLEN);
	delArr2(VX,nxn*XLEN);
	delArr2(VY,nxn*XLEN);
	delArr2(VZ,nxn*XLEN);
	delArr2(TXX,nxn*XLEN);
	delArr2(TXY,nxn*XLEN);
	delArr2(TXZ,nxn*XLEN);
	delArr2(TYY,nxn*XLEN);
	delArr2(TYZ,nxn*XLEN);
	delArr2(TZZ,nxn*XLEN);
	delArr2(TPAR,nxn*XLEN);
	delArr2(TPER1,nxn*XLEN);
	delArr2(TPER2,nxn*XLEN);
	delArr2(NE,nxn*XLEN);
	delArr2(NI,nxn*XLEN);
	delete[] temp_storageX;
	delete[] temp_storageY;
	delete[] temp_storageZ;
	
	
	return(0);
}


int readsettings(){
// Open the  settings file
file_id = H5Fopen("settings.hdf", H5F_ACC_RDWR, H5P_DEFAULT);
if (file_id < 0){
	cout << "couldn't open file: settings.hdf" << endl;

	return -1;
}
else
{
// First read the topology
int nproc;
dataset_id = H5Dopen2(file_id, "/topology/Nprocs", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nproc);
status = H5Dclose(dataset_id);

dataset_id = H5Dopen2(file_id, "/topology/XLEN", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&XLEN);
status = H5Dclose(dataset_id);

dataset_id = H5Dopen2(file_id, "/topology/YLEN", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&YLEN);
status = H5Dclose(dataset_id);

dataset_id = H5Dopen2(file_id, "/topology/ZLEN", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ZLEN);
status = H5Dclose(dataset_id);

// read Lx
dataset_id = H5Dopen2(file_id, "/collective/Lx", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lx);
status = H5Dclose(dataset_id);
// read Ly
dataset_id = H5Dopen2(file_id, "/collective/Ly", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Ly);
status = H5Dclose(dataset_id);
// read Lz
dataset_id = H5Dopen2(file_id, "/collective/Lz", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lz);
status = H5Dclose(dataset_id);
// read nxc
dataset_id = H5Dopen2(file_id, "/collective/Nxc", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nxc);
status = H5Dclose(dataset_id);
// read nyc
dataset_id = H5Dopen2(file_id, "/collective/Nyc", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nyc);
status = H5Dclose(dataset_id);
// read nyc
dataset_id = H5Dopen2(file_id, "/collective/Nzc", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nzc);
status = H5Dclose(dataset_id);
// read ns
dataset_id = H5Dopen2(file_id, "/collective/Ns", H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ns);
// read qom
qom = new double[ns];
stringstream specie;
string temp;
for (int is=0; is<ns; is++){
specie.clear();
specie.str("");
specie << is;
temp = "/collective/species_"+specie.str()+"/qom";
dataset_id = H5Dopen2(file_id, temp.c_str(), H5P_DEFAULT);
status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&qom[is]);
}
// at this point you can close settings
status = H5Fclose(file_id);
return 0;
}
}


void readvect(string campo, string vectname, double** EX, double** EY,double** EZ) {

	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    for (int it=0; it < nlevels; it++){
    cycle = it * DeltaT + InitLevel;
	  cout << "READING VECTOR FROM HDF5 FILES  Time Level="<<cycle << endl;

    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	int j=jproc;
	int k=kproc;
    for (int i=0; i < XLEN;i++){
	    cout << "i="<<i << " j="<<j<< " k="<<k << endl;
        proc= i*YLEN*ZLEN+j*ZLEN+k;
		stringstream ss;
		ss << proc;
		cout << "ss="<<ss.str() << endl;
		temp = "proc" + ss.str() + ".hdf";
		proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		// read data
	    cout << "file = " << temp << endl;
		temp = campo+vectname+"x/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
		temp = campo+vectname+"y/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageY); status = H5Dclose(dataset_id);
		temp = campo+vectname+"z/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageZ); status = H5Dclose(dataset_id);
		int node=0;
		for (int ii=0; ii < (nxn+1);ii++)
			for (int jj=0; jj < (nyn+1);jj++)
				for (int kk=0; kk < (nzn+1);kk++){
					if (jj== 0 && ii!= nxn && kk== 0){
						EX[ii + nxn*i][it] = temp_storageX[node];
						EY[ii + nxn*i][it] = temp_storageY[node];
						EZ[ii + nxn*i][it] = temp_storageZ[node];
					}
					node++;
				}
		// close the file
		H5Fclose(proc_file_id);
		// go to other proc
	}
    }

}

void addreadvect(string campo, string vectname, double** EX, double** EY,double** EZ) {

	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    for (int it=0; it < nlevels; it++){
    cycle = it * DeltaT + InitLevel;
	  cout << "READING VECTOR FROM HDF5 FILES  Time Level="<<cycle << endl;

    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	int j=jproc;
	int k=kproc;
    for (int i=0; i < XLEN;i++){
	    cout << "i="<<i << " j="<<j<< " k="<<k << endl;
        proc= i*YLEN*ZLEN+j*ZLEN+k;
		stringstream ss;
		ss << proc;
		cout << "ss="<<ss.str() << endl;
		temp = "proc" + ss.str() + ".hdf";
		proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		// read data
	    cout << "file = " << temp << endl;
		temp = campo+vectname+"x/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
		temp = campo+vectname+"y/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageY); status = H5Dclose(dataset_id);
		temp = campo+vectname+"z/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageZ); status = H5Dclose(dataset_id);
		int node=0;
		for (int ii=0; ii < (nxn+1);ii++)
			for (int jj=0; jj < (nyn+1);jj++)
				for (int kk=0; kk < (nzn+1);kk++){
					if (jj== 0 && ii!= nxn && kk== 0){
						EX[ii + nxn*i][it] += temp_storageX[node];
						EY[ii + nxn*i][it] += temp_storageY[node];
						EZ[ii + nxn*i][it] += temp_storageZ[node];
					}
					node++;
				}
		// close the file
		H5Fclose(proc_file_id);
		// go to other proc
	}
    }

}


void readscalar(string campo, string scalarname, double** EX) {

	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    for (int it=0; it < nlevels; it++){
    cycle = it * DeltaT + InitLevel;
	  cout << "READING VECTOR FROM HDF5 FILES  Time Level="<<cycle << endl;

    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	int j=jproc;
	int k=kproc;
    for (int i=0; i < XLEN;i++){
	    cout << "i="<<i << " j="<<j<< " k="<<k << endl;
	    proc= i*YLEN*ZLEN+j*ZLEN+k;
		stringstream ss;
		ss << proc;
		cout << "ss="<<ss.str() << endl;
		temp = "proc" + ss.str() + ".hdf";
		proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		// read data
	    cout << "file = " << temp << endl;
		temp = campo+scalarname+"/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);

		int node=0;
		for (int ii=0; ii < (nxn+1);ii++)
			for (int jj=0; jj < (nyn+1);jj++)
				for (int kk=0; kk < (nzn+1);kk++){
					if (jj== 0 && ii!= nxn && kk== 0){
						EX[ii + nxn*i][it] = temp_storageX[node];
					}
					node++;
				}
		// close the file
		H5Fclose(proc_file_id);
		// go to other proc
	}
    }

}

void addreadscalar(string campo, string scalarname, double** EX) {

	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    for (int it=0; it < nlevels; it++){
    cycle = it * DeltaT + InitLevel;
	  cout << "READING VECTOR FROM HDF5 FILES  Time Level="<<cycle << endl;

    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	int j=jproc;
	int k=kproc;
    for (int i=0; i < XLEN;i++){
	    cout << "i="<<i << " j="<<j<< " k="<<k << endl;
	    proc= i*YLEN*ZLEN+j*ZLEN+k;
		stringstream ss;
		ss << proc;
		cout << "ss="<<ss.str() << endl;
		temp = "proc" + ss.str() + ".hdf";
		proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		// read data
	    cout << "file = " << temp << endl;
		temp = campo+scalarname+"/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);

		int node=0;
		for (int ii=0; ii < (nxn+1);ii++)
			for (int jj=0; jj < (nyn+1);jj++)
				for (int kk=0; kk < (nzn+1);kk++){
					if (jj== 0 && ii!= nxn && kk== 0){
						EX[ii + nxn*i][it] += temp_storageX[node];
					}
					node++;
				}
		// close the file
		H5Fclose(proc_file_id);
		// go to other proc
	}
    }

}

void readpotential(string campo, string scalarname, double** EX) {

	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    for (int it=0; it < nlevels; it++){
    cycle = it * DeltaT + InitLevel;
	  cout << "READING VECTOR FROM HDF5 FILES  Time Level="<<cycle << endl;

    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	int j=jproc;
	int k=kproc;
    for (int i=0; i < XLEN;i++){
	    cout << "i="<<i << " j="<<j<< " k="<<k << endl;
	    proc= i*YLEN*ZLEN+j*ZLEN+k;
		stringstream ss;
		ss << proc;
		cout << "ss="<<ss.str() << endl;
		temp = "proc" + ss.str() + ".hdf";
		proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		// read data
	    cout << "file = " << temp << endl;
		temp = campo+scalarname+"/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);

		int node=0;
		for (int ii=0; ii < (nxn);ii++)
			for (int jj=0; jj < (nyn);jj++)
				for (int kk=0; kk < (nzn);kk++){
					if (ii== 0 && kk== 0){
						EX[ii + nxn*i][it] = temp_storageX[node];
					}
					node++;
				}
		// close the file
		H5Fclose(proc_file_id);
		// go to other proc
	}
    }

}



void extract_pressure(double qom, double** BX, double** BY, double** BZ,
		double** VX, double** VY, double** VZ,
		double** N, double** pXX, double** pXY,
		double** pXZ, double** pYY, double** pYZ, double** pZZ,
		double** pPAR, double** pPER1, double** pPER2) {
for (int it=0; it < nlevels;it++)
    for (int jj=0; jj < nxn*XLEN;jj++){
    	if(N[jj][it]!=0.0)
    		pXX[jj][it] = pXX[jj][it]-VX[jj][it]*VX[jj][it]/N[jj][it];
	    pXX[jj][it] = pXX[jj][it] / qom;
    	if(N[jj][it]!=0.0)
	    pXY[jj][it] = pXY[jj][it]-VX[jj][it]*VY[jj][it]/N[jj][it];
	    pXY[jj][it] = pXY[jj][it] / qom;
    	if(N[jj][it]!=0.0)
	    pXZ[jj][it] = pXZ[jj][it]-VX[jj][it]*VZ[jj][it]/N[jj][it];
	    pXZ[jj][it] = pXZ[jj][it] / qom;
    	if(N[jj][it]!=0.0)
	    pYY[jj][it] = pYY[jj][it]-VY[jj][it]*VY[jj][it]/N[jj][it];
	    pYY[jj][it] = pYY[jj][it] / qom;
    	if(N[jj][it]!=0.0)
	    pYZ[jj][it] = pYZ[jj][it]-VY[jj][it]*VZ[jj][it]/N[jj][it];
	    pYZ[jj][it] = pYZ[jj][it] / qom;
    	if(N[jj][it]!=0.0)
	    pZZ[jj][it] = pZZ[jj][it]-VZ[jj][it]*VZ[jj][it]/N[jj][it];
	    pZZ[jj][it] = pZZ[jj][it] / qom;
	    
	double b2D = 1e-10 + BX[jj][it]*BX[jj][it] + BY[jj][it]*BY[jj][it];
	double b = b2D + BZ[jj][it]*BZ[jj][it];
	double perp2x = BZ[jj][it]*BX[jj][it] / sqrt(b*b2D);
	double perp2y = BZ[jj][it]*BY[jj][it] / sqrt(b*b2D);
	double perp2z = -sqrt(b2D/b);
				  
	pPAR[jj][it] = BX[jj][it]*pXX[jj][it]*BX[jj][it] + 
	                2*BX[jj][it]*pXY[jj][it]*BY[jj][it] + 
	                2*BX[jj][it]*pXZ[jj][it]*BZ[jj][it];
	pPAR[jj][it]+= BY[jj][it]*pYY[jj][it]*BY[jj][it] + 
					2*BY[jj][it]*pYZ[jj][it]*BZ[jj][it];
	pPAR[jj][it]+= BZ[jj][it]*pZZ[jj][it]*BZ[jj][it];
	
	pPAR[jj][it] = pPAR[jj][it]/b;
	
	pPER1[jj][it] = BY[jj][it]*pXX[jj][it]*BY[jj][it] - 
					2*BY[jj][it]*pXY[jj][it]*BX[jj][it];
	pPER1[jj][it]+= BX[jj][it]*pYY[jj][it]*BX[jj][it];

	pPER1[jj][it] = pPER1[jj][it]/b2D;

	pPER2[jj][it] = perp2x*pXX[jj][it]*perp2x + 2*perp2x*pXY[jj][it]*perp2y + 2*perp2x*pXZ[jj][it]*perp2z;
	pPER2[jj][it]+= perp2y*pYY[jj][it]*perp2y + 2*perp2y*pYZ[jj][it]*perp2z;
	pPER2[jj][it]+= perp2z*pZZ[jj][it]*perp2z;  

	    		    
	}
}


				
