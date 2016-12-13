/***************************************************************************
 convHDF5.cpp  -  Convert program to open Parsek Output
 -------------------
 begin                : Jun 2008
 copyright            : (C) 2004 by Stefano Markidis, Giovanni Lapenta
 ************************************************************************** */

#include "hdf5.h"
#include "Alloc.h"
#include "math.h"


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

void writeVTKvect(string vectname, string addname, double **EX, double **EY,double **EZ);

void writeVTKscalar(string scalarname, string addname, double **f);

void writeVTKscalar_species(string scalarname,  double **f1, double **f2);

void writeVTKtensor(string tensorname, string addname, double** pXX, double** pXY,
		double** pXZ, double** pYY, double** pYZ, double** pZZ,		
		double** pPAR, double** pPER1, double** pPER2);

void extract_pressure(double qom, 
		double** BX, double** BY, double** BZ,
		double** VX, double** VY, double** VZ,
		double** N, double** pXX, double** pXY,
		double** pXZ, double** pYY, double** pYZ, double** pZZ,
		double** pPAR, double** pPER1, double** pPER2);

int nxn, nyn, nzn;
int nxc, nyc, nzc;
int XLEN, YLEN, ZLEN;
// x position in processor topology. it will read the first cell in x in that processor
int iproc;
// z position in processor topology. it will read the first cell in x in that processor
int kproc;

int nproc;
int ns;
double* qom;

double Lx, Ly, Lz;

int MaxLevel;
int InitLevel;
int DeltaT;
int nlevels;

double *temp_storageX;
double *temp_storageY;
double *temp_storageZ;

// hdf stuff
hid_t    file_id;
hid_t    dataset_id;
herr_t   status;


int main (int argc, char **argv) {
	// cycle we want to open
	sscanf(argv[3],"%d",&MaxLevel);
	sscanf(argv[2],"%d",&DeltaT);
	sscanf(argv[1],"%d",&InitLevel);
	sscanf(argv[4],"%d",&iproc);
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

	double** BX = newArr2(double,nyn*YLEN,nlevels);
	double** BY = newArr2(double,nyn*YLEN,nlevels);
	double** BZ = newArr2(double,nyn*YLEN,nlevels);
	double** VX = newArr2(double,nyn*YLEN,nlevels);
	double** VY = newArr2(double,nyn*YLEN,nlevels);
	double** VZ = newArr2(double,nyn*YLEN,nlevels);
	double** TXX= newArr2(double,nyn*YLEN,nlevels);
	double** TXY = newArr2(double,nyn*YLEN,nlevels);
	double** TXZ = newArr2(double,nyn*YLEN,nlevels);
	double** TYY = newArr2(double,nyn*YLEN,nlevels);
	double** TYZ = newArr2(double,nyn*YLEN,nlevels);
	double** TZZ = newArr2(double,nyn*YLEN,nlevels);
	double** TPAR= newArr2(double,nyn*YLEN,nlevels);
	double** TPER1 = newArr2(double,nyn*YLEN,nlevels);
	double** TPER2 = newArr2(double,nyn*YLEN,nlevels);
	double** NE = newArr2(double,nyn*YLEN,nlevels);
	double** NI = newArr2(double,nyn*YLEN,nlevels);
	
	//Electric field
    readvect("/fields/","E", VX, VY, VZ);
    writeVTKvect("E", "", VX, VY, VZ);

	//Magnetic field
    readvect("/fields/","B", BX, BY, BZ);
    writeVTKvect("B", "", BX, BY, BZ);

	//Rho by species
    readscalar("/moments/species_0/","rho",  NE);
    readscalar("/moments/species_1/","rho",  NI);
    if (ns >2) addreadscalar("/moments/species_2/","rho",  NE);
    if (ns >3) addreadscalar("/moments/species_3/","rho",  NI);
    writeVTKscalar_species("rho", NE, NI);
   // writeVTKscalar("rho", "e", VX);

	//Currents species0
    readvect("/moments/species_0/","J",  VX, VY, VZ);
    if (ns >2) addreadvect("/moments/species_2/","J",  VX, VY, VZ);
    writeVTKvect("J", "e", VX, VY, VZ);
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
    writeVTKtensor("P", "e", TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2);

    //Currents species1
    readvect("/moments/species_1/","J",  VX, VY, VZ);
    if (ns >3) addreadvect("/moments/species_3/","J",  VX, VY, VZ);
    writeVTKvect("J", "i", VX, VY, VZ);
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
    writeVTKtensor("P", "i", TXX, TXY, TXZ, TYY, TYZ, TZZ, TPAR, TPER1, TPER2);
    

    //Potential
/*    readpotential("/potentials/","phi",  NI);
    writeVTKscalar("phi","", NI);
*/

	delArr2(BX,nyn*YLEN);
	delArr2(BY,nyn*YLEN);
	delArr2(BZ,nyn*YLEN);
	delArr2(VX,nyn*YLEN);
	delArr2(VY,nyn*YLEN);
	delArr2(VZ,nyn*YLEN);
	delArr2(TXX,nyn*YLEN);
	delArr2(TXY,nyn*YLEN);
	delArr2(TXZ,nyn*YLEN);
	delArr2(TYY,nyn*YLEN);
	delArr2(TYZ,nyn*YLEN);
	delArr2(TZZ,nyn*YLEN);
	delArr2(TPAR,nyn*YLEN);
	delArr2(TPER1,nyn*YLEN);
	delArr2(TPER2,nyn*YLEN);
	delArr2(NE,nyn*YLEN);
	delArr2(NI,nyn*YLEN);
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
dataset_id = H5Dopen(file_id, "/topology/Nprocs");
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nproc);
status = H5Dclose(dataset_id);

dataset_id = H5Dopen(file_id, "/topology/XLEN");
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&XLEN);
status = H5Dclose(dataset_id);

dataset_id = H5Dopen(file_id, "/topology/YLEN");
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&YLEN);
status = H5Dclose(dataset_id);

dataset_id = H5Dopen(file_id, "/topology/ZLEN");
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ZLEN);
status = H5Dclose(dataset_id);

// read Lx
dataset_id = H5Dopen(file_id, "/collective/Lx");
status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lx);
status = H5Dclose(dataset_id);
// read Ly
dataset_id = H5Dopen(file_id, "/collective/Ly");
status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Ly);
status = H5Dclose(dataset_id);
// read Lz
dataset_id = H5Dopen(file_id, "/collective/Lz");
status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lz);
status = H5Dclose(dataset_id);
// read nxc
dataset_id = H5Dopen(file_id, "/collective/Nxc");
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nxc);
status = H5Dclose(dataset_id);
// read nyc
dataset_id = H5Dopen(file_id, "/collective/Nyc");
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nyc);
status = H5Dclose(dataset_id);
// read nyc
dataset_id = H5Dopen(file_id, "/collective/Nzc");
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nzc);
status = H5Dclose(dataset_id);
// read ns
dataset_id = H5Dopen(file_id, "/collective/Ns");
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
dataset_id = H5Dopen(file_id, temp.c_str());
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
	int i=iproc;
	int k=kproc;
    for (int j=0; j < YLEN;j++){
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
		dataset_id = H5Dopen(proc_file_id,temp.c_str());
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
		temp = campo+vectname+"y/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen(proc_file_id,temp.c_str());
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageY); status = H5Dclose(dataset_id);
		temp = campo+vectname+"z/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen(proc_file_id,temp.c_str());
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageZ); status = H5Dclose(dataset_id);
		int node=0;
		for (int ii=0; ii < (nxn+1);ii++)
			for (int jj=0; jj < (nyn+1);jj++)
				for (int kk=0; kk < (nzn+1);kk++){
					if (ii== 0 && jj!= nyn && kk== 0){
						EX[jj + nyn*j][it] = temp_storageX[node];
						EY[jj + nyn*j][it] = temp_storageY[node];
						EZ[jj + nyn*j][it] = temp_storageZ[node];
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
	int i=iproc;
	int k=kproc;
    for (int j=0; j < YLEN;j++){
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
		dataset_id = H5Dopen(proc_file_id,temp.c_str());
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
		temp = campo+vectname+"y/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen(proc_file_id,temp.c_str());
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageY); status = H5Dclose(dataset_id);
		temp = campo+vectname+"z/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen(proc_file_id,temp.c_str());
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageZ); status = H5Dclose(dataset_id);
		int node=0;
		for (int ii=0; ii < (nxn+1);ii++)
			for (int jj=0; jj < (nyn+1);jj++)
				for (int kk=0; kk < (nzn+1);kk++){
					if (ii== 0 && jj!= nyn && kk== 0){
						EX[jj + nyn*j][it] += temp_storageX[node];
						EY[jj + nyn*j][it] += temp_storageY[node];
						EZ[jj + nyn*j][it] += temp_storageZ[node];
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
	int i=iproc;
	int k=kproc;
    for (int j=0; j < YLEN;j++){
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
		dataset_id = H5Dopen(proc_file_id,temp.c_str());
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);

		int node=0;
		for (int ii=0; ii < (nxn+1);ii++)
			for (int jj=0; jj < (nyn+1);jj++)
				for (int kk=0; kk < (nzn+1);kk++){
					if (ii== 0 && jj!= nyn && kk== 0){
						EX[jj + nyn*j][it] = temp_storageX[node];
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
	int i=iproc;
	int k=kproc;
    for (int j=0; j < YLEN;j++){
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
		dataset_id = H5Dopen(proc_file_id,temp.c_str());
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);

		int node=0;
		for (int ii=0; ii < (nxn+1);ii++)
			for (int jj=0; jj < (nyn+1);jj++)
				for (int kk=0; kk < (nzn+1);kk++){
					if (ii== 0 && jj!= nyn && kk== 0){
						EX[jj + nyn*j][it] += temp_storageX[node];
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
	int i=iproc;
	int k=kproc;
    for (int j=0; j < YLEN;j++){
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
		dataset_id = H5Dopen(proc_file_id,temp.c_str());
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);

		int node=0;
		for (int ii=0; ii < (nxn);ii++)
			for (int jj=0; jj < (nyn);jj++)
				for (int kk=0; kk < (nzn);kk++){
					if (ii== 0 && kk== 0){
						EX[jj + nyn*j][it] = temp_storageX[node];
					}
					node++;
				}
		// close the file
		H5Fclose(proc_file_id);
		// go to other proc
	}
    }

}


void writeVTKvect(string vectname, string addname, double** EX, double** EY,double** EZ) {
	double dx = Lx/nxc;
	double dy = Ly/nyc;
	double dz = Lz/nzc;

	string temp;
	temp = vectname +addname +"_yt";
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for" << vectname+addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << vectname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nyn*YLEN << " " << nlevels << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << dy << " " << dy*nyn*YLEN/nlevels << " " << 1.0 << endl;
	my_fileE << "POINT_DATA " << nyn*YLEN*nlevels << endl;
	my_fileE << "VECTORS " << vectname+addname << " float" << endl;
cout << "WRITING VECTOR " << vectname +addname<<" TO VTK FILE" << endl;
for (int it=0; it < nlevels;it++)
	for (int jj=0; jj < nyn*YLEN;jj++){
		my_fileE << EX[jj][it] << " " << EY[jj][it] << " " << EZ[jj][it] << endl;
			// my_file << temp_storageX[node] << " " << temp_storageY[node] << " " << temp_storageZ[node] << endl;
		}
my_fileE.close();
}



void writeVTKscalar(string scalarname, string addname, double** EX) {
	double dx = Lx/nxc;
	double dy = Ly/nyc;
	double dz = Lz/nzc;

	string temp;
	temp = scalarname +addname +"_yt";
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname+addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nyn*YLEN << " " << nlevels << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << dy << " " << dy*nyn*YLEN/nlevels << " " << dz << endl;
	my_fileE << "POINT_DATA " << nyn*YLEN*nlevels << endl;
	my_fileE << "SCALARS " << scalarname+addname << " float" << endl;
	  my_fileE << "LOOKUP_TABLE default" << endl;
cout << "WRITING SCALAR " << scalarname +addname<<" TO VTK FILE" << endl;
for (int it=0; it < nlevels;it++)
	for (int jj=0; jj < nyn*YLEN;jj++){
		my_fileE << EX[jj][it] << endl;
		}
my_fileE.close();
}


void writeVTKscalar_species(string scalarname, double** EX, double** EY) {
	double dx = Lx/nxc;
	double dy = Ly/nyc;
	double dz = Lz/nzc;

	string temp;
	temp = scalarname +"_yt";
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nyn*YLEN << " " << nlevels << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << dy << " " << dy*nyn*YLEN/nlevels << " " << 1.0 << endl;
	my_fileE << "POINT_DATA " << nyn*YLEN*nlevels << endl;
	my_fileE << "SCALARS " << scalarname << "0 float" << endl;
	  my_fileE << "LOOKUP_TABLE default" << endl;
cout << "WRITING SCALAR " << scalarname <<"0 TO VTK FILE" << endl;
for (int it=0; it < nlevels;it++)
	for (int jj=0; jj < nyn*YLEN;jj++){
		my_fileE << EX[jj][it] << endl;
		}
my_fileE << "SCALARS " << scalarname << "1 float" << endl;
  my_fileE << "LOOKUP_TABLE default" << endl;
cout << "WRITING SCALAR " << scalarname <<"1 TO VTK FILE" << endl;
for (int it=0; it < nlevels;it++)
	for (int jj=0; jj < nyn*YLEN;jj++){
	    my_fileE << EY[jj][it] << endl;
	}
my_fileE.close();
}

void writeVTKtensor(string tensorname, string addname, double** EXX, double** EXY,
		double** EXZ, double** EYY, double** EYZ, double** EZZ,
		double** EPAR, double** EPER1, double** EPER2) {
	double dx = Lx/nxc;
	double dy = Ly/nyc;
	double dz = Lz/nzc;

	string temp;
	temp = tensorname + addname +"_yt";
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << tensorname +addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << tensorname +addname << " Tensor from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nyn*YLEN << " " << nlevels << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << dy << " " << dy*nyn*YLEN/nlevels << " " << 1.0 << endl;
	my_fileE << "POINT_DATA " << nyn*YLEN*nlevels << endl;

    cout << "WRITING Tensor " << tensorname +addname<<" TO VTK FILE" << endl;

	my_fileE << "SCALARS " << tensorname+addname << "xx float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
    for (int it=0; it < nlevels;it++)
	    for (int jj=0; jj < nyn*YLEN;jj++){
		    my_fileE << EXX[jj][it] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "xy float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
    for (int it=0; it < nlevels;it++)
	    for (int jj=0; jj < nyn*YLEN;jj++){
		    my_fileE << EXY[jj][it] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "xz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
    for (int it=0; it < nlevels;it++)
	    for (int jj=0; jj < nyn*YLEN;jj++){
		    my_fileE << EXZ[jj][it] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "yy float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
    for (int it=0; it < nlevels;it++)
	    for (int jj=0; jj < nyn*YLEN;jj++){
		    my_fileE << EYY[jj][it] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "yz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
    for (int it=0; it < nlevels;it++)
	    for (int jj=0; jj < nyn*YLEN;jj++){
		    my_fileE << EYZ[jj][it] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "zz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
    for (int it=0; it < nlevels;it++)
	    for (int jj=0; jj < nyn*YLEN;jj++){
		    my_fileE << EZZ[jj][it] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "par float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
    for (int it=0; it < nlevels;it++)
	    for (int jj=0; jj < nyn*YLEN;jj++){
			my_fileE << EPAR[jj][it] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "perp1 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
    for (int it=0; it < nlevels;it++)
	    for (int jj=0; jj < nyn*YLEN;jj++){
		    my_fileE << EPER1[jj][it] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "perp2 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
    for (int it=0; it < nlevels;it++)
	    for (int jj=0; jj < nyn*YLEN;jj++){
		    my_fileE << EPER2[jj][it] << endl;
		}
my_fileE.close();
}

void extract_pressure(double qom, double** BX, double** BY, double** BZ,
		double** VX, double** VY, double** VZ,
		double** N, double** pXX, double** pXY,
		double** pXZ, double** pYY, double** pYZ, double** pZZ,
		double** pPAR, double** pPER1, double** pPER2) {
for (int it=0; it < nlevels;it++)
    for (int jj=0; jj < nyn*YLEN;jj++){
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


				