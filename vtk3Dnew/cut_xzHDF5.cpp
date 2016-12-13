/***************************************************************************
 convHDF5.cpp  -  Convert program to open Parsek Output
 -------------------
 begin                : Jun 2008
 copyright            : (C) 2004 by Stefano Markidis, Giovanni Lapenta
 ************************************************** */

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

void readvect(int it, string campo, string vectname, double **EX, double **EY,double **EZ);

void readscalar(int it, string campo, string scalarname, double **EX);

void readpotential(int it, string campo, string scalarname, double **POT);

void writeVTKvect(int it, string vectname, string addname, double **EX, double **EY,double **EZ);

void writeVTKscalar(int it, string scalarname, string addname, double **f);

void writeVTKscalar(int it, string scalarname, string addname, double **f);

void writeVTKscalar_species(int it, string scalarname,  double **f1, double **f2);

void writeVTKscalar_species(int it, string scalarname,  double **f1, double **f2, double **f3);

void writeVTKtensor(int it, string tensorname, string addname, double** EXX, double** EXY,
		double** EXZ, double** EYY, double** EYZ, double** EZZ);

void extract_pressure(double qom, double** VX, double** VY, double** VZ,
		double** N, double** pXX, double** pXY,
		double** pXZ, double** pYY, double** pYZ, double** pZZ);

void cross(double** AX, double** AY, double** AZ,
		double** BX, double** BY, double** BZ,
		double** CX, double** CY, double** CZ);

int nxn, nyn, nzn;
int nxc, nyc, nzc;
int XLEN, YLEN, ZLEN;
int nproc;
int ns;
double* qom;

double Lx, Ly, Lz;

int InitLevel;
int MaxLevel;
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
	sscanf(argv[1],"%d",&InitLevel);
	sscanf(argv[2],"%d",&MaxLevel);
	sscanf(argv[3],"%d",&DeltaT);
	nlevels = MaxLevel/DeltaT;

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

	double** EX = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** EY = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** EZ = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** BX = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** BY = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** BZ = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** VX = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** VY = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** VZ = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** TXX = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** TXY = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** TXZ = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** TYY = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** TYZ = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** TZZ = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** NE = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** NI = newArr(double,nxn*XLEN,nzn*ZLEN);
	double** NH = newArr(double,nxn*XLEN,nzn*ZLEN);
	
	//Electric field
	for (int it=InitLevel/DeltaT; it < nlevels+1;it++){
    readvect(it, "/fields/","E", EX, EY, EZ);
    writeVTKvect(it, "E", "", EX, EY, EZ);

	//Magnetic field
    readvect(it,"/fields/","B", BX, BY, BZ);
    writeVTKvect(it,"B", "", BX, BY, BZ);

    cross(EX, EY, EZ, BX, BY, BZ, VX, VY, VZ);
    writeVTKvect(it,"S", "", VX, VY, VZ);

	//Rho by species
    readscalar(it,"/moments/species_0/","rho",  NE);
    readscalar(it,"/moments/species_1/","rho",  NI);
    readscalar(it,"/moments/species_2/","rho",  NH);
    writeVTKscalar_species(it,"rho", NE, NI, NH);
   // writeVTKscalar("rho", "e", VX);

	//Currents species0
    readvect(it,"/moments/species_0/","J",  VX, VY, VZ);
    writeVTKvect(it,"J", "e", VX, VY, VZ);
    //Pressure tensor species 0
    readscalar(it,"/moments/species_0/","pXX",  TXX);
    readscalar(it,"/moments/species_0/","pXY",  TXY);
    readscalar(it,"/moments/species_0/","pXZ",  TXZ);
    readscalar(it,"/moments/species_0/","pYY",  TYY);
    readscalar(it,"/moments/species_0/","pYZ",  TYZ);
    readscalar(it,"/moments/species_0/","pZZ",  TZZ);
    extract_pressure(qom[0], VX, VY, VZ, NE, TXX, TXY, TXZ, TYY, TYZ, TZZ);
    writeVTKtensor(it, "P", "e", TXX, TXY, TXZ, TYY, TYZ, TZZ);

    //Currents species1
    readvect(it,"/moments/species_1/","J",  VX, VY, VZ);
    writeVTKvect(it,"J", "i", VX, VY, VZ);
    //Pressure tensor species 1
    readscalar(it,"/moments/species_1/","pXX",  TXX);
    readscalar(it,"/moments/species_1/","pXY",  TXY);
    readscalar(it,"/moments/species_1/","pXZ",  TXZ);
    readscalar(it,"/moments/species_1/","pYY",  TYY);
    readscalar(it,"/moments/species_1/","pYZ",  TYZ);
    readscalar(it,"/moments/species_1/","pZZ",  TZZ);
    extract_pressure(qom[1], VX, VY, VZ, NE, TXX, TXY, TXZ, TYY, TYZ, TZZ);
    writeVTKtensor(it, "P", "i", TXX, TXY, TXZ, TYY, TYZ, TZZ);

    //Potential
    readpotential(it,"/potentials/","phi",  NI);
    writeVTKscalar(it,"phi","", NI);
}
	delArr(EX,nxn*XLEN);
	delArr(EY,nxn*XLEN);
	delArr(EZ,nxn*XLEN);
	delArr(BX,nxn*XLEN);
	delArr(BY,nxn*XLEN);
	delArr(BZ,nxn*XLEN);
	delArr(VX,nxn*XLEN);
	delArr(VY,nxn*XLEN);
	delArr(VZ,nxn*XLEN);
	delArr(TXX,nxn*XLEN);
	delArr(TXY,nxn*XLEN);
	delArr(TXZ,nxn*XLEN);
	delArr(TYY,nxn*XLEN);
	delArr(TYZ,nxn*XLEN);
	delArr(TZZ,nxn*XLEN);
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


void readvect(int it, string campo, string vectname, double** EX, double** EY,double** EZ) {

	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    cycle = it * DeltaT;
	  cout << "READING VECTOR FROM HDF5 FILES  Time Level="<<cycle << endl;

    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	int j=YLEN/2;
	for (int i=0; i < XLEN;i++)
    for (int k=0; k < ZLEN;k++){
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
					if (ii!=nxn && jj== 0 && kk!= nzn){
						EX[ii + nxn*i][kk + nzn*k] = temp_storageX[node];
						EY[ii + nxn*i][kk + nzn*k] = temp_storageY[node];
						EZ[ii + nxn*i][kk + nzn*k] = temp_storageZ[node];
					}
					node++;
				}
		// close the file
		H5Fclose(proc_file_id);
		// go to other proc
	}

}


void readscalar(int it, string campo, string scalarname, double** EX) {

	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    cycle = it * DeltaT;
	  cout << "READING VECTOR FROM HDF5 FILES  Time Level="<<cycle << endl;

    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	int j = YLEN/2;
	for (int i=0; i < XLEN;i++)
    for (int k=0; k < ZLEN;k++){
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
					if (ii!=nxn && jj== 0 && kk!= nzn){
						EX[ii + nxn*i][kk + nzn*k] = temp_storageX[node];
					}
					node++;
				}
		// close the file
		H5Fclose(proc_file_id);
		// go to other proc
    }
}

void readpotential(int it, string campo, string scalarname, double** EX) {

	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    cycle = it * DeltaT;
	  cout << "READING VECTOR FROM HDF5 FILES  Time Level="<<cycle << endl;

    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	int j = YLEN/2;
	for (int i=0; i < XLEN;i++)
    for (int k=0; k < ZLEN;k++){
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
					if (jj== 0) EX[ii + nxn*i][kk + nzn*k] = temp_storageX[node];
					node++;
				}
		// close the file
		H5Fclose(proc_file_id);
		// go to other proc
	}
}

void writeVTKvect(int it, string vectname, string addname, double** EX, double** EY,double** EZ) {
	double dx = Lx/nxc;
	double dy = Ly/nyc;
	double dz = Lz/nzc;

    int cycle = it * DeltaT;
	stringstream stringcycle;
	stringcycle << cycle;
	string temp;
	temp = vectname +addname +"_cut_xz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for" << vectname+addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << vectname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nzn*ZLEN << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << dy << " " << dy << " " << dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nzn*ZLEN << endl;
	my_fileE << "VECTORS " << vectname+addname << " float" << endl;
cout << "WRITING VECTOR " << vectname +addname<<" TO VTK FILE" << endl;
for (int kk=0; kk < nzn*ZLEN;kk++)
	   for (int ii=0; ii < nxn*XLEN;ii++){
		my_fileE << EX[ii][kk] << " " << EY[ii][kk] << " " << EZ[ii][kk] << endl;
			// my_file << temp_storageX[node] << " " << temp_storageY[node] << " " << temp_storageZ[node] << endl;
		}
my_fileE.close();
}


void writeVTKscalar(int it, string scalarname, string addname, double** EX) {
	double dx = Lx/nxc;
	double dy = Ly/nyc;
	double dz = Lz/nzc;

    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +addname +"_cut_xz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname+addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nzn*ZLEN << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << dy << " " << dy << " " << dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nzn*ZLEN << endl;
	my_fileE << "SCALARS " << scalarname+addname << " float" << endl;
	  my_fileE << "LOOKUP_TABLE default" << endl;
cout << "WRITING SCALAR " << scalarname +addname<<" TO VTK FILE" << endl;
for (int kk=0; kk < nzn*ZLEN;kk++)
	   for (int ii=0; ii < nxn*XLEN;ii++){
		my_fileE << EX[ii][kk] << endl;
		}
my_fileE.close();
}


void writeVTKscalar_species(int it, string scalarname, double** EX, double** EY) {
	double dx = Lx/nxc;
	double dy = Ly/nyc;
	double dz = Lz/nzc;

    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_cut_xz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nzn*ZLEN << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << dy << " " << dy << " " << dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nzn*ZLEN << endl;
	my_fileE << "SCALARS " << scalarname << "0 float" << endl;
	  my_fileE << "LOOKUP_TABLE default" << endl;
cout << "WRITING SCALAR " << scalarname <<"0 TO VTK FILE" << endl;
for (int kk=0; kk < nzn*ZLEN;kk++)
	   for (int ii=0; ii < nxn*XLEN;ii++){
		my_fileE << EX[ii][kk] << endl;
		}
my_fileE << "SCALARS " << scalarname << "1 float" << endl;
  my_fileE << "LOOKUP_TABLE default" << endl;
cout << "WRITING SCALAR " << scalarname <<"1 TO VTK FILE" << endl;
for (int kk=0; kk < nzn*ZLEN;kk++)
	   for (int ii=0; ii < nxn*XLEN;ii++){
	    my_fileE << EY[ii][kk] << endl;
	}
my_fileE.close();
}

void writeVTKscalar_species(int it, string scalarname, double** EX, double** EY, double** EZ) {
	double dx = Lx/nxc;
	double dy = Ly/nyc;
	double dz = Lz/nzc;

    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_cut_xz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nzn*ZLEN << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << dy << " " << dy << " " << dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nzn*ZLEN << endl;
	my_fileE << "SCALARS " << scalarname << "0 float" << endl;
	  my_fileE << "LOOKUP_TABLE default" << endl;
cout << "WRITING SCALAR " << scalarname <<"0 TO VTK FILE" << endl;
for (int kk=0; kk < nzn*ZLEN;kk++)
	   for (int ii=0; ii < nxn*XLEN;ii++){
		my_fileE << EX[ii][kk] << endl;
		}
my_fileE << "SCALARS " << scalarname << "1 float" << endl;
  my_fileE << "LOOKUP_TABLE default" << endl;
cout << "WRITING SCALAR " << scalarname <<"1 TO VTK FILE" << endl;
for (int kk=0; kk < nzn*ZLEN;kk++)
	   for (int ii=0; ii < nxn*XLEN;ii++){
	    my_fileE << EY[ii][kk] << endl;
	}
my_fileE << "SCALARS " << scalarname << "2 float" << endl;
  my_fileE << "LOOKUP_TABLE default" << endl;
cout << "WRITING SCALAR " << scalarname <<"2 TO VTK FILE" << endl;
for (int kk=0; kk < nzn*ZLEN;kk++)
	   for (int ii=0; ii < nxn*XLEN;ii++){
	    my_fileE << EZ[ii][kk] << endl;
	}
my_fileE.close();
}


void writeVTKtensor(int it, string tensorname, string addname, double** EXX, double** EXY,
		double** EXZ, double** EYY, double** EYZ, double** EZZ) {
	double dx = Lx/nxc;
	double dy = Ly/nyc;
	double dz = Lz/nzc;

    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = tensorname + addname +"_cut_xz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << tensorname +addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << tensorname +addname << " Tensor from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nzn*ZLEN << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << dy << " " << dy << " " << dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nzn*ZLEN << endl;

    cout << "WRITING Tensor " << tensorname +addname<<" TO VTK FILE" << endl;

	my_fileE << "SCALARS " << tensorname+addname << "xx float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		   for (int ii=0; ii < nxn*XLEN;ii++){
		    my_fileE << EXX[ii][kk] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "xy float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		   for (int ii=0; ii < nxn*XLEN;ii++){
		    my_fileE << EXY[ii][kk] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "xz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		   for (int ii=0; ii < nxn*XLEN;ii++){
		    my_fileE << EXZ[ii][kk] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "yy float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		   for (int ii=0; ii < nxn*XLEN;ii++){
		    my_fileE << EYY[ii][kk] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "yz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		   for (int ii=0; ii < nxn*XLEN;ii++){
		    my_fileE << EYZ[ii][kk] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "zz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		   for (int ii=0; ii < nxn*XLEN;ii++){
		    my_fileE << EZZ[ii][kk] << endl;
		}
my_fileE.close();
}

void extract_pressure(double qom, double** VX, double** VY, double** VZ,
		double** N, double** pXX, double** pXY,
		double** pXZ, double** pYY, double** pYZ, double** pZZ) {
	for (int kk=0; kk < nzn*ZLEN;kk++)
		   for (int ii=0; ii < nxn*XLEN;ii++){
    	if(N[ii][kk]!=0.0)
    		pXX[ii][kk] = pXX[ii][kk]-VX[ii][kk]*VX[ii][kk]/N[ii][kk];
	    pXX[ii][kk] = pXX[ii][kk] / qom;
    	if(N[ii][kk]!=0.0)
	    pXY[ii][kk] = pXY[ii][kk]-VX[ii][kk]*VY[ii][kk]/N[ii][kk];
	    pXY[ii][kk] = pXY[ii][kk] / qom;
    	if(N[ii][kk]!=0.0)
	    pXZ[ii][kk] = pXZ[ii][kk]-VX[ii][kk]*VZ[ii][kk]/N[ii][kk];
	    pXZ[ii][kk] = pXZ[ii][kk] / qom;
    	if(N[ii][kk]!=0.0)
	    pYY[ii][kk] = pYY[ii][kk]-VY[ii][kk]*VY[ii][kk]/N[ii][kk];
	    pYY[ii][kk] = pYY[ii][kk] / qom;
    	if(N[ii][kk]!=0.0)
	    pYZ[ii][kk] = pYZ[ii][kk]-VY[ii][kk]*VZ[ii][kk]/N[ii][kk];
	    pYZ[ii][kk] = pYZ[ii][kk] / qom;
    	if(N[ii][kk]!=0.0)
	    pZZ[ii][kk] = pZZ[ii][kk]-VZ[ii][kk]*VZ[ii][kk]/N[ii][kk];
	    pZZ[ii][kk] = pZZ[ii][kk] / qom;
	}
}


void cross(double** AX, double** AY, double** AZ,
		   double** BX, double** BY, double** BZ,
		   double** CX, double** CY, double** CZ) {
	for (int kk=0; kk < nzn*ZLEN;kk++)
		   for (int ii=0; ii < nxn*XLEN;ii++){
    		CX[ii][kk] = AY[ii][kk] * BZ[ii][kk] - AZ[ii][kk] * BY[ii][kk];
    	    CY[ii][kk] = AZ[ii][kk] * BX[ii][kk] - AX[ii][kk] * BZ[ii][kk];
    	    CZ[ii][kk] = AX[ii][kk] * BY[ii][kk] - AY[ii][kk] * BX[ii][kk];

	}
}


