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

void readpart(int is, string varname, double *VX);

void writeTXTpart(int is, double *Q, double *PX, double *PY, double *PZ, double *VX, double *VY, double *VZ);

int nxn, nyn, nzn;
int nxc, nyc, nzc;
int XLEN, YLEN, ZLEN;
int nproc;
int ns;
long* Np;
long* Nptot;
double* qom;

double Lx, Ly, Lz;

int MaxLevel;
int DeltaT;
int nlevels;

double *temp_storageX;


// hdf stuff
hid_t    file_id;
hid_t    dataset_id;
herr_t   status;


int main (int argc, char **argv) {
	// cycle we want to open

    int out;
	out = readsettings();
	if(out<0)
		return -1;

	nxn = nxc/XLEN;
	nyn = nyc/YLEN;
	nzn = nzc/ZLEN;


    for(int is=0; is < ns; is++){

	temp_storageX = new double[Np[is]];

	double* Q = new double[Np[is]];
	double* PX0 = new double[Np[is]];
	double* PY0 = new double[Np[is]];
	double* PZ0 = new double[Np[is]];
	double* VX0 = new double[Np[is]];
	double* VY0 = new double[Np[is]];
	double* VZ0 = new double[Np[is]];
	//Particles

    readpart(is, "x", PX0);
    readpart(is, "y", PY0);
    readpart(is, "z", PZ0);
    readpart(is, "u", VX0);
    readpart(is, "v", VY0);
    readpart(is, "w", VZ0);
    readpart(is, "q", Q);
    writeTXTpart(is, Q, PX0, PY0, PZ0, VX0, VY0, VZ0);

    delete[] PX0;
    delete[] PY0;
    delete[] PZ0;
    delete[] VX0;
    delete[] VY0;
    delete[] VZ0;
	delete[] temp_storageX;

}
	
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
Np = new long[ns];
Nptot = new long[ns];
stringstream specie;
string temp;
for (int is=0; is<ns; is++){
specie.clear();
specie.str("");
specie << is;
temp = "/collective/species_"+specie.str()+"/qom";
dataset_id = H5Dopen(file_id, temp.c_str());
status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&qom[is]);
temp = "/collective/species_"+specie.str()+"/Np";
dataset_id = H5Dopen(file_id, temp.c_str());
status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Np[is]);
cout << "Number of Particles" <<Np[0]<<"  "<<Np[1]<< endl;
}
// at this point you can close settings
status = H5Fclose(file_id);
return 0;
}
}

void readpart(int is, string varname,double* VX) {

	hid_t proc_file_id;
    hid_t       filespace;
    hsize_t     dims[2];
    int rank;
	stringstream cc;
    int cycle = 0;
	  cout << "READING PARTICLES FROM HDF5 FILES  " << endl;

    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	stringstream specie;
	specie.clear();
	specie.str("");
	specie << is;
	int proc;
	int part=0;
	for (int i=0; i < XLEN;i++)
    for (int j=0; j < YLEN;j++)
	for (int k=0; k < ZLEN;k++)
	{
	    cout << "i="<<i << " j="<<j<< " k="<<k << endl;
        proc= i*YLEN*ZLEN+j*ZLEN+k;
		stringstream ss;
		ss << proc;
		cout << "ss="<<ss.str() << endl;
		temp = "restart" + ss.str() + ".hdf";
		proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		// read data
	    cout << "file = " << temp << endl;
		temp = "/particles/species_"+specie.str()+"/"+varname+"/cycle_"+ cc.str();
	    cout << "dataset = " << temp << endl;
		dataset_id = H5Dopen(proc_file_id,temp.c_str());
		filespace = H5Dget_space(dataset_id);    /* Get filespace handle first. */
		rank      = H5Sget_simple_extent_ndims(filespace);
		status    = H5Sget_simple_extent_dims(filespace, dims, NULL);
        cout << "dimensions " << dims[0]<<endl;
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
		for (int ip=0; ip < dims[0];ip++){
			VX[part] = temp_storageX[ip];
					part++;
				}
		// close the file
		H5Fclose(proc_file_id);
		Nptot[is]=part;
		// go to other proc
	}
}


void writeTXTpart(int is, double* Q, double* PX, double* PY, double* PZ, double* VX, double* VY, double* VZ) {


	stringstream stringspec;
	stringspec << is;
	string temp;
	temp = "species" +stringspec.str();
	temp += ".txt";
	cout << "Writing file: " << temp << endl;
	ofstream my_file(temp.c_str());
	cout << "writing particles NpToT="<<Nptot[is]  << endl;
for (int part=0; part < Nptot[is]; part++){
		my_file << Q[part]  << " "<< PX[part] << " "<< PY[part]<< " " << PZ[part] << " "<< VX[part]  << " " << VY[part]<< " " << VZ[part] <<endl;
					}
my_file.close();
}

