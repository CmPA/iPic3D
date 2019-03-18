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


int readsettings();

void readvect(int it, string campo, string vectname, double ***EX, double ***EY,double ***EZ);

void addreadvect(int it, string campo, string vectname, double ***EX, double ***EY,double ***EZ);

void readscalar(int it, string campo, string scalarname, double ***EX);

void addreadscalar(int it, string campo, string scalarname, double ***EX);

void readpotential(int it, string campo, string scalarname, double ***POT);

// Write 3D VTK Files

void writeVTKvect(int it, string vectname, string addname, double ***EX, double ***EY,double ***EZ);

void writeVTKvect_binary(int it, string vectname, string addname, double ***EX, double ***EY,double ***EZ);

void writeVTKscalar(int it, string scalarname, string addname, double ***f);

void writeVTKscalar_binary(int it, string scalarname, string addname, double ***f);

void writeVTKscalar_species(int it, string scalarname,  double ***f1, double ***f2);

void writeVTKscalar_species_binary(int it, string scalarname,  double ***f1, double ***f2);

void writeVTKscalar_species(int it, string scalarname,  double ***f1, double ***f2, double ***f3);

void writeVTKscalar_species_binary(int it, string scalarname,  double ***f1, double ***f2, double ***f3);

void writeVTKtensor(int it, string tensorname, string addname, double*** EXX, double*** EXY,
					double*** EXZ, double*** EYY, double*** EYZ, double*** EZZ, 
					double*** pPAR, double*** pPER1, double*** pPER2, double*** EPS);

void writeVTKtensor_binary(int it, string tensorname, string addname, double*** EXX, double*** EXY,
					double*** EXZ, double*** EYY, double*** EYZ, double*** EZZ,
					double*** pPAR, double*** pPER1, double*** pPER2, double*** EPS);

void readVTKvect(string filename, double ***BX, double ***BY,double ***BZ);

void readVTKpreamble(string filename);


// Write 2D VTK Files

void writeVTKvect(int it, string vectname, string addname, double **EX, double **EY,double **EZ);

void writeVTKscalar(int it, string scalarname, string addname, double **f);

void writeVTKscalar(int it, string scalarname, string addname, double **f);

void writeVTKscalar_species(int it, string scalarname,  double **f1, double **f2);

void writeVTKscalar_species(int it, string scalarname,  double **f1, double **f2, double **f3);

void writeVTKtensor(int it, string tensorname, string addname, double** EXX, double** EXY,
					double** EXZ, double** EYY, double** EYZ, double** EZZ, 
					double** pPAR, double** pPER1, double** pPER2, double** EPS);
// Function needed fro switching from little-Endian to big-Endian needed by paraview
float ReverseFloat( const float inFloat );

//function to rotate indexes
int irot(int i);

// Various useful variables about the dataset

int MaxLevel;
int DeltaT;
int InitT;
int Species;
int nlevels;
int initlevel;
int NdimCode = 3;
int xshift = 0;
int yshift = 0;
int zshift = 0;

double *temp_storageX;
double *temp_storageY;
double *temp_storageZ;

// hdf stuff
hid_t    file_id;
hid_t    dataset_id;
herr_t   status;

string grid_str;

int readsettings(){
	string temp;
	// Open the  settings file
	temp = "settings" + grid_str  + ".hdf";
	file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
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
		
		if(NdimCode==2){
			ZLEN=1;
		}
		else
		{
		dataset_id = H5Dopen2(file_id, "/topology/ZLEN", H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ZLEN);
		status = H5Dclose(dataset_id);
		}
		
		// read Lx
		dataset_id = H5Dopen2(file_id, "/collective/Lx", H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lx);
		status = H5Dclose(dataset_id);
		// read Ly
		dataset_id = H5Dopen2(file_id, "/collective/Ly", H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Ly);
		status = H5Dclose(dataset_id);
		// read Lz
		if(NdimCode==2){
			Lz=0;
		}
		else
		{
		dataset_id = H5Dopen2(file_id, "/collective/Lz", H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lz);
		status = H5Dclose(dataset_id);
		}
		// read Dx
		dataset_id = H5Dopen2(file_id, "/collective/Dx", H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Dx);
		status = H5Dclose(dataset_id);
		// read Dy
		dataset_id = H5Dopen2(file_id, "/collective/Dy", H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Dy);
		status = H5Dclose(dataset_id);
		// read Dz
		if(NdimCode == 2)
		{
			Dz=Dy;
		}
		else
		{
		dataset_id = H5Dopen2(file_id, "/collective/Dz", H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Dz);
		status = H5Dclose(dataset_id);
		}
		
		// read nxc
		dataset_id = H5Dopen2(file_id, "/collective/Nxc", H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nxc);
		status = H5Dclose(dataset_id);
		// read nyc
		dataset_id = H5Dopen2(file_id, "/collective/Nyc", H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nyc);
		status = H5Dclose(dataset_id);
		// read nyc
		if(NdimCode == 2)
		{
			nzc=1;
		}
		else
		{
		dataset_id = H5Dopen2(file_id, "/collective/Nzc", H5P_DEFAULT);
		status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nzc);
		status = H5Dclose(dataset_id);
		}

		nxn = nxc/XLEN;
		nyn = nyc/YLEN;
		nzn = nzc/ZLEN;
		nodes_z=nzn+1;
		if(NdimCode==2) nodes_z=nzn;

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


void readvect(int it, string campo, string vectname, double*** EX, double*** EY,double*** EZ) {
	
	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    cycle = it * DeltaT;
	cout << "READING VECTOR " << vectname <<" FROM HDF5 FILES  Time Level="<<cycle << endl;
	
    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	for (int i=0; i < XLEN;i++)
		for (int j=0; j < YLEN;j++)
			for (int k=0; k < ZLEN;k++){
				//cout << "i="<<i << " j="<<j<< " k="<<k << endl;
				proc= i*YLEN*ZLEN+j*ZLEN+k;
				stringstream ss;
				ss << proc;
				//cout << "ss="<<ss.str() + grid_str  << endl;
				temp = "proc" + ss.str() + grid_str  + ".hdf";
				proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
				// read data
				//cout << "file = " << temp << endl;
				temp = campo+vectname+"x/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
				temp = campo+vectname+"y/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageY); status = H5Dclose(dataset_id);
				temp = campo+vectname+"z/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageZ); status = H5Dclose(dataset_id);
				int node=0;
				for (int ii=0; ii < (nxn+1);ii++)
					for (int jj=0; jj < (nyn+1);jj++)
						for (int kk=0; kk < nodes_z;kk++){
							if (ii!=nxn && jj!= nyn && kk!= nzn){
								EX[irot(ii + nxn*i)][jj + nyn*j][kk + nzn*k] = temp_storageX[node];
								EY[irot(ii + nxn*i)][jj + nyn*j][kk + nzn*k] = temp_storageY[node];
								EZ[irot(ii + nxn*i)][jj + nyn*j][kk + nzn*k] = temp_storageZ[node];
							}
							node++;
						}
				// close the file
				H5Fclose(proc_file_id);
				// go to other proc
			}
	
}

void addreadvect(int it, string campo, string vectname, double*** EX, double*** EY,double*** EZ) {
	
	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    cycle = it * DeltaT;
	cout << "READING VECTOR "<< vectname<<" FROM HDF5 FILES  Time Level="<<cycle << endl;
	
    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	for (int i=0; i < XLEN;i++)
		for (int j=0; j < YLEN;j++)
			for (int k=0; k < ZLEN;k++){
				//cout << "i="<<i << " j="<<j<< " k="<<k << endl;
				proc= i*YLEN*ZLEN+j*ZLEN+k;
				stringstream ss;
				ss << proc;
				//cout << "ss="<<ss.str() + grid_str  << endl;
				temp = "proc" + ss.str() + grid_str  + ".hdf";
				proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
				// read data
				//cout << "file = " << temp << endl;
				temp = campo+vectname+"x/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
				temp = campo+vectname+"y/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageY); status = H5Dclose(dataset_id);
				temp = campo+vectname+"z/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageZ); status = H5Dclose(dataset_id);
				int node=0;
				for (int ii=0; ii < (nxn+1);ii++)
					for (int jj=0; jj < (nyn+1);jj++)
						for (int kk=0; kk < nodes_z;kk++){
							if (ii!=nxn && jj!= nyn && kk!= nzn){
								EX[irot(ii + nxn*i)][jj + nyn*j][kk + nzn*k] += temp_storageX[node];
								EY[irot(ii + nxn*i)][jj + nyn*j][kk + nzn*k] += temp_storageY[node];
								EZ[irot(ii + nxn*i)][jj + nyn*j][kk + nzn*k] += temp_storageZ[node];
							}
							node++;
						}
				// close the file
				H5Fclose(proc_file_id);
				// go to other proc
			}
	
}

void readscalar(int it, string campo, string scalarname, double*** EX) {
	
	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    cycle = it * DeltaT;
	cout << "READING SCALAR " << scalarname <<" FROM HDF5 FILES  Time Level="<<cycle << endl;
	
    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	int i=XLEN/2;
	int k=ZLEN/2;
	for (int i=0; i < XLEN;i++)
		for (int j=0; j < YLEN;j++)
			for (int k=0; k < ZLEN;k++){
				//cout << "i="<<i << " j="<<j<< " k="<<k << endl;
				proc= i*YLEN*ZLEN+j*ZLEN+k;
				stringstream ss;
				ss << proc;
				//cout << "ss="<<ss.str() + grid_str  << endl;
				temp = "proc" + ss.str() + grid_str  + ".hdf";
				proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
				// read data
				//cout << "file = " << temp << endl;
				temp = campo+scalarname+"/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
				
				int node=0;
				for (int ii=0; ii < (nxn+1);ii++)
					for (int jj=0; jj < (nyn+1);jj++)
						for (int kk=0; kk < nodes_z;kk++){
							if (ii!=nxn && jj!= nyn && kk!= nzn){
								EX[irot(ii + nxn*i)][jj + nyn*j][kk + nzn*k] = temp_storageX[node];
							}
							node++;
						}
				// close the file
				H5Fclose(proc_file_id);
				// go to other proc
			}
}
void addreadscalar(int it, string campo, string scalarname, double*** EX) {
	
	hid_t proc_file_id;
	int cycle;
	stringstream cc;
    cycle = it * DeltaT;
	cout << "READING SCALAR " << scalarname <<" FROM HDF5 FILES  Time Level="<<cycle << endl;
	
    cc.clear();
    cc.str("");
	cc << cycle;
	string temp;
	int proc;
	int i=XLEN/2;
	int k=ZLEN/2;
	for (int i=0; i < XLEN;i++)
		for (int j=0; j < YLEN;j++)
			for (int k=0; k < ZLEN;k++){
				//cout << "i="<<i << " j="<<j<< " k="<<k << endl;
				proc= i*YLEN*ZLEN+j*ZLEN+k;
				stringstream ss;
				ss << proc;
				//cout << "ss="<<ss.str() + grid_str  << endl;
				temp = "proc" + ss.str() + grid_str  + ".hdf";
				proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
				// read data
				//cout << "file = " << temp << endl;
				temp = campo+scalarname+"/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
				
				int node=0;
				for (int ii=0; ii < (nxn+1);ii++)
					for (int jj=0; jj < (nyn+1);jj++)
						for (int kk=0; kk < nodes_z;kk++){
							if (ii!=nxn && jj!= nyn && kk!= nzn){
								EX[irot(ii + nxn*i)][jj + nyn*j][kk + nzn*k] += temp_storageX[node];
							}
							node++;
						}
				// close the file
				H5Fclose(proc_file_id);
				// go to other proc
			}
}

void readpotential(int it, string campo, string scalarname, double*** EX) {
	
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
	for (int i=0; i < XLEN;i++)
		for (int j=0; j < YLEN;j++)
			for (int k=0; k < ZLEN;k++){
				//cout << "i="<<i << " j="<<j<< " k="<<k << endl;
				proc= i*YLEN*ZLEN+j*ZLEN+k;
				stringstream ss;
				ss << proc;
				//cout << "ss="<<ss.str() + grid_str  << endl;
				temp = "proc" + ss.str() + grid_str  + ".hdf";
				proc_file_id = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
				// read data
				//cout << "file = " << temp << endl;
				temp = campo+scalarname+"/cycle_"+ cc.str();
				//cout << "dataset = " << temp << endl;
				dataset_id = H5Dopen2(proc_file_id,temp.c_str(), H5P_DEFAULT);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storageX); status = H5Dclose(dataset_id);
				
				int node=0;
				for (int ii=0; ii < (nxn);ii++)
					for (int jj=0; jj < (nyn);jj++)
						for (int kk=0; kk < (nzn);kk++){
							EX[irot(ii + nxn*i)][jj + nyn*j][kk + nzn*k] = temp_storageX[node];
							node++;
						}
				// close the file
				H5Fclose(proc_file_id);
				// go to other proc
			}
}

void writeVTKvect(int it, string vectname, string addname, double*** EX, double*** EY,double*** EZ) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;
	
    int cycle = it * DeltaT;
	stringstream stringcycle;
	stringcycle << cycle;
	string temp;
	temp = vectname +addname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for" << vectname+addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << vectname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_fileE << "VECTORS " << vectname+addname << " float" << endl;
	cout << "WRITING VECTOR " << vectname +addname<<" TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
 				float fbx = EX[ii][jj][kk];
 				float fby = EY[ii][jj][kk];
 				float fbz = EZ[ii][jj][kk];
//if(abs(fbx)<1e-10) fbx=0.0;
//if(abs(fby)<1e-10) fby=0.0;
//if(abs(fbz)<1e-10) fbz=0.0;
//				my_fileE << EX[ii][jj][kk] << " " << EY[ii][jj][kk] << " " << EZ[ii][jj][kk] << endl;
				my_fileE << fbx << " " << fby << " " << fbz << endl;
				// my_file << temp_storageX[node] << " " << temp_storageY[node] << " " << temp_storageZ[node] << endl;
			}
	my_fileE.close();
}


void writeVTKvect_binary(int it, string vectname, string addname, double*** EX, double*** EY,double*** EZ) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;

    int cycle = it * DeltaT;
	stringstream stringcycle;
	stringcycle << cycle;
	string temp;
	temp = vectname +addname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for" << vectname+addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << vectname << " Field from Parsek" << endl;
	my_fileE << "BINARY" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_fileE << "VECTORS " << vectname+addname << " float" << endl;
	cout << "WRITING VECTOR " << vectname +addname<<" TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
 				float fbx = ReverseFloat(EX[ii][jj][kk]);
 				float fby = ReverseFloat(EY[ii][jj][kk]);
 				float fbz = ReverseFloat(EZ[ii][jj][kk]);
 				my_fileE.write((char*)&fbx,sizeof(float));
 				my_fileE.write((char*)&fby,sizeof(float));
 				my_fileE.write((char*)&fbz,sizeof(float));
			}
	my_fileE.close();
}


void writeVTKscalar(int it, string scalarname, string addname, double*** EX) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;
	
    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +addname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname+addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_fileE << "SCALARS " << scalarname+addname << " float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname +addname<<" TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EX[ii][jj][kk] << endl;
			}
	my_fileE.close();
}

void writeVTKscalar_binary(int it, string scalarname, string addname, double*** EX) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;

    int cycle = it * DeltaT;
    //float fl[nxn*XLEN*nyn*YLEN*nzn*ZLEN];
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +addname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str(),std::ofstream::binary);
	cout << "writing to file mesh points for " << scalarname+addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "BINARY" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_fileE << "SCALARS " << scalarname+addname << " float" << endl;
	my_fileE << "LOOKUP_TABLE default" <<"\n"  ;
	cout << "WRITING SCALAR " << scalarname +addname<<" TO VTK FILE" << endl;
	//int i=0;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EX[ii][jj][kk]);
				//i++;
				//my_fileE << fl << endl;
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	//my_fileE.write((char*)&EX[0][0][0],nxn*XLEN*nyn*YLEN*nzn*ZLEN);
	my_fileE.close();
}

float ReverseFloat( const float inFloat )
{
   float retVal;
   char *floatToConvert = ( char* ) & inFloat;
   char *returnFloat = ( char* ) & retVal;

   // swap the bytes into a temporary buffer
   returnFloat[0] = floatToConvert[3];
   returnFloat[1] = floatToConvert[2];
   returnFloat[2] = floatToConvert[1];
   returnFloat[3] = floatToConvert[0];

   return retVal;
}

void writeVTKscalar_species(int it, string scalarname, double*** EX, double*** EY) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;
	
    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_fileE << "SCALARS " << scalarname << "0 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"0 TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EX[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << scalarname << "1 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"1 TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EY[ii][jj][kk] << endl;
			}
	my_fileE.close();
}


void writeVTKscalar_species_binary(int it, string scalarname, double*** EX, double*** EY) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;

    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "BINARY" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_fileE << "SCALARS " << scalarname << "0 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"0 TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EX[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << scalarname << "1 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"1 TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EY[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE.close();
}


void writeVTKscalar_species(int it, string scalarname, double*** EX, double*** EY, double*** EZ) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;
	
    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_fileE << "SCALARS " << scalarname << "0 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"0 TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EX[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << scalarname << "1 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"1 TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EY[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << scalarname << "2 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"2 TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EZ[ii][jj][kk] << endl;
			}
	my_fileE.close();
}


void writeVTKscalar_species_binary(int it, string scalarname, double*** EX, double*** EY, double*** EZ) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;

    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "BINARY" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	my_fileE << "SCALARS " << scalarname << "0 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"0 TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EX[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << scalarname << "1 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"1 TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EY[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << scalarname << "2 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"2 TO VTK FILE" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EZ[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE.close();
}



void writeVTKtensor(int it, string tensorname, string addname, double*** EXX, double*** EXY,
					double*** EXZ, double*** EYY, double*** EYZ, double*** EZZ,
					double*** pPAR, double*** pPER1, double*** pPER2, double*** EPS) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;
	
    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = tensorname + addname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << tensorname +addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << tensorname +addname << " Tensor from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;
	
    cout << "WRITING Tensor " << tensorname +addname<<" TO VTK FILE" << endl;
	
	my_fileE << "SCALARS " << tensorname+addname << "xx float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EXX[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << tensorname+addname << "xy float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EXY[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << tensorname+addname << "xz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EXZ[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << tensorname+addname << "yy float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EYY[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << tensorname+addname << "yz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EYZ[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << tensorname+addname << "zz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EZZ[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << tensorname+addname << "par float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << pPAR[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << tensorname+addname << "per1 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << pPER1[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << tensorname+addname << "per2 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << pPER2[ii][jj][kk] << endl;
			}
	my_fileE << "SCALARS " << tensorname+addname << "eps float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				my_fileE << EPS[ii][jj][kk] << endl;
			}
	my_fileE.close();
}


void writeVTKtensor_binary(int it, string tensorname, string addname, double*** EXX, double*** EXY,
					double*** EXZ, double*** EYY, double*** EYZ, double*** EZZ,
					double*** pPAR, double*** pPER1, double*** pPER2, double*** EPS) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;

    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = tensorname + addname +"_xyz_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << tensorname +addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << tensorname +addname << " Tensor from Parsek" << endl;
	my_fileE << "BINARY" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << nzn*ZLEN << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN*nzn*ZLEN << endl;

    cout << "WRITING Tensor " << tensorname +addname<<" TO VTK FILE" << endl;

	my_fileE << "SCALARS " << tensorname+addname << "xx float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EXX[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << tensorname+addname << "xy float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EXY[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << tensorname+addname << "xz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EXZ[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << tensorname+addname << "yy float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EYY[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << tensorname+addname << "yz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EYZ[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << tensorname+addname << "zz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EZZ[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << tensorname+addname << "par float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(pPAR[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << tensorname+addname << "per1 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(pPER1[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << tensorname+addname << "per2 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(pPER2[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE << "SCALARS " << tensorname+addname << "eps float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	for (int kk=0; kk < nzn*ZLEN;kk++)
		for (int jj=0; jj < nyn*YLEN;jj++)
			for (int ii=0; ii < nxn*XLEN;ii++){
				float fl=ReverseFloat(EPS[ii][jj][kk]);
				my_fileE.write((char*)&fl,sizeof(fl));
			}
	my_fileE.close();
}

void writeVTKvect(int it, string vectname, string addname, double** EX, double** EY,double** EZ) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;
	
    int cycle = it * DeltaT;
	stringstream stringcycle;
	stringcycle << cycle;
	string temp;
	temp = vectname +addname +"_AVG_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for" << vectname+addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << vectname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN << endl;
	my_fileE << "VECTORS " << vectname+addname << " float" << endl;
	cout << "WRITING VECTOR " << vectname +addname<<" TO VTK FILE" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EX[ii][jj] << " " << EY[ii][jj] << " " << EZ[ii][jj] << endl;
			// my_file << temp_storageX[node] << " " << temp_storageY[node] << " " << temp_storageZ[node] << endl;
		}
	my_fileE.close();
}


void writeVTKscalar(int it, string scalarname, string addname, double** EX) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;
	
    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +addname +"_AVG_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname+addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN << endl;
	my_fileE << "SCALARS " << scalarname+addname << " float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname +addname<<" TO VTK FILE" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EX[ii][jj] << endl;
		}
	my_fileE.close();
}


void writeVTKscalar_species(int it, string scalarname, double** EX, double** EY) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;
	
    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_AVG_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN << endl;
	my_fileE << "SCALARS " << scalarname << "0 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"0 TO VTK FILE" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EX[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << scalarname << "1 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"1 TO VTK FILE" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EY[ii][jj] << endl;
		}
	my_fileE.close();
}

void writeVTKscalar_species(int it, string scalarname, double** EX, double** EY, double** EZ) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;
	
    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = scalarname +"_AVG_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << scalarname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << scalarname << " Field from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN << endl;
	my_fileE << "SCALARS " << scalarname << "0 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"0 TO VTK FILE" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EX[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << scalarname << "1 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"1 TO VTK FILE" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EY[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << scalarname << "2 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	cout << "WRITING SCALAR " << scalarname <<"2 TO VTK FILE" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EZ[ii][jj] << endl;
		}
	my_fileE.close();
}



void writeVTKtensor(int it, string tensorname, string addname, double** EXX, double** EXY,
					double** EXZ, double** EYY, double** EYZ, double** EZZ,
					double** pPAR, double** pPER1, double** pPER2, double** EPS) {
//	double dx = Lx/nxc;
//	double dy = Ly/nyc;
//	double dz = Lz/nzc;
	
    int cycle = it * DeltaT;
	string temp;
	stringstream stringcycle;
	stringcycle << cycle;
	temp = tensorname + addname +"_AVG_cycle" +stringcycle.str();
	temp += ".vtk";
	cout << "Writing file: " << temp << endl;
	ofstream my_fileE(temp.c_str());
	cout << "writing to file mesh points for " << tensorname +addname << endl;
	my_fileE << "# vtk DataFile Version 1.0" << endl;
	my_fileE << tensorname +addname << " Tensor from Parsek" << endl;
	my_fileE << "ASCII" << endl;
	my_fileE << "DATASET STRUCTURED_POINTS" << endl;
	my_fileE << "DIMENSIONS " << nxn*XLEN << " " << nyn*YLEN << " " << 1 << endl;
	my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
	my_fileE << "SPACING " << Dx << " " << Dy << " " << Dz << endl;
	my_fileE << "POINT_DATA " << nxn*XLEN*nyn*YLEN << endl;
	
    cout << "WRITING Tensor " << tensorname +addname<<" TO VTK FILE" << endl;
	
	my_fileE << "SCALARS " << tensorname+addname << "xx float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EXX[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "xy float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EXY[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "xz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EXZ[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "yy float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EYY[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "yz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EYZ[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "zz float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EZZ[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "par float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << pPAR[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "per1 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << pPER1[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "per2 float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << pPER2[ii][jj] << endl;
		}
	my_fileE << "SCALARS " << tensorname+addname << "eps float" << endl;
	my_fileE << "LOOKUP_TABLE default" << endl;
	
	for (int jj=0; jj < nyn*YLEN;jj++)
		for (int ii=0; ii < nxn*XLEN;ii++){
			my_fileE << EPS[ii][jj] << endl;
		}
	my_fileE.close();
}

void readVTKvect(string filename, double*** BX, double*** BY,double*** BZ) {
	string line;
	string s;
	ifstream myfile(filename.c_str());
	
	
	getline (myfile,line);
	getline (myfile,line);
	myfile >> s; // read ASCII
	myfile >> s; // read dataset
	myfile >> s; // read strucutred_points
	myfile >> s; //read dimensions
	
	myfile >> nxc;
	myfile >> nyc; 
	myfile >> nzc;
	
//	cout <<"nxc= " << nxc << " nyc=" << nyc << " nzc=" <<nzc << endl; 
	myfile >> s; //read origin
	double ox;
	double oy;
	double oz;
	
	myfile >> ox; 
	myfile >> oy; 
	myfile >> oz; 
	myfile >> s; //read spacing
	
	double dx, dy, dz;

	myfile >> dx; 
	myfile >> dy; 
	myfile >> dz; 
//	cout <<"dx= " << dx << " dy=" << dy << " dz=" <<dz << endl; 
	
	myfile >> s; //read point_data
	
	int npoints;
	myfile >> npoints;
	myfile >> s; //read vectros
	myfile >> s; //read b
	myfile >> s; //read double
	cout << s<< endl;
		
	for(int i = 0; i < nxc; i++) 
		for(int j = 0; j < nyc; j++)
			for(int k = 0; k < nzc; k++){
				myfile >> BX[i][j][k];
				myfile >> BY[i][j][k];
				myfile >> BZ[i][j][k];
				//cout <<"bx= " << bx << " by=" << by << " bz=" <<bz << endl; 
			}
    myfile.close();	
}	
void readVTKpreamble(string filename) {
	string line;
	string s;
	ifstream myfile(filename.c_str());
	
	
	getline (myfile,line);
	getline (myfile,line);
	myfile >> s; // read ASCII
	myfile >> s; // read dataset
	myfile >> s; // read strucutred_points
	myfile >> s; //read dimensions
	
	myfile >> nxc;
	myfile >> nyc; 
	myfile >> nzc;
	
	//	cout <<"nxc= " << nxc << " nyc=" << nyc << " nzc=" <<nzc << endl; 
	myfile >> s; //read origin
	double ox;
	double oy;
	double oz;
	
	myfile >> ox; 
	myfile >> oy; 
	myfile >> oz; 
	myfile >> s; //read spacing
	
	double dx, dy, dz;
	myfile >> dx; 
	myfile >> dy; 
	myfile >> dz; 
	//	cout <<"dx= " << dx << " dy=" << dy << " dz=" <<dz << endl; 
	
	myfile >> s; //read point_data
	
	int npoints;
	myfile >> npoints;
    myfile.close();	
}	

int irot(int i){
int ii = (i+ xshift) % (nxn*XLEN);
if (ii<0) ii = (nxn*XLEN) + ii;;
return ii;
}


