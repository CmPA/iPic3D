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



int main(int argc, char **argv) {
  // cycle we want to open
  int n_cycle;
  sscanf(argv[1], "%d", &n_cycle);
  // hdf stuff 
  hid_t file_id;
  hid_t dataset_id;
  herr_t status;
  // Open the settings file 
  file_id = H5Fopen("settings.hdf", H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0) {
    cout << "couldn't open file: settings.hdf" << endl;
    return -1;
  }
  // First read the topology
  int nproc;
  dataset_id = H5Dopen(file_id, "/topology/Nprocs");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nproc);
  status = H5Dclose(dataset_id);
  int XLEN;
  dataset_id = H5Dopen(file_id, "/topology/XLEN");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &XLEN);
  status = H5Dclose(dataset_id);
  int YLEN;
  dataset_id = H5Dopen(file_id, "/topology/YLEN");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &YLEN);
  status = H5Dclose(dataset_id);
  int ZLEN;
  dataset_id = H5Dopen(file_id, "/topology/ZLEN");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ZLEN);
  status = H5Dclose(dataset_id);

  // read Lx 
  double Lx;
  dataset_id = H5Dopen(file_id, "/collective/Lx");
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lx);
  status = H5Dclose(dataset_id);
  // read Ly
  double Ly;
  dataset_id = H5Dopen(file_id, "/collective/Ly");
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Ly);
  status = H5Dclose(dataset_id);
  // read Lz
  double Lz;
  dataset_id = H5Dopen(file_id, "/collective/Lz");
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Lz);
  status = H5Dclose(dataset_id);
  // read nxc
  int nxc;
  dataset_id = H5Dopen(file_id, "/collective/Nxc");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nxc);
  status = H5Dclose(dataset_id);
  // read nyc
  int nyc;
  dataset_id = H5Dopen(file_id, "/collective/Nyc");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nyc);
  status = H5Dclose(dataset_id);
  // read nyc
  int nzc;
  dataset_id = H5Dopen(file_id, "/collective/Nzc");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nzc);
  status = H5Dclose(dataset_id);
  // read ns
  int ns;
  dataset_id = H5Dopen(file_id, "/collective/Ns");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ns);
  // at this point you can close settings
  status = H5Fclose(file_id);
  // prepare to read the proc files
  hid_t *proc_file_id = new hid_t[nproc];
  string temp;
  int *cartesian_cor = new int[3];
  int mappa[140][80][40];
  for (int i = 0; i < nproc; i++) {
    stringstream ss;
    ss << i;
    temp = "proc" + ss.str() + ".hdf";
    // temp = "proc" << ss.c_str() <<".hdf";
    proc_file_id[i] = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (proc_file_id[i] < 0) {
      cout << "couldn't open file:  " << temp << endl;
      return -1;
    }
    // read the position in the topology
    dataset_id = H5Dopen(proc_file_id[i], "/topology/cartesian_coord");
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, cartesian_cor);
    mappa[cartesian_cor[0]][cartesian_cor[1]][cartesian_cor[2]] = i;
    cout << "file" << i << " in topology[" << cartesian_cor[0] << "][" << cartesian_cor[1] << "][" << cartesian_cor[2] << "]" << endl;
    status = H5Dclose(dataset_id);
    H5Fclose(proc_file_id[i]);

  }
  // open the output file
  stringstream cc;
  cc << n_cycle;
  // prepare the file
  // int nxn = nxc/XLEN + 1;
  // int nyn = nyc/YLEN + 1;
  // int nzn = nzc/ZLEN + 1;
  int nxn = nxc / XLEN;
  int nyn = nyc / YLEN;
  int nzn = nzc / ZLEN;
  double dx = Lx / nxc;
  double dy = Ly / nyc;
  double dz = Lz / nzc;
  // E
  temp = "E_cycle" + cc.str();
  temp += ".vtk";
  cout << "Preparing file: " << temp << endl;
  ofstream my_fileE(temp.c_str());
  // phi
  temp = "phi_cycle" + cc.str();
  temp += ".vtk";
  cout << "Preparing file: " << temp << endl;
  ofstream my_filePhi(temp.c_str());
  // rhoi
  temp = "rhoi_cycle" + cc.str();
  temp += ".vtk";
  cout << "Preparing file: " << temp << endl;
  ofstream my_fileRhoi(temp.c_str());
  // rhoe
  temp = "rhoe_cycle" + cc.str();
  temp += ".vtk";
  cout << "Preparing file: " << temp << endl;
  ofstream my_fileRhoe(temp.c_str());
  // E
  cout << "writing to file mesh points for E" << endl;
  my_fileE << "# vtk DataFile Version 1.0" << endl;
  my_fileE << "Electric Field from Parsek" << endl;
  my_fileE << "ASCII" << endl;
  my_fileE << "DATASET STRUCTURED_POINTS" << endl;
  my_fileE << "DIMENSIONS " << nxn * XLEN << " " << nyn * YLEN << " " << nzn * ZLEN << endl;
  my_fileE << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
  my_fileE << "SPACING " << dx << " " << dy << " " << dz << endl;
  my_fileE << "POINT_DATA " << nxn * nyn * nzn * nproc << endl;
  my_fileE << "VECTORS E float" << endl;
  // phi
  cout << "writing to file mesh points for Phi" << endl;
  my_filePhi << "# vtk DataFile Version 1.0" << endl;
  my_filePhi << "Potential from Parsek" << endl;
  my_filePhi << "ASCII" << endl;
  my_filePhi << "DATASET STRUCTURED_POINTS" << endl;
  my_filePhi << "DIMENSIONS " << nxn * XLEN << " " << nyn * YLEN << " " << nzn * ZLEN << endl;
  my_filePhi << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
  my_filePhi << "SPACING " << dx << " " << dy << " " << dz << endl;
  my_filePhi << "POINT_DATA " << nxn * nyn * nzn * nproc << endl;
  my_filePhi << "SCALARS Phi float" << endl;
  my_filePhi << "LOOKUP_TABLE default" << endl;
  // rhoi
  cout << "writing to file mesh points for Rhoi" << endl;
  my_fileRhoi << "# vtk DataFile Version 1.0" << endl;
  my_fileRhoi << "Ion Density from Parsek" << endl;
  my_fileRhoi << "ASCII" << endl;
  my_fileRhoi << "DATASET STRUCTURED_POINTS" << endl;
  my_fileRhoi << "DIMENSIONS " << nxn * XLEN << " " << nyn * YLEN << " " << nzn * ZLEN << endl;
  my_fileRhoi << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
  my_fileRhoi << "SPACING " << dx << " " << dy << " " << dz << endl;
  my_fileRhoi << "POINT_DATA " << nxn * nyn * nzn * nproc << endl;
  my_fileRhoi << "SCALARS Rhoi float" << endl;
  my_fileRhoi << "LOOKUP_TABLE default" << endl;
  // rhoe
  cout << "writing to file mesh points for Rhoe" << endl;
  my_fileRhoe << "# vtk DataFile Version 1.0" << endl;
  my_fileRhoe << "Ion Density from Parsek" << endl;
  my_fileRhoe << "ASCII" << endl;
  my_fileRhoe << "DATASET STRUCTURED_POINTS" << endl;
  my_fileRhoe << "DIMENSIONS " << nxn * XLEN << " " << nyn * YLEN << " " << nzn * ZLEN << endl;
  my_fileRhoe << "ORIGIN " << 0.0 << " " << 0.0 << " " << 0.0 << endl;
  my_fileRhoe << "SPACING " << dx << " " << dy << " " << dz << endl;
  my_fileRhoe << "POINT_DATA " << nxn * nyn * nzn * nproc << endl;
  my_fileRhoe << "SCALARS Rhoe float" << endl;
  my_fileRhoe << "LOOKUP_TABLE default" << endl;
  cout << "READING VECTOR FROM HDF5 FILES" << endl;
  double *temp_storageX = new double[(nxn + 1) * (nyn + 1) * (nzn + 1)];
  double *temp_storageY = new double[(nxn + 1) * (nyn + 1) * (nzn + 1)];
  double *temp_storageZ = new double[(nxn + 1) * (nyn + 1) * (nzn + 1)];
  double *temp_storagePHI = new double[(nxn) * (nyn) * (nzn)];
  double *temp_storageRHO = new double[(nxn + 1) * (nyn + 1) * (nzn + 1)];

  double ***EX = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);
  double ***EY = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);
  double ***EZ = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);
  double ***PHI = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);
  double ***RHO = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);



  // write Electric field
  int proc = 0;
  for (int i = 0; i < XLEN; i++)
    for (int j = 0; j < YLEN; j++)
      for (int k = 0; k < ZLEN; k++) {
        stringstream ss;
        ss << mappa[i][j][k];
        temp = "proc" + ss.str() + ".hdf";
        proc_file_id[mappa[i][j][k]] = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        // read data

        temp = "/fields/Ex/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageX);
        status = H5Dclose(dataset_id);
        temp = "/fields/Ey/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageY);
        status = H5Dclose(dataset_id);
        temp = "/fields/Ez/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageZ);
        status = H5Dclose(dataset_id);
        int node = 0;
        for (int ii = 0; ii < (nxn + 1); ii++)
          for (int jj = 0; jj < (nyn + 1); jj++)
            for (int kk = 0; kk < (nzn + 1); kk++) {
              if (ii != nxn && jj != nyn && kk != nzn) {
                EX[ii + nxn * i][jj + nyn * j][kk + nzn * k] = temp_storageX[node];
                EY[ii + nxn * i][jj + nyn * j][kk + nzn * k] = temp_storageY[node];
                EZ[ii + nxn * i][jj + nyn * j][kk + nzn * k] = temp_storageZ[node];
              }
              node++;
            }
        // close the file
        H5Fclose(proc_file_id[mappa[i][j][k]]);
        // go to other proc
        proc++;
      }
  cout << "WRITING VECTOR E TO VTK FILE" << endl;
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++) {
        my_fileE << EX[ii][jj][kk] << " " << EY[ii][jj][kk] << " " << EZ[ii][jj][kk] << endl;
        // my_file << temp_storageX[node] << " " << temp_storageY[node] << " " << temp_storageZ[node] << endl;
      }
  my_fileE.close();

  // write PHI
  proc = 0;
  for (int i = 0; i < XLEN; i++)
    for (int j = 0; j < YLEN; j++)
      for (int k = 0; k < ZLEN; k++) {
        stringstream ss;
        ss << mappa[i][j][k];
        temp = "proc" + ss.str() + ".hdf";
        proc_file_id[mappa[i][j][k]] = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        // read data

        temp = "/potentials/phi/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storagePHI);
        status = H5Dclose(dataset_id);
        int node = 0;
        for (int ii = 0; ii < (nxn); ii++)
          for (int jj = 0; jj < (nyn); jj++)
            for (int kk = 0; kk < (nzn); kk++) {
              if (ii != nxn && jj != nyn && kk != nzn) {
                PHI[ii + nxn * i][jj + nyn * j][kk + nzn * k] = temp_storagePHI[node];
              }
              node++;
            }
        // close the file
        H5Fclose(proc_file_id[mappa[i][j][k]]);
        // go to other proc
        proc++;
      }
  cout << "WRITING SCALAR PHI TO VTK FILE" << endl;
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++) {
        my_filePhi << PHI[ii][jj][kk] << endl;

      }
  my_filePhi.close();

  // write Rhoi
  proc = 0;
  for (int i = 0; i < XLEN; i++)
    for (int j = 0; j < YLEN; j++)
      for (int k = 0; k < ZLEN; k++) {
        stringstream ss;
        ss << mappa[i][j][k];
        temp = "proc" + ss.str() + ".hdf";
        proc_file_id[mappa[i][j][k]] = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        // read data

        temp = "/moments/species_1/rho/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageRHO);
        status = H5Dclose(dataset_id);
        int node = 0;
        for (int ii = 0; ii < (nxn + 1); ii++)
          for (int jj = 0; jj < (nyn + 1); jj++)
            for (int kk = 0; kk < (nzn + 1); kk++) {
              if (ii != nxn && jj != nyn && kk != nzn) {
                RHO[ii + nxn * i][jj + nyn * j][kk + nzn * k] = temp_storageRHO[node];
              }
              node++;
            }
        // close the file
        H5Fclose(proc_file_id[mappa[i][j][k]]);
        // go to other proc
        proc++;
      }
  cout << "WRITING SCALAR RHOi TO VTK FILE" << endl;
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++) {
        my_fileRhoi << RHO[ii][jj][kk] << endl;

      }
  my_fileRhoi.close();

  // write Rhoe
  proc = 0;
  for (int i = 0; i < XLEN; i++)
    for (int j = 0; j < YLEN; j++)
      for (int k = 0; k < ZLEN; k++) {
        stringstream ss;
        ss << mappa[i][j][k];
        temp = "proc" + ss.str() + ".hdf";
        proc_file_id[mappa[i][j][k]] = H5Fopen(temp.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        // read data

        temp = "/moments/species_0/rho/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageRHO);
        status = H5Dclose(dataset_id);
        int node = 0;
        for (int ii = 0; ii < (nxn + 1); ii++)
          for (int jj = 0; jj < (nyn + 1); jj++)
            for (int kk = 0; kk < (nzn + 1); kk++) {
              if (ii != nxn && jj != nyn && kk != nzn) {
                RHO[ii + nxn * i][jj + nyn * j][kk + nzn * k] = temp_storageRHO[node];
              }
              node++;
            }
        // close the file
        H5Fclose(proc_file_id[mappa[i][j][k]]);
        // go to other proc
        proc++;
      }
  cout << "WRITING SCALAR RHOi TO VTK FILE" << endl;
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++) {
        my_fileRhoe << RHO[ii][jj][kk] << endl;

      }
  my_fileRhoe.close();

  delete[]proc_file_id;
  delete[]temp_storageX;
  delete[]temp_storageY;
  delete[]temp_storageZ;
  delArr3(EX, nxn * XLEN, nyn * YLEN);
  delArr3(EY, nxn * XLEN, nyn * YLEN);
  delArr3(EZ, nxn * XLEN, nyn * YLEN);



  return (0);
}
