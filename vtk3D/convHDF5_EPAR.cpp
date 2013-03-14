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
  status = H5Dclose(file_id);
  // prepare to read the proc files
  hid_t *proc_file_id = new hid_t[nproc];
  string temp;
  int *cartesian_cor = new int[3];
  int mappa[40][40][40];
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
  // B
  temp = "B_cycle" + cc.str();
  temp += ".vtk";
  cout << "Preparing file: " << temp << endl;
  ofstream my_file(temp.c_str());
  // E
  temp = "E_cycle" + cc.str();
  temp += ".vtk";
  cout << "Preparing file: " << temp << endl;
  ofstream my_fileE(temp.c_str());
  // Epar
  temp = "Epar_cycle" + cc.str();
  temp += ".vtk";
  cout << "Preparing file: " << temp << endl;
  ofstream my_fileEpar(temp.c_str());
  // Eper
  temp = "Eper_cycle" + cc.str();
  temp += ".vtk";
  cout << "Preparing file: " << temp << endl;
  ofstream my_fileEper(temp.c_str());
  // B
  my_file << "# vtk DataFile Version 1.0" << endl;
  my_file << "Magnetic Field from Parsek" << endl;
  my_file << "ASCII" << endl;
  my_file << "DATASET STRUCTURED_GRID" << endl;
  my_file << "DIMENSIONS " << nxn * XLEN << " " << nyn * YLEN << " " << nzn * ZLEN << endl;
  my_file << "POINTS " << nxn * nyn * nzn * nproc << " float" << endl;
  cout << "writing to file mesh points for B" << endl;
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++)
        my_file << ii * dx << " " << jj * dy << " " << kk * dz << endl;
  my_file << endl;
  my_file << "POINT_DATA " << nxn * nyn * nzn * nproc << endl;
  my_file << "VECTORS B float" << endl;

  // E
  cout << "writing to file mesh points for E" << endl;
  my_fileE << "# vtk DataFile Version 1.0" << endl;
  my_fileE << "Electric Field from Parsek" << endl;
  my_fileE << "ASCII" << endl;
  my_fileE << "DATASET STRUCTURED_GRID" << endl;
  my_fileE << "DIMENSIONS " << nxn * XLEN << " " << nyn * YLEN << " " << nzn * ZLEN << endl;
  my_fileE << "POINTS " << nxn * nyn * nzn * nproc << " float" << endl;
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++)
        my_fileE << ii * dx << " " << jj * dy << " " << kk * dz << endl;

  my_fileE << endl;
  my_fileE << "POINT_DATA " << nxn * nyn * nzn * nproc << endl;
  my_fileE << "VECTORS E float" << endl;
  // Epar
  cout << "writing to file mesh points for Epar" << endl;
  my_fileEpar << "# vtk DataFile Version 1.0" << endl;
  my_fileEpar << "Parallel Electric Field from Parsek" << endl;
  my_fileEpar << "ASCII" << endl;
  my_fileEpar << "DATASET STRUCTURED_GRID" << endl;
  my_fileEpar << "DIMENSIONS " << nxn * XLEN << " " << nyn * YLEN << " " << nzn * ZLEN << endl;
  my_fileEpar << "POINTS " << nxn * nyn * nzn * nproc << " float" << endl;
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++)
        my_fileEpar << ii * dx << " " << jj * dy << " " << kk * dz << endl;

  my_fileEpar << endl;
  my_fileEpar << "POINT_DATA " << nxn * nyn * nzn * nproc << endl;
  my_fileEpar << "SCALARS Epar float" << endl;
  my_fileEpar << "LOOKUP_TABLE default" << endl;
  // Eper
  cout << "writing to file mesh points for Eper" << endl;
  my_fileEper << "# vtk DataFile Version 1.0" << endl;
  my_fileEper << "Parallel Electric Field from Parsek" << endl;
  my_fileEper << "ASCII" << endl;
  my_fileEper << "DATASET STRUCTURED_GRID" << endl;
  my_fileEper << "DIMENSIONS " << nxn * XLEN << " " << nyn * YLEN << " " << nzn * ZLEN << endl;
  my_fileEper << "POINTS " << nxn * nyn * nzn * nproc << " float" << endl;
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++)
        my_fileEper << ii * dx << " " << jj * dy << " " << kk * dz << endl;

  my_fileEper << endl;
  my_fileEper << "POINT_DATA " << nxn * nyn * nzn * nproc << endl;
  my_fileEper << "SCALARS Eper float" << endl;
  my_fileEper << "LOOKUP_TABLE default" << endl;

  // my_file << "SCALARS Bx float" << endl;
  // my_file << "LOOKUP_TABLE default" << endl;
  cout << "READING VECTOR FROM HDF5 FILES" << endl;
  double *temp_storageX = new double[(nxn + 1) * (nyn + 1) * (nzn + 1)];
  double *temp_storageY = new double[(nxn + 1) * (nyn + 1) * (nzn + 1)];
  double *temp_storageZ = new double[(nxn + 1) * (nyn + 1) * (nzn + 1)];
  double ***BX = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);
  double ***BY = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);
  double ***BZ = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);
  double ***EX = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);
  double ***EY = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);
  double ***EZ = newArr3(double, nxn * XLEN, nyn * YLEN, nzn * ZLEN);
  double Bmod;
  double Epar;
  double Eper;
  double Eperx;
  double Epery;
  double Eperz;

  int proc = 0;
  for (int i = 0; i < XLEN; i++)
    for (int j = 0; j < YLEN; j++)
      for (int k = 0; k < ZLEN; k++) {
        temp = "/fields/Bx/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageX);
        temp = "/fields/By/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageY);
        temp = "/fields/Bz/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageZ);
        int node = 0;
        for (int ii = 0; ii < (nxn + 1); ii++)
          for (int jj = 0; jj < (nyn + 1); jj++)
            for (int kk = 0; kk < (nzn + 1); kk++) {
              if (ii != nxn && jj != nyn && kk != nzn) {
                BX[ii + nxn * i][jj + nyn * j][kk + nzn * k] = temp_storageX[node];
                BY[ii + nxn * i][jj + nyn * j][kk + nzn * k] = temp_storageY[node];
                BZ[ii + nxn * i][jj + nyn * j][kk + nzn * k] = temp_storageZ[node];
              }
              node++;
            }
        proc++;
      }
  cout << "WRITING VECTOR B TO VTK FILE" << endl;
  // write to disc
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++) {
        my_file << BX[ii][jj][kk] << " " << BY[ii][jj][kk] << " " << BZ[ii][jj][kk] << endl;
        // my_file << temp_storageX[node] << " " << temp_storageY[node] << " " << temp_storageZ[node] << endl;
      }
  // write Electric field
  proc = 0;
  for (int i = 0; i < XLEN; i++)
    for (int j = 0; j < YLEN; j++)
      for (int k = 0; k < ZLEN; k++) {
        temp = "/fields/Ex/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageX);
        temp = "/fields/Ey/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageY);
        temp = "/fields/Ez/cycle_" + cc.str();
        dataset_id = H5Dopen(proc_file_id[mappa[i][j][k]], temp.c_str());
        status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_storageZ);
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
        proc++;
      }
  cout << "WRITING VECTOR E TO VTK FILE" << endl;
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++) {
        my_fileE << EX[ii][jj][kk] << " " << EY[ii][jj][kk] << " " << EZ[ii][jj][kk] << endl;
        // my_file << temp_storageX[node] << " " << temp_storageY[node] << " " << temp_storageZ[node] << endl;
      }
  cout << "WRITING VECTOR Epar TO VTK FILE" << endl;
  for (int kk = 0; kk < nzn * ZLEN; kk++)
    for (int jj = 0; jj < nyn * YLEN; jj++)
      for (int ii = 0; ii < nxn * XLEN; ii++) {
        Bmod = sqrt(BX[ii][jj][kk] * BX[ii][jj][kk] + BY[ii][jj][kk] * BY[ii][jj][kk] + BZ[ii][jj][kk] * BZ[ii][jj][kk]);
        Epar = (EX[ii][jj][kk] * BX[ii][jj][kk] + EY[ii][jj][kk] * BY[ii][jj][kk] + EZ[ii][jj][kk] * BZ[ii][jj][kk]) / Bmod;
        Eperx = EX[ii][jj][kk] - Epar * BX[ii][jj][kk] / Bmod;
        Epery = EY[ii][jj][kk] - Epar * BY[ii][jj][kk] / Bmod;
        Eperz = EZ[ii][jj][kk] - Epar * BZ[ii][jj][kk] / Bmod;
        Eper = sqrt(Eperx * Eperx + Epery * Epery + Eperz * Eperz);
        my_fileEpar << Epar << endl;
        my_fileEper << Eper << endl;

      }

  delete[]proc_file_id;
  delete[]temp_storageX;
  delete[]temp_storageY;
  delete[]temp_storageZ;
  delArr3(BX, nxn * XLEN, nyn * YLEN);
  delArr3(BY, nxn * XLEN, nyn * YLEN);
  delArr3(BZ, nxn * XLEN, nyn * YLEN);
  delArr3(EX, nxn * XLEN, nyn * YLEN);
  delArr3(EY, nxn * XLEN, nyn * YLEN);
  delArr3(EZ, nxn * XLEN, nyn * YLEN);
  my_file.close();
  my_fileE.close();
  my_fileEpar.close();
  my_fileEper.close();
  return (0);
}
