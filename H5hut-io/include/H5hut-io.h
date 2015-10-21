#ifndef __H5HUTIO_H__
#define __H5HUTIO_H__

#include <mpi.h>
#include <H5hut.h>
#include <stdio.h>

#include <iostream>
#include <iomanip>
#include <ostream>
#include <string>
#include <sstream>

class H5output{
  public:
    void SetNameCycle(std::string name, int c);

    void OpenFieldsFile(std::string dtype, int nspec, int ntx, int nty, int ntz, int *coord, int *pdims, MPI_Comm CART_COMM);
    void WriteFields(double ***field, std::string fname, int nx, int ny, int nz, int rank=-1);
    void CloseFieldsFile();

    void OpenPartclFile(int nspec, MPI_Comm CART_COMM);
    void WriteParticles(int ispec, long long np, double *q, double *x, double *y, double *z, double *u, double *v, double *w, MPI_Comm CART_COMM);
    void ClosePartclFile();
    
  private:
    int         cycle;
    std::string basename;
    h5_file_t   *partfile;
    h5_file_t   *fldsfile;
};

class H5input{
  public:
    void SetNameCycle(std::string name, int rc);

    // Field read functions:
    void OpenFieldsFile(std::string dtype, int nspec, int ntx, int nty, int ntz, int *coord, int *pdims, MPI_Comm CART_COMM);
    void ReadFields(double ***field, std::string fname, int nx, int ny, int nz, int rank=-1);
    void CloseFieldsFile();

    // Particle read functions:
    void      OpenPartclFile(int ns, int rank, int nproc, MPI_Comm CART_COMM);
    void      ReadParticles(int rank, int nproc, int ispec, int *pdims, double *L, MPI_Comm CART_COMM);
    void      LoadParticles(long long sizevec, int rank, int nproc, int ispec, int *pdims, double *L, MPI_Comm CART_COMM, double* q,
                            double *x, double *y, double *z,
                            double *u, double *v, double *w);
    void      ClosePartclFile();
    long long get_nops() {return f_nop;};

  private:
    int         nspec;
    int         recycle;
    std::string basename;

    // Field input:
    h5_file_t   *fldsfile;

    // Particle input PHDF5:
    unsigned int *inproc;
    double       *i_q;
    double       *i_x;
    double       *i_y;
    double       *i_z;
    double       *i_u;
    double       *i_v;
    double       *i_w;
    long long    *s_nop;
    long long    *r_nop;
    long long    *r_beg;
    long long    *nops_beg;
    long long    *nops;
    long long    *h5npart;
    long long    f_nop;
    long long    S_MAX_NOP;
    long long    R_MAX_NOP;
    hid_t        partfile;
    hid_t        pclgroup;

    // Particle read private functions:
    void ReadPartDataset(hid_t group, std::string dtset, long long nops, long long nops_beg, double *arr);
    void FillPartVectors(long long sizevec, int rank, int jproc, int ispec, long long r_nop, long long r_beg, double* r_buffer,
                         double *q, double *x, double *y, double *z, double *u, double *v, double *w);
    void ExchangeParticles(long long sizevec, int nproc, int myrank, int ispec, long long nop, int *pdims, double *L, MPI_Comm CART_COMM,
                           double *q,
                           double *x, double *y, double *z,
                           double *u, double *v, double *w );
    void SortParticles(int nproc, int myrank, int ispec, int ndim, long long nop, int *pdims, double *L, MPI_Comm CART_COMM);
};

#endif
