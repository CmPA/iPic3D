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

class H5hutpart{
  public:
    H5hutpart();
    ~H5hutpart();

    void init(long long, double);
    void memalloc(long long);
    void info();
    long long getNP(){return np;};
    double getQ(int i){return q[i];};

    double *getQ(){return q;};
    double *getX(){return x;};
    double *getY(){return y;};
    double *getZ(){return z;};
    double *getU(){return u;};
    double *getV(){return v;};
    double *getW(){return w;};

  private:
    long long np;
    double *q;
    double *x;
    double *y;
    double *z;
    double *u;
    double *v;
    double *w;
    bool    allocated;
};

class H5output{
  public:
    void SetBaseName(std::string name);

    void OpenFieldsCellFile(int nspec, int ntcx, int ntcy, int ntcz, int cycle, int *coord, int *dimns, MPI_Comm CART_COMM);
    void WriteFieldsCell(double ***field, std::string fname, int nxc, int nyc, int nzc);
    void CloseFieldsCellFile();

    void OpenPartclFile(int nspec, int cycle, MPI_Comm CART_COMM);
    void WriteParticles(int ispec, long long np, double *q, double *x, double *y, double *z, double *u, double *v, double *w, MPI_Comm CART_COMM);
    void ClosePartclFile();
    
  private:
    std::string basename;
    h5_file_t *partfile;
    h5_file_t *fldsfile;
};

class H5input{
  public:
    void OpenPartclFile(int ns, std::string fn);
    void ReadParticles(int ndim, int rank, int nproc, int *dimns, double *L, MPI_Comm CART_COMM);
    void ClosePartclFile();

  private:
    H5hutpart   *part;
    int         nspec;
    std::string filename;

    void LoadParticles(int ndim, int rank, int nproc, int *dimns, double *L, MPI_Comm CART_COMM);
    void InitParticles(int ndim, int rank, MPI_Comm CART_COMM);

    void LoadLocalParticles(long long *np,
         double *q_loc,
         double *x_loc,
         double *y_loc,
         double *z_loc,
         double *u_loc,
         double *v_loc,
         double *w_loc);
    
    void FindLocalParticles(int nproc, int ndim, h5_int64_t *h5npart, long long nop, int *dimns, double *L, MPI_Comm CART_COMM,
         double *q,
         double *x, double *y, double *z,
         double *u, double *v, double *w);
    
};

#endif
