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

    void init(long long, double, int, int, int);
    void SetNp(long long n){np = n;};

    void memalloc(long long);
    void memallocX() { x = new double[np];};
    void memallocY() { y = new double[np];};
    void memallocZ() { z = new double[np];};
    void memallocU() { u = new double[np];};
    void memallocV() { v = new double[np];};
    void memallocW() { w = new double[np];};
    void memallocQ() { q = new double[np];};

    void freeX(){delete [] x;};
    void freeY(){delete [] y;};
    void freeZ(){delete [] z;};
    void freeU(){delete [] u;};
    void freeV(){delete [] v;};
    void freeW(){delete [] w;};
    void freeQ(){delete [] q;};
    void info();
    long long getNP(){return np;};
    double getQ(int i){return q[i];};
    double getX(int i){return x[i];};
    double getY(int i){return y[i];};
    double getZ(int i){return z[i];};
    double getU(int i){return u[i];};
    double getV(int i){return v[i];};
    double getW(int i){return w[i];};

    double *getQ(){return q;};
    double *getX(){return x;};
    double *getY(){return y;};
    double *getZ(){return z;};
    double *getU(){return u;};
    double *getV(){return v;};
    double *getW(){return w;};

    long long & getNref(){return np;};
    double *& getQref(){return q;};
    double *& getXref(){return x;};
    double *& getYref(){return y;};
    double *& getZref(){return z;};
    double *& getUref(){return u;};
    double *& getVref(){return v;};
    double *& getWref(){return w;};

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
    void SetNameCycle(std::string name, int c);

    void OpenFieldsFile(std::string dtype, int nspec, int ntx, int nty, int ntz, int *coord, int *dimns, MPI_Comm CART_COMM);
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

    void OpenFieldsFile(std::string dtype, int nspec, int ntx, int nty, int ntz, int *coord, int *dimns, MPI_Comm CART_COMM);
    void ReadFields(double ***field, std::string fname, int nx, int ny, int nz, int rank=-1);
    void CloseFieldsFile();

    void OpenPartclFile(int ns);
    void ReadParticles(int rank, int nproc, int *dimns, double *L, MPI_Comm CART_COMM);
    long long GetNp(int s) {return part[s].getNP();};
    void DumpPartclX(double *& tgt, int s);
    void DumpPartclY(double *& tgt, int s);
    void DumpPartclZ(double *& tgt, int s);
    void DumpPartclU(double *& tgt, int s);
    void DumpPartclV(double *& tgt, int s);
    void DumpPartclW(double *& tgt, int s);
    void DumpPartclQ(double *& tgt, int s);
    void ClosePartclFile();

  private:
    H5hutpart   *part;
    int         nspec;
    int         recycle;
    std::string basename;
    h5_file_t   *partfile;
    h5_file_t   *fldsfile;

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
