#ifndef Moments_H
#define Moments_H
#include "Alloc.h"

// class to accumulate node-centered species moments
// 
class Moments {
  private:
    array_ref3_double rho;

    /** current density, defined on nodes */
    array_ref3_double Jx;
    array_ref3_double Jy;
    array_ref3_double Jz;

    /** pressure tensor components, defined on nodes */
    array_ref3_double pXX;
    array_ref3_double pXY;
    array_ref3_double pXZ;
    array_ref3_double pYY;
    array_ref3_double pYZ;
    array_ref3_double pZZ;
    int nx;
    int ny;
    int nz;
  public:
    // get accessors (read access)
    int get_nx() const { return nx; }
    int get_ny() const { return ny; }
    int get_nz() const { return nz; }
    double get_rho(int i, int j, int k) const { return rho.get(i,j,k); }
    double get_Jx (int i, int j, int k) const { return Jx .get(i,j,k); }
    double get_Jy (int i, int j, int k) const { return Jy .get(i,j,k); }
    double get_Jz (int i, int j, int k) const { return Jz .get(i,j,k); }
    double get_pXX(int i, int j, int k) const { return pXX.get(i,j,k); }
    double get_pXY(int i, int j, int k) const { return pXY.get(i,j,k); }
    double get_pXZ(int i, int j, int k) const { return pXZ.get(i,j,k); }
    double get_pYY(int i, int j, int k) const { return pYY.get(i,j,k); }
    double get_pYZ(int i, int j, int k) const { return pYZ.get(i,j,k); }
    double get_pZZ(int i, int j, int k) const { return pZZ.get(i,j,k); }
    // fetch accessors (write access)
    array_ref3_double& fetch_rho() { return rho; }
    array_ref3_double& fetch_Jx () { return Jx ; }
    array_ref3_double& fetch_Jy () { return Jy ; }
    array_ref3_double& fetch_Jz () { return Jz ; }
    array_ref3_double& fetch_Pxx() { return pXX; }
    array_ref3_double& fetch_Pxy() { return pXY; }
    array_ref3_double& fetch_Pxz() { return pXZ; }
    array_ref3_double& fetch_Pyy() { return pYY; }
    array_ref3_double& fetch_Pyz() { return pYZ; }
    array_ref3_double& fetch_Pzz() { return pZZ; }
  public:
    Moments(int nxn, int nyn, int nzn) :
      nx(nxn),
      ny(nyn),
      nz(nzn),
      rho (nxn, nyn, nzn),
      Jx  (nxn, nyn, nzn),
      Jy  (nxn, nyn, nzn),
      Jz  (nxn, nyn, nzn),
      pXX (nxn, nyn, nzn),
      pXY (nxn, nyn, nzn),
      pXZ (nxn, nyn, nzn),
      pYY (nxn, nyn, nzn),
      pYZ (nxn, nyn, nzn),
      pZZ (nxn, nyn, nzn)
    {
    };
    ~Moments(){};
    void set_to_zero();
};

#endif
