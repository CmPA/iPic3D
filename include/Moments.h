#ifndef Moments_H
#define Moments_H
#include "Alloc.h"

// class to accumulate node-centered species moments
// 
class Moments {
  private:
    doubleArr3 rho;

    /** current density, defined on nodes */
    doubleArr3 Jx;
    doubleArr3 Jy;
    doubleArr3 Jz;

    /** pressure tensor components, defined on nodes */
    doubleArr3 pXX;
    doubleArr3 pXY;
    doubleArr3 pXZ;
    doubleArr3 pYY;
    doubleArr3 pYZ;
    doubleArr3 pZZ;
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
    doubleArr3& fetch_rho() { return rho; }
    doubleArr3& fetch_Jx () { return Jx ; }
    doubleArr3& fetch_Jy () { return Jy ; }
    doubleArr3& fetch_Jz () { return Jz ; }
    doubleArr3& fetch_Pxx() { return pXX; }
    doubleArr3& fetch_Pxy() { return pXY; }
    doubleArr3& fetch_Pxz() { return pXZ; }
    doubleArr3& fetch_Pyy() { return pYY; }
    doubleArr3& fetch_Pyz() { return pYZ; }
    doubleArr3& fetch_Pzz() { return pZZ; }
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
