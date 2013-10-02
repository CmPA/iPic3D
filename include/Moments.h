#ifndef Moments_H
#define Moments_H
#include "Alloc.h"

class Moments10
{
  private:
    arr4_double arr;
    int nx;
    int ny;
    int nz;
  public:
    void set_to_zero();

    // fetch accessors (write access)
    arr4_double fetch_arr() { return arr; }

    Moments10(int nxn, int nyn, int nzn) :
      nx(nxn),
      ny(nyn),
      nz(nzn),
      arr (nxn, nyn, nzn,10)
    {
    };
    ~Moments10(){};
};

// class to accumulate node-centered species moments
// 
#ifdef TENMOMENTS
class TenMoments {
  private:
    arr3_double rho;

    /** current density, defined on nodes */
    arr3_double Jx;
    arr3_double Jy;
    arr3_double Jz;

    /** pressure tensor components, defined on nodes */
    arr3_double pXX;
    arr3_double pXY;
    arr3_double pXZ;
    arr3_double pYY;
    arr3_double pYZ;
    arr3_double pZZ;
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
    arr3_double fetch_rho() { return rho; }
    arr3_double fetch_Jx () { return Jx ; }
    arr3_double fetch_Jy () { return Jy ; }
    arr3_double fetch_Jz () { return Jz ; }
    arr3_double fetch_Pxx() { return pXX; }
    arr3_double fetch_Pxy() { return pXY; }
    arr3_double fetch_Pxz() { return pXZ; }
    arr3_double fetch_Pyy() { return pYY; }
    arr3_double fetch_Pyz() { return pYZ; }
    arr3_double fetch_Pzz() { return pZZ; }
  public:
    TenMoments(int nxn, int nyn, int nzn) :
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
    ~TenMoments(){};
    void set_to_zero();
};
#endif // TENMOMENTS

#endif
