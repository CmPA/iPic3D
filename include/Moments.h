#ifndef Moments_H
#define Moments_H

// class to accumulate node-centered species moments
// 
class Moments {
  private:
    double ***rho;

    /** current density, defined on nodes */
    double ***Jx;
    double ***Jy;
    double ***Jz;

    /** pressure tensor components, defined on nodes */
    double ***pXX;
    double ***pXY;
    double ***pXZ;
    double ***pYY;
    double ***pYZ;
    double ***pZZ;
    int nx;
    int ny;
    int nz;
  public:
    // get accessors (read access)
    int get_nx() const { return nx; }
    int get_ny() const { return ny; }
    int get_nz() const { return nz; }
    double get_rho(int i, int j, int k) const { return rho[i][j][k]; }
    double get_Jx (int i, int j, int k) const { return Jx [i][j][k]; }
    double get_Jy (int i, int j, int k) const { return Jy [i][j][k]; }
    double get_Jz (int i, int j, int k) const { return Jz [i][j][k]; }
    double get_pXX(int i, int j, int k) const { return pXX[i][j][k]; }
    double get_pXY(int i, int j, int k) const { return pXY[i][j][k]; }
    double get_pXZ(int i, int j, int k) const { return pXZ[i][j][k]; }
    double get_pYY(int i, int j, int k) const { return pYY[i][j][k]; }
    double get_pYZ(int i, int j, int k) const { return pYZ[i][j][k]; }
    double get_pZZ(int i, int j, int k) const { return pZZ[i][j][k]; }
    // fetch accessors (write access)
    double*** fetch_rho() { return rho; }
    double*** fetch_Jx () { return Jx ; }
    double*** fetch_Jy () { return Jy ; }
    double*** fetch_Jz () { return Jz ; }
    double*** fetch_Pxx() { return pXX; }
    double*** fetch_Pxy() { return pXY; }
    double*** fetch_Pxz() { return pXZ; }
    double*** fetch_Pyy() { return pYY; }
    double*** fetch_Pyz() { return pYZ; }
    double*** fetch_Pzz() { return pZZ; }
  public:
    Moments() {
    };
    Moments(int nx_, int ny_, int nz_){
      init(nx_,ny_,nz_);
    }
    void init(int nx_, int ny_, int nz_);
    ~Moments();
    void set_to_zero();
};

#endif
