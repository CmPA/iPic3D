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
    void addRho(double weight[][2][2], int X, int Y, int Z);
    void addJx(double weight[][2][2], int X, int Y, int Z);
    void addJy(double weight[][2][2], int X, int Y, int Z);
    void addJz(double weight[][2][2], int X, int Y, int Z);

    void addPxx(double weight[][2][2], int X, int Y, int Z);
    void addPxy(double weight[][2][2], int X, int Y, int Z);
    void addPxz(double weight[][2][2], int X, int Y, int Z);
    void addPyy(double weight[][2][2], int X, int Y, int Z);
    void addPyz(double weight[][2][2], int X, int Y, int Z);
    void addPzz(double weight[][2][2], int X, int Y, int Z);
};

/** add an amount of charge density to charge density field at node X,Y */
inline void Moments::addRho(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        rho[X - i][Y - j][Z - k] += weight[i][j][k];
      }
}
/** add an amount of charge density to current density - direction X to current density field on the node*/
inline void Moments::addJx(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        Jx[X - i][Y - j][Z - k] += weight[i][j][k];
      }
}
/** add an amount of current density - direction Y to current density field on the node */
inline void Moments::addJy(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        Jy[X - i][Y - j][Z - k] += weight[i][j][k];
      }
}
/** add an amount of current density - direction Z to current density field on the node */
inline void Moments::addJz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        Jz[X - i][Y - j][Z - k] += weight[i][j][k];
      }
}
/** add an amount of pressure density - direction XX to current density field on the node */
inline void Moments::addPxx(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        pXX[X - i][Y - j][Z - k] += weight[i][j][k];
      }
}
/** add an amount of pressure density - direction XY to current density field on the node*/
inline void Moments::addPxy(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        pXY[X - i][Y - j][Z - k] += weight[i][j][k];
      }
}
/** add an amount of pressure density - direction XZ to current density field on the node */
inline void Moments::addPxz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        pXZ[X - i][Y - j][Z - k] += weight[i][j][k];
      }
}
/** add an amount of pressure density - direction YY to current density field on the node*/
inline void Moments::addPyy(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        pYY[X - i][Y - j][Z - k] += weight[i][j][k];
      }
}
/** add an amount of pressure density - direction YZ to current density field on the node */
inline void Moments::addPyz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        pYZ[X - i][Y - j][Z - k] += weight[i][j][k];
      }
}
/** add an amount of pressure density - direction ZZ to current density field on the node */
inline void Moments::addPzz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        pZZ[X - i][Y - j][Z - k] += weight[i][j][k];
      }
}

#endif
