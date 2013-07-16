#ifndef Moments_H
#define Moments_H

// class to accumulate node-centered species moments
// 
class Moments
{
private:
  double invVOL;
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
  int get_nx()const{return nx;}
  int get_ny()const{return ny;}
  int get_nz()const{return nz;}
  double get_invVOL()const{return invVOL;}
  double get_rho(int i,int j,int k)const{return rho[i][j][k];}
  double get_Jx (int i,int j,int k)const{return Jx [i][j][k];}
  double get_Jy (int i,int j,int k)const{return Jy [i][j][k];}
  double get_Jz (int i,int j,int k)const{return Jz [i][j][k];}
  double get_pXX(int i,int j,int k)const{return pXX[i][j][k];}
  double get_pXY(int i,int j,int k)const{return pXY[i][j][k];}
  double get_pXZ(int i,int j,int k)const{return pXZ[i][j][k];}
  double get_pYY(int i,int j,int k)const{return pYY[i][j][k];}
  double get_pYZ(int i,int j,int k)const{return pYZ[i][j][k];}
  double get_pZZ(int i,int j,int k)const{return pZZ[i][j][k];}
public:
  Moments(){};
  Moments(int nx_, int ny_, int nz_, double invVOL_){
    init(nx_,ny_,nz_,invVOL_);
  }
  void init(int nx_, int ny_, int nz_, double invVOL_);
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

// construct empty instance (not zeroed)
inline void Moments::init(int nx_, int ny_, int nz_,double invVOL_)
{
  nx = nx_;
  ny = ny_;
  nz = nz_;
  invVOL = invVOL_;
  rho = newArr3(double, nx, ny, nz);
  Jx  = newArr3(double, nx, ny, nz);
  Jy  = newArr3(double, nx, ny, nz);
  Jz  = newArr3(double, nx, ny, nz);
  pXX = newArr3(double, nx, ny, nz);
  pXY = newArr3(double, nx, ny, nz);
  pXZ = newArr3(double, nx, ny, nz);
  pYY = newArr3(double, nx, ny, nz);
  pYZ = newArr3(double, nx, ny, nz);
  pZZ = newArr3(double, nx, ny, nz);
}

inline Moments::~Moments() {
  // nodes and species
  delArr3(rho, nx, ny);
  delArr3(Jx,  nx, ny);
  delArr3(Jy,  nx, ny);
  delArr3(Jz,  nx, ny);
  delArr3(pXX, nx, ny);
  delArr3(pXY, nx, ny);
  delArr3(pXZ, nx, ny);
  delArr3(pYY, nx, ny);
  delArr3(pYZ, nx, ny);
  delArr3(pZZ, nx, ny);
}

inline void Moments::set_to_zero() {
  //#pragma omp parallel for collapse(1)
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++) {
        rho[i][j][k] = 0.0;
        Jx [i][j][k] = 0.0;
        Jy [i][j][k] = 0.0;
        Jz [i][j][k] = 0.0;
        pXX[i][j][k] = 0.0;
        pXY[i][j][k] = 0.0;
        pXZ[i][j][k] = 0.0;
        pYY[i][j][k] = 0.0;
        pYZ[i][j][k] = 0.0;
        pZZ[i][j][k] = 0.0;
      }
}

/** add an amount of charge density to charge density field at node X,Y */
inline void Moments::addRho(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++) {
        const double temp = weight[i][j][k] * invVOL;
        rho[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of charge density to current density - direction X to current density field on the node*/
inline void Moments::addJx(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        const double temp = weight[i][j][k] * invVOL;
        Jx[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of current density - direction Y to current density field on the node */
inline void Moments::addJy(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        const double temp = weight[i][j][k] * invVOL;
        Jy[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of current density - direction Z to current density field on the node */
inline void Moments::addJz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        const double temp = weight[i][j][k] * invVOL;
        Jz[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction XX to current density field on the node */
inline void Moments::addPxx(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        const double temp = weight[i][j][k] * invVOL;
        pXX[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction XY to current density field on the node*/
inline void Moments::addPxy(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        const double temp = weight[i][j][k] * invVOL;
        pXY[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction XZ to current density field on the node */
inline void Moments::addPxz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        const double temp = weight[i][j][k] * invVOL;
        pXZ[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction YY to current density field on the node*/
inline void Moments::addPyy(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        const double temp = weight[i][j][k] * invVOL;
        pYY[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction YZ to current density field on the node */
inline void Moments::addPyz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        const double temp = weight[i][j][k] * invVOL;
        pYZ[X - i][Y - j][Z - k] += temp;
      }
}
/** add an amount of pressure density - direction ZZ to current density field on the node */
inline void Moments::addPzz(double weight[][2][2], int X, int Y, int Z) {
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++){
        const double temp = weight[i][j][k] * invVOL;
        pZZ[X - i][Y - j][Z - k] += temp;
      }
}

#endif
