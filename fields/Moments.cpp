#include "Moments.h"
#include "Alloc.h"

// construct empty instance (not zeroed)
void Moments::init(int nx_, int ny_, int nz_)
{
  nx = nx_;
  ny = ny_;
  nz = nz_;
  rho = newArr3(double, nx, ny, nz);
  Jx = newArr3(double, nx, ny, nz);
  Jy = newArr3(double, nx, ny, nz);
  Jz = newArr3(double, nx, ny, nz);
  pXX = newArr3(double, nx, ny, nz);
  pXY = newArr3(double, nx, ny, nz);
  pXZ = newArr3(double, nx, ny, nz);
  pYY = newArr3(double, nx, ny, nz);
  pYZ = newArr3(double, nx, ny, nz);
  pZZ = newArr3(double, nx, ny, nz);
}

Moments::~Moments() {
  // nodes and species
  delArr3(rho, nx, ny);
  delArr3(Jx, nx, ny);
  delArr3(Jy, nx, ny);
  delArr3(Jz, nx, ny);
  delArr3(pXX, nx, ny);
  delArr3(pXY, nx, ny);
  delArr3(pXZ, nx, ny);
  delArr3(pYY, nx, ny);
  delArr3(pYZ, nx, ny);
  delArr3(pZZ, nx, ny);
}

void Moments::set_to_zero() {
  // #pragma omp parallel for collapse(1)
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++) {
        rho[i][j][k] = 0.0;
        Jx[i][j][k] = 0.0;
        Jy[i][j][k] = 0.0;
        Jz[i][j][k] = 0.0;
        pXX[i][j][k] = 0.0;
        pXY[i][j][k] = 0.0;
        pXZ[i][j][k] = 0.0;
        pYY[i][j][k] = 0.0;
        pYZ[i][j][k] = 0.0;
        pZZ[i][j][k] = 0.0;
      }
}

