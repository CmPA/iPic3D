#include "Moments.h"
#include "Alloc.h"

void Moments10::set_to_zero()
{
  #pragma omp parallel for collapse(4)
  for (register int i = 0; i < nx; i++)
  for (register int j = 0; j < ny; j++)
  for (register int k = 0; k < nz; k++)
  for (register int m = 0; m < 10; m++)
  {
    arr[i][j][k][m] = 0.0;
  }
}

void Moments::set_to_zero() {
  #pragma omp parallel for collapse(3)
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

