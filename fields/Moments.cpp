#include "Moments.h"
#include "Alloc.h"

void Moments10::set_to_zero()
{
  arr.setall(0);
  //#pragma omp parallel for collapse(4)
  //for (register int i = 0; i < nx; i++)
  //for (register int j = 0; j < ny; j++)
  //for (register int k = 0; k < nz; k++)
  //for (register int m = 0; m < 10; m++)
  //{
  //  arr[i][j][k][m] = 0.0;
  //}
}

