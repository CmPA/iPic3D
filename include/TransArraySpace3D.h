/***************************************************************************
  TransArraySpace.h  -  
  -------------------
begin                : May 2008
copyright            : (C) KUL
developers           : Stefano Markidis, Giovanni Lapenta

 ***************************************************************************/

#ifndef TransArraySpace3D_H
#define TransArraySpace3D_H

/** method to convert a 1D field in a 3D field not considering guard cells*/
inline void solver2phys(array_ref3_double& vectPhys, double *vectSolver, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++)
        vectPhys[i][j][k] = *vectSolver++;

}
/** method to convert a 1D field in a 3D field not considering guard cells*/
inline void solver2phys(array_ref3_double& vectPhys1, array_ref3_double& vectPhys2, array_ref3_double& vectPhys3, double *vectSolver, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++) {
        vectPhys1[i][j][k] = *vectSolver++;
        vectPhys2[i][j][k] = *vectSolver++;
        vectPhys3[i][j][k] = *vectSolver++;
      }
}
/** method to convert a 3D field in a 1D field not considering guard cells*/
inline void phys2solver(double *vectSolver, const array_ref3_double& vectPhys, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++)
        *vectSolver++ = vectPhys.get(i,j,k);
}
/** method to convert a 3D field in a 1D field not considering guard cells*/
inline void phys2solver(double *vectSolver, const array_ref3_double& vectPhys1, const array_ref3_double& vectPhys2, const array_ref3_double& vectPhys3, int nx, int ny, int nz) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      for (register int k = 1; k < nz - 1; k++) {
        *vectSolver++ = vectPhys1.get(i,j,k);
        *vectSolver++ = vectPhys2.get(i,j,k);
        *vectSolver++ = vectPhys3.get(i,j,k);
      }
}
#endif
