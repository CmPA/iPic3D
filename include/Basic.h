/*******************************************************************************************
  Basic.h  -  Basic operations 
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/
#ifndef Basic_H
#define Basic_H
#include "arraysfwd.h"
#include <math.h>

/**
 *  
 * Basic operations defined. This library provides methods to calculate:
 *
 * - dot product of two vectors
 * - square norm of a vector
 * - norm of a vector
 * - difference of two vector
 * - scalar-vector product
 * - vector1 = vector1 + alfa*vector2
 * - vector1 = beta*vector1 + alfa*vector2
 * - opposite of a vector
 *
 * 
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */


/** method to calculate the parallel dot product with vect1, vect2 having the ghost cells*/
double dotP(double *vect1, double *vect2, int n);
/** method to calculate dot product */
double dot(double *vect1, double *vect2, int n);
/** method to calculate the square norm of a vector */
double norm2(double **vect, int nx, int ny);
/** method to calculate the square norm of a vector */
double norm2(const arr3_double vect, int nx, int ny);
/** method to calculate the square norm of a vector */
double norm2(double *vect, int nx);
/** method to calculate the parallel dot product */
double norm2P(const arr3_double vect, int nx, int ny, int nz);
/** method to calculate the parallel norm of a vector on different processors with the ghost cell */
double norm2P(double *vect, int n);
/** method to calculate the parallel norm of a vector on different processors with the gost cell*/
double normP(double *vect, int n);
/** method to calculate the difference of two vectors*/
void sub(double *res, double *vect1, double *vect2, int n);
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(double *vect1, double *vect2, int n);
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(arr3_double vect1, const arr3_double vect2, int nx, int ny, int nz);
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(arr3_double vect1, const arr3_double vect2, int nx, int ny);
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(arr3_double vect1, const arr4_double vect2, int nx, int ny, int nz, int ns);
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
void sum(arr3_double vect1, const arr4_double vect2, int nx, int ny, int ns);
/** method to calculate the subtraction of two vectors vector1 = vector1 - vector2*/
void sub(arr3_double vect1, const arr3_double vect2, int nx, int ny, int nz);
/** method to calculate the subtraction of two vectors vector1 = vector1 - vector2*/
void sub(arr3_double vect1, const arr3_double vect2, int nx, int ny);
/** method to sum 4 vectors vector1 = alfa*vector1 + beta*vector2 + gamma*vector3 + delta*vector4 */
void sum4(arr3_double vect1, double alfa, const arr3_double vect2, double beta, const arr3_double vect3, double gamma, const arr3_double vect4, double delta, const arr3_double vect5, int nx, int ny, int nz);
/** method to calculate the scalar-vector product */
void scale(double *vect, double alfa, int n);
/** method to calculate the scalar-vector product */
void scale(arr3_double vect, double alfa, int nx, int ny);
/** method to calculate the scalar-vector product */
void scale(arr3_double vect, double alfa, int nx, int ny, int nz);
///** method to calculate the scalar product */
//inline void scale(double vect[][2][2], double alfa, int nx, int ny, int nz) {
//  for (int i = 0; i < nx; i++)
//    for (int j = 0; j < ny; j++)
//      for (int k = 0; k < nz; k++)
//        vect[i][j][k] *= alfa;
//}
/** method to calculate the scalar-vector product */
void scale(arr3_double vect1, const arr3_double vect2, double alfa, int nx, int ny, int nz);
/** method to calculate the scalar-vector product */
void scale(arr3_double vect1, const arr3_double vect2, double alfa, int nx, int ny);
/** method to calculate the scalar-vector product */
void scale(double *vect1, double *vect2, double alfa, int n);
/** method to calculate vector1 = vector1 + alfa*vector2   */
void addscale(double alfa, arr3_double vect1, const arr3_double vect2, int nx, int ny, int nz);
/** add scale for weights */
void addscale(double alfa, double vect1[][2][2], double vect2[][2][2], int nx, int ny, int nz);
/** method to calculate vector1 = vector1 + alfa*vector2   */
void addscale(double alfa, arr3_double vect1, const arr3_double vect2, int nx, int ny);
/** method to calculate vector1 = vector1 + alfa*vector2   */
void addscale(double alfa, double *vect1, double *vect2, int n);
/** method to calculate vector1 = beta*vector1 + alfa*vector2   */
void addscale(double alfa, double beta, double *vect1, double *vect2, int n);
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
void addscale(double alfa, double beta, arr3_double vect1, const arr3_double vect2, int nx, int ny, int nz);
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
void addscale(double alfa, double beta, arr3_double vect1, const arr3_double vect2, int nx, int ny);
/** method to calculate vector1 = alfa*vector2 + beta*vector3 */
void scaleandsum(arr3_double vect1, double alfa, double beta, const arr3_double vect2, const arr3_double vect3, int nx, int ny, int nz);
/** method to calculate vector1 = alfa*vector2 + beta*vector3 with vector2 depending on species*/
void scaleandsum(arr3_double vect1, double alfa, double beta, const arr4_double vect2, const arr3_double vect3, int ns, int nx, int ny, int nz);
/** method to calculate vector1 = alfa*vector2*vector3 with vector2 depending on species*/
void prod(arr3_double vect1, double alfa, const arr4_double vect2, int ns, const arr3_double vect3, int nx, int ny, int nz);
/** method to calculate vect1 = vect2/alfa */
void div(arr3_double vect1, double alfa, const arr3_double vect2, int nx, int ny, int nz);
void prod6(arr3_double vect1, const arr3_double vect2, const arr3_double vect3, const arr3_double vect4, const arr3_double vect5, const arr3_double vect6, const arr3_double vect7, int nx, int ny, int nz);
/** method used for calculating PI */
void proddiv(arr3_double vect1, const arr3_double vect2, double alfa, const arr3_double vect3, const arr3_double vect4, const arr3_double vect5, const arr3_double vect6, double beta, const arr3_double vect7, const arr3_double vect8, double gamma, const arr3_double vect9, int nx, int ny, int nz);
/** method to calculate the opposite of a vector */
void neg(arr3_double vect, int nx, int ny, int nz);
/** method to calculate the opposite of a vector */
void neg(arr3_double vect, int nx, int ny);
/** method to calculate the opposite of a vector */
void neg(arr3_double vect, int nx);
/** method to calculate the opposite of a vector */
void neg(double *vect, int n);
/** method to set equal two vectors */
void eq(arr3_double vect1, const arr3_double vect2, int nx, int ny, int nz);
/** method to set equal two vectors */
void eq(arr3_double vect1, const arr3_double vect2, int nx, int ny);
/** method to set equal two vectors */
void eq(arr4_double vect1, const arr3_double vect2, int nx, int ny, int is);
/** method to set equal two vectors */
void eq(arr4_double vect1, const arr3_double vect2, int nx, int ny, int nz, int is);
inline void eq(double *vect1, double *vect2, int n){
  for (register int i = 0; i < n; i++)
    vect1[i] = vect2[i];
}
/** method to set a vector to a Value */
void eqValue(double value, arr3_double vect, int nx, int ny, int nz);
//void eqValue(double value, double vect[][2][2], int nx, int ny, int nz);
/** method to set a vector to a Value */
void eqValue(double value, arr3_double vect, int nx, int ny);
/** method to set a vector to a Value */
void eqValue(double value, arr3_double vect, int nx);
/** method to set a vector to a Value */
void eqValue(double value, double *vect, int n);
/** method to put a column in a matrix 2D */
void putColumn(double **Matrix, double *vect, int column, int n);
/** method to get a column in a matrix 2D */
void getColumn(double *vect, double **Matrix, int column, int n);

/** method to get rid of the ghost cells */
inline void getRidGhost(double **out, double **in, int nx, int ny);

/** RIFAI QUESTA PARTE questo e' la tomba delle performance*/
inline void MODULO(double *x, double L) {
  *x = *x - floor(*x / L) * L;

}
/** method to calculate the epsilon machine */
inline double eps() {
  double eps;
  int i = 1;
  double num = 1;
  double newsum = 1;
  double oldsum = 1;
  while (true) {
    num = num / (2 * i);
    newsum += num;
    if (newsum == oldsum)
      break;
    oldsum = newsum;
    i++;
  }
  eps = num * 2;
  return (eps);
}
/** method to calculate cross product of two vectors C= A x B */
inline void cross_product(double a1, double a2, double a3, double b1, double b2, double b3, double *c){
  c[0] = a2 * b3 - a3 * b2;
  c[1] = a3 * b1 - a1 * b3;
  c[2] = a1 * b2 - a2 * b1;
}

void loopX(double *b, double z, double x, double y, double a, double zc, double xc, double yc, double m);
void loopY(double *b, double y, double z, double x, double a, double yc, double zc, double xc, double m);
void loopZ(double *b, double x, double y, double z, double a, double xc, double yc, double zc, double m);

#endif
