/*******************************************************************************************
  Basic.h  -  Basic operations 
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/
#ifndef Basic_H
#define Basic_H

#include <iostream>
#include <math.h>

#include "MPIdata.h"
#include "EllipticF.h"
#include "Alloc.h"

using std::cout;
using std::endl;


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
inline double dotP(double *vect1, double *vect2, int n) {
  double result = 0;
  double local_result = 0;
  for (register int i = 0; i < n; i++)
    local_result += vect1[i] * vect2[i];
  MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (result);

}
/** method to calculate dot product */
inline double dot(double *vect1, double *vect2, int n) {
  double result = 0;
  for (int i = 0; i < n; i++)
    result += vect1[i] * vect2[i];
  return (result);
}
/** method to calculate the square norm of a vector */
inline double norm2(double **vect, int nx, int ny) {
  double result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      result += vect[i][j] * vect[i][j];
  return (result);
}
/** method to calculate the square norm of a vector */
inline double norm2(const array_ref3_double& vect, int nx, int ny) {
  double result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      result += vect.get(i,j,0) * vect.get(i,j,0);
  return (result);
}
/** method to calculate the square norm of a vector */
inline double norm2(double *vect, int nx) {
  double result = 0;
  for (int i = 0; i < nx; i++)
    result += vect[i] * vect[i];
  return (result);
}



/** method to calculate the parallel dot product */
inline double norm2P(const array_ref3_double& vect, int nx, int ny, int nz) {
  double result = 0;
  double local_result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        local_result += vect.get(i,j,k) * vect.get(i,j,k);

  MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (result);
}
/** method to calculate the parallel norm of a vector on different processors with the ghost cell */
inline double norm2P(double *vect, int n) {
  double result = 0;
  double local_result = 0;
  for (int i = 0; i < n; i++)
    local_result += vect[i] * vect[i];
  MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return (result);
}
/** method to calculate the parallel norm of a vector on different processors with the gost cell*/
inline double normP(double *vect, int n) {
  double result = 0.0;
  double local_result = 0.0;
  for (register int i = 0; i < n; i++)
    local_result += vect[i] * vect[i];


  MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return (sqrt(result));

}
/** method to calculate the difference of two vectors*/
inline void sub(double *res, double *vect1, double *vect2, int n) {
  for (register int i = 0; i < n; i++)
    res[i] = vect1[i] - vect2[i];
}
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double *vect1, double *vect2, int n) {
  for (register int i = 0; i < n; i++)
    vect1[i] += vect2[i];


}
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(array_ref3_double& vect1, const array_ref3_double& vect2, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) += vect2.get(i,j,k);
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(array_ref3_double& vect1, const array_ref3_double& vect2, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) += vect2.get(i,j,0);
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(array_ref3_double& vect1, const array_ref4_double& vect2, int nx, int ny, int nz, int ns) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) += vect2.get(ns,i,j,k);
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(array_ref3_double& vect1, const array_ref4_double& vect2, int nx, int ny, int ns) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) += vect2.get(ns,i,j,0);
}
/** method to calculate the subtraction of two vectors vector1 = vector1 - vector2*/
inline void sub(array_ref3_double& vect1, const array_ref3_double& vect2, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) -= vect2.get(i,j,k);
}

/** method to calculate the subtraction of two vectors vector1 = vector1 - vector2*/
inline void sub(array_ref3_double& vect1, const array_ref3_double& vect2, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) -= vect2.get(i,j,0);
}


/** method to sum 4 vectors vector1 = alfa*vector1 + beta*vector2 + gamma*vector3 + delta*vector4 */
inline void sum4(array_ref3_double& vect1, double alfa, const array_ref3_double& vect2, double beta, const array_ref3_double& vect3, double gamma, const array_ref3_double& vect4, double delta, const array_ref3_double& vect5, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = alfa * (vect2.get(i,j,k) + beta * vect3.get(i,j,k) + gamma * vect4.get(i,j,k) + delta * vect5.get(i,j,k));

}
/** method to calculate the scalar-vector product */
inline void scale(double *vect, double alfa, int n) {
  for (register int i = 0; i < n; i++)
    vect[i] *= alfa;
}

/** method to calculate the scalar-vector product */
inline void scale(array_ref3_double& vect, double alfa, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect.fetch(i,j,0) *= alfa;
}


/** method to calculate the scalar-vector product */
inline void scale(array_ref3_double& vect, double alfa, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect.fetch(i,j,k) *= alfa;
}
/** method to calculate the scalar product */
inline void scale(double vect[][2][2], double alfa, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect[i][j][k] *= alfa;
}
/** method to calculate the scalar-vector product */
inline void scale(array_ref3_double& vect1, const array_ref3_double& vect2, double alfa, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = vect2.get(i,j,k) * alfa;
}

/** method to calculate the scalar-vector product */
inline void scale(array_ref3_double& vect1, const array_ref3_double& vect2, double alfa, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) = vect2.get(i,j,0) * alfa;
}

/** method to calculate the scalar-vector product */
inline void scale(double *vect1, double *vect2, double alfa, int n) {
  for (register int i = 0; i < n; i++)
    vect1[i] = vect2[i] * alfa;
}

/** method to calculate vector1 = vector1 + alfa*vector2   */
inline void addscale(double alfa, array_ref3_double& vect1, const array_ref3_double& vect2, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = vect1.get(i,j,k) + alfa * vect2.get(i,j,k);
}
/** add scale for weights */
inline void addscale(double alfa, double vect1[][2][2], double vect2[][2][2], int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect1[i][j][k] + alfa * vect2[i][j][k];

}
/** method to calculate vector1 = vector1 + alfa*vector2   */
inline void addscale(double alfa, array_ref3_double& vect1, const array_ref3_double& vect2, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) += alfa * vect2.get(i,j,0);
}
/** method to calculate vector1 = vector1 + alfa*vector2   */
inline void addscale(double alfa, double *vect1, double *vect2, int n) {
  for (register int i = 0; i < n; i++)
    vect1[i] += alfa * vect2[i];

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2   */
inline void addscale(double alfa, double beta, double *vect1, double *vect2, int n) {
  for (register int i = 0; i < n; i++)
    vect1[i] = vect1[i] * beta + alfa * vect2[i];

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
inline void addscale(double alfa, double beta, array_ref3_double& vect1, const array_ref3_double& vect2, int nx, int ny, int nz) {

  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++) {
        vect1.fetch(i,j,k) = beta * vect1.get(i,j,k) + alfa * vect2.get(i,j,k);
      }

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
inline void addscale(double alfa, double beta, array_ref3_double& vect1, const array_ref3_double& vect2, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) = beta * vect1.get(i,j,0) + alfa * vect2.get(i,j,0);

}


/** method to calculate vector1 = alfa*vector2 + beta*vector3 */
inline void scaleandsum(array_ref3_double& vect1, double alfa, double beta, const array_ref3_double& vect2, const array_ref3_double& vect3, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = alfa * vect2.get(i,j,k) + beta * vect3.get(i,j,k);
}
/** method to calculate vector1 = alfa*vector2 + beta*vector3 with vector2 depending on species*/
inline void scaleandsum(array_ref3_double& vect1, double alfa, double beta, const array_ref4_double& vect2, const array_ref3_double& vect3, int ns, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = alfa * vect2.get(ns,i,j,k) + beta * vect3.get(i,j,k);
}
/** method to calculate vector1 = alfa*vector2*vector3 with vector2 depending on species*/
inline void prod(array_ref3_double& vect1, double alfa, const array_ref4_double& vect2, int ns, const array_ref3_double& vect3, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = alfa * vect2.get(ns,i,j,k) * vect3.get(i,j,k);

}
/** method to calculate vect1 = vect2/alfa */
inline void div(array_ref3_double& vect1, double alfa, const array_ref3_double& vect2, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = vect2.get(i,j,k) / alfa;

}
inline void prod6(array_ref3_double& vect1, const array_ref3_double& vect2, const array_ref3_double& vect3, const array_ref3_double& vect4, const array_ref3_double& vect5, const array_ref3_double& vect6, const array_ref3_double& vect7, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = vect2.get(i,j,k) * vect3.get(i,j,k) + vect4.get(i,j,k) * vect5.get(i,j,k) + vect6.get(i,j,k) * vect7.get(i,j,k);
}
/** method used for calculating PI */
inline void proddiv(array_ref3_double& vect1, const array_ref3_double& vect2, double alfa, const array_ref3_double& vect3, const array_ref3_double& vect4, const array_ref3_double& vect5, const array_ref3_double& vect6, double beta, const array_ref3_double& vect7, const array_ref3_double& vect8, double gamma, const array_ref3_double& vect9, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = (vect2.get(i,j,k) + alfa * (vect3.get(i,j,k) * vect4.get(i,j,k) - vect5.get(i,j,k) * vect6.get(i,j,k)) + beta * vect7.get(i,j,k) * vect8.get(i,j,k)) / (1 + gamma * vect9.get(i,j,k));

  // questo mi convince veramente poco!!!!!!!!!!!!!! CAZZO!!!!!!!!!!!!!!!!!!
  // ***vect1++ = (***vect2++ + alfa*((***vect3++)*(***vect4++) - (***vect5++)*(***vect6++)) + beta*(***vect7++)*(***vect8++))/(1+gamma*(***vect9++));
}
/** method to calculate the opposite of a vector */
inline void neg(array_ref3_double& vect, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect.fetch(i,j,k) = -vect.get(i,j,k);
}

/** method to calculate the opposite of a vector */
inline void neg(array_ref3_double& vect, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect.fetch(i,j,0) = -vect.get(i,j,0);
}
/** method to calculate the opposite of a vector */
inline void neg(array_ref3_double& vect, int nx) {
  for (register int i = 0; i < nx; i++)
    vect.fetch(i,0,0) = -vect.get(i,0,0);
}
/** method to calculate the opposite of a vector */
inline void neg(double *vect, int n) {
  for (register int i = 0; i < n; i++)
    vect[i] = -vect[i];


}
/** method to set equal two vectors */
inline void eq(array_ref3_double& vect1, const array_ref3_double& vect2, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(i,j,k) = vect2.get(i,j,k);

}
/** method to set equal two vectors */
inline void eq(array_ref3_double& vect1, const array_ref3_double& vect2, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(i,j,0) = vect2.get(i,j,0);

}

/** method to set equal two vectors */
inline void eq(array_ref4_double& vect1, const array_ref3_double& vect2, int nx, int ny, int is) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect1.fetch(is,i,j,0) = vect2.get(i,j,0);

}
/** method to set equal two vectors */
inline void eq(array_ref4_double& vect1, const array_ref3_double& vect2, int nx, int ny, int nz, int is) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect1.fetch(is,i,j,k) = vect2.get(i,j,k);

}

inline void eq(double *vect1, double *vect2, int n) {
  for (register int i = 0; i < n; i++)
    vect1[i] = vect2[i];
}
/** method to set a vector to a Value */
inline void eqValue(double value, array_ref3_double& vect, int nx, int ny, int nz) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      for (register int k = 0; k < nz; k++)
        vect.fetch(i,j,k) = value;

}
inline void eqValue(double value, double vect[][2][2], int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect[i][j][k] = value;

}
/** method to set a vector to a Value */
inline void eqValue(double value, array_ref3_double& vect, int nx, int ny) {
  for (register int i = 0; i < nx; i++)
    for (register int j = 0; j < ny; j++)
      vect.fetch(i,j,0) = value;

}
/** method to set a vector to a Value */
inline void eqValue(double value, array_ref3_double& vect, int nx) {
  for (register int i = 0; i < nx; i++)
    vect.fetch(i,0,0) = value;

}
/** method to set a vector to a Value */
inline void eqValue(double value, double *vect, int n) {
  for (register int i = 0; i < n; i++)
    vect[i] = value;
}
/** method to put a column in a matrix 2D */
inline void putColumn(double **Matrix, double *vect, int column, int n) {
  for (int i = 0; i < n; i++)
    Matrix[i][column] = vect[i];

}
/** method to get a column in a matrix 2D */
inline void getColumn(double *vect, double **Matrix, int column, int n) {
  for (int i = 0; i < n; i++)
    vect[i] = Matrix[i][column];
}
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
/** method to get rid of the ghost cells */
inline void getRidGhost(double **out, double **in, int nx, int ny) {
  for (register int i = 1; i < nx - 1; i++)
    for (register int j = 1; j < ny - 1; j++)
      out[i - 1][j - 1] = in[i][j];
}
/** method to calculate cross product of two vectors C= A x B */
inline void cross_product(double a1, double a2, double a3, double b1, double b2, double b3, double *c){
  c[0] = a2 * b3 - a3 * b2;
  c[1] = a3 * b1 - a1 * b3;
  c[2] = a1 * b2 - a2 * b1;
}

inline void loopX(double *b, double z, double x, double y, double a, double zc, double xc, double yc, double m){

  double r = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc));
  double theta = acos((z-zc+1e-10)/(r+1e-10));
  double phi = atan2(y-yc,x-xc);
  //double Rho = r * sin(theta);
  double Rho = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));

  double Alpha = Rho/a;
  double Beta = (z-zc)/a;
  double Gamma = (z-zc+1e-10)/(Rho+1e-10);

  double Q = ((1 + Alpha)*(1 + Alpha) + Beta*Beta);
  double k = sqrt(4*Alpha/Q);
  double B0 = m / (2*a); //m * (C_LIGHT * MU0)/(2*a*a*a*M_PI);

  int err = 0;

  double Bz = B0*(EllipticE(k,err)*(1-Alpha*Alpha-Beta*Beta)/(Q-4*Alpha)+EllipticF(k,err))/(M_PI*sqrt(Q));
  double BRho = B0*Gamma*(EllipticE(k,err)*(1+Alpha*Alpha+Beta*Beta)/(Q-4*Alpha)-EllipticF(k,err))/(M_PI*sqrt(Q));

  if (err)
    cout << "Err came back :" << err << endl;

  if ( isnan(BRho) )
    BRho = 0;
  if ( isnan(Bz) )
    Bz = 0;

  double Bx = BRho * cos(phi);
  double By = BRho * sin(phi);

  //for debugging
  /*cout << "\n\nAt (" << x << "," << y << "," << z << "), the field is :" << endl;
    cout << "Bx: " << Bx << " T" << endl;
    cout << "By: " << By << " T" << endl;
    cout << "Bz: " << Bz << " T" << endl;
    cout << "BRho: " << BRho << " T" << endl;*/

  b[1] = Bx;
  b[2] = By;
  b[0] = Bz;
}

inline void loopY(double *b, double y, double z, double x, double a, double yc, double zc, double xc, double m){

  double r = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc));
  double theta = acos((z-zc+1e-10)/(r+1e-10));
  double phi = atan2(y-yc,x-xc);
  //double Rho = r * sin(theta);
  double Rho = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));

  double Alpha = Rho/a;
  double Beta = (z-zc)/a;
  double Gamma = (z-zc+1e-10)/(Rho+1e-10);

  double Q = ((1 + Alpha)*(1 + Alpha) + Beta*Beta);
  double k = sqrt(4*Alpha/Q);
  double B0 = m / (2*a); //m * (C_LIGHT * MU0)/(2*a*a*a*M_PI);

  int err = 0;

  double Bz = B0*(EllipticE(k,err)*(1-Alpha*Alpha-Beta*Beta)/(Q-4*Alpha)+EllipticF(k,err))/(M_PI*sqrt(Q));
  double BRho = B0*Gamma*(EllipticE(k,err)*(1+Alpha*Alpha+Beta*Beta)/(Q-4*Alpha)-EllipticF(k,err))/(M_PI*sqrt(Q));

  if (err)
    cout << "Err came back :" << err << endl;

  if ( isnan(BRho) )
    BRho = 0;
  if ( isnan(Bz) )
    Bz = 0;

  double Bx = BRho * cos(phi);
  double By = BRho * sin(phi);

  //for debugging
  /*cout << "\n\nAt (" << x << "," << y << "," << z << "), the field is :" << endl;
    cout << "Bx: " << Bx << " T" << endl;
    cout << "By: " << By << " T" << endl;
    cout << "Bz: " << Bz << " T" << endl;
    cout << "BRho: " << BRho << " T" << endl;*/

  b[2] = Bx;
  b[0] = By;
  b[1] = Bz;
}

inline void loopZ(double *b, double x, double y, double z, double a, double xc, double yc, double zc, double m){

  double r = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc));
  double theta = acos((z-zc+1e-10)/(r+1e-10));
  double phi = atan2(y-yc,x-xc);

  double Rho = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));

  double Alpha = Rho/a;
  double Beta = (z-zc)/a;
  double Gamma = (z-zc+1e-10)/(Rho+1e-10);

  double Q = ((1 + Alpha)*(1 + Alpha) + Beta*Beta);
  double k = sqrt(4*Alpha/Q);
  double B0 = m / (2*a); //m * (C_LIGHT * MU0)/(2*a*a*a*M_PI);

  int err = 0;

  double Bz = B0*(EllipticE(k,err)*(1-Alpha*Alpha-Beta*Beta)/(Q-4*Alpha)+EllipticF(k,err))/(M_PI*sqrt(Q));
  double BRho = B0*Gamma*(EllipticE(k,err)*(1+Alpha*Alpha+Beta*Beta)/(Q-4*Alpha)-EllipticF(k,err))/(M_PI*sqrt(Q));

  if (err)
    cout << "Err came back :" << err << endl;

  if ( isnan(BRho) )
    BRho = 0;
  if ( isnan(Bz) )
    Bz = 0;

  double Bx = BRho * cos(phi);
  double By = BRho * sin(phi);

  b[0] = Bx;
  b[1] = By;
  b[2] = Bz;
}

#endif
