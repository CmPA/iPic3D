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
  for (int i = 0; i < n; i++)
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
inline double norm2(double ***vect, int nx, int ny) {
  double result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      result += vect[i][j][0] * vect[i][j][0];
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
inline double norm2P(double ***vect, int nx, int ny, int nz) {
  double result = 0;
  double local_result = 0;
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        local_result += vect[i][j][k] * vect[i][j][k];

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
  for (int i = 0; i < n; i++)
    local_result += vect[i] * vect[i];


  MPI_Allreduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return (sqrt(result));

}
/** method to calculate the difference of two vectors*/
inline void sub(double *res, double *vect1, double *vect2, int n) {
  for (int i = 0; i < n; i++)
    res[i] = vect1[i] - vect2[i];
}
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double *vect1, double *vect2, int n) {
  for (int i = 0; i < n; i++)
    vect1[i] += vect2[i];


}
/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double ***vect1, double ***vect2, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] += vect2[i][j][k];
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double ***vect1, double ***vect2, int nx, int ny) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect1[i][j][0] += vect2[i][j][0];
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double ***vect1, double ****vect2, int nx, int ny, int nz, int ns) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] += vect2[ns][i][j][k];
}

/** method to calculate the sum of two vectors vector1 = vector1 + vector2*/
inline void sum(double ***vect1, double ****vect2, int nx, int ny, int ns) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect1[i][j][0] += vect2[ns][i][j][0];
}
/** method to calculate the subtraction of two vectors vector1 = vector1 - vector2*/
inline void sub(double ***vect1, double ***vect2, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] -= vect2[i][j][k];


}

/** method to calculate the subtraction of two vectors vector1 = vector1 - vector2*/
inline void sub(double ***vect1, double ***vect2, int nx, int ny) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect1[i][j][0] -= vect2[i][j][0];


}


/** method to sum 4 vectors vector1 = alfa*vector1 + beta*vector2 + gamma*vector3 + delta*vector4 */
inline void sum4(double ***vect1, double alfa, double ***vect2, double beta, double ***vect3, double gamma, double ***vect4, double delta, double ***vect5, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = alfa * (vect2[i][j][k] + beta * vect3[i][j][k] + gamma * vect4[i][j][k] + delta * vect5[i][j][k]);

}
/** method to calculate the scalar-vector product */
inline void scale(double *vect, double alfa, int n) {
  for (int i = 0; i < n; i++)
    vect[i] *= alfa;
}

/** method to calculate the scalar-vector product */
inline void scale(double ***vect, double alfa, int nx, int ny) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect[i][j][0] *= alfa;
}


/** method to calculate the scalar-vector product */
inline void scale(double ***vect, double alfa, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect[i][j][k] *= alfa;
}
/** method to calculate the scalar product */
inline void scale(double vect[][2][2], double alfa, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect[i][j][k] *= alfa;
}
/** method to calculate the scalar-vector product */
inline void scale(double ***vect1, double ***vect2, double alfa, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect2[i][j][k] * alfa;
}

/** method to calculate the scalar-vector product */
inline void scale(double ***vect1, double ***vect2, double alfa, int nx, int ny) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect1[i][j][0] = vect2[i][j][0] * alfa;
}

/** method to calculate the scalar-vector product */
inline void scale(double *vect1, double *vect2, double alfa, int n) {
  for (int i = 0; i < n; i++)
    vect1[i] = vect2[i] * alfa;
}

/** method to calculate vector1 = vector1 + alfa*vector2   */
inline void addscale(double alfa, double ***vect1, double ***vect2, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect1[i][j][k] + alfa * vect2[i][j][k];
}
/** add scale for weights */
inline void addscale(double alfa, double vect1[][2][2], double vect2[][2][2], int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect1[i][j][k] + alfa * vect2[i][j][k];

}
/** method to calculate vector1 = vector1 + alfa*vector2   */
inline void addscale(double alfa, double ***vect1, double ***vect2, int nx, int ny) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect1[i][j][0] += alfa * vect2[i][j][0];
}
/** method to calculate vector1 = vector1 + alfa*vector2   */
inline void addscale(double alfa, double *vect1, double *vect2, int n) {
  for (int i = 0; i < n; i++)
    vect1[i] += alfa * vect2[i];

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2   */
inline void addscale(double alfa, double beta, double *vect1, double *vect2, int n) {
  for (int i = 0; i < n; i++)
    vect1[i] = vect1[i] * beta + alfa * vect2[i];

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
inline void addscale(double alfa, double beta, double ***vect1, double ***vect2, int nx, int ny, int nz) {

  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++) {
        vect1[i][j][k] = beta * vect1[i][j][k] + alfa * vect2[i][j][k];
      }

}
/** method to calculate vector1 = beta*vector1 + alfa*vector2 */
inline void addscale(double alfa, double beta, double ***vect1, double ***vect2, int nx, int ny) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect1[i][j][0] = beta * vect1[i][j][0] + alfa * vect2[i][j][0];

}

/** method to calculate v1 = v1 + alpha * v2 * v3 */
inline void sumscalprod(double ***vect1, double alfa, double ***vect2, double ***vect3, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] += alfa * vect2[i][j][k] * vect3[i][j][k];
}

/** method to calculate vector1 = alfa*vector2 + beta*vector3 */
inline void scaleandsum(double ***vect1, double alfa, double beta, double ***vect2, double ***vect3, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = alfa * vect2[i][j][k] + beta * vect3[i][j][k];
}
/** method to calculate vector1 = alfa*vector2 + beta*vector3 with vector2 depending on species*/
inline void scaleandsum(double ***vect1, double alfa, double beta, double ****vect2, double ***vect3, int ns, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = alfa * vect2[ns][i][j][k] + beta * vect3[i][j][k];
}
/** method to calculate vector1 = alfa*vector2*vector3 with vector2 depending on species*/
inline void prod(double ***vect1, double alfa, double ****vect2, int ns, double ***vect3, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = alfa * vect2[ns][i][j][k] * vect3[i][j][k];

}
/** method to calculate vect1 = vect2/alfa */
inline void div(double ***vect1, double alfa, double ***vect2, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect2[i][j][k] / alfa;

}
inline void prod6(double ***vect1, double ***vect2, double ***vect3, double ***vect4, double ***vect5, double ***vect6, double ***vect7, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect2[i][j][k] * vect3[i][j][k] + vect4[i][j][k] * vect5[i][j][k] + vect6[i][j][k] * vect7[i][j][k];
}
/** method used for calculating PI */
inline void proddiv(double ***vect1, double ***vect2, double alfa, double ***vect3, double ***vect4, double ***vect5, double ***vect6, double beta, double ***vect7, double ***vect8, double gamma, double ***vect9, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = (vect2[i][j][k] + alfa * (vect3[i][j][k] * vect4[i][j][k] - vect5[i][j][k] * vect6[i][j][k]) + beta * vect7[i][j][k] * vect8[i][j][k]) / (1 + gamma * vect9[i][j][k]);

  // questo mi convince veramente poco!!!!!!!!!!!!!! CAZZO!!!!!!!!!!!!!!!!!!
  // ***vect1++ = (***vect2++ + alfa*((***vect3++)*(***vect4++) - (***vect5++)*(***vect6++)) + beta*(***vect7++)*(***vect8++))/(1+gamma*(***vect9++));
}
/** method to calculate the opposite of a vector */
inline void neg(double ***vect, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect[i][j][k] = -vect[i][j][k];


}

/** method to calculate the opposite of a vector */
inline void neg(double ***vect, int nx, int ny) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect[i][j][0] = -vect[i][j][0];
}
/** method to calculate the opposite of a vector */
inline void neg(double ***vect, int nx) {
  for (int i = 0; i < nx; i++)
    vect[i][0][0] = -vect[i][0][0];
}
/** method to calculate the opposite of a vector */
inline void neg(double *vect, int n) {
  for (int i = 0; i < n; i++)
    vect[i] = -vect[i];


}
/** method to set equal two vectors */
inline void eq(double ***vect1, double ***vect2, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[i][j][k] = vect2[i][j][k];

}
/** method to set equal two vectors */
inline void eq(double ***vect1, double ***vect2, int nx, int ny) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect1[i][j][0] = vect2[i][j][0];

}

/** method to set equal two vectors */
inline void eq(double ****vect1, double ***vect2, int nx, int ny, int is) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect1[is][i][j][0] = vect2[i][j][0];

}
/** method to set equal two vectors */
inline void eq(double ****vect1, double ***vect2, int nx, int ny, int nz, int is) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect1[is][i][j][k] = vect2[i][j][k];

}

inline void eq(double *vect1, double *vect2, int n) {
  for (int i = 0; i < n; i++)
    vect1[i] = vect2[i];
}
/** method to set a vector to a Value */
inline void eqValue(double value, double ****vect, int ns, int nx, int ny, int nz) {

  for (int s = 0; s < ns; s++)
    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++)
          vect[s][i][j][k] = value;

}
inline void eqValue(double value, double ****vect, int ns, int nx, int ny, int nz, int is) {

    for (int i = 0; i < nx; i++)
      for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++)
          vect[is][i][j][k] = value;

}
inline void eqValue(double value, double ***vect, int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect[i][j][k] = value;

}
inline void eqValue(double value, double vect[][2][2], int nx, int ny, int nz) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      for (int k = 0; k < nz; k++)
        vect[i][j][k] = value;

}
/** method to set a vector to a Value */
inline void eqValue(double value, double ***vect, int nx, int ny) {
  for (int i = 0; i < nx; i++)
    for (int j = 0; j < ny; j++)
      vect[i][j][0] = value;

}
/** method to set a vector to a Value */
inline void eqValue(double value, double ***vect, int nx) {
  for (int i = 0; i < nx; i++)
    vect[i][0][0] = value;

}
/** method to set a vector to a Value */
inline void eqValue(double value, double *vect, int n) {
  for (int i = 0; i < n; i++)
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
/** method to calculate vector1 = vectro1 * vector2 + (1-vector2) * vector3 */
inline void weight_avg(double ***vect1, double ***vect2,  double ***vect3, int nx, int ny, int nz){
   for (int i=0; i< nx; i++)
    for (int j=0; j< ny; j++)
      for (int k=0; k< nz; k++){
	double teta = fabs(vect2[i][j][k]);
	if(teta>1.0) teta=1.0;

         vect1[i][j][k] = vect1[i][j][k]*(1.0 - pow(teta,2.0)) + vect3[i][j][k]* pow(teta,2.0);
}
}
/** method to overwrite vector1 = vector2 where vector3 is above threshold */
inline void weight_threshold(double ***vect1, double ***vect2, double ***vect3, double gate,  int nx, int ny, int nz){
   for (int i=0; i< nx; i++)
    for (int j=0; j< ny; j++)
      for (int k=0; k< nz; k++){

    	  if(vect3[i][j][k]>gate) {
    		  vect1[i][j][k] = vect2[i][j][k];
    	  }
      }
}
/** method to calculate tapering vector1 = vectro1 * (1-vector2) */
inline void weight_tapering(double ***vect1, double ***vect2,   int nx, int ny, int nz){
   for (int i=0; i< nx; i++)
    for (int j=0; j< ny; j++)
      for (int k=0; k< nz; k++){
	double teta = fabs(vect2[i][j][k]);
	if(teta>1.0) teta=1.0;

         vect1[i][j][k] = vect1[i][j][k]*(1.0 - pow(teta,2.0)) ;
}
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
  for (int i = 1; i < nx - 1; i++)
    for (int j = 1; j < ny - 1; j++)
      out[i - 1][j - 1] = in[i][j];
}
/** method to calculate cross product of two vectors C= A x B */
inline void cross_product(double a1, double a2, double a3, double b1, double b2, double b3, double *c){
  c[0] = a2 * b3 - a3 * b2;
  c[1] = a3 * b1 - a1 * b3;
  c[2] = a1 * b2 - a2 * b1;
}

/** method to calculate loop normal to z */
inline void loopapproxZ(double *b, double x, double y, double z, double a, double xc, double yc, double zc, double m){

  double r = sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc)+(z-zc)*(z-zc));
  double theta = acos((z-zc+1e-10)/(r+1e-10));
  double phi = atan2(y-yc,x-xc);

  double br = m*cos(theta)*(2.0*a*a+2*r*r+a*r*sin(theta))/pow(1e-10+a*a+r*r+2*a*r*sin(theta),2.5);
  double bt = -m*sin(theta)*(2.0*a*a-r*r+a*r*sin(theta))/pow(1e-10+a*a+r*r+2*a*r*sin(theta),2.5);

  b[0] = br* cos(phi)*sin(theta)+bt*cos(phi)*cos(theta);
  b[1] = br* sin(phi)*sin(theta)+bt*sin(phi)*cos(theta);
  b[2] = br* cos(theta)-bt*sin(theta);
  //cout << bx << "   " << x << "   " << xc << "   " << m << endl;
}

// This function is identically loopZ
inline void loopWork(double *b, double x, double y, double z, double a, double xc, double yc, double zc, double m){

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
	//double B0 = m / (2*a); //m * (C_LIGHT * MU0)/(2*a*a*a*M_PI);

	int err = 0;

	double central_value = (EllipticE(0.0,err)+EllipticF(0.0,err))/M_PI;
	//cout << central_value << "   " << xc << "   " << yc << "   " <<zc <<endl;
	double B0 = m  /central_value;

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

   b[0] = Bx;
   b[1] = By;
   b[2] = Bz;
 }

// Here, we will break up the calculation into a superposition of current loops
inline void loopWorker(double *b, double x, double y, double z, double a, double xc, double yc, double zc, double m, char *printed){

	int nLoops = 1;
	double DeltaA[] = { 0 };
	double DeltaZ[] = { 0 };
	/*int nLoops = 351;
	double DeltaA[] = {-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.01444,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.013267692,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.012095385,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.010923077,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.009750769,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.008578462,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.007406154,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.006233846,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.005061538,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.003889231,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.002716923,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.001544615,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308,-0.000372308};
	double DeltaZ[] = {-0.02,-0.018768485,-0.01753697,-0.016305455,-0.015073939,-0.013842424,-0.012610909,-0.011379394,-0.010147879,-0.008916364,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.018768485,-0.01753697,-0.016305455,-0.015073939,-0.013842424,-0.012610909,-0.011379394,-0.010147879,-0.008916364,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.01753697,-0.016305455,-0.015073939,-0.013842424,-0.012610909,-0.011379394,-0.010147879,-0.008916364,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.016305455,-0.015073939,-0.013842424,-0.012610909,-0.011379394,-0.010147879,-0.008916364,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.015073939,-0.013842424,-0.012610909,-0.011379394,-0.010147879,-0.008916364,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.013842424,-0.012610909,-0.011379394,-0.010147879,-0.008916364,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.012610909,-0.011379394,-0.010147879,-0.008916364,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.011379394,-0.010147879,-0.008916364,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.010147879,-0.008916364,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.008916364,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.007684848,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.006453333,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485,-0.005221818,-0.003990303,-0.002758788,-0.001527273,-0.000295758,0.000935758,0.002167273,0.003398788,0.004630303,0.005861818,0.007093333,0.008324848,0.009556364,0.010787879,0.012019394,0.013250909,0.014482424,0.015713939,0.016945455,0.01817697,0.019408485};
	*/double mFactor[nLoops];

	double PhysicalL = 0.25;
	double CodeL = 100;

	//Initialize the loop-dependant variables
	for (int i = 0; i < nLoops; i++){
		// If we're under the center of the machine, then we need to flip
		//  the sign on the DeltaZ
		if (z > zc)
			DeltaZ[i] *= -1.0;

		DeltaZ[i] *= CodeL/PhysicalL;
		DeltaA[i] *= CodeL/PhysicalL;

		mFactor[i] = 1.0/nLoops;
	}


	// Now, do the calculate the loop fields nLoops times and add them together
	for (int i = 0; i < nLoops; i++){
		double buffer[3] = {0,0,0};

		// Lets only print once per coil, but print info about each loop
		if (*printed < 2 && x < 1.0 && y < 1.0 && z < 1.0){
			cout << "Loop " << i+1 << " centered at (" << xc;
			cout << "," << yc << "," << zc+DeltaZ[i] << ") with ";
			cout << m*mFactor[i] << " Amp-Turns and radius " << a+DeltaA[i] << endl;

		}

		loopWork(buffer, x, y, z, a+DeltaA[i], xc, yc, zc+DeltaZ[i], m*mFactor[i]);

		b[0] += buffer[0];
		b[1] += buffer[1];
		b[2] += buffer[2];
	}

	(*printed)++;
}

inline void loopY(double *b, double y, double z, double x, double a, double yc, double zc, double xc, double m){

	double buffer[3] = {0,0,0};
	static char printed = 0;

	if (printed < 2 && x < 1.0 && y < 1.0 && z < 1.0)
		cout << "LoopY: Location (y, z, x)\n";

	loopWorker(buffer, x, y, z, a, xc, yc, zc, m, &printed);


	b[2] = buffer[0];
	b[0] = buffer[1];
	b[1] = buffer[2];
 }
inline void loopX(double *b, double z, double x, double y, double a, double zc, double xc, double yc, double m){

	double buffer[3] = {0,0,0};
	static char printed = 0;

	if (printed < 2 && x < 1.0 && y < 1.0 && z < 1.0)
		cout << "LoopX: Location (z, x, y)\n";


	loopWorker(buffer, x, y, z, a, xc, yc, zc, m, &printed);

	b[1] = buffer[0];
	b[2] = buffer[1];
	b[0] = buffer[2];
 }
inline void loopZ(double *b, double x, double y, double z, double a, double xc, double yc, double zc, double m){

	double buffer[3] = {0,0,0};
	static char printed = 0;

	if (printed < 2 && x < 1.0 && y < 1.0 && z < 1.0)
		cout << "LoopZ: Location (x,y,z)\n";

	loopWorker(buffer, x, y, z, a, xc, yc, zc, m, &printed);

	b[0] = buffer[0];
	b[1] = buffer[1];
	b[2] = buffer[2];
}



#endif

