/*******************************************************************************************
  CG.h  -  Conjugate Gradient Solver
  -------------------
developers: Stefano Markidis, Giovanni Lapenta
 ********************************************************************************************/


#ifndef CG_H
#define CG_H

#include <iostream>

#include "../mathlib/Basic.h"
#include "../utility/TransArraySpace3D.h"

class EMfields3D;
typedef EMfields3D Field;
typedef void (Field::*FIELD_IMAGE) (double *, double *, Grid *, VirtualTopology3D *);
typedef void (*GENERIC_IMAGE) (double *, double *, Grid *, VirtualTopology3D *);
using std::cout;
using std::cerr;
using std::endl;


/**
 * 
 * Electrostatic Field with 3 components(x,y,z) defined on a 2-D grid
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
 * @version 2.0
 *
 */


inline bool CG(double *xkrylov, int xkrylovlen, double *b, int maxit, double tol, FIELD_IMAGE FunctionImage, Grid * grid, VirtualTopology3D * vct, Field * field) {
  // allocate residual, image, p, b, calculated on central points
  double *r = new double[xkrylovlen];
  double *v = new double[xkrylovlen];
  double *z = new double[xkrylovlen];
  double *im = new double[xkrylovlen];
  double c, t, d, initial_error;
  eqValue(0.0, r, xkrylovlen);
  eqValue(0.0, v, xkrylovlen);
  eqValue(0.0, z, xkrylovlen);
  eqValue(0.0, im, xkrylovlen);

  int i = 0;
  bool CONVERGED = false;
  bool CGVERBOSE = false;
  // initial guess for x: all the components are equal to 0
  eqValue(0, xkrylov, xkrylovlen);
  // Compute r = b -Ax
  (field->*FunctionImage) (im, xkrylov, grid, vct);
  sub(r, b, im, xkrylovlen);
  // v = r
  eq(v, r, xkrylovlen);
  c = dotP(r, r, xkrylovlen);
  initial_error = sqrt(c);
  if (vct->getCartesian_rank() == 0)
    cout << "CG Initial error: " << initial_error << endl;
  if (initial_error < 1E-16)
    return (true);
  while (i < maxit) {
    (field->*FunctionImage) (z, v, grid, vct);
    t = c / dotP(v, z, xkrylovlen);
    // x(i+1) = x + t*v
    addscale(t, xkrylov, v, xkrylovlen);
    // r(i+1) = r - t*z
    addscale(-t, r, z, xkrylovlen);
    d = dotP(r, r, xkrylovlen);
    if (CGVERBOSE && vct->getCartesian_rank() == 0)
      cout << "Iteration # " << i << " - norm of residual relative to initial error " << sqrt(d) / initial_error << endl;
    if (sqrt(d) < tol * initial_error) {
      // if (d < tol){

      if (vct->getCartesian_rank() == 0)
        cout << "CG converged at iteration # " << i << " with error " << sqrt(d) << endl;
      CONVERGED = true;
      break;
    }
    else if (sqrt(d) > 10E8 * initial_error) {
      if (vct->getCartesian_rank() == 0) {
        cerr << "CG not converging" << endl;
        cerr << "CG stopped" << endl;

      }
      CONVERGED = false;
      return (CONVERGED);
      break;

    }

    // calculate the new v
    addscale(1, d / c, v, r, xkrylovlen);
    c = d;
    i++;

  }
  if (i == maxit)
    cout << "CG not converged after " << maxit << " iterations" << endl;
  // deallocate
  delete[]r;
  delete[]im;
  delete[]v;
  delete[]z;
  return (CONVERGED);

}

inline bool CG(double *xkrylov, int xkrylovlen, double *b, int maxit, double tol, GENERIC_IMAGE FunctionImage, Grid * grid, VirtualTopology3D * vct) {
  // allocate residual, image, p, b, calculated on central points
  double *r = new double[xkrylovlen];
  double *v = new double[xkrylovlen];
  double *z = new double[xkrylovlen];
  double *im = new double[xkrylovlen];
  double c, t, d, initial_error;
  bool CONVERGED = false;
  int i = 0;
  bool CGVERBOSE = true;
  // initial guess for x: all the components are equal to 0
  eqValue(0, xkrylov, xkrylovlen);
  // Compute r = b -Ax
  (*FunctionImage) (im, xkrylov, grid, vct);
  sub(r, b, im, xkrylovlen);
  // v = r
  eq(v, r, xkrylovlen);
  c = dotP(r, r, xkrylovlen);
  initial_error = sqrt(c);
  if (vct->getCartesian_rank() == 0)
    cout << "Initial error: " << initial_error << endl;
  // for i=0,1,..., until convergence
  if (CGVERBOSE && vct->getCartesian_rank() == 0) {
    cout << "------------------------------------" << endl;
    cout << "-               CG                 -" << endl;
    cout << "------------------------------------" << endl;
    cout << endl;
  }
  while (i < maxit) {
    (*FunctionImage) (z, v, grid, vct);

    t = c / dotP(v, z, xkrylovlen);
    // x(i+1) = x + t*v
    addscale(t, xkrylov, v, xkrylovlen);
    // r(i+1) = r - t*z
    addscale(-t, r, z, xkrylovlen);
    d = dotP(r, r, xkrylovlen);

    if (CGVERBOSE && vct->getCartesian_rank() == 0)
      cout << "Iteration # " << i << " - norm of residual relative to initial error " << sqrt(d) / initial_error << endl;
    if (sqrt(d) < tol * initial_error) {
      if (vct->getCartesian_rank() == 0)
        cout << "CG converged at iteration # " << i << endl;
      CONVERGED = true;
      break;
    }
    else if (sqrt(d) > 10E8 * initial_error) {
      if (vct->getCartesian_rank() == 0) {
        cerr << "CG not converging" << endl;
        cerr << "CG stopped" << endl;
      }
      CONVERGED = false;
      break;
    }
    addscale(1, d / c, v, r, xkrylovlen);
    c = d;
    i++;

  }
  // deallocate
  delete[]r;
  delete[]im;
  delete[]v;
  delete[]z;
  return (CONVERGED);


}
#endif
