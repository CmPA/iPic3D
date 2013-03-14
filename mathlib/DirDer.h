/***************************************************************************
  DirDer.h  - Finite Difference Directional Derivative
  -------------------
begin             : Wed Jul 14 2004
copyright         : (C) 2004 Los Alamos National Laboratory
developers        : Stefano Markidis, Giovanni Lapenta
email             : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/


#ifndef DirDer_H
#define DirDer_H

#include <math.h>
#include "Basic.h"
using std::cout;
using std::cerr;
using std::endl;

typedef void (Field::*fieldFUNCTION) (double *, double *, Grid *, VirtualTopology *);
typedef void (*genericFUNCTION) (double *, double *, Grid *, VirtualTopology *);

/**
 *  DirDer  - Finite Difference Directional Derivative
 *  Appximate f'(x)w
 *
 * @date Fri Jun 4 2004
 * @par Copyright:
 * (C) 2004 Los Alamos National Laboratory
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 1.0
 */
inline void DirDer(double *dirder, double *x, int xlen, double *w, double *f0, fieldFUNCTION FunctionImage, Grid * grid, VirtualTopology * vct, Field * field) {
  double epsnew = 1E-7;
  double norm;
  norm = normP(w, xlen);
  if (norm > 0) {

    double *del = new double[xlen];
    double *f1 = new double[xlen];
    epsnew /= norm;
    norm = normP(x, xlen);
    if (norm > 0)
      epsnew *= norm;
    for (int i = 0; i < xlen; i++)
      *del++ = *x++ + epsnew * (*w++);
    (field->*FunctionImage) (f1, del, grid, vct);
    sub(dirder, f1, f0, xlen);
    scale(dirder, 1.0 / epsnew, xlen);
  }
  else {
    eqValue(0.0, dirder, xlen);

  }

}
inline void DirDer(double *dirder, double *x, int xlen, double *w, double *f0, genericFUNCTION FunctionImage, Grid * grid, VirtualTopology * vct) {
  double epsnew = 1E-7;
  double norm;
  norm = normP(w, xlen);
  if (norm == 0) {
    eqValue(0.0, dirder, xlen);
  }
  else {
    double *del = new double[xlen];
    double *f1 = new double[xlen];
    epsnew /= norm;
    norm = normP(x, xlen);
    if (norm > 0)
      epsnew *= norm;
    for (int i = 0; i < xlen; i++)
      *del++ = *x++ + epsnew * (*w++);
    (FunctionImage) (f1, del, grid, vct);
    sub(dirder, f1, f0, xlen);
    scale(dirder, 1.0 / epsnew, xlen);
    delete[]del;
    delete[]f1;
  }
}
#endif
