#ifndef GMRES_new2_H
#define GMRES_new2_H

#include <iostream>
#include <math.h>

#include "Alloc.h"
#include "Basic.h"
#include "EMfields3D.h"
#include "VirtualTopology3D.h"

using std::cout;
using std::cerr;
using std::endl;

class EMfields3D;
typedef void (EMfields3D::*FIELD_IMAGE) (double *, double *, Grid *, VirtualTopology3D *);
typedef void (EMfields3D::*FIELD_IMAGE_ECSIM) (double *, double *, Grid *, VirtualTopology3D *,Collective *col );
typedef void (*GENERIC_IMAGE) (double *, double *, Grid *, VirtualTopology3D *);

typedef void (EMfields3D::*FIELD_IMAGE_2) (double *, double *, Grid *, VirtualTopology3D *, int, int, int);

void GMRES(FIELD_IMAGE FunctionImage, double *xkrylov, int xkrylovlen, double *b, int m, int max_iter, double tol, Grid * grid, VirtualTopology3D * vct, EMfields3D * field);
void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn);

void GMRES_ECSIM(FIELD_IMAGE_ECSIM FunctionImage, double *xkrylov, int xkrylovlen, double *b, int m, int max_iter, double tol, Grid * grid, VirtualTopology3D * vct, Collective *col, EMfields3D * field);

/** to use in Poisson correction per face **/
void GMRES_FACE(FIELD_IMAGE_2 FunctionImage, double *xkrylov, int xkrylovlen, double *b, int m, int max_iter, double tol, Grid * grid, VirtualTopology3D * vct, EMfields3D * field, MPI_Comm COMM, int nxc, int nyc, int nzc);

#endif
