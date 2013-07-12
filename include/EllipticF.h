/*******************************************************************************************
EllipticF.h  -  Elliptic Functions Tool
                            -------------------
developers: Found on the web
********************************************************************************************/

#ifndef ELLIPTIC_H
#define ELLIPTIC_H

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

double drf(double x, double y, double z, int *piErr);
double drd(double x, double y, double z, int *piErr);

#define EllipticF(k,ierr) drf(0.0,1.0-pow(k,2),1.0,&ierr)

#define EllipticE(k,ierr) (drf(0.0,1.0-pow(k,2),1.0,&ierr)-(pow(k,2)/3.0)*drd(0.0,1.0-pow(k,2),1.0,&ierr))

#endif
