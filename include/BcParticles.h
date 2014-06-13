/***************************************************************************
  BcParticles.h  -  Library to manage boundary conditions for particles
  -------------------
begin                : Fri Jun 4 2004
copyright            : (C) 2004 Los Alamos National Laboratory
developers           : Stefano Markidis, Giovanni Lapenta
email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef BcParticles_H
#define BcParticles_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "VirtualTopology3D.h"
#include "Basic.h"

void BCpart(double *x, double *u, double Lx, double ut, int bcFaceXright, int bcFaceXleft);
void BCpart(double *x, double *u, double *v, double *w, double Lx, double ut, double vt, double wt, int bcFaceXright, int bcFaceXleft);
void BCpart(double *x, double *y, double *z, double *u, double *v, double *w, double Lx, double Ly, double Lz, double ut, double vt, double wt, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, int bcFaceZright, int bcFaceZleft, VirtualTopology3D * vct);

#endif
