
#ifndef PETSCSOLVER_H
#define PETSCSOLVER_H

#include <petscdmda.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>
#include "Collective.h"
#include "Grid.h"
#include "EMfields3D.h"
#include "Particles3D.h"
#include "VirtualTopology3D.h"

typedef struct {
    Grid *grid;
    VirtualTopology3D *vct;
    EMfields3D *EMf;
    Collective *col;
    Particles3D *part;
  } CtxSolver;

class PetscSolver {
public:
  PetscSolver(EMfields3D *EMf, Grid* grid, VirtualTopology3D *vct, Collective *col, Particles3D* part);
  ~PetscSolver();

  void solve();

private:
  // Dimension of the system
  int dim;
  // Global vectors
  Vec globX, globR;
  // Jacobian Matrix 
  Mat J;
  // Solver
  KSP ksp;
  SNES snes;
  // Preconditioner
  PC pc;
  // Context for the A*x product
  CtxSolver ctx;
  // Tolerances
  double tolRel, tolAbs;
  // Maximum number of iterations
  int iter_max;
};

#endif

extern PetscErrorCode FormFunction(SNES snes, Vec x, Vec r, void* ctx);
extern PetscErrorCode FormJacobian(SNES snes,Vec x,Mat jac,Mat B,void *dummy);

