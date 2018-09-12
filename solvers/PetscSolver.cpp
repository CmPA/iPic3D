#include "PetscSolver.h"


/******************************************************************************/
/*   Field Solver using PETSc.                                                */
/*                                                                            */
/* Matrix-free method                                                         */
/* It uses the original routines MaxwellImage and MaxwellSources              */
/* in EMfields3D.                                                             */
/*                                                                            */
/* Diego Gonzalez                                                             */
/* May 2016                                                                   */
/*                                                                            */
/******************************************************************************/

PetscSolver::PetscSolver(EMfields3D *EMf, Grid* grid, VirtualTopology3D *vct, Collective *col, Particles3D* part)
{
  
  int na = 2;
  char** args = new char*[na];
  for (int i=0; i<na; i++) args[i] = new char[128];


  sprintf(args[0], "");
  sprintf(args[1], "-snes_mf");

  PetscInitialize(&na,&args , NULL, NULL);
  //PetscInitialize(NULL, NULL , NULL, NULL);
  //
  ctx.EMf  = EMf;
  ctx.grid = grid;
  ctx.vct  = vct;
  ctx.col  = col;
  ctx.part = part;

  int nxn = grid->getNXN();
  int nyn = grid->getNYN();
  int nzn = grid->getNZN();
  dim = 3*(nxn-2)*(nyn-2)*(nzn-2);

  tolRel = col->getGMREStol();
  tolAbs = 1e-14;
  iter_max = 5000;

  SNESCreate(PETSC_COMM_WORLD, &snes);
  SNESSetType(snes, "qn");
  SNESSetTolerances(snes,tolAbs,tolAbs,tolAbs,5000,5000);
  
  // Create the vectors
  VecCreate(PETSC_COMM_WORLD, &globX);
  VecSetSizes(globX, dim,  PETSC_DECIDE);
  VecSetFromOptions(globX);
  VecDuplicate(globX, &globR);

  MatCreate(PETSC_COMM_WORLD,&J);
  MatSetSizes(J,dim, dim, PETSC_DECIDE,PETSC_DECIDE);
  MatSetFromOptions(J);
  MatSetUp(J);

  SNESSetFunction(snes,globR,FormFunction, &ctx);
  //SNESSetJacobian(snes,J,J,FormJacobian,NULL);

  // Preconditioning
  SNESGetKSP(snes,&ksp);
  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCNONE);
  KSPSetTolerances(ksp,1.e-14,1.e-14,PETSC_DEFAULT,50);

  SNESSetFromOptions(snes);

  
  // Preconditioner
  //KSPGetPC(ksp, &pc);
  //PCSetType(pc, PCJACOBI);
  
  // Tolerance and number of iterations
  //KSPSetTolerances(ksp,tolRel,tolAbs,10000,iter_max);

  //KSPSetFromOptions(ksp);

  //finish the assemble of A and globX
  //MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  //MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);
  //


}

PetscSolver::~PetscSolver()
{
  VecDestroy(&globX);
  VecDestroy(&globR);
  //MatDestroy(&A);
  //KSPDestroy(&ksp);
  PetscFinalize();
}

void PetscSolver::solve()
{
  if (ctx.vct->getCartesian_rank() == 0)
    cout << "*** MAXWELL SOLVER (PETSc) ***" << endl;

  ctx.EMf->MaxwellSource(ctx.grid, ctx.vct, ctx.col);

  ctx.EMf->startEcalc(ctx.grid, ctx.vct, ctx.col);
  SNESSolve(snes,NULL,globX);
  double *d_x;
  VecGetArray(globX, &d_x);
  ctx.EMf->endEcalc(d_x, ctx.grid, ctx.vct, ctx.col);
  VecRestoreArray(globX, &d_x);

  //PetscInt iter;
  //PetscScalar residual;
  SNESView(snes,PETSC_VIEWER_STDOUT_WORLD);
  //KSPGetIterationNumber(ksp,&iter);
  //KSPGetResidualNorm(ksp, &residual);
  //PetscPrintf(PETSC_COMM_WORLD,"Solver converges at iteration %d with error %10.3e\n"
                              //,iter,residual);
  //Vec f;
  //VecView(globX,PETSC_VIEWER_STDOUT_WORLD);
  //SNESGetFunction(snes,&f,0,0);
  //VecView(globR,PETSC_VIEWER_STDOUT_WORLD);
  
}


/*****************************/
/*  FUNCTION TO COMPUTE A*x  */
/*****************************/
/* This function does exaclty the same than MaxwellImage(...) from EMfields3D */
PetscErrorCode FormFunction(SNES snes, Vec x, Vec r, void* _ctx) 
{
  CtxSolver* ctx = (CtxSolver*) _ctx;

  double *d_r, *d_x;
  VecGetArray(r, &d_r);
  VecGetArray(x, &d_x);
  ctx->EMf->MaxwellImage(d_r, d_x, ctx->grid, ctx->vct, ctx->part); 

  VecRestoreArray(x, &d_x);
  VecRestoreArray(r, &d_r);

  return 0;
}

PetscErrorCode FormJacobian(SNES snes,Vec x,Mat jac,Mat B,void *dummy)
{

  MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
  if (jac != B) {
    MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);
  }
  return 0;
}
