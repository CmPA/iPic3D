#include "PetscSolver.h"


#ifdef __PETSC_SOLVER__

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

PetscSolver::PetscSolver(EMfields3D *EMf, Grid* grid, VirtualTopology3D *vct, Collective *col)
{
  PetscInitialize(NULL, NULL, NULL, NULL);

  ctx.EMf  = EMf;
  ctx.grid = grid;
  ctx.vct = vct;
  ctx.col = col;

  int nxn = grid->getNXN();
  int nyn = grid->getNYN();
  int nzn = grid->getNZN();
  dim = 3*(nxn-2)*(nyn-2)*(nzn-2);

  double tol = col->getGMREStol();

  // Create the coefficient matrix
  MatCreateShell(PETSC_COMM_WORLD ,dim ,dim, PETSC_DETERMINE, PETSC_DETERMINE, &ctx, &A);
  MatShellSetOperation(A, MATOP_MULT, (void(*)(void))MyMatMult);
  MatSetFromOptions(A);  /*must be called*/

  // Create the vectors
  VecCreate(PETSC_COMM_WORLD, &globX);
  VecSetSizes(globX, dim,  PETSC_DECIDE);
  VecSetFromOptions(globX);
  VecDuplicate(globX, &globB);

  // Create the Krylov subspace solver
  KSPCreate(PETSC_COMM_WORLD, &ksp);
  KSPSetType(ksp, KSPGMRES);
  KSPSetOperators(ksp, A, A);
  
  // Preconditioner
  //KSPGetPC(ksp, &pc);
  //PCSetType(pc, PCJACOBI);
  
  // Tolerance and number of iterations
  KSPSetTolerances(ksp,1e-20,tol,10000,2000);

  KSPSetFromOptions(ksp);

  //finish the assemble of A and globX
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (A, MAT_FINAL_ASSEMBLY);

  VecAssemblyBegin(globX);
  VecAssemblyEnd  (globX);
}

PetscSolver::~PetscSolver()
{
  VecDestroy(&globX);
  VecDestroy(&globB);
  MatDestroy(&A);
  KSPDestroy(&ksp);
  PetscFinalize();
}

void PetscSolver::solveE()
{
  if (ctx.vct->getCartesian_rank() == 0)
    cout << "*** MAXWELL SOLVER (PETSc) ***" << endl;

  double *d_b;
  VecGetArray(globB, &d_b);
  ctx.EMf->MaxwellSource(d_b, ctx.grid, ctx.vct, ctx.col);
  VecRestoreArray(globB, &d_b);

  /*finish the assemble of X and B, no modification allowed*/
  VecAssemblyBegin(globB);
  VecAssemblyEnd  (globB);

  ctx.EMf->startEcalc(ctx.grid, ctx.vct, ctx.col);
  KSPSolve(ksp,globB,globX);
  double *d_x;
  VecGetArray(globX, &d_x);
  ctx.EMf->endEcalc(d_x, ctx.grid, ctx.vct, ctx.col);
  VecRestoreArray(globX, &d_x);

  PetscInt iter;
  PetscScalar residual;
  KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
  KSPGetIterationNumber(ksp,&iter);
  KSPGetResidualNorm(ksp, &residual);
  PetscPrintf(PETSC_COMM_WORLD,"Solver converges at iteration %d with error %10.3e\n"
                              ,iter,residual);

}


/*****************************/
/*  FUNCTION TO COMPUTE A*x  */
/*****************************/
/* This function does exaclty the same than MaxwellImage(...) from EMfields3D */
PetscErrorCode MyMatMult(Mat A,Vec x,Vec y) 
{
  CtxSolver* ctx;
  MatShellGetContext(A, &ctx);
  

  double *d_y, *d_x;
  VecGetArray(y, &d_y);
  VecGetArray(x, &d_x);
  ctx->EMf->MaxwellImage(d_y, d_x, ctx->grid, ctx->vct); 

  VecRestoreArray(x, &d_x);
  VecRestoreArray(y, &d_y);

  return 0;
}

#endif
