
#include <iomanip>
#include "iPic3D.h"
#include "MyClock.h"

MyClock *clocks;

using namespace iPic3D;

int main(int argc, char **argv) {

  iPic3D::c_Solver KCode;
  bool b_err = false;

  
  /* ------------------------------ */
  /* 0- Initialize the solver class */
  /* ------------------------------ */

  KCode.Init(argc, argv);
    
  KCode.InjectBoundaryParticles();
  // without mlmd ops
  KCode.GatherMoments_Init();

  // added to compare with ECSIM        
  KCode.WriteOutput(KCode.FirstCycle());
  KCode.WriteConserved(KCode.FirstCycle());
  /* ------------ */
  /* 1- Main loop */
  /* ------------ */

  // otherwise restart skips a cycle (the +1 in the for)
  if (KCode.FirstCycle()>0){
    KCode.DecrementFirstCycle();
  }
  
  for (int i = KCode.FirstCycle()+1; i <= KCode.LastCycle()+1; i++) {
    /*! mlmd: KCode.get_myrank() is on the local grid communicator */
    if (KCode.get_myrank() == 0) cout << " ======= Grid " << KCode.get_numGrid()  <<"  Cycle " << i << " ======= " << endl;

    /* ----------------------------------------------------- */
    /* 2- Calculate fields and move particles                */
    /*    Exit if there is a memory issue with the particles */
    /* ----------------------------------------------------- */

    KCode.UpdateCycleInfo(i, KCode.FirstCycle());

    KCode.CalculateField();


    /* this is the normal and tested */
    /*
    b_err = KCode.ParticlesMover(i); 

    // with mlmd ops
    if (!b_err) KCode.GatherMoments();
    */
    /* end this is the normal and tested  */


    KCode.Mover_GatherMoments_Interleaved(i, KCode.FirstCycle());

    if (!b_err) KCode.CalculateBField(i); 

    if ( b_err) i = KCode.LastCycle() + 1;

    /* --------------- */
    /* 3- Output files */
    /* --------------- */


    KCode.WriteOutput(i);
    KCode.WriteConserved(i);
    KCode.WriteRestart(i);

    
  }

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  KCode.Finalize();

  return 0;
}
