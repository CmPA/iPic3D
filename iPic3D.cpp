
#include <iomanip>
#include "iPic3D.h"

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

  for (int i = KCode.FirstCycle()+1; i <= KCode.LastCycle(); i++) {
    /*! mlmd: KCode.get_myrank() is on the local grid communicator */
    if (KCode.get_myrank() == 0) cout << " ======= Grid " << KCode.get_numGrid()  <<"  Cycle " << i << " ======= " << endl;

    /* ----------------------------------------------------- */
    /* 2- Calculate fields and move particles                */
    /*    Exit if there is a memory issue with the particles */
    /* ----------------------------------------------------- */

    KCode.UpdateCycleInfo(i);
    KCode.CalculateField();

    b_err = KCode.ParticlesMover(i); 

    // with mlmd ops
    if (!b_err) KCode.GatherMoments();

    if (!b_err) KCode.CalculateBField(i); 

    if ( b_err) i = KCode.LastCycle() + 1;

    /* --------------- */
    /* 3- Output files */
    /* --------------- */

    
    KCode.WriteOutput(i);
    KCode.WriteConserved(i);
    KCode.WriteRestart(i);
    
  }

  KCode.Finalize();

  return 0;
}
