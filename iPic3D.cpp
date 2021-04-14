
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
  KCode.GatherMoments();
  KCode.WriteOutput(KCode.FirstCycle());
  

  KCode.WriteConserved(0);
  KCode.WriteChunks(0);
  /* ------------ */
  /* 1- Main loop */
  /* ------------ */
  for (int i = KCode.FirstCycle()+1; i <= KCode.LastCycle(); i++) {

    if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

    /* ----------------------------------------------------- */
    /* 2- Calculate fields and move particles                */
    /*    Exit if there is a memory issue with the particles */
    /* ----------------------------------------------------- */

    KCode.UpdateCycleInfo(i);
    KCode.CalculateField(i);

    //cout << "*** out of field ***" << endl;

    b_err = KCode.ParticlesMover();

    if (!b_err) KCode.CalculateBField();
    if (!b_err) KCode.GatherMoments();
    if ( b_err) i = KCode.LastCycle() + 1;

    /* --------------- */
    /* 3- Output files */
    /* --------------- */
    KCode.WriteOutput(i);
    KCode.WriteConserved(i);
    KCode.WriteRestart(i);

    KCode.WriteChunks(i);
  }

  KCode.Finalize();

  return 0;
}
