
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
  KCode.WriteOutput(0);
  KCode.WriteConserved(0);

  /* ------------ */
  /* 1- Main loop */
  /* ------------ */

  for (int i = KCode.FirstCycle()+1; i <= KCode.LastCycle()+1; i++) {

    if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

    // Push particle position half a step
    KCode.ParticlesMover();

    // Calculate E and B fields
    KCode.CalculateField();

    // Final position/velocity update, moment gathering for output
    KCode.FinalMomentumUpdate();
    KCode.ParticlesMover();

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
