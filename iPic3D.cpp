
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

  /* ------------ */
  /* 1- Main loop */
  /* ------------ */

  for (int i = KCode.FirstCycle(); i <= KCode.LastCycle(); i++) {

    if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

    // Update particle position
    KCode.ParticlesMover();

    // Calculate E and B fields
    KCode.CalculateField();

    // Update the momentum and velocity
    KCode.FinalMomentumUpdate();

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
