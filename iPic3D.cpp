
#include "iPic3D.h"

using namespace iPic3D;

int main(int argc, char **argv) {

  iPic3D::c_Solver KCode;
  bool b_err = false;

  KCode.Init(argc, argv);

  for (int i = KCode.FirstCycle(); i < KCode.LastCycle(); i++) {

    if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

    KCode.WriteOutput(i);

    if (!b_err) {
      KCode.CalculateField();
      b_err = KCode.ParticlesMover();
    }
    else {
      i = KCode.LastCycle() + 1;
    }

    KCode.WriteConserved(i);
    KCode.WriteRestart(i);

  }

  KCode.Finalize();

  return 0;
}
