
#include <iomanip>
#include "iPic3D.h"

#include "MonteCarlo.h"

int main(int argc, char **argv) {

  c_Solver KCode;
  bool b_err = false;

  /* ------------------------------ */
  /* 0- Initialize the solver class */
  /* ------------------------------ */

  KCode.Init(argc, argv);
  KCode.InjectBoundaryParticles();
  KCode.GatherMoments(-1);

  MonteCarlo *MC= new MonteCarlo(&KCode);

  /* ------------ */
  /* 1- Main loop */
  /* ------------ */

  for (int i = KCode.FirstCycle(); i <= KCode.LastCycle(); i++) {

    if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

    /* ----------------------------------------------------- */
    /* 2- Calculate fields and move particles                */
    /*    Exit if there is a memory issue with the particles */
    /* ----------------------------------------------------- */

    KCode.UpdateCycleInfo(i);
    KCode.CalculateField();

    b_err = KCode.ParticlesMover();

    if (MC->getMonteCarloPlugIn()== "yes"){
      MC->SelectCollidingParticles(&KCode);
      MC->MCWrapper(&KCode, i);
    }

    if (!b_err) KCode.CalculateBField();
    if (!b_err) KCode.GatherMoments(i);
    if ( b_err) i = KCode.LastCycle() + 1;

    /* --------------- */
    /* 3- Output files */
    /* --------------- */

    // MonteCarlo diagnostics - do this BEFORE printing                    
    if (MC->getMonteCarloPlugIn()== "yes"){
      KCode.WriteCollisionDiagnostics(MC, i);
    }

    KCode.WriteOutput(i, MC);
    KCode.WriteConserved(i);
    KCode.WriteRestart(i);

  }

  KCode.Finalize();

  return 0;
}
