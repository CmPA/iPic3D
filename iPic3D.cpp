
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

//  clocks->start(0); // Activating this clock causes SegFault...
  KCode.Init(argc, argv);
//  clocks->stop(0);
//  KCode.InjectBoundaryParticles();
//  KCode.GatherMoments();
  KCode.WriteOutput(KCode.FirstCycle());
  KCode.WriteConserved(KCode.FirstCycle());

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

    KCode.PushParticles();

    clocks->start(1);
    if (!b_err) KCode.GatherMoments();
    clocks->stop(1);

    clocks->start(2);
    KCode.CalculateField();
    clocks->stop(2);

    clocks->start(3);
    b_err = KCode.ParticlesMover();
    clocks->stop(3);

    if (!b_err) KCode.CalculateBField();
    if ( b_err) i = KCode.LastCycle() + 1;

    /* --------------- */
    /* 3- Output files */
    /* --------------- */

    clocks->start(4);
    KCode.WriteOutput(i);
    KCode.WriteConserved(i);
    KCode.WriteRestart(i);
    clocks->stop(4);

//    if (i == 0 || i%(10)==0) {
//      if (KCode.get_myrank() == 0) {
//        std::cout << "################################################" << std::endl
//                  << "Initialization                : " << clocks->get(0) << " s" << std::endl
//                  << "Moments Gathering             : " << clocks->get(1) << " s" << std::endl
//                  << "Field Calculation             : " << clocks->get(2) << " s" << std::endl
//                  << "Particle Mover                : " << clocks->get(3) << " s" << std::endl
//                  << "Writing                       : " << clocks->get(4) << " s" << std::endl
//                  << "**************************************************" << std::endl;
//        }
//      }
  }

  KCode.Finalize();

  if (KCode.get_myrank() == 0) {
    std::cout << "################################################" << std::endl
              << "Initialization                : " << clocks->get(0) << " s" << std::endl
              << "Moments Gathering             : " << clocks->get(1) << " s" << std::endl
              << "Field Calculation             : " << clocks->get(2) << " s" << std::endl
              << "Particle Mover                : " << clocks->get(3) << " s" << std::endl
              << "Writing                       : " << clocks->get(4) << " s" << std::endl
              << "----------------------------------------------------------" << std::endl
              << "Total                         : " << clocks->get(5)<< " s" << std::endl;
  }
  return 0;
}
