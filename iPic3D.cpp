
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

  for (int i = KCode.FirstCycle()+1; i <= KCode.LastCycle(); i++) {
    /*! mlmd: KCode.get_myrank() is on the local grid communicator */
    if (KCode.get_myrank() == 0) cout << " ======= Grid " << KCode.get_numGrid()  <<"  Cycle " << i << " ======= " << endl;

    /* ----------------------------------------------------- */
    /* 2- Calculate fields and move particles                */
    /*    Exit if there is a memory issue with the particles */
    /* ----------------------------------------------------- */

    KCode.UpdateCycleInfo(i);

    clocks->start(1);
    KCode.CalculateField();
    clocks->stop(1);

    clocks->start(2);
    b_err = KCode.ParticlesMover(i); 
    clocks->stop(2);

    // with mlmd ops
    clocks->start(3);
    if (!b_err) KCode.GatherMoments();
    clocks->stop(3);

    clocks->start(4);
    if (!b_err) KCode.CalculateBField(i); 
    clocks->stop(4);

    if ( b_err) i = KCode.LastCycle() + 1;

    /* --------------- */
    /* 3- Output files */
    /* --------------- */

    clocks->start(5);
    KCode.WriteOutput(i);
    KCode.WriteConserved(i);
    KCode.WriteRestart(i);
    clocks->stop(5);
    
  }

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  KCode.Finalize();
  // Profiling
#ifdef __PROFILING__ // add more
  if (myrank == 0) {
    cout << "##############################################" <<endl
	 << "Initialization                      : " << clocks->get(0) << endl
	 << "Calculate Field                     : " << clocks->get(1) << endl
	 << "Particle Mover                      : " << clocks->get(2) << endl
	 << " --- moving particle                : " << clocks->get(6) << endl
	 << " --- repopulating particle          : " << clocks->get(7) << endl
	 << "     --- communicate after mover    : " << clocks->get(8) << endl
         << "     --- send PBC                   : " << clocks->get(9) << endl 
	 << "         --- particle packing       : " << clocks->get(12) << endl
	 << "         --- resizing ops           : " << clocks->get(13) << endl
	 << "         --- sending particles      : " << clocks->get(14) << endl
	 << "     --- receive PBC                : " << clocks->get(10) << endl
	 << "         --- resizing ops           : " << clocks->get(15) << endl 
	 << "         --- receiving ops          : " << clocks->get(16) << endl
	 << "         --- apply PBC              : " << clocks->get(17) << endl
	 << "     --- parent/ child barrier      : " << clocks->get(11) << endl
	 << "Moment Gathering                    : " << clocks->get(3) << endl
	 << "B calculation                       : " << clocks->get(4) << endl
	 << "Output                              : " << clocks->get(5) << endl
	 << "##############################################" <<endl;
  }
  #else
  if (myrank == 0) {
    cout << "##############################################" <<endl
	 << "Initialization             : " << clocks->get(0) << endl
	 << "Calculate Field            : " << clocks->get(1) << endl
	 << "Particle Mover             : " << clocks->get(2) << endl
	 << "Moment Gathering           : " << clocks->get(3) << endl
	 << "B calculation              : " << clocks->get(4) << endl
	 << "Output                     : " << clocks->get(5) << endl
	 << "##############################################" <<endl;
  }
  #endif


  return 0;
}
