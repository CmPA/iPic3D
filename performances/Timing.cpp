
#include "Timing.h"
#include "ipicdefs.h"

/**
 * 
 * series of methods for timing and profiling PARSEK 
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */

/** default constructor */
Timing::Timing() {
}

/** constructor with the initialization of the log file */
Timing::Timing(int my_rank) {
  rank_id = my_rank;
  // initialize the logger
  // MPE_Init_log();
  // start the timer
  startTiming();
  // get event ID, from MPE
  // event1a = MPE_Log_get_event_number();
  // event1b = MPE_Log_get_event_number();
  // event2a = MPE_Log_get_event_number();
  // event2b = MPE_Log_get_event_number();
  // event3a = MPE_Log_get_event_number();
  // event3b = MPE_Log_get_event_number();
  // described the events
  // if (my_rank==0){
  // MPE_Describe_state(event1a,event1b,"Mover","red"); // the mover is red in the visualizer
  // MPE_Describe_state(event2a,event2b,"Field","blue"); // the mover is blue in the visualizer
  // MPE_Describe_state(event3a,event3b,"Interp P->G","yellow"); // the interpolation particle->Grid is yellow in the visualizer
  // }
  former_MPI_Barrier(MPI_COMM_WORLD);
  // start the log
  // MPE_Start_log();

}

/** start the timer */
void Timing::startTiming() {
  ttick = MPI_Wtick();
  former_MPI_Barrier(MPI_COMM_WORLD);
  tstart = MPI_Wtime();
}
/** stop the timer */
void Timing::stopTiming() {
  former_MPI_Barrier(MPI_COMM_WORLD);
  tend = MPI_Wtime();
  texecution = tend - tstart;
  if (rank_id == 0) {
    cout << endl;
    cout << endl;
    cout << "*** SIMULATION ENDED SUCESSFULLY ***" << endl;
    cout << " PARSEK Simulation Time: " << texecution << " sec" << " (" << texecution / 3600 << " hours)" << endl;
    cout << "***" << endl;
    cout << endl;
  }
  // close the log file
  // MPE_Finish_log("PARSEK_LOG");
}

/** start timing the mover */
void Timing::start_mover() {
  // MPE_Log_event(event1a,0,"start mover");
}
/** stop timing the mover */
void Timing::stop_mover() {
  // MPE_Log_event(event1b,0,"end mover");
}
/** start timing the field solver */
void Timing::start_field() {
  // MPE_Log_event(event2a,0,"start Field solver");
}
/** stop timing the field solver */
void Timing::stop_field() {
  // MPE_Log_event(event2b,0,"stop Field solver");
}
/** start timing the interpolation Particle -> Grid */
void Timing::start_interpP2G() {
  // MPE_Log_event(event3a,0,"start interpolation");
}
/** stop timing the interpolation Particle -> Grid*/
void Timing::stop_interpP2G() {
  // MPE_Log_event(event3b,0,"stop interpolation");
}
/** get the elapsed time from start_timng and stop_timing */
double Timing::getExecutionTime() {
  return (texecution);
}
/** print to screen the elapsed time */
void Timing::Print() {
  cout << "Execution Time: " << texecution << " sec" << " (" << texecution / 3600 << " hours)" << endl;
}
/** print to screen the elapsed time from t_start to the call to print function*/
void Timing::Print_OnAir() {
  cout << "Execution Time: " << MPI_Wtime() - tstart << " sec" << endl;
}
