/*******************************************************************************************
  Timing.h  -  series of methods for timing and profiling PARSEK 
  -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
 ********************************************************************************************/

#ifndef TIMING_H
#define TIMING_H

#include "MPIdata.h"

#include <iostream>

using std::cout;
using std::endl;
/**
 * 
 * series of methods for timing and profiling PARSEK 
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */
class Timing {
public:
  /** default constructor */
  Timing();
  /** default constructor */
  Timing(int my_rank);
  /** start timing the mover */
  void start_mover();
  /** stop timing the mover */
  void stop_mover();
  /** start timing the field solver */
  void start_field();
  /** stop timing the field solver */
  void stop_field();
  /** start timing the interpolation Particle -> Grid */
  void start_interpP2G();
  /** stop timing the interpolation Particle -> Grid */
  void stop_interpP2G();
  /** start the timer */
  void startTiming();
  /** stop the timer */
  void stopTiming();
  /** get the execution time */
  double getExecutionTime();
  /** print the execution time */
  void Print();
  /** print the execution time without stopping the clock*/
  void Print_OnAir();

private:
  /** rank of the processor */
  int rank_id;
  /** wall-clock time when startTiming method is called */
  double tstart;
  /** wall-clock time when endTiming method is called */
  double tend;
  /** execution time */
  double texecution;
  /** time precision */
  double ttick;
  /** events to write the logging: particle mover, field solver. it can be extended with other events */
  int event1a, event1b, event2a, event2b, event3a, event3b, event4a, event4b;

};

#endif
