
#include <mpi.h>
#include <stdarg.h>
#include "TimeTasks.h"
#include "asserts.h"
#include "MPIdata.h" // for get_rank

/** implementation of declarations in utility/TimeTasks.h **/

TimeTasks timeTasks;

void TimeTasks::resetCycle()
{
  for(int e=0;e<LAST;e++)
  {
    //compute[e]=0.;
    start_times[e]=0.;
    task_duration[e]=0.;
    communicate[e]=0.;
  }
  active_task=NONE;
  active_mode=COMPUTATION;
  t_start_communicate = 0.;
}
void TimeTasks::start(int taskid)
{
  assert_eq(active_task+1,taskid);
  active_task = taskid;
  double now = MPI_Wtime();
  start_times[active_task] = now;
}
void TimeTasks::end(int taskid)
{
  assert_eq(taskid,active_task);
  double now = MPI_Wtime();
  task_duration[active_task] = now - start_times[active_task];
  compute[active_task] = task_duration[active_task]-communicate[active_task];
}
void TimeTasks::start_communicate()
{
  if(!active_task) return;
  assert_eq(active_mode,COMPUTATION);
  t_start_communicate = MPI_Wtime();
  active_mode=COMMUNICATION;
}
void TimeTasks::addto_communicate()
{
  if(!active_task) return;
  assert_eq(active_mode,COMMUNICATION);
  assert_ne(t_start_communicate,0.);
  communicate[active_task] += MPI_Wtime()-t_start_communicate;
  t_start_communicate = 0.;
  active_mode=COMPUTATION;
}
#define TIMING_PREFIX "| "
void TimeTasks::print_cycle_times()
{
  // we could report average for all processes
  if(!MPIdata::get_rank())
  {
    fflush(stdout);
    fprintf(stdout,"=== times for cycle for rank %d === \n",
      MPIdata::get_rank());
    fprintf(stdout, TIMING_PREFIX
      "moms flds pcls Bfld cycl\n");
    fprintf(stdout, TIMING_PREFIX
      "%4.2f "
      "%4.2f "
      "%4.2f "
      "%4.2f "
      "%4.2f (total time)\n",
      get_time(TimeTasks::MOMENTS),
      get_time(TimeTasks::FIELDS),
      get_time(TimeTasks::PARTICLES),
      get_time(TimeTasks::BFIELD),
      get_time()
      );
    fprintf(stdout, TIMING_PREFIX
      "%4.2f "
      "%4.2f "
      "%4.2f "
      "%4.2f "
      "%4.2f (communication)\n",
      get_communicate(TimeTasks::MOMENTS),
      get_communicate(TimeTasks::FIELDS),
      get_communicate(TimeTasks::PARTICLES),
      get_communicate(TimeTasks::BFIELD),
      get_communicate()
      );
    fprintf(stdout, TIMING_PREFIX
      "%4.2f "
      "%4.2f "
      "%4.2f "
      "%4.2f "
      "%4.2f (computation)\n",
      get_compute(TimeTasks::MOMENTS),
      get_compute(TimeTasks::FIELDS),
      get_compute(TimeTasks::PARTICLES),
      get_compute(TimeTasks::BFIELD),
      get_compute()
      );
    //fprintf(stdout, TIMING_PREFIX
    //  "MOMS comm  FLDS comm  PCLS comm  CYCL comm\n");
    //fprintf(stdout, TIMING_PREFIX
    //  "%4.2f "
    //  "%4.2f  "
    //  "%4.2f "
    //  "%4.2f  "
    //  "%4.2f "
    //  "%4.2f  "
    //  "%4.2f "
    //  "%4.2f\n",
    //  get_time(TimeTasks::MOMENTS),
    //  get_communicate(TimeTasks::MOMENTS),
    //  get_time(TimeTasks::FIELDS),
    //  get_communicate(TimeTasks::FIELDS),
    //  get_time(TimeTasks::PARTICLES),
    //  get_communicate(TimeTasks::PARTICLES),
    //  get_time(),
    //  get_communicate()
    //  );
    fflush(stdout);
  }
}
