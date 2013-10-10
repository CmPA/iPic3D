
#include <mpi.h>
#include <stdarg.h>
#include "TimeTasks.h"
#include "asserts.h"
#include "MPIdata.h" // for get_rank
#include "debug.h"

/** implementation of declarations in utility/TimeTasks.h **/

TimeTasks timeTasks;

static const char *taskNames[] = // order must agree with Tasks in TimeTasks.h
{
  "none",
  "moments",
  "fields",
  "particles",
  "last",
  "bfield",
  "moment_accumulation",
  "moment_reduction",
  "number_of_tasks"
};

const char* TimeTasks::get_taskname(int arg)
{
  assert_le(arg,NUMBER_OF_TASKS);
  return taskNames[arg];
}

void TimeTasks::resetCycle()
{
  for(int e=0;e<NUMBER_OF_TASKS;e++)
  {
    task_duration[e]=0.;
    compute[e]=0.;
    communicate[e]=0.;
    active[e]=false;
    stack_depth[e]=0;
    start_times[e]=0.;
  }
  active_task=NONE;
  communicating=false;
}
void TimeTasks::start_main_task(TimeTasks::Tasks taskid)
{
  assert(is_exclusive(taskid));
  assert_ne(active_task, taskid);
  active_task = taskid;
  assert(!active[taskid]);
  active[taskid]=true;
  //if(!MPIdata::get_rank())
  //dprintf("starting task %s at time %24.16e\n", get_taskname(taskid), MPI_Wtime());
}
void TimeTasks::start_task(TimeTasks::Tasks taskid)
{
  assert(!is_exclusive(taskid));
  assert(!active[taskid]);
  active[taskid]=true;
  //dprintf("starting task %s at time %24.16e\n", get_taskname(taskid), MPI_Wtime());
}
// have to manage the task stack explicitly
void TimeTasks::start_task(TimeTasks::Tasks taskid, double start_time)
{
  if(stack_depth[taskid]==0)
  {
    start_times[taskid]=start_time;
    start_task(taskid);
  }
  stack_depth[taskid]++;
  //dprintf("starting task %s at time %24.16e\n", get_taskname(taskid), start_time);
}
void TimeTasks::end_main_task(TimeTasks::Tasks taskid, double start_time)
{
  end_task(taskid, start_time);
  active_task = NONE;
}
void TimeTasks::end_task(TimeTasks::Tasks taskid, double start_time)
{
  assert(active[taskid]);
  double now = MPI_Wtime();
  // compute time spent on task
  task_duration[taskid] += now - start_time;
  active[taskid] = false;
}
// have to manage the task stack explicitly
void TimeTasks::end_task(TimeTasks::Tasks taskid)
{
  stack_depth[taskid]--;
  assert_ge(stack_depth[taskid],0);
  if(stack_depth[taskid]==0)
  {
    end_task(taskid, start_times[taskid]);
  }
}
void TimeTasks::end_communicating(double start_time)
{
  //if(!active_task) return;
  assert(active_task);
  assert(communicating);
  double additional_communication_time = MPI_Wtime()-start_time;
  //dprint(additional_communication_time);
  communicate[active_task] += additional_communication_time;
  communicating=false;
}
#define TIMING_PREFIX "| "
void TimeTasks::print_cycle_times(int cycle)
{
  // calculate portion of time spent computing
  //
  for(int e=NONE+1; e<NUMBER_OF_TASKS; e++)
  {
    compute[e] = task_duration[e]-communicate[e];
  }

  FILE* file = stdout;
  // we could report average for all processes
  if(!MPIdata::get_rank())
  {
    fflush(file);
    fprintf(file,"=== times for cycle %d for rank %d === \n",
      cycle,
      MPIdata::get_rank());
    fprintf(file, TIMING_PREFIX "total  comput commun task\n");
    for(int e=NONE+1; e<LAST; e++)
    {
      fprintf(file, TIMING_PREFIX "%6.3f %6.3f %6.3f %s\n",
      get_time(e),
      get_compute(e),
      get_communicate(e),
      get_taskname(e));
    }
    fprintf(file, TIMING_PREFIX "%6.3f %6.3f %6.3f %s\n",
      get_time(),
      get_compute(),
      get_communicate(),
      "[total times]");

    fprintf(file, TIMING_PREFIX "time  subtask\n");
    for(int e=LAST+1; e<NUMBER_OF_TASKS; e++)
    {
      assert_eq(stack_depth[e],0);
      fprintf(file, TIMING_PREFIX "%5.3f %s\n",
      get_time(e),
      get_taskname(e));
    }
    
    fflush(file);
  }
}

TimeTasks_caller_to_set_communication_mode_for_scope::
TimeTasks_caller_to_set_communication_mode_for_scope()
{
  already_communicating = timeTasks.get_communicating();
  if(!already_communicating)
  {
    start_time = MPI_Wtime();
    timeTasks.set_communicating(true);
  }
}
TimeTasks_caller_to_set_communication_mode_for_scope::
~TimeTasks_caller_to_set_communication_mode_for_scope()
{
  if(!already_communicating)
  {
    timeTasks.end_communicating(start_time);
  }
}
