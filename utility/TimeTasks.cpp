
#include <mpi.h>
#include <stdarg.h>
#include "TimeTasks.h"
#include "asserts.h"
#include "MPIdata.h" // for get_rank
#include "parallel.h"
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
  "moment_pcl_sorting",
  "moment_accumulation",
  "moment_reduction",
  "mover_pcl_sorting",
  "mover_pcl_moving",
  "transpose_pcls_to_AoS",
  "transpose_pcls_to_SoA",
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
  if(!is_output_thread()) return;
  assert(is_exclusive(taskid));
  assert_ne(active_task, taskid);
  active_task = taskid;
  assert(!active[taskid]);
  active[taskid]=true;
}
void TimeTasks::start_task(TimeTasks::Tasks taskid)
{
  if(!is_output_thread()) return;
  assert(!is_exclusive(taskid));
  assert(!active[taskid]);
  active[taskid]=true;
}
// have to manage the task stack explicitly
void TimeTasks::start_task(TimeTasks::Tasks taskid, double start_time)
{
  if(!is_output_thread()) return;
  if(stack_depth[taskid]==0)
  {
    start_times[taskid]=start_time;
    start_task(taskid);
  }
  stack_depth[taskid]++;
}
void TimeTasks::end_main_task(TimeTasks::Tasks taskid, double start_time)
{
  if(!is_output_thread()) return;
  end_task(taskid, start_time);
  active_task = NONE;
}
void TimeTasks::end_task(TimeTasks::Tasks taskid, double start_time)
{
  if(!is_output_thread()) return;
  assert(active[taskid]);
  double now = MPI_Wtime();
  // compute time spent on task
  task_duration[taskid] += now - start_time;
  active[taskid] = false;
}
// have to manage the task stack explicitly
void TimeTasks::end_task(TimeTasks::Tasks taskid)
{
  if(!is_output_thread()) return;
  stack_depth[taskid]--;
  assert_ge(stack_depth[taskid],0);
  if(stack_depth[taskid]==0)
  {
    end_task(taskid, start_times[taskid]);
  }
}
void TimeTasks::end_communicating(double start_time)
{
  if(!is_output_thread()) return;
  assert(active_task);
  assert(communicating);
  double additional_communication_time = MPI_Wtime()-start_time;
  communicate[active_task] += additional_communication_time;
  communicating=false;
}
#define TIMING_PREFIX "| "
void TimeTasks::print_cycle_times(int cycle)
{
  if(!is_output_thread()) return;
  FILE* file = stdout;
  // we could report average for all processes
  //if(!MPIdata::get_rank())
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
      communicate[e],
      get_taskname(e));
    }

    // report total times
    //
    // get total time spent on exclusive tasks
    //
    double total_task_duration = 0.;
    for (int i = NONE + 1; i < LAST; i++) {
      total_task_duration += task_duration[i];
    }
    // get total time spent in exclusive tasks spent communicating
    //
    double total_communicate = 0.;
    for (int i = NONE + 1; i < LAST; i++) {
      total_communicate += communicate[i];
    }
    const double total_computing_time = total_task_duration - total_communicate;
    fprintf(file, TIMING_PREFIX "%6.3f %6.3f %6.3f %s\n",
      total_task_duration,
      total_computing_time,
      total_communicate,
      "[total times]");

    fprintf(file, TIMING_PREFIX "time   subtask\n");
    for(int e=LAST+1; e<NUMBER_OF_TASKS; e++)
    {
      // do not show tasks that are not executed
      double elapsed_time = get_time(e);
      if(!elapsed_time)
        continue;

      assert_eq(stack_depth[e],0);
      fprintf(file, TIMING_PREFIX "%6.3f %s\n",
      elapsed_time,
      get_taskname(e));
    }
    
    fflush(file);
  }
}

// The following three methods provide for a hack by which
// the timeTasks copies of all threads are averaged.
// 
void TimeTasks::operator/=(int num)
{
  assert(false); // this method is not in use.
  for(int e=NONE+1;e<NUMBER_OF_TASKS;e++)
  {
    task_duration[e]/=num;
    start_times[e]/=num;
    communicate[e]/=num;
  }
}
void TimeTasks::operator+=(const TimeTasks& arg)
{
  assert(false); // this method is not in use.
  active_task = arg.active_task;
  communicating = arg.communicating;
  for(int e=NONE+1;e<NUMBER_OF_TASKS;e++)
  {
    active[e] = arg.active[e];
    task_duration[e]+=arg.task_duration[e];
    stack_depth[e] = arg.stack_depth[e];
    start_times[e]+=arg.start_times[e];
    communicate[e]+=arg.communicate[e];
  }
}
void TimeTasks::operator=(const TimeTasks& arg)
{
  assert(false); // this method is not in use.
  active_task = arg.active_task;
  communicating = arg.communicating;
  for(int e=NONE+1;e<NUMBER_OF_TASKS;e++)
  {
    active[e] = arg.active[e];
    task_duration[e]=arg.task_duration[e];
    stack_depth[e] = arg.stack_depth[e];
    start_times[e]=arg.start_times[e];
    communicate[e]=arg.communicate[e];
  }
}

TimeTasks_caller_to_set_main_task_for_scope::
TimeTasks_caller_to_set_main_task_for_scope(TimeTasks::Tasks _task) :
  task(_task)
{
  if(!is_output_thread()) return;
  start_time = MPI_Wtime();
  timeTasks.start_main_task(task);
}
TimeTasks_caller_to_set_main_task_for_scope::
~TimeTasks_caller_to_set_main_task_for_scope()
{
  if(!is_output_thread()) return;
  timeTasks.end_main_task(task, start_time);
}

TimeTasks_caller_to_set_task_for_scope::
TimeTasks_caller_to_set_task_for_scope(TimeTasks::Tasks _task)
{
  if(!is_output_thread()) return;
  task = _task;
  already_active = timeTasks.is_active(task);
  if(!already_active)
  {
    start_time = MPI_Wtime();
    timeTasks.start_task(task);
  }
}
TimeTasks_caller_to_set_task_for_scope::
~TimeTasks_caller_to_set_task_for_scope()
{
  if(!is_output_thread()) return;
  if(already_active)
  {
    assert(timeTasks.is_active(task));
  }
  else
  {
    timeTasks.end_task(task, start_time);
  }
}

TimeTasks_caller_to_set_communication_mode_for_scope::
TimeTasks_caller_to_set_communication_mode_for_scope()
{
  if(!is_output_thread()) return;
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
  if(!is_output_thread()) return;
  if(!already_communicating)
  {
    timeTasks.end_communicating(start_time);
  }
}
