#ifndef __TimeTasks_H__
#define __TimeTasks_H__
#include "assert.h"

/* Avoid direct use of this class.
   Instead, use and add to the macros at the bottom
   so that we can redefine the macros when desired
   (e.g. defining them to the empty string to
   remove performance penalty).
 */

class TimeTasks
{
 public:

  // legitimate active subcycle values
  //
  // timeTasks_set_task(0) is a no-op, so
  // MOMENT_REDUCTION=0
  // would prevent monitoring of this task.
  //
  enum Tasks // order must agree with taskNames in TimeTasks.cpp
  {
    NONE = 0,
    MOMENTS,
    FIELDS,
    PARTICLES,
    LAST, // no more exclusive tasks
    BFIELD,
    MOMENT_PCL_SORTING,
    MOMENT_ACCUMULATION,
    MOMENT_REDUCTION,
    MOVER_PCL_SORTING,
    MOVER_PCL_MOVING,
    TRANSPOSE_PCLS_TO_AOS,
    TRANSPOSE_PCLS_TO_SOA,
    NUMBER_OF_TASKS // this line should be last
  };

 private:
  //enum Modes // for exclusive tasks
  //{
  //  COMPUTATION = 0,
  //  COMMUNICATION,
  //};

 public: // methods

  TimeTasks() {
    resetCycle();
  }

  // monitoring
  //
  void resetCycle();
  //
  // hack to support averaging timeTasks copies of all threads.
  //
  void operator+=(const TimeTasks& arg);
  void operator/=(int num);
  void operator=(const TimeTasks& arg);
  //
  // provide start_time on ending call
  //
  void end_communicating(double start_time);
  void start_main_task(TimeTasks::Tasks taskid);
  void end_main_task(TimeTasks::Tasks taskid, double start_time);
  void start_task(TimeTasks::Tasks taskid);
  void end_task(TimeTasks::Tasks taskid, double start_time);
  //
  // provide start_time at starting call
  //
  void start_task(TimeTasks::Tasks taskid, double start_time);
  void end_task(TimeTasks::Tasks taskid);

  // accessors
  //
  bool is_active(Tasks taskid){ return active[taskid]; }
  bool get_communicating() { return communicating; }
  void set_communicating(bool val) { communicating = val; }
  int get_stack_depth(TimeTasks::Tasks taskid) { return stack_depth[taskid]; }

  // reporting
  //
  void print_cycle_times(int cycle);

 private:

  // is task exclusive?
  bool is_exclusive(Tasks taskid) { return (taskid < LAST); }

  // reporting
  //
  double get_time(int arg) {
    return task_duration[arg];
  }
  double get_communicate(int arg) {
    return communicate[arg];
  }
  double get_compute(int arg) {
    return get_time(arg) - get_communicate(arg);
  }
  const char* get_taskname(int arg);

 private:
  int active_task;
  bool active[NUMBER_OF_TASKS];
  bool communicating;
  double task_duration[NUMBER_OF_TASKS];
  double communicate[NUMBER_OF_TASKS];
  int stack_depth[NUMBER_OF_TASKS];
  double start_times[NUMBER_OF_TASKS];
};

extern TimeTasks timeTasks;

// construct an anonymous instance of TimeTasksCaller
class TimeTasks_caller_to_set_main_task_for_scope
{
  double start_time;
  TimeTasks::Tasks task;
 public:
  TimeTasks_caller_to_set_main_task_for_scope(TimeTasks::Tasks _task);
  ~TimeTasks_caller_to_set_main_task_for_scope();
};

class TimeTasks_caller_to_set_task_for_scope
{
  bool already_active;
  double start_time;
  TimeTasks::Tasks task;
 public:
  TimeTasks_caller_to_set_task_for_scope(TimeTasks::Tasks _task);
  ~TimeTasks_caller_to_set_task_for_scope();
};

class TimeTasks_caller_to_set_communication_mode_for_scope
{
 private:
  bool already_communicating;
  double start_time;
 public:
  TimeTasks_caller_to_set_communication_mode_for_scope();
  ~TimeTasks_caller_to_set_communication_mode_for_scope();
};

// These macros could be changed to provide file and line number
//
// We need to create nonanonymous instances so that the destructor
// will not be called until the end of the scope, so we use the preprocessor
// to generate unique names of nonanonymous instances.
//
#define timeTasks_set_main_task(task) \
  TimeTasks_caller_to_set_main_task_for_scope myFunnyInstance(task);
#define timeTasks_set_task(task) \
  TimeTasks_caller_to_set_task_for_scope myFunnyName##__func__##__LINE__(task);
#define timeTasks_set_communicating() \
  TimeTasks_caller_to_set_communication_mode_for_scope myFunnyCommunicationInstance;
//
// The scoping trick does not work if the timeTasks call needs to be conditional,
// so we also provide the ability to explicitly begin and end.
#define timeTasks_begin_task(task) if(task) timeTasks.start_task(task, MPI_Wtime());
#define timeTasks_end_task(task) if(task) timeTasks.end_task(task);
//

#endif
