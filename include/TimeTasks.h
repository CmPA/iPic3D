#ifndef __TimeTasks_H__
#define __TimeTasks_H__

class TimeTasks {

public:
  // legitimate active subcycle values
  enum Tasks {
    NONE = 0,
    MOMENTS,
    FIELDS,
    PARTICLES,
    BFIELD,
    LAST,
  };
  enum Modes {
    COMPUTATION = 0,
    COMMUNICATION,
  };

public:
  void setActiveTask(int arg) {
    active_task = arg;
  } void setActiveMode(int in) {
    t_start_communicate = in;
  }
  void resetCycle();
  void start(int taskid);
  void end(int taskid);
  void start_communicate();
  void addto_communicate();
  void print_cycle_times();
  TimeTasks() {
    resetCycle();
  }
  double get_time(int arg) {
    return task_duration[arg];
  }
  double get_communicate(int arg) {
    return communicate[arg];
  }

  double get_communicate() {
    double total = 0.;
    for (int i = NONE + 1; i < LAST; i++) {
      total += communicate[i];
    }
    return total;
  }

  double get_time() {
    double total = 0.;
    for (int i = NONE + 1; i < LAST; i++) {
      total += task_duration[i];
    }
    return total;
  }

  double get_compute(int arg) {
    return get_time(arg) - get_communicate(arg);
  }
  double get_compute() {
    return get_time() - get_communicate();
  }

private:
  int active_task;
  int active_mode;
  double t_start_communicate;
  double start_times[LAST];
  double task_duration[LAST];
  double communicate[LAST];
  double compute[LAST];

};

extern TimeTasks timeTasks;

#endif
