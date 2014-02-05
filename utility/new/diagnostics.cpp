
/** implementation of declarations in utility/debug.h **/

#include <stdarg.h>
#include "TimeTasks.h"
#include "debug.h"
#include "asserts.h"
#include "../mpidata/MPIdata.h" // for rank

#define implement_dprintvar_fileLine(code,type) \
  void dprintvar_fileLine(const char* func, const char* file, int line, \
    const char* name, type val) \
  { \
    dfprintf_fileLine(stderr,func,file,line, \
      code " == %s",val,name); \
  }
implement_dprintvar_fileLine("%s",const char*);
implement_dprintvar_fileLine("%d",int);
//implement_dprintvar_fileLine("%24.16e",double);
implement_dprintvar_fileLine("%f",double);

void dfprintf_fileLine(FILE* fptr, const char *func, const char *file, int line_number,
  const char *format, ...)
{
  fflush(fptr);
  va_list args;
  va_start(args, format);
  fprintf(fptr, "(%d) DEBUG %s(), %s:%d: ",
    get_rank(), func,
    file, // my_basename(file),
    line_number);
  /* print out remainder of message */
  vfprintf(fptr, format, args);
  va_end(args);
  fprintf(fptr,"\n");
  fflush(fptr);
}

int get_rank() { return mpi->rank; }

/** implementation of declarations in utility/assert.h **/

// so that we can print doubles to desired precision
//
void assert_error(const char* file, int line, const char* func,
  const char* op, const char* lhs_str, const char* rhs_str,
  double lhs, double rhs)
{
  fprintf(stderr,"ERROR in file %s, line %d, function %s"
      "\n\tassertion failed: %s %s %s, i.e., %24.16e %s %24.16e\n",
    file, line, func, lhs_str, op, rhs_str, lhs, op, rhs);
  abort();
}

#define implement_assert_errmsg(t1,t2) \
  void assert_error(const char* file, int line, const char* func, \
    const char* op, const char* lhs_str, const char* rhs_str, \
    t1 lhs, t2 rhs) \
  { \
    std::cerr<< "ERROR in file " << file << ", line " << line  \
      << ", function " << func  \
      <<"\n\tassertion failed: " << lhs_str << op << rhs_str \
      << ", i.e., " << lhs << op << rhs << endl; \
      abort(); \
  }

implement_assert_errmsg(int,int);
implement_assert_errmsg(const char*,const char*);
implement_assert_errmsg(const string&,const string&);

/** implementation of declarations in utility/TimeTasks.h **/

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
  if(!get_rank())
  {
    fflush(stdout);
    fprintf(stdout,"=== timing information for cycle=== \n");
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
    fprintf(stdout, TIMING_PREFIX
      "MOMS comm  FLDS comm  PCLS comm  CYCL comm\n");
    fprintf(stdout, TIMING_PREFIX
      "%4.2f "
      "%4.2f  "
      "%4.2f "
      "%4.2f  "
      "%4.2f "
      "%4.2f  "
      "%4.2f "
      "%4.2f\n",
      get_time(TimeTasks::MOMENTS),
      get_communicate(TimeTasks::MOMENTS),
      get_time(TimeTasks::FIELDS),
      get_communicate(TimeTasks::FIELDS),
      get_time(TimeTasks::PARTICLES),
      get_communicate(TimeTasks::PARTICLES),
      get_time(),
      get_communicate()
      );
    fflush(stdout);
  }
}
