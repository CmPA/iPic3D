
#include <stdarg.h>
#include "asserts.h"
#include "debug.h"
#include "../mpidata/MPIdata.h" // for rank

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

/** implementation of declarations in utility/debug.h **/

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

