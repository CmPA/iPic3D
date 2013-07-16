
#include "MPIdata.h" // for get_rank
#include "debug.h"

#define implement_dprintvar_fileLine(code,type) \
  void dprintvar_fileLine(const char* func, const char* file, int line, \
    const char* name, type val) \
  { \
    dfprintf_fileLine(stdout,func,file,line, code " == %s",val,name); \
  }

implement_dprintvar_fileLine("%s", const char *);
implement_dprintvar_fileLine("%d", int);
implement_dprintvar_fileLine("%f", double);

void dfprintf_fileLine(FILE * fptr, const char *func, const char *file, int line_number, const char *format, ...) {
  fflush(fptr);
  va_list args;
  va_start(args, format);
  fprintf(fptr, "(%d) DEBUG %s(), %s:%d: ",
    MPIdata::get_rank(),
    func, file, // my_basename(file),
    line_number);
  /* print out remainder of message */
  vfprintf(fptr, format, args);
  va_end(args);
  fprintf(fptr, "\n");
  fflush(fptr);
}
