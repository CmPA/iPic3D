 
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include "errors.h"
//#include "MPIdata.h" // for rank

/** implementation of declarations in errors.h **/

void errmsg_printf_fileLine(const char *func, const char *file, int line_number,
  const char *format, ...)
{
  FILE* fptr = stdout;
  fflush(fptr);
  va_list args;
  va_start(args, format);
  fprintf(fptr, "ERROR in function %s, file %s, line %d: \n\t",
    func, file, line_number);
  /* print out remainder of message */
  vfprintf(fptr, format, args);
  va_end(args);
  // append terminating newline so user does not have to do it
  fprintf(fptr, "\n");
  fflush(fptr);

  abort();
}

#include <iostream>
using namespace std;
#define implement_invalid_value_error(t1) \
  void invalid_value_error_fileLine(const char* file, int line, const char* func, \
    const char* type, const char* expr, t1 val) \
  { \
    std::cerr<< "ERROR in file " << file << ", line " << line  \
      << ", function " << func  \
      <<"\n\t" << type << " value: " << expr << " = " << val << endl; \
      abort(); \
  }

implement_invalid_value_error(double);
implement_invalid_value_error(int);
implement_invalid_value_error(const char*);

