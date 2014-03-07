 
#ifndef NO_MPI
  #include "MPIdata.h" // for get_rank
#endif
#include "ompdefs.h" // for omp_get_thread_num
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include "errors.h"

/** implementation of declarations in errors.h **/

void fprintf_fileLine(FILE * fptr, const char *type, const char *func, const char *file, int line_number, const char *format, ...);

// This is not thread-safe.
//void errmsg_printf_fileLine(const char *func, const char *file, int line_number,
//  const char *format, ...)
//{
//  FILE* fptr = stdout;
//  fflush(fptr);
//  va_list args;
//  va_start(args, format);
//  fprintf(fptr, "ERROR in function %s, file %s, line %d: \n\t",
//    func, file, line_number);
//  /* print out remainder of message */
//  vfprintf(fptr, format, args);
//  va_end(args);
//  // append terminating newline so user does not have to do it
//  fprintf(fptr, "\n");
//  fflush(fptr);
//
//  abort();
//}

// lazy implementation using streams class
//
using namespace std;
#define implement_invalid_value_error(t1) \
  void invalid_value_error_fileLine(const char* file, int line, const char* func, \
    const char* type, const char* expr, t1 val) \
  { \
    /* To be thread-safe, write the message to a string and \
     * then print it out as an atomic operation. */ \
    std::stringstream ss; \
    ss << "(" << MPIdata::get_rank() << "." << omp_get_thread_num() << ") " \
      << "ERROR in file " << file << ", line " << line  \
      << ", function " << func  \
      <<"\n\t" << type << " value: " << expr << " = " << val << endl; \
    fflush(stdout); \
      { fprintf(stdout,ss.str().c_str()); } \
    fflush(stdout); \
    abort(); \
  }

implement_invalid_value_error(double);
implement_invalid_value_error(int);
implement_invalid_value_error(const char*);

/*! a more verbose version of fprintf_fileLine for use in 
 * warnings and error messages */
void eprintf_fileLine(FILE * fptr, const char *type,
  const char *func, const char *file, int line_number,
  const char *format, ...) 
{
  // writing directly to fptr would avoid limiting the length
  // of the output string, but by first writing to a string
  // we achieve thread safety.
  //
  // write the message to a string.
  //
  const int maxchars = 1024;
  char error_msg[maxchars+2];
  // identify the process and thread
  char process_thread_str[50];
  #ifndef NO_MPI
    #ifdef _OPENMP
      snprintf(process_thread_str, 50, ", process %d, thread %d",
        MPIdata::get_rank(), omp_get_thread_num());
    #else
      snprintf(process_thread_str, 50, ", process %d",
        MPIdata::get_rank());
    #endif
  #else
    #ifdef _OPENMP
      snprintf(process_thread_str, 50, ", thread %d",
        omp_get_thread_num());
    #else
      sprintf(process_thread_str, "");
    #endif
  #endif
  char *sptr = error_msg;
  int chars_so_far=0;
  va_list args;
  va_start(args, format);
  chars_so_far = snprintf(sptr, maxchars,
    "%s in method %s(), file %s, line %d%s:\n\t",
    type,
    func, file, // my_basename(file),
    line_number, process_thread_str);
  /* print out remainder of message */
  chars_so_far += vsnprintf(sptr+chars_so_far, maxchars-chars_so_far, format, args);
  va_end(args);
  sprintf(sptr+chars_so_far, "\n");

  // print the message
  fflush(fptr);
    // #pragma omp critical // need this?
    { fprintf(fptr,error_msg); }
  fflush(fptr);
  abort();
}
