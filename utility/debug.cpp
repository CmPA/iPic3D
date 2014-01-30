
#ifndef NO_MPI
  #include "MPIdata.h" // for get_rank
#endif
#include "ompdefs.h" // for omp_get_thread_num
#include "debug.h"
#include "parallel.h" // temporary

#define implement_dprintvar_fileLine(code,type) \
  void printvar_fileLine(const char* func, const char* file, int line, \
    const char* name, type val) \
  { \
    fprintf_fileLine(stdout,"DEBUG", func,file,line, code " == %s",val,name); \
  }

implement_dprintvar_fileLine("%s", const char *);
implement_dprintvar_fileLine("%d", int);
implement_dprintvar_fileLine("%e", double);
implement_dprintvar_fileLine("%p", const void *);

// void dfprintf_fileLine(FILE * fptr, const char *func, const char *file, int line_number, const char *format, ...)
// {
//   // writing directly to fptr would avoid limiting the length
//   // of the output string, but by first writing to a string
//   // we achieve thread safety.
//   //
//   // write the message to a string.
//   //
//   const int maxchars = 1024;
//   char error_msg[maxchars+2];
//   // identify the process and thread
//   char process_thread_str[20];
//   #ifndef NO_MPI
//     #ifdef _OPENMP
//       snprintf(process_thread_str, 20, "(%d.%d) ",
//         MPIdata::get_rank(), omp_get_thread_num());
//     #else
//       snprintf(process_thread_str, 20, "(%d)",
//         MPIdata::get_rank());
//     #endif
//   #else
//     #ifdef _OPENMP
//       snprintf(process_thread_str, 20, "(.%d) ",
//         omp_get_thread_num());
//     #else
//       snprintf(process_thread_str, 20, "");
//     #endif
//   #endif
//   char *sptr = error_msg;
//   int chars_so_far=0;
//   va_list args;
//   va_start(args, format);
//   chars_so_far = snprintf(sptr, maxchars,
//     "%sDEBUG %s(), %s:%d: ",
//     process_thread_str,
//     func, file, // my_basename(file),
//     line_number);
//   /* print out remainder of message */
//   chars_so_far += vsnprintf(sptr+chars_so_far, maxchars-chars_so_far, format, args);
//   va_end(args);
//   sprintf(sptr+chars_so_far, "\n");
// 
//   // print the message
//   fflush(fptr);
//     fprintf(fptr,error_msg);
//   fflush(fptr);
// }

void fprintf_fileLine(FILE * fptr,
  const char *type, const char *func, const char *file, int line_number,
  const char *format, ...)
{
  if(!is_output_thread()) return; // temporary

  // writing directly to fptr would avoid limiting the length
  // of the output string, but by first writing to a string
  // we achieve thread safety.
  //
  // write the message to a string.
  //
  const int maxchars = 1024;
  char error_msg[maxchars+2];
  // identify the process and thread
  char process_thread_str[20];
  #ifndef NO_MPI
    #ifdef _OPENMP
      snprintf(process_thread_str, 20, "(%d.%d) ",
        MPIdata::get_rank(), omp_get_thread_num());
    #else
      snprintf(process_thread_str, 20, "(%d)",
        MPIdata::get_rank());
    #endif
  #else
    #ifdef _OPENMP
      snprintf(process_thread_str, 20, "(.%d) ",
        omp_get_thread_num());
    #else
      snprintf(process_thread_str, 20, "");
    #endif
  #endif
  char *sptr = error_msg;
  int chars_so_far=0;
  va_list args;
  va_start(args, format);
  chars_so_far = snprintf(sptr, maxchars,
    "%s%s %s(), %s:%d: ",
    process_thread_str,
    type,
    func, file, // my_basename(file),
    line_number);
  /* print out remainder of message */
  chars_so_far += vsnprintf(sptr+chars_so_far, maxchars-chars_so_far, format, args);
  va_end(args);
  sprintf(sptr+chars_so_far, "\n");

  // print the message
  fflush(fptr);
    fprintf(fptr,error_msg);
  fflush(fptr);
}

