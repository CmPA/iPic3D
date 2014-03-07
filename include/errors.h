#ifndef ipic_errors_H
#define ipic_errors_H

//void errmsg_printf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
//void eprintf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
//void Wprintf_fileLine(const char *func, const char *file, int line_number, const char *format, ...);
void eprintf_fileLine(FILE * fptr, const char *type,
  const char *func, const char *file, int line_number,
  const char *format, ...);

#define eprintf(args...) \
  eprintf_fileLine(stdout,"ERROR",__func__, __FILE__, __LINE__, ## args);
#define error_printf(args...) \
  eprintf_fileLine(stdout,"ERROR",__func__, __FILE__, __LINE__, ## args);
//#define eprintf(args...) \
//  eprintf_fileLine("ERROR",__func__, __FILE__, __LINE__, ## args);
#define warning_printf(args...) \
  fprintf_fileLine(stdout,"WARNING",__func__, __FILE__, __LINE__, ## args);
#define declare_invalid_value_error(t1) \
  void invalid_value_error_fileLine(const char* file, int line, const char* func, \
    const char* type, const char* expr, t1 val);
declare_invalid_value_error(double);
declare_invalid_value_error(int);
declare_invalid_value_error(const char*);
#define unsupported_value_error(val) invalid_value_error_fileLine( \
  __FILE__, __LINE__, __func__, "unsupported", #val, val);
#define invalid_value_error(val) invalid_value_error_fileLine( \
  __FILE__, __LINE__, __func__, "invalid", #val, val);

#endif
