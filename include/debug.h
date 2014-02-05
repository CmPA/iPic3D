/*******************************************************************************************
  debug.h  -  definitions for debug and diagnostics (written by Alec Johnson)
 ********************************************************************************************/
#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <cstdarg>
#include <cstdio>

#include "errors.h"

void fprintf_fileLine(FILE * fptr, const char *type, const char *func,
  const char *file, int line_number, const char *format, ...);

#define dprintf(args...) fprintf_fileLine(stdout, "DEBUG", __func__, __FILE__, __LINE__,## args)
#define dprint(var) printvar_fileLine(__func__, __FILE__,__LINE__,#var,var);
#define dprint0(var) dprint(var)
#define declare_dprintvar_fileLine(type) \
void printvar_fileLine(const char*,const char*,int,const char*,type);

declare_dprintvar_fileLine(int);
declare_dprintvar_fileLine(double);
declare_dprintvar_fileLine(const char *);
declare_dprintvar_fileLine(const void *);

#endif
