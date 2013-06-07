/*******************************************************************************************
  debug.h  -  definitions for debug and diagnostics (written by Alec Johnson)
 ********************************************************************************************/
#ifndef utility_debug_H
#define utility_debug_H

void dfprintf_fileLine(FILE* fptr, const char *func, const char *file, int line_number,
  const char *format, ...);
#define dprintf(args...) dfprintf_fileLine(stderr, __func__, __FILE__, __LINE__,## args)
#define dprint(var) dprintvar_fileLine(__func__,__FILE__,__LINE__,#var,var);
#define dprint0(var) if(!get_rank()) dprint(var)
#define declare_dprintvar_fileLine(type) \
  void dprintvar_fileLine(const char*,const char*,int,const char*,type);
declare_dprintvar_fileLine(int);
declare_dprintvar_fileLine(double);
declare_dprintvar_fileLine(const char*);

int get_rank();
#endif
