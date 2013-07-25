
#include <iostream>
#include "asserts.h"

void assert_error(const char *file, int line, const char *func, const char *op, const char *lhs_str, const char *rhs_str, double lhs, double rhs) {
  fprintf(stdout, "ERROR in file %s, line %d, function %s" "\n\tassertion failed: %s %s %s, i.e., %24.16e %s %24.16e\n", file, line, func, lhs_str, op, rhs_str, lhs, op, rhs);
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
      << ", i.e., " << lhs << op << rhs << std::endl; \
      abort(); \
  }

implement_assert_errmsg(size_t, size_t);
implement_assert_errmsg(int, size_t);
implement_assert_errmsg(size_t, int);
implement_assert_errmsg(int, int);
implement_assert_errmsg(const char *, const char *);
