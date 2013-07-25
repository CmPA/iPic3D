#ifndef __ASSERTS_H__
#define __ASSERTS_H__

#include <cstdlib>
#include <cstdio>

// asserts higher than MAX_ASSERT_LEVEL are compiled out
// (to avoid slowing execution)

// uncomment next line to make this variable global
// #undef MAX_ASSERT_LEVEL
#ifndef MAX_ASSERT_LEVEL
#define MAX_ASSERT_LEVEL 2
#endif

#ifdef NDEBUG                   // completely turns off assert statements

// if debug is turned off create vacuous versions
// of everything that user might use
// (a list of what we currently actually support;
// other stuff below is experimental)

#define assert1(e) ((void)0)
#define assert2(e) ((void)0)
#define assert3(e) ((void)0)
#define assert_printf(args...) ((void)0)
#define assert_printf1(args...) ((void)0)
#define assert_printf2(args...) ((void)0)
#define assert_printf3(args...) ((void)0)
#define assert_almost_eq(a) ((void)0)
#define assert_eq(a,b) ((void)0)
#define assert_ne(a,b) ((void)0)
#define assert_gt(a,b) ((void)0)
#define assert_lt(a,b) ((void)0)
#define assert_ge(a,b) ((void)0)
#define assert_le(a,b) ((void)0)
#define assert_isnum(a) ((void)0)

#else // ifndef NDEBUG

// override system assert.h
// #define assert_fileLine(e, file, line) \
// ((void)printf ("%s:%u: failed assertion `%s'\n", file, line, e), abort())
// void eprintf_fileLine(const char *func, const char *file, int line_number,
// const char *format, ...);

#define dassert_fileLine(e, file, line, func) \
  (void)(printf("ERROR: file %s, line %d, function %s:\n\tfailed assertion: (%s)\n", file, line, func,e),abort())
#define dassert_printf_fileLine(e, file, line, func, args...) \
  (void)(printf("ERROR: file %s, line %d, function %s:\n\tfailed assertion: (%s)\n\t", file, line, func,e), printf(args), printf("\n"), abort())

// comment out the next line if __builtin_expect causes problems
#define USE_GCC_OPTIMIZATION
#ifndef USE_GCC_OPTIMIZATION

// override system assert.h
// #define assert(e) \
// ((void) ((e) ? (void)0 : assert_fileLine (#e, __FILE__, __LINE__)))

#define dassert_(e)  \
  ((void) ((e) ? (void)0 : dassert_fileLine(#e, __FILE__, __LINE__, __func__)))
#define dassert_printf_(e, args...)  \
  ((void) ((e) ? (void)0 : dassert_printf_fileLine(#e, __FILE__, __LINE__, __func__,##args)))
#else // ifdef USE_GCC_OPTIMIZATION
// optimized version of preceding
// #define assert(e) \
// (__builtin_expect(!(e), 0) ? assert_fileLine (#e, __FILE__, __LINE__) : (void)0)
#define dassert_(e)  \
  (__builtin_expect(!(e), 0) ? dassert_fileLine (#e, __FILE__, __LINE__, __func__) : (void)0)
#define dassert_printf(e, args...)  \
  (__builtin_expect(!(e), 0) ? dassert_printf_fileLine (#e, __FILE__, __LINE__, __func__,##args) : (void)0)
#endif // USE_GCC_OPTIMIZATION

#if(MAX_ASSERT_LEVEL>=1)
#define assert1 dassert_
#define assert_printf1 dassert_printf
#else
#define assert1(e)
#define assert_printf1(args...)
#endif
#if(MAX_ASSERT_LEVEL>=2)
#define assert2 dassert_
// #define assert dassert_
#define assert_printf2 dassert_printf
#define assert_printf  dassert_printf
#else
#define assert2(e)
// #define assert(e)
#define assert_printf2(args...)
#define assert_printf(args...)
#endif
#if(MAX_ASSERT_LEVEL>=3)
#define assert3 dassert_
#define assert_printf3 dassert_printf
#else
#define assert3(e)
#define assert_printf3(args...)
#endif

// asserting specific relationships

// void assert_error_double(const char* file, int line, const char* func,
// const char* op, const char* lhs_str, const char* rhs_str,
// double lhs, double rhs);

#define declare_assert_errmsg(t1,t2) \
  void assert_error(const char* file, int line, const char* func, \
      const char* op, const char* lhs_str, const char* rhs_str, \
      t1 lhs, t2 rhs);
declare_assert_errmsg(double, double);
declare_assert_errmsg(size_t, size_t);
declare_assert_errmsg(int, size_t);
declare_assert_errmsg(size_t, int);
declare_assert_errmsg(int, int);
declare_assert_errmsg(const char *, const char *);
// put in assert_string.h:
// #include "assert.h"
// #include<string>

extern "C" {
  int fcmp(double x1, double x2, double epsilon);
}
#ifndef USE_GCC_OPTIMIZATION
#define builtin_expect(a,b) (a)
#else
#define builtin_expect(a,b) __builtin_expect(a,b)
#endif
#define assert_not_almost_eq(lhs,rhs) \
  (fcmp(lhs, rhs, 1e-14) \
   ? (void)0 \
   : assert_error(__FILE__, __LINE__, __func__, " !=~= ", #lhs, #rhs, lhs, rhs))
#define assert_almost_eq(lhs,rhs) \
  (builtin_expect(fcmp(lhs, rhs, 1e-10),0) \
   ? assert_error(__FILE__, __LINE__, __func__, " =~= ", #lhs, #rhs, lhs, rhs) \
   : (void)0)
#define assert_divides(lhs,rhs) \
  (builtin_expect(rhs%lhs,0) \
   ? assert_error(__FILE__, __LINE__, __func__, "(divides)", #lhs, #rhs, lhs, rhs) \
   : (void)0)
#define assert_streq(lhs,rhs) \
  (builtin_expect(strcmp(lhs,rhs),0) \
   ? assert_error(__FILE__, __LINE__, __func__, "==", #lhs, #rhs, lhs, rhs) \
   : (void)0)
#define assert_op(op,lhs,rhs) \
  (builtin_expect(!(lhs op rhs),0) \
   ? assert_error(__FILE__, __LINE__, __func__, #op, #lhs, #rhs, lhs, rhs) \
   : (void)0)
// these implementations are much faster than if using // std::isnan(a) and std::isinf(a)// #include <float.h> // need for DBL_MAX if invoking assert_isfinite// check that the number is between -inf and inf// (so that this will work even in case -ffast-math is set)
#define assert_isfinite(a) \
  (builtin_expect(!(a>=-DBL_MAX && a<=DBL_MAX),0) ? \
   (void)(printf("ERROR: file %s, line %d, function %s:\n\t: %s = %24.16e " \
       "is not finite\n", __FILE__, __LINE__, __func__,#a,a),abort()) : (void) 0)
// // if -funsafe-math-optimizations is set this will fail to detect// nan for old gcc compilers. We could work around this by// implementing assert_isnum in an object module compiled without// -funsafe-math-optimizations, but this would incur unacceptable// function call overhead. So the user will only get checks// against nan for modules not compiled with this option.// 
#define report_nan(a) \
  (builtin_expect(!(a==a),0) ? \
   (printf("error: file %s, line %d, function %s:\n\t: %s = %24.16e " \
           "is nan\n", __FILE__, __LINE__, __func__,#a,a),1) : 0)
#define assert_isnum(a) \
  (builtin_expect(!(a==a),0) ? \
   (void)(printf("ERROR: file %s, line %d, function %s:\n\t: %s = %24.16e " \
       "is nan\n", __FILE__, __LINE__, __func__,#a,a),abort()) : (void) 0)
#define assert_eq(a,b) assert_op(==,a,b);
#define assert_ne(a,b) assert_op(!=,a,b);
#define assert_gt(a,b) assert_op(>,a,b);
#define assert_lt(a,b) assert_op(<,a,b);
#define assert_ge(a,b) assert_op(>=,a,b);
#define assert_le(a,b) assert_op(<=,a,b);
#endif                          // NDEBUG
#endif                          // __ASSERTS_H__
