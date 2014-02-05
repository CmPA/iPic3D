
#ifndef NO_MPI
  #include "MPIdata.h" // for get_rank
#endif
#include "ompdefs.h" // for omp_get_thread_num
#include <iostream>
#include "asserts.h"

void assert_error(const char *file, int line, const char *func, const char *op, const char *lhs_str, const char *rhs_str, double lhs, double rhs) {

  eprintf_fileLine(stdout, "ERROR", func,file,line,
    "\n\tassertion failed: %s %s %s, i.e., %24.16e %s %24.16e\n", lhs_str, op, rhs_str, lhs, op, rhs);
  abort();
}

#ifndef NO_MPI
  #ifdef _OPENMP
    #define process_string \
      std::cerr << "(" << MPIdata::get_rank() << "." <<  omp_get_thread_num() << ")";
  #else
    #define process_string \
      std::cerr << "(" << MPIdata::get_rank() << ")";
  #endif
#else
  #ifdef _OPENMP
    #define process_string \
      std::cerr << "(." << omp_get_thread_num() << ")";
  #else
    #define process_string 
  #endif
#endif

 #define implement_assert_errmsg(t1,t2) \
   void assert_error(const char* file, int line, const char* func, \
     const char* op, const char* lhs_str, const char* rhs_str, \
     t1 lhs, t2 rhs) \
   { \
     process_string \
     std::cerr << " ERROR in file " << file << ", line " << line  \
       << ", function " << func  \
       <<"\n\tassertion failed: " << lhs_str << op << rhs_str \
       << ", i.e., " << lhs << op << rhs << std::endl; \
       abort(); \
   }

implement_assert_errmsg(size_t, size_t);
implement_assert_errmsg(int, size_t);
implement_assert_errmsg(size_t, int);
implement_assert_errmsg(int, int);
implement_assert_errmsg(long long, long long);
implement_assert_errmsg(const char *, const char *);

/*
 fcmp
 Copyright (c) 1998-2000 Theodore C. Belding
 University of Michigan Center for the Study of Complex Systems
 <mailto:Ted.Belding@umich.edu>
 <http://www-personal.umich.edu/~streak/>		

 This file is part of the fcmp distribution. fcmp is free software;
 you can redistribute and modify it under the terms of the GNU Library
 General Public License (LGPL), version 2 or later.  This software
 comes with absolutely no warranty. See the file COPYING for details
 and terms of copying.

 File: fcmp.h 

 Description:
 
 Knuth's floating point comparison operators, from:
 Knuth, D. E. (1998). The Art of Computer Programming.
 Volume 2: Seminumerical Algorithms. 3rd ed. Addison-Wesley.
 Section 4.2.2, p. 233. ISBN 0-201-89684-2.

 Input parameters:
 x1, x2: numbers to be compared
 epsilon: determines tolerance

 epsilon should be carefully chosen based on the machine's precision,
 the observed magnitude of error, the desired precision, and the
 magnitude of the numbers to be compared. See the fcmp README file for
 more information.

 This routine may be used for both single-precision (float) and
 double-precision (double) floating-point numbers.
 
 Returns:
 -1 if x1 < x2
  0 if x1 == x2
  1 if x1 > x2
*/

/*
 fcmp
 Copyright (c) 1998-2000 Theodore C. Belding
 University of Michigan Center for the Study of Complex Systems
 <mailto:Ted.Belding@umich.edu>
 <http://www-personal.umich.edu/~streak/>		

 This file is part of the fcmp distribution. fcmp is free software;
 you can redistribute and modify it under the terms of the GNU Library
 General Public License (LGPL), version 2 or later.  This software
 comes with absolutely no warranty. See the file COPYING for details
 and terms of copying.

 File: fcmp.c

 Description: see fcmp.h and README files.
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>

int fcmp(double x1, double x2, double epsilon)
{
  double diff = x1-x2;
  if(diff>epsilon) return 1;
  if(diff<-epsilon) return -1;
  return 0;

  // the code below was failing for some reason. -eaj

  int exponent;
  double delta;
  double difference;
  
  /* Get exponent(max(fabs(x1), fabs(x2))) and store it in exponent. */

  /* If neither x1 nor x2 is 0, */
  /* this is equivalent to max(exponent(x1), exponent(x2)). */

  /* If either x1 or x2 is 0, its exponent returned by frexp would be 0, */
  /* which is much larger than the exponents of numbers close to 0 in */
  /* magnitude. But the exponent of 0 should be less than any number */
  /* whose magnitude is greater than 0. */
  
  /* So we only want to set exponent to 0 if both x1 and */
  /* x2 are 0. Hence, the following works for all x1 and x2. */

  frexp(fabs(x1) > fabs(x2) ? x1 : x2, &exponent);

  /* Do the comparison. */

  /* delta = epsilon * pow(2, exponent) */

  /* Form a neighborhood around x2 of size delta in either direction. */
  /* If x1 is within this delta neighborhood of x2, x1 == x2. */
  /* Otherwise x1 > x2 or x1 < x2, depending on which side of */
  /* the neighborhood x1 is on. */
  
  delta = ldexp(epsilon, exponent); 
  
  difference = x1 - x2;

  if (difference > delta)
    return 1; /* x1 > x2 */
  else if (difference < -delta) 
    return -1;  /* x1 < x2 */
  else /* -delta <= difference <= delta */
    return 0;  /* x1 == x2 */
}

