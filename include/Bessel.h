/***************************************************************************
  Bessel.h  -  
  -------------------
begin                : Tue Jan 30 2007

 ***************************************************************************/

#ifndef Bessel_H
#define Bessel_H

#include <iostream>
#include <math.h>
// #include "gsl/gsl_sf_bessel.h"


using std::cout;
using std::endl;

/*! \brief Calculate Bessel functions of the first kind J_n \param lambda_s Argument of function \param nmax Sequence J_n(x) for n=0 ... nmax+1 is calculated \param bessel_Jn_array Array of values for J_n corresponding to n=0 ... nmax+1 \param bessel_Jn_prime_array Array of values for J_n_prime ((ie derivate of J_n) for n=0 ... nmax \note J_n array goes to n=nmax + 1, but the J_n_prime array goes to nmax only. \note This function uses the corresponding call in the GSL. */
void calc_bessel_Jn_seq(double lambda, int nmax, double bessel_Jn_array[], double bessel_Jn_prime_array[]) {
  int gsl_fn_result = 0.0;
  // CHANGE IT HERE: PUT HERE THE BESSEL FUNCTION
  // gsl_fn_result = gsl_sf_bessel_Jn_array( 0, nmax+1, lambda, bessel_Jn_array );

  if (gsl_fn_result)
    cout << "Error in calc_bessel_Jn_seq !!!!" << endl;

  bessel_Jn_prime_array[0] = -bessel_Jn_array[1];

  for (int i = 1; i <= nmax; ++i)
    bessel_Jn_prime_array[i] = 0.5 * (bessel_Jn_array[i - 1] - bessel_Jn_array[i + 1]);

}


#endif
