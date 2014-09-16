//______________________________________________________________________________
//
// TVpIntegral performs numerical integration of a member function.

#include <cmath>
#include <iostream>
#include "TVpIntegral.h"

ClassImp(TVpIntegral)

//______________________________________________________________________________
Double_t TVpIntegral::Integral(Int_t selector, Double_t par, Double_t a, Double_t b)
{
  // Integrate from a to b.  Use the extended trapezoidal rule, see
  // e.g. Numerical Recipes in C.
  //
  // x_i = a + i * h,  i = 0, ..., n
  // x_0 = a,  x_n = b
  // I = h * (0.5*f_0 + f_1 + ... + f_{n-1} + 0.5*f_n)

  Double_t h = (b - a) / fNdiv;
  Double_t sum = 0.5 * (EvaluateIntegrand(selector, par, a) 
			+ EvaluateIntegrand(selector, par, b));
  for (Int_t i = 1; i < fNdiv; i++)
    {
      Double_t x = a + i * h;
      sum += EvaluateIntegrand(selector, par, x);
    }
  return h * sum;
}
