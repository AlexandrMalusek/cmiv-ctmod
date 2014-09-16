#ifndef TVpIntegral_h
#define TVpIntegral_h

#include "TObject.h"

class TVpIntegral
{
 public:
  Int_t  fNdiv;  // Number of divisions

  TVpIntegral() {fNdiv = 1024;};
  virtual ~TVpIntegral() {};
  Double_t         Integral(Int_t selector, Double_t par, Double_t a, Double_t b);
  virtual Double_t EvaluateIntegrand(Int_t selector, Double_t par, Double_t x) = 0;
  
  ClassDef(TVpIntegral,0)   // Integrator implementation
};

#endif  // TVpIntegral_h
