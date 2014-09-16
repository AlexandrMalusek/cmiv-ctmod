#ifndef TVpRandom3_h
#define TVpRandom3_h

#include "TObject.h"
#include "TVpRandom.h"

class TVpRandom3 : public TVpRandom
{
 private:
  UInt_t   fMt[624];
  Int_t    fCount624;
  
 public:
  TVpRandom3(UInt_t seed=65539);
  virtual ~TVpRandom3();
  virtual  Double_t  Rand(Int_t i=0);
  virtual  void      SetSeed(UInt_t seed=0);
  
  ClassDef(TVpRandom3,1) // Random number generator, Mersenne twister, period 2**19937-1
};

#endif // TVpRandom3_h
