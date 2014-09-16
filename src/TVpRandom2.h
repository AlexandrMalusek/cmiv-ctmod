#ifndef TVpRandom2_h
#define TVpRandom2_h

#include "TObject.h"
#include "TVpRandom.h"

class TVpRandom2 : public TVpRandom
{
 protected:
  Double_t   fSeed1;  //Random number generator seed 1
  Double_t   fSeed2;  //Random number generator seed 2
  
 public:
  TVpRandom2(UInt_t seed=65539);
  virtual ~TVpRandom2();
  virtual  void     GetSeed2(UInt_t &seed1, UInt_t &seed2);
  virtual  Double_t Rand(Int_t i=0);
  virtual  void     SetSeed(UInt_t seed=0);
  virtual  void     SetSeed2(UInt_t seed1, UInt_t seed2);

  ClassDef(TVpRandom2,1)  // Random number generator, periodicity > 10**14
};

#endif // TVpRandom2_h
