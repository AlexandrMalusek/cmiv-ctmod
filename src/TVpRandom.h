#ifndef TVpRandom_h
#define TVpRandom_h

#include "TObject.h"

class TVpRandom
{
 protected:
  UInt_t   fSeed;  //Random number generator seed
  
 public:
  TVpRandom(UInt_t seed=65539);
  virtual ~TVpRandom();
  virtual  UInt_t   GetSeed() {return fSeed;}
  virtual  void     SetSeed(UInt_t seed=65539);
  virtual  Double_t Rand(Int_t i=0);

  ClassDef(TVpRandom,1)  // Random number generator, congruential
};

#endif // TVpRandom_h
