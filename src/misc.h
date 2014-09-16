#ifndef misc_h
#define misc_h

#include <stdlib.h>
#include "TVpRandom3.h"
#include "TObject.h"

// Miscelaneous functions

Double_t getRand();
Double_t getRand(Double_t min, Double_t max);

inline Double_t getRandTriang()
{
	double r1=getRand();
	double r2=getRand();
	return r1 > r2 ? r1 : r2;
}

extern TVpRandom3 rng;

inline Double_t getRand()
{
  return rng.Rand();
}

#endif // misc_h
