#include "misc.h"

TVpRandom3 rng;

Double_t getRand(Double_t min, Double_t max)
{
  return min + (max-min) * getRand();
}
