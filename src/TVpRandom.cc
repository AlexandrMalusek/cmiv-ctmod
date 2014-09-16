//______________________________________________________________________________
//
// Random number generator taken from ROOT.

#include <time.h>
#include "TVpRandom.h"

ClassImp(TVpRandom)

//______________________________________________________________________________
TVpRandom::TVpRandom(UInt_t seed)
{
  // Default constructor

  SetSeed(seed);
}

//______________________________________________________________________________
TVpRandom::~TVpRandom()
{
  // Default destructor

}

//______________________________________________________________________________
Double_t TVpRandom::Rand(Int_t i)
{
  //  Machine independent random number generator.
  //  Produces uniformly-distributed floating points between 0 and 1.
  //  Identical sequence on all machines of >= 32 bits.
  //  Periodicity = 10**8
  //  Universal version (Fred james 1985).
  //  generates a number in [0,1]

   const Float_t kCONS = 4.6566128730774E-10;
   const Int_t kMASK31 = 2147483647;

   fSeed *= 69069;
   // keep only lower 31 bits
   fSeed &= kMASK31;
   // Set lower 8 bits to zero to assure exact float
   Int_t jy = (fSeed/256)*256;
   Float_t random = kCONS*jy;
   return Double_t(random);
}

//______________________________________________________________________________
void TVpRandom::SetSeed(UInt_t seed)
{
  //  Set the random generator seed
  //  if seed is zero, the seed is set to the current  machine clock
  //  Note that the machine clock is returned with a precision of 1 second.
  //  If one calls SetSeed(0) within a loop and the loop time is less than 1s, 
  //  all generated numbers will be identical!
  
  if( seed == 0 )
    {
      time_t curtime;      // Set 'random' seed number  if seed=0
      time(&curtime);      // Get current time in fSeed.
      fSeed = (UInt_t)curtime;
    }
  else
    fSeed = seed;
}
