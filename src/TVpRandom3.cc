//______________________________________________________________________________
//
// Random number generator taken from ROOT.
//
// Random number generator class based on
//   M. Matsumoto and T. Nishimura,
//   Mersenne Twister: A 623-diminsionally equidistributed
//   uniform pseudorandom number generator
//   ACM Transactions on Modeling and Computer Simulation,
//   Vol. 8, No. 1, January 1998, pp 3--30.
//
// For more information see the Mersenne Twister homepage
//   http://www.math.keio.ac.jp/~matumoto/emt.html
//
// Advantage: large period 2**19937-1
//            relativly fast
//              (only two times slower than TRandom, but
//               two times faster than TRandom2)
// Drawback:  a relative large internal state of 624 integers
//
//
// Aug.99 ROOT implementation based on CLHEP by P.Malzacher
//
// the original code contains the following copyright notice:
// /* This library is free software; you can redistribute it and/or   */
// /* modify it under the terms of the GNU Library General Public     */
// /* License as published by the Free Software Foundation; either    */
// /* version 2 of the License, or (at your option) any later         */
// /* version.                                                        */
// /* This library is distributed in the hope that it will be useful, */
// /* but WITHOUT ANY WARRANTY; without even the implied warranty of  */
// /* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.            */
// /* See the GNU Library General Public License for more details.    */
// /* You should have received a copy of the GNU Library General      */
// /* Public License along with this library; if not, write to the    */
// /* Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA   */
// /* 02111-1307  USA                                                 */
// /* Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.       */
// /* When you use this, send an email to: matumoto@math.keio.ac.jp   */
// /* with an appropriate reference to your work.                     */

#include "TVpRandom3.h"

ClassImp(TVpRandom3)

//______________________________________________________________________________
TVpRandom3::TVpRandom3(UInt_t seed)
{
  // Default constructor

  SetSeed(seed);
}

//______________________________________________________________________________
TVpRandom3::~TVpRandom3()
{
  // Default destructor
}

//______________________________________________________________________________
Double_t TVpRandom3::Rand(Int_t)
{
  //  Machine independent random number generator.
  //  Produces uniformly-distributed floating points between 0 and 1.
  //  Method: Mersenne Twistor
  
  UInt_t y;
  
  const Int_t  kM = 397;
  const Int_t  kN = 624;
  const UInt_t kTemperingMaskB =  0x9d2c5680;
  const UInt_t kTemperingMaskC =  0xefc60000;
  const UInt_t kUpperMask =       0x80000000;
  const UInt_t kLowerMask =       0x7fffffff;
  const UInt_t kMatrixA =         0x9908b0df;
  
  //const UInt_t kMAG01[2] = {0x0,kMatrixA};
  
  if (fCount624 >= kN)
    {
      register Int_t i;
      
      for (i=0; i < kN-kM; i++)
	{
	  y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
	  fMt[i] = fMt[i+kM] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
	  //fMt[i] = fMt[i+kM] ^ (y >> 1) ^ kMAG01[y&0x1];
	}
      
      for (   ; i < kN-1    ; i++)
	{
	  y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
	  fMt[i] = fMt[i+kM-kN] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
	  //fMt[i] = fMt[i+kM-kN] ^ (y >> 1) ^ kMAG01[y&0x1];
	}
      
      y = (fMt[kN-1] & kUpperMask) | (fMt[0] & kLowerMask);
      fMt[kN-1] = fMt[kM-1] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
      //fMt[kN-1] = fMt[kM-1] ^ (y >> 1) ^ kMAG01[y & 0x1];
      fCount624 = 0;
    }

  y = fMt[fCount624++];
  y ^=  (y >> 11);
  y ^= ((y << 7 ) & kTemperingMaskB );
  y ^= ((y << 15) & kTemperingMaskC );
  y ^=  (y >> 18);
  
  return ( (Double_t) y * 2.3283064365386963e-10); // * Power(2,-32)
}

//______________________________________________________________________________
void TVpRandom3::SetSeed(UInt_t seed)
{
  //  Set the random generator sequence
  
  TVpRandom::SetSeed(seed);
  fCount624 = 624;
  fMt[0] = fSeed;
  Int_t i;
  for(i=1; i<624; i++) {
    fMt[i] = (69069 * fMt[i-1]) & 0xffffffff;
  }
}
