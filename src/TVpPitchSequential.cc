//______________________________________________________________________________
//
// The TVpPitchSequential class inmplements the sequential mode of a tomograph
// source and detector trajectory.

#include <math.h>
#include "TVpPitchSequential.h"
#include "TVpMatrix3x3.h"

ClassImp(TVpPitchSequential)

//______________________________________________________________________________
TVpPitchSequential::TVpPitchSequential(Double_t pitchSize)
{
  // Constructor with full initialization.
  //
  // Input parameters:
  // - pitchSize - pitch size
  
  fPitchSize = pitchSize;
}

//______________________________________________________________________________
void TVpPitchSequential::GetActiveTransformation
( Double_t rotationAngle, Int_t sliceIndex, TVpMatrix3x3& rotMatOld, TVpVector3& traVecOld,
  TVpMatrix3x3& rotMatNew, TVpVector3& traVecNew) const
{
  // Get active transformation

  TVpMatrix3x3 rotMat = rotationMatrix(TVpVector3(1,0,0), rotationAngle);
  TVpVector3 traVec = (sliceIndex * fPitchSize) * TVpVector3(1,0,0);
  traVecNew = rotMat * traVecOld + traVec;
  rotMatNew = rotMat * rotMatOld;
}

//______________________________________________________________________________
void TVpPitchSequential::PrintStatus(std::ostream &out) const
{
  // Print the object status

  out << "<TVpPitchSequential>\n"
      << "$Id: TVpPitchSequential.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "</TVpPitchSequential>\n";
}
