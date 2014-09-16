//______________________________________________________________________________
//
// PitchHelical describes helical scanning mode of a CT scanner.  It provides
// transformation matrices which are used to move the X-ray source and the
// detector array from their default positions.

#include <math.h>
#include "TVpPitchHelical.h"
#include "TVpMatrix3x3.h"

ClassImp(TVpPitchHelical)

//______________________________________________________________________________
TVpPitchHelical::TVpPitchHelical(Double_t pitchSize)
{
  // Constructor with full initialization.
  //
  // Input parameters:
  // - pitchSize - pitch size
  
  fPitchSize = pitchSize;
}

//______________________________________________________________________________
void TVpPitchHelical::GetActiveTransformation
( Double_t rotationAngle, Int_t sliceIndex, TVpMatrix3x3& rotMatOld, TVpVector3& traVecOld,
  TVpMatrix3x3& rotMatNew, TVpVector3& traVecNew) const
{
  // Get active transformation

  if (sliceIndex != 0)
    std::cerr << "Warning: TVpPitchHelical::GetActiveTransformation: "
	      << "Any non-zero sliceIndex ignored.\n";
  
  TVpMatrix3x3 rotMat = rotationMatrix(TVpVector3(1,0,0), rotationAngle);
  TVpVector3 traVec = (fPitchSize * rotationAngle / (2 * M_PI)) * TVpVector3(1,0,0);
  traVecNew = rotMat * traVecOld + traVec;
  rotMatNew = rotMat * rotMatOld;
}

//______________________________________________________________________________
void TVpPitchHelical::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  // Input parameters:
  // - out - output stream (default = cout)

  out << "<TVpPitchHelical>\n"
      << "$Id: TVpPitchHelical.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "</TVpPitchHelical>\n";
}
