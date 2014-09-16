#ifndef TVpPitchHelical_h
#define TVpPitchHelical_h

#include "TObject.h"
#include <iosfwd>
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"
#include "TVpPitch.h"

class TVpPitchHelical : public TVpPitch
{
 public:
  Double_t     fPitchSize; // The size of the pitch
  
  TVpPitchHelical();
  TVpPitchHelical(Double_t pitchSize);
  virtual ~TVpPitchHelical();

  void GetActiveTransformation
    (Double_t rotationAngle, Int_t sliceIndex,
     TVpMatrix3x3& rotMatOld, TVpVector3& traVecOld,
     TVpMatrix3x3& rotMatNew, TVpVector3& traVecNew) const;
  void     PrintStatus(std::ostream &out = std::cout) const;

  ClassDef(TVpPitchHelical,1) // Helical trajectory of CT scanner components
};

//______________________________________________________________________________
inline TVpPitchHelical::TVpPitchHelical()
{
  // Default constructor
}

//______________________________________________________________________________
inline TVpPitchHelical::~TVpPitchHelical()
{
  // Destructor
}

#endif
