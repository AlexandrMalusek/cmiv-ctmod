#ifndef TVpPitchSequential_h
#define TVpPitchSequential_h

#include "TObject.h"
#include <iosfwd>
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"
#include "TVpPitch.h"

class TVpPitchSequential : public TVpPitch
{
 public:
  Double_t     fPitchSize; // The size of the pitch
  
  TVpPitchSequential();
  TVpPitchSequential(Double_t pitchSize);
  virtual ~TVpPitchSequential();

  void GetActiveTransformation
    (Double_t rotationAngle, Int_t sliceIndex,
     TVpMatrix3x3& rotMatOld, TVpVector3& traVecOld,
     TVpMatrix3x3& rotMatNew, TVpVector3& traVecNew) const;
  void     PrintStatus(std::ostream &out = std::cout) const;

  ClassDef(TVpPitchSequential,1) // Sequential mode trajectory of CT scanner components
};

//______________________________________________________________________________
inline TVpPitchSequential::TVpPitchSequential()
{
  // Default constructor
}

//______________________________________________________________________________
inline TVpPitchSequential::~TVpPitchSequential()
{
  // Destructor
}

#endif
