#ifndef TVpPitch_h
#define TVpPitch_h

#include "TObject.h"
#include <iosfwd>
#include "TVpVector3.h"
#include "TVpMatrix3x3.h"

class TVpPitch
{
 public:
  TVpPitch();
  virtual ~TVpPitch();

  virtual void GetActiveTransformation
    (Double_t rotationAngle, Int_t sliceIndex,
     TVpMatrix3x3& rotMatOld, TVpVector3& traVecOld,
     TVpMatrix3x3& rotMatNew, TVpVector3& traVecNew) const = 0;
  virtual void PrintStatus(std::ostream &out = std::cout) const = 0;

  ClassDef(TVpPitch,1) // The trajectory of tomograph components, abstract class.
};

//______________________________________________________________________________
inline TVpPitch::TVpPitch()
{
  // Default constructor
}

//______________________________________________________________________________
inline TVpPitch::~TVpPitch()
{
  // Destructor
}

#endif
