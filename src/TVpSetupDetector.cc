#include <math.h>
#include "TVpSetupDetector.h"

ClassImp(TVpSetupDetector)

//______________________________________________________________________________
TVpSetupDetector::TVpSetupDetector()
{
  // Default constructor

  fSource = 0;
  fGeometry = 0;
}

//______________________________________________________________________________
TVpSetupDetector::TVpSetupDetector(TVpSource *source, TVpGeometry *geometry)
{
  // Constructor with full initialization.

  fSource = source;
  fGeometry = geometry;
}

//______________________________________________________________________________
TVpSetupDetector::~TVpSetupDetector()
{
  // Destructor
}


//______________________________________________________________________________
void TVpSetupDetector::SetPosition(TVpVector3& helTraVec, TVpMatrix3x3& helRotMat)
{
  // Set the source and detector position.

  // Set the source position
  //fSource->SetDefaultPosition();
  //fSource->Move(helRotMat*fSrcRotMat, helTraVec+fSrcTraVec);
}



