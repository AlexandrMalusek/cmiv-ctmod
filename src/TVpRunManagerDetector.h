#ifndef TVpRunManagerDetector_h
#define TVpRunManagerDetector_h

#include "TObject.h"
#include "TVpRunManager.h"
#include "TVpGeometry.h"
#include "TVpSetupDetector.h"

class TVpRunManagerDetector : public TVpRunManager
{
 private:
  TVpGeometry      *fGeometryPtr;              //! the same as fSetup->fGeometry
  
 public:
  TVpSetupDetector *fSetupDetectorPtr;         //! Detector setup
  Double_t          fPhotonWeight;             //! Weight cutoff

  TVpRunManagerDetector();
  TVpRunManagerDetector(TVpSetupDetector *setupDetectorPtr);
  virtual ~TVpRunManagerDetector();

  Int_t     PhotonLife(TVpPhoton *photonPtr);

 public:
  ClassDef(TVpRunManagerDetector,1) // Run manager, detector simulation setup
};

#endif  // TVpRunManagerDetector_h
