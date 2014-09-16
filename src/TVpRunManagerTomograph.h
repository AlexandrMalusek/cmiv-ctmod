#ifndef TVpRunManagerTomograph_h
#define TVpRunManagerTomograph_h

#include "TObject.h"
#include "TVpRunManager.h"
#include "TVpGeometry.h"
#include "TVpSetupTomograph.h"

class TVpRunManagerTomograph : public TVpRunManager
{
 private:
  TVpGeometry        *fGeometryPtr;            //! the same as fTomograph->fGeometry

 public:
  TVpSetupTomograph  *fSetupTomographPtr;      //! Tomograph setup

  TVpRunManagerTomograph();
  TVpRunManagerTomograph(TVpSetupTomograph *setupTomographPtr);
  virtual ~TVpRunManagerTomograph();

  Int_t     PhotonLife(TVpPhoton *photonPtr);
  void      Run(Long_t numOfHistories);
  void      FillFomEvent(UInt_t time) const;

 public:
  ClassDef(TVpRunManagerTomograph,1) // Run manager, tomograph simulation setup
};

#endif  // TVpRunManagerTomograph_h
