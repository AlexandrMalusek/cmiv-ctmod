#ifndef TVpVolumeIntegrator_h
#define TVpVolumeIntegrator_h

#include "TVpVector3.h"

class TVpVolumeIntegrator
{
 public:
  Int_t fNumOfVolumeElements;

  TVpVolumeIntegrator();
  virtual ~TVpVolumeIntegrator();
  void GetVolumeElement(Int_t index, TVpVector3 &positionU, Double_t &volume);
  inline Int_t GetNumOfVolumeElements();
  
  ClassDef(TVpVolumeIntegrator,1) // Volume integrator
};

//______________________________________________________________________________
inline Int_t TVpVolumeIntegrator::GetNumOfVolumeElements()
{
  return fNumOfVolumeElements;
}

#endif  // TVpVolumeIntegrator_h
