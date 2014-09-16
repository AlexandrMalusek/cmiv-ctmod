#ifndef TVpConstant_h
#define TVpConstant_h

#include "TObject.h"

class TVpConstant
{
 public:
  static Double_t hc;                           // = h * c  [keV cm]
  static Double_t fElectronRestEnergy;          // [keV]
  static Double_t fClassicalElectronRadius;     // [cm]
  static Double_t fClassicalElectronRadius22;   // [cm^2]

  static Double_t fBarn;                        // barn to cm^2 conversion factor
  static Double_t fMicroMeter;                  // micrometer to cm conversion

  virtual ~TVpConstant() {};

  ClassDef(TVpConstant,0) // Physical constants
};

#endif  // TVpConstant_h
