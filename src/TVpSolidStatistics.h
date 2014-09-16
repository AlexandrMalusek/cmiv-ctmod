#ifndef TVpSolidStatistics_h
#define TVpSolidStatistics_h

#include "TObject.h"
#include <iostream>

class TVpSolidStatistics
{
 public:
  Double_t     fEnergyImpartedPh;    //  Energy imparted by photoeffect
  Double_t     fEnergyImpartedIn;    //  Energy imparted by incoherent (Compton) scattering
  Double_t     fEnergyImpartedTr;    //  Energy imparted by the track end
  ULong_t      fNumOfInterPh;        //  Number of photoeffects
  ULong_t      fNumOfInterIn;        //  Number of incoherent scatterings
  ULong_t      fNumOfInterCo;        //  Number of coherent scatterings
  ULong_t      fNumOfInterTr;        //  Number of track ends
  Double_t     fNumOfParticlesIn;    //  Number of particles entering the solid
  Double_t     fNumOfParticlesOut;   //  Number of particles leaving the solid

  TVpSolidStatistics();
  virtual ~TVpSolidStatistics();
  void     PrintStatistics(std::ostream &out = std::cout, Long_t numOfHist = 1) const;
  Double_t GetEnergyImparted(Long_t numOfHist = 1) const;

  ClassDef(TVpSolidStatistics,1) // Solid specific MC simulation statistics
};

#endif  // TVpSolidStatistics_h
