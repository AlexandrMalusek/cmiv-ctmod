#ifndef TVpSpectrum_h
#define TVpSpectrum_h

#include "TObject.h"
#include "TH1.h"
#include <iostream>

class TVpSpectrum
{
 protected:
  Double_t *fDiscDist;           //! Distribution function

  Double_t GetContRandEnergy() const;
  Double_t GetDiscRandEnergy() const;

 public:
  Char_t   *fName;               //! Descriptive name of the spectrum
  Int_t     fNumOfChannels;      //  Number of cont.+disc. channels
  Int_t     fNumOfContChannels;  //  Number of continuous channels
  Int_t     fNumOfDiscChannels;  //  Number of discrete channels
  Double_t  fContFraction;       //  Weight fraction of cont. channels
  Double_t  fContMaximum;        //  Maximum weight for cont. channels 
  Double_t *fContChannelWeight;  //! Array of continuous channel weights
  Double_t *fDiscChannelEnergy;  //! Array of discrete channel energies in keV
  Double_t *fDiscChannelWeight;  //! Array of discrete channel weights

  TVpSpectrum();
  TVpSpectrum(const Char_t *fileName);
  TVpSpectrum(Double_t energy);
  virtual ~TVpSpectrum();

  Int_t    ReadSpeFile(const Char_t *fileName);
  void     Normalize();
  void     GetEnergyAndWeight(Int_t channel, Double_t& energy, Double_t& weight) const;
  void     GetContEnergyAndWeight(Int_t channel, Double_t& energy, Double_t& weight) const;
  void     GetDiscEnergyAndWeight(Int_t channel, Double_t& energy, Double_t& weight) const;
  Int_t    GetNumOfChannels() const { return fNumOfChannels; }
  Int_t    GetNumOfContChannels() const { return fNumOfContChannels; }
  Int_t    GetNumOfDiscChannels() const { return fNumOfDiscChannels; }
  Double_t GetContWeight(Double_t energy) const { return fContChannelWeight[(int) energy]; }
  Double_t GetRandEnergy() const;
  Double_t GetPhi() const;
  Double_t GetPsi() const;
  void     PrintStatus(std::ostream &out = std::cout) const;

  TH1F   *GetHist();
  TH1F   *GetRandHist(Int_t nEvent);

  ClassDef(TVpSpectrum,1) // Energy spectrum of a source (discrete lines and continuum)
};

#endif // TVpSpectrum_h
