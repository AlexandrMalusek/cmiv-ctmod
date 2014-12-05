#ifndef TVpPointDetector_h
#define TVpPointDetector_h

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TVpVector3.h"
#include "TVpMatrix3x3.h"
#include "TVpParticle.h"
#include "TVpPhoton.h"
#include "TVpGeometry.h"
#include "TVpDetectorResponse.h"
#include "dedEvent.h"

#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

extern DedEvent dedEvent;

class TVpPointDetector : public TVpObjectLocation
{
 private:
  Double_t    fAc2e;             // linear transform, energy = fAc2e * channel + fBc2e
  Double_t    fBc2e;             // linear transform, energy = fAc2e * channel + fBc2e
  Double_t    fAe2c;             // linear transform, channel = fAe2c * energy + fBe2c
  Double_t    fBe2c;             // linear transform, channel = fAe2c * energy + fBe2c
  Int_t       fUnderflowChannel; // fUnderflowChannel = fNumOfChannels + 0
  Int_t       fOverflowChannel;  // fOverflowChannel = fNumOfChannels + 1

  inline void      UpdateCounters(Int_t channel, Double_t value, ULong_t history);

 public:
  enum EScoredQuantity {kPlanarEnergyFluence = 0, kEnergyFluence = 1,
			kDetectorResponse = 2, kAirKerma = 3,
			kFluence = 4, kPlanarFluence = 5};
  
  Int_t       fNumOfChannels;      // Number of channels
  Double_t    fMinEnergy;          // Lower energy cutoff
  Double_t    fMaxEnergy;          // Upper energy cutoff
  Double_t   *fChannelContent;     //! Array of channel sum
  Double_t   *fChannelContent2;    //! Array of channel sum of squares
  Double_t   *fChannelCounter;     //! Array of channel counters, see class description
  ULong_t    *fChannelLastHistory; //! Array of last histories which contributed
  ULong_t     fNumOfHistories;     // Number of histories
  TVpDetectorResponse *fDetectorResponse; //! The detector response function class
  EScoredQuantity fScoredQuantity; // Scored quantity
  TTree      *fDedEventTreePtr;    //! Pointer to the DED tree
  TFile      *fDedEventFilePtr;    //! Pointer to the DED file
  
  TVpPointDetector();
  virtual ~TVpPointDetector();

  inline Int_t     GetChannel(Double_t energy) const;
  inline Double_t  GetEnergy(Int_t channel) const;
  inline void      GetPosition(TVpVector3 *position) const;
  void             SetNumOfChannels(Int_t numOfChannels,
				    Double_t minEnergy, Double_t maxEnergy);
  inline void      SetDetectorResponse(TVpDetectorResponse *detectorResponse);
  inline void      SetScoredQuantity(EScoredQuantity scoredQuantity);
  inline void      RegisterScoredQuantity(TVpParticle *particle);
  void             RegisterVirtualPhoton(TVpGeometry *geometry, TVpPhoton *photonPtr);
  void             RegisterParticle(TVpParticle *particle);
  void             MultiplyChannelContentBy(Double_t factor);
  Double_t         GetChannelMeanEstimate(Int_t channel) const;
  Double_t         GetChannelVarianceEstimate(Int_t channel) const;
  Double_t         GetChannelSigmaEstimate(Int_t channel) const;
  Double_t         GetMeanEstimate() const;
  Double_t         GetVarianceEstimate() const;
  Double_t         GetSigmaEstimate() const;
  inline ULong_t   GetNumberOfHistories() const;
  Int_t            WritePDS(FILE *fp, Int_t withErrors = 0);
  Int_t            ReadPDS(FILE *fp, Int_t withErrors = 0, ULong_t numOfHistories = 1);
  void             Zero();
  void             UpdateCountersAtEndOfHistory(ULong_t numOfHistories);
  virtual void     PrintStatus(std::ostream &out = std::cout) const;
  Double_t         GetError(Int_t channel) const;

  enum EDrawOption {kDrawPoint = 1, kDrawAxes = 2};
  void             Draw(Int_t option = kDrawPoint) const;
  void             DedEventInitialize(const Char_t *fileName);
  void             DedEventClose();
  TH1F            *GetHistEnergySpectrum() const;

  ClassDef(TVpPointDetector,1)   // Point detector
};

//______________________________________________________________________________
inline Int_t TVpPointDetector::GetChannel(Double_t energy) const
{
  // Return the channel number corresponding to a given photon's energy.  If
  // energy is out of the range then return either fUnderflowChannel or
  // fOverflowChannel.  Linear interpolation is used.

  if (energy < fMinEnergy)
    return fUnderflowChannel;
  if (energy > fMaxEnergy)
    return fOverflowChannel;
  return (int) (fAe2c * energy + fBe2c);
}

//______________________________________________________________________________
inline Double_t TVpPointDetector::GetEnergy(Int_t channel) const
{
  // Return the energy corresponding to a channel number.

  if (channel < 0 || channel > fNumOfChannels)
    {
      fprintf(stderr, "TVpPointDetector::GetEnergy: channel %d out of range\n", channel);
      exit(1);
    }
  
  return (double)(fAc2e * channel + fBc2e);
}

//______________________________________________________________________________
inline void TVpPointDetector::GetPosition(TVpVector3 *position) const
{
  *position = fTraVecL2u;
}

//______________________________________________________________________________
inline void TVpPointDetector::SetDetectorResponse(TVpDetectorResponse *detectorResponse)
{
  // Set the detector response function.  Also set the scored quantity.

  fDetectorResponse = detectorResponse;
  fScoredQuantity = kDetectorResponse;
}

//______________________________________________________________________________
inline void TVpPointDetector::SetScoredQuantity(EScoredQuantity scoredQuantity)
{
  // Set the scored quantity.

  fScoredQuantity = scoredQuantity;
}

//______________________________________________________________________________
inline ULong_t TVpPointDetector::GetNumberOfHistories() const
{
  // Return the number of histories
  
  return fNumOfHistories;
}

//______________________________________________________________________________
inline void TVpPointDetector::UpdateCounters(Int_t channel, Double_t value, ULong_t history)
{
  // Update counters.  This routine is used by RegisterScoredQuantity().

  if (history != fChannelLastHistory[channel])
    {  // New history - update counters
      fChannelLastHistory[channel] = history;
      Double_t x = fChannelCounter[channel];
      fChannelContent[channel] += x;
      fChannelContent2[channel] += x * x;
      fChannelCounter[channel] = value;
    }
  else
    fChannelCounter[channel] += value;

  if (fDedEventTreePtr != 0)
    {
      dedEvent.c = value;
      fDedEventTreePtr->Fill();
    }
}

//______________________________________________________________________________
inline void TVpPointDetector::RegisterScoredQuantity(TVpParticle *particle)
{
  // Register the scored quantity.

  // Convert particle's direction from LCS of the solid to UCS.
  TVpVector3 particleDirU = particle->fSolid->DirLocToUni(particle->fDirection);

  // Convert the particle's direction from UCS to LCS of the point detector.
  TVpVector3 particleDirL = DirUniToLoc(particleDirU);

  Double_t xi = std::abs(particleDirU.GetZ());  // cosine of the incidence angle
  Double_t energy = particle->GetEnergy();
  Int_t channel = GetChannel(energy);
  ULong_t history = particle->fHistoryNumber;
  Double_t value;

  switch (fScoredQuantity) 
    {
    case kPlanarEnergyFluence:
      value = particle->GetWeight() * xi * energy;
      break;
    case kEnergyFluence:
      value = particle->GetWeight() * energy;
      break;
    case kDetectorResponse:
      value = particle->GetWeight() * xi * energy
	* fDetectorResponse->GetResponse(energy, particleDirL);
      break;
    case kFluence:
      value = particle->GetWeight();
      break;
    case kPlanarFluence:
      value = particle->GetWeight() * xi;
      break;
    case kAirKerma:
      value = particle->GetWeight() * energy 
	* fDetectorResponse->GetResponse(energy, particleDirL);
      break;
    }
  
  if (fDedEventTreePtr != 0)
    dedEvent.a = xi;

  UpdateCounters(channel, value, history);
}

#endif // TVpPointDetector
