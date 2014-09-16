//______________________________________________________________________________
//
// A point detector uses the collision density estimator variance reduction
// technique to score one of the following quantities:
//
// - fluence (kFluence)
// - plane fluence (kPlanarFluence)
// - energy fluence (kEnergyFluence)
// - plane energy fluence (kPlanarEnergyFluence)
// - air kerma (kAirKerma)
// - energy imparted per unit surface area (kDetectorResponse)
//
// The selection is made by SetScoredQuantity().  Scored data are histogrammed
// by photon energy in an internal array.
//
// *******
// Method:
// *******
//
// In a Monte Carlo simulation, each coherent or incoherent scattering
// contributes to the scored quantity via RegisterVirtualPhoton.  At the
// interaction point, a virtual photon heading towards the point detector is
// created and transported to the point detector.
//
// A point detector is located in the point fPosition and the normal to its
// surface (a unit vector) is fDirectionZ. It has fNumOfChannels energy
// channels, and one underflow and one overflow channel. Relation between the
// photon energy and the channel number is shown in the following figure:
//
// fMinEnergy                  fMaxEnergy
// |                                    |
// | 0 | 1 | 2 | ... | fNumOfChannels-1 | fNumOfChannels | fNumOfChannels+1 |
//                                        underflow        overflow
//
// If the energy of the scored photon is lower than fMinEnergy then the
// contribution is scored in the underflow channel.  If the energy of the
// scored photon is greater than fMaxEnergy then the contribution is scored to
// the overflow channel.
//
// Calculation of error: fChannelCounter accumulates contributions from a
// single history.  If a new history has been started then its value is added
// to fChannelContent, its square is added to fChannelContent2, and
// fChannelCounter is zeroed.  UpdateCountersAtEndOfHistory() adds the content
// of fChannelCounter to fChannelContent, its square to fChannelContent2 and
// finally zeroes fChannelCounter.  This function is used when mean and
// variance estimates are needed in the middle and at the end of the run.
//
// Details:
//
// The default scored quantity is plane energy fluence.  The function
// SetDetectorResponse() switches to energy imparted per unit surface area.
// If air kerma is scored then TVpDetectorResponse must contain values of
// (muEn/rho)_air.
//
//
// *******************************
// DED (Detector Event Data) files
// *******************************
//
// The DED file structure is defined in dedEvent.h as follows:
//
// typedef struct {
//  Int_t   i;    //  Interaction type
//  Float_t x0;   //  x[0]
//  Float_t x1;   //  x[1]
//  Float_t x2;   //  x[2]
//  Float_t u0;   //  u[0]
//  Float_t u1;   //  u[1]
//  Float_t u2;   //  u[2]
//  Float_t e;    //  Energy
//  Float_t w;    //  Weight
//  Float_t c;    //  Contribution
//  Float_t a;    //  xi = cos(theta), where theta is the incidence angle
//  Int_t s;      //  Solid index
//  Int_t t;      //  Solid SubIndex (e.g. tissue number)
//  Int_t h;      //  History number
// } DedEvent;
//
// Interaction types are defined in TVpParticle.h:
// 3 - kComptonScattering
// 4 - kCoherentScattering
// Other interactions do not contribute to the point detector.

#include <math.h>
#include <stdio.h>
#include "TFile.h"
#include "TVpPointDetector.h"

ClassImp(TVpPointDetector)

DedEvent dedEvent;

//______________________________________________________________________________
TVpPointDetector::TVpPointDetector()
{
  // Default constructor. The number of channels is set to 0.  All
  // related data members (fMinEnergy, fAc2e, fDetectorResponse, ...)
  // have no meaning and are set to 0.

  fNumOfChannels = 0;
  fUnderflowChannel = fOverflowChannel = 0;
  fMinEnergy = fMaxEnergy = 0.0;
  fChannelContent = fChannelContent2 = fChannelCounter = 0;
  fChannelLastHistory = 0;
  fAc2e = fBc2e = fAe2c = fBe2c = 0.0;
  fDetectorResponse = 0;
  fNumOfHistories = 0;
  fScoredQuantity = kPlanarEnergyFluence;
  fDedEventTreePtr = 0;
  fDedEventFilePtr = 0;
}

//______________________________________________________________________________
TVpPointDetector::~TVpPointDetector()
{
  // Destructor. Delete all allocated arrays, do not delete fDetectorResponse
  // as it may be shared by more point detectors.

  delete [] fChannelContent;
  delete [] fChannelContent2;
  delete [] fChannelCounter;
  delete [] fChannelLastHistory;
}

//______________________________________________________________________________
void TVpPointDetector::UpdateCountersAtEndOfHistory(ULong_t numOfHistories)
{
  // Update counters at the end of a run.
  
  fNumOfHistories = numOfHistories;
  for (Int_t i = 0; i < fNumOfChannels+2; i++)
    {
      Double_t x = fChannelCounter[i];
      fChannelContent[i] += x;
      fChannelContent2[i] += x * x;
      fChannelLastHistory[i] = numOfHistories;
      fChannelCounter[i] = 0.0;
    }
}

//______________________________________________________________________________
Double_t TVpPointDetector::GetError(Int_t channel) const
{
  // Return the 3sigma error for a given channel
  
  Double_t error = 3 * sqrt((fChannelContent2[channel] - fChannelContent[channel]
			     * fChannelContent[channel]
			     / fNumOfHistories) / (fNumOfHistories - 1));
  return error;
}

//______________________________________________________________________________
void TVpPointDetector::SetNumOfChannels(Int_t numOfChannels,
					Double_t minEnergy, Double_t maxEnergy)
{
  // Allocate the channel array and initialize it to 0. Deallocate
  // previously allocated channel array and calculate energy-channel
  // transformation coefficients.

  fMinEnergy = minEnergy;
  fMaxEnergy = maxEnergy;
  fNumOfChannels = numOfChannels;
  fUnderflowChannel = fNumOfChannels;
  fOverflowChannel = fNumOfChannels + 1;
  
  delete [] fChannelContent;
  delete [] fChannelContent2;
  delete [] fChannelCounter;
  delete [] fChannelLastHistory;
  fChannelContent = new Double_t[fNumOfChannels + 2];    // + underflow and overflow channels
  fChannelContent2 = new Double_t[fNumOfChannels + 2];   // + underflow and overflow channels
  fChannelCounter = new Double_t[fNumOfChannels + 2];    // + underflow and overflow channels
  fChannelLastHistory = new ULong_t[fNumOfChannels + 2]; // + underflow and overflow channels

  // Initialize the array
  for (Int_t i = 0; i < fNumOfChannels + 2; i++)
    {
      fChannelContent[i] = fChannelContent2[i] = fChannelCounter[i] = 0.0;
      fChannelLastHistory[i] = 0;
    }

  // Calculate linear transformation coeficients
  fAc2e = (fMaxEnergy - fMinEnergy) / fNumOfChannels;
  fBc2e = fMinEnergy;
  fAe2c = 1.0 / fAc2e;
  fBe2c = fAe2c * fBc2e;
}

//______________________________________________________________________________
Double_t TVpPointDetector::GetChannelMeanEstimate(Int_t channel) const
{
  // Return the estimate of the mean value (the average).

  return fChannelContent[channel] / fNumOfHistories;
}

//______________________________________________________________________________
Double_t TVpPointDetector::GetChannelVarianceEstimate(Int_t channel) const
{
  // Return the estimate of the variance of the mean.

  Double_t average = GetChannelMeanEstimate(channel);
  return (fChannelContent2[channel] - fNumOfHistories * average * average)
    / (fNumOfHistories * (fNumOfHistories - 1.0));
}

//______________________________________________________________________________
Double_t TVpPointDetector::GetChannelSigmaEstimate(Int_t channel) const
{
  // Return the estimate of the sigma of the mean.

  return sqrt(GetChannelVarianceEstimate(channel));
}

//______________________________________________________________________________
Double_t TVpPointDetector::GetMeanEstimate() const
{
  // Get mean estimate.

  Double_t sum = 0.0;
  for (Int_t i = 0; i < fNumOfChannels+2; i++)
    sum += fChannelContent[i];
  return sum / fNumOfHistories;
}

//______________________________________________________________________________
Double_t TVpPointDetector::GetVarianceEstimate() const
{
  // Get variance estimate of the mean.

  Double_t average = GetMeanEstimate();
  Double_t sum2 = 0.0;
  for (Int_t i = 0; i < fNumOfChannels+2; i++)
    sum2 += fChannelContent2[i];
  return (sum2 - fNumOfHistories * average * average) /
    (fNumOfHistories * (fNumOfHistories - 1.0));
}

//______________________________________________________________________________
Double_t TVpPointDetector::GetSigmaEstimate() const
{
  // Get sigma estimate of the mean.

  return sqrt(GetVarianceEstimate());
}

//______________________________________________________________________________
Int_t TVpPointDetector::WritePDS(FILE *fp, Int_t withErrors)
{
  // Write the PDS file (Point Detector Spectrum).
  //
  // withErrors = 0 do not write errors (for backward compatibility)
  // withErrors = 1 write absolute errors
  
  Int_t count = 0;

  fprintf(fp, "\n");
  for (Int_t i = 0; i < fNumOfChannels+2; i++)
    {
      if (withErrors == 0)
	{
	  fprintf(fp, "%14e", GetChannelMeanEstimate(i));
	  if (++count % 5 == 0)
	    fprintf(fp, "\n");
	}
      else
	{
	  fprintf(fp, "%14e%11.3e", GetChannelMeanEstimate(i), GetChannelSigmaEstimate(i));
	  if (++count % 3 == 0)
	    fprintf(fp, "\n"); 
	}
    }
  fprintf(fp, "\n");
  
  return 0;
}

//______________________________________________________________________________
Int_t TVpPointDetector::ReadPDS(FILE *fp, Int_t withErrors, ULong_t numOfHistories)
{
  // Read the PDS file (Point Detector Spectrum)
  //
  // withErrors = 0 do not read errors (for backward compatibility)
  // withErrors = 1 read absolute errors of the average
  //
  // if withErrors is specified then numOfHistories must be specified too.
  // Otherwise fChannelContent2[] cannot be calculated.

  Double_t average, sigmaAvg;
  fNumOfHistories = numOfHistories;

  // Read average and sigmaAvg
  for (Int_t i = 0; i < fNumOfChannels+2; i++)
    {
      if (withErrors == 0)
	{
	  if (fscanf(fp, "%lf", &average) != 1)
	    goto error;
	  sigmaAvg = 0.0;
	}
      else
	{
	  if (fscanf(fp, "%lf%lf", &average, &sigmaAvg) != 2)
	    goto error;
	}
      fChannelContent[i] = fNumOfHistories * average;
      fChannelContent2[i] = sigmaAvg*sigmaAvg * fNumOfHistories * (fNumOfHistories - 1.0) +
	fNumOfHistories * average * average;
    }
  return 0;

 error:
  fprintf(stderr, "TVpPointDetector::ReadPDS: Unable to read ChannelContent.\n");
  return 1;
}

//______________________________________________________________________________
void TVpPointDetector::RegisterParticle(TVpParticle *particlePtr)
{
  // Register contribution from the particle.  This function is used by
  // TVpSetupTomograph::AnalyticProjection().  It does not transport the
  // particle through the geometry.  MC simulation uses
  // RegisterVirtualPhoton().

  if (fDedEventTreePtr != 0)
    { // Initialize the event data.  They will be written to a ROOT file by
      // the routine UpdateCounters() since event.c is still not known.
      TVpVector3 pos = particlePtr->fSolid->PosLocToUni(particlePtr->fPosition);
      TVpVector3 dir = particlePtr->fSolid->DirLocToUni(particlePtr->fDirection);
      dedEvent.i = particlePtr->fLastInteraction;
      dedEvent.x0 = pos.GetX();
      dedEvent.x1 = pos.GetY();
      dedEvent.x2 = pos.GetZ();
      dedEvent.u0 = dir.GetX();
      dedEvent.u1 = dir.GetY();
      dedEvent.u2 = dir.GetZ();
      dedEvent.e = particlePtr->fEnergy;
      dedEvent.w = particlePtr->fWeight;
      dedEvent.s = particlePtr->fSolid->GetIndex();
      dedEvent.t = particlePtr->fSolid->GetSubIndex(particlePtr->GetPos());
      dedEvent.h = particlePtr->fHistoryNumber;
    }

  RegisterScoredQuantity(particlePtr);
}

//______________________________________________________________________________
void TVpPointDetector::RegisterVirtualPhoton(TVpGeometry *geometry,
					     TVpPhoton *photonPtr)
{
  // Register contribution from the particle.  This function generates a
  // virtual photon heading towards the point detector, transports it through
  // the geometry and scores its contribution.
  //
  // If the point detector is associated with a DED file then a DED record is
  // written.

  if (fDedEventTreePtr != 0)
    { // Initialize the event data.  They will be written to a ROOT file by
      // the routine UpdateCounters() since event.c is still not known.
      TVpVector3 posU = photonPtr->fSolid->PosLocToUni(photonPtr->fPosition);
      TVpVector3 dirU = photonPtr->fSolid->DirLocToUni(photonPtr->fDirection);
      dedEvent.i = photonPtr->fLastInteraction;
      dedEvent.x0 = posU.GetX();
      dedEvent.x1 = posU.GetY();
      dedEvent.x2 = posU.GetZ();
      dedEvent.u0 = dirU.GetX();
      dedEvent.u1 = dirU.GetY();
      dedEvent.u2 = dirU.GetZ();
      dedEvent.e = photonPtr->fEnergy;
      dedEvent.w = photonPtr->fWeight;
      dedEvent.s = photonPtr->fSolid->GetIndex();
      dedEvent.t = photonPtr->fSolid->GetSubIndex(photonPtr->GetPos());
      dedEvent.h = photonPtr->fHistoryNumber;
    }

  // Generate a virtual particle
  TVpPhoton photonVir = *photonPtr;

  // Convert point detector position to local coordinates of the
  // solid where the particle is
  TVpVector3 pdPositionL = photonPtr->fSolid->PosUniToLoc(fTraVecL2u);
  TVpVector3 translationL = pdPositionL - photonPtr->fPosition;
  Double_t distance = norm(translationL);
  TVpVector3 dirToPdL = (1.0/distance) * translationL;
  photonVir.VirtualScatterToDirection(dirToPdL);
  Double_t opticalPath = geometry->GetOpticalPath(photonVir, distance);
  Double_t weight = photonVir.GetWeight() *
    exp(-opticalPath) / (distance * distance);
  photonVir.SetWeight(weight);
  
  // Register the virual particle
  RegisterScoredQuantity(&photonVir);
}

//______________________________________________________________________________
void TVpPointDetector::MultiplyChannelContentBy(Double_t factor)
{
  // Multiply the content of all channels including the owerflow and underflow
  // ones by a number.  Also multiply the second moment channels.

  Double_t factor2 = factor * factor;
  for (Int_t i = 0; i < fNumOfChannels + 2; i++)
    {
      fChannelContent[i] *= factor;
      fChannelContent2[i] *= factor2;
      fChannelCounter[i] *= factor;
    }
}

//______________________________________________________________________________
void TVpPointDetector::Zero()
{
  // Zero the content of all channels including the owerflow and underflow
  // ones.

  for (Int_t i = 0; i < fNumOfChannels + 2; i++)
    {
      fChannelContent[i] = fChannelContent2[i] = fChannelCounter[i] = 0.0;
      fChannelLastHistory[i] = 0;
    }
}

//______________________________________________________________________________
void TVpPointDetector::PrintStatus(std::ostream &out) const
{
  // Print the object status
  //
  // Input parameters:
  // - out - output stream
  //
  // Example:
  // root [] TVpPointDetector *pdPtr = new TVpPointDetector();
  // root [] pdPtr->PrintStatus();
  // <TVpPointDetector>
  // Position:    0.00000000e+00   0.00000000e+00   0.00000000e+00
  // Rotation matrix:
  //    1.00000000e+00   0.00000000e+00   0.00000000e+00
  //    0.00000000e+00   1.00000000e+00   0.00000000e+00
  //    0.00000000e+00   0.00000000e+00   1.00000000e+00
  // </TVpPointDetector>

  out << "<TVpPointDetector>\n"
      << "$Id: TVpPointDetector.cc 62 2009-06-27 10:54:08Z malusek $\n";
  PrintLocation(out);
  out << "</TVpPointDetector>\n";
}


#include "TPolyMarker3D.h"
  
//______________________________________________________________________________
void TVpPointDetector::Draw(Int_t option) const
{
  // Draw a dot
  
  TPolyMarker3D *pos;
  if (option & kDrawAxes)
    DrawLcsAxes(1, 1, 1);
  if (option & kDrawPoint)
    {
      pos= new TPolyMarker3D (1, 8,"");
      pos->SetMarkerColor(1);
      pos->SetMarkerStyle(8);
      pos->SetMarkerSize(0.1);
      pos->SetPoint(0, fTraVecL2u.GetX(), fTraVecL2u.GetY(), fTraVecL2u.GetZ());
      pos->Draw();
    }
}

//______________________________________________________________________________
void TVpPointDetector::DedEventInitialize(Char_t *fileName)
{
  // Create a ROOT tree.  The branch address is initialized to an external
  // object.

  fDedEventFilePtr = new TFile(fileName, "RECREATE");
  fDedEventTreePtr = new TTree("Tded", "DedEvent");
  fDedEventTreePtr->Branch("DedEvent", &dedEvent, "i/I:x0/F:x1:x2:u0:u1:u2:e:w:c:a:s/I:t:h");
  fDedEventTreePtr->Print();
}

//______________________________________________________________________________
void TVpPointDetector::DedEventClose()
{
  fDedEventFilePtr->Write();
  fDedEventTreePtr->Print();
  fDedEventFilePtr->Close();
}

//______________________________________________________________________________
TH1F *TVpPointDetector::GetHistEnergySpectrum() const
{
  // Return histogram containing energy spectrum of the scored quantity.
  //
  // errors will be added later

  TH1F *h = new TH1F("h", "energy spectrum", fNumOfChannels, fMinEnergy, fMaxEnergy);
  for (Int_t i = 0; i < fNumOfChannels; i++)
    {
      h->SetBinContent(i+1, fChannelContent[i]);
    }
  return h;
}
