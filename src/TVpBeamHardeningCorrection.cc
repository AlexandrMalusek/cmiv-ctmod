//______________________________________________________________________________
//
// TVpBeamHardeningCorrection
//
// 
// -1        0         1                  n         -2
// -----|---------|-------.........|------------|-------------
//     x[0]      x[1]            x[n-2]       x[n-1]
//
// x[n]

#include <cmath>
#include <iostream>
#include "TVpBeamHardeningCorrection.h"
#include "TVpMath.h"

ClassImp(TVpBeamHardeningCorrection)

//______________________________________________________________________________
TVpBeamHardeningCorrection::TVpBeamHardeningCorrection(TVpMaterial *materialPtr,
	              TVpSpectrum *spectrumPtr, TVpDetectorResponse *detResponsePtr,
						       Int_t dim, Double_t maxThickness)
{
  // Constructor.
  //
  // Input parameters:
  // - materialPtr - material for which the correction is performed (mostly water)
  // - spectrumPtr - spectrum of photons
  // - detResponsePtr - detector response 
  // - dim - dimension of the grid of points (array)
  // - maxThickness - max thickness in cm for which the grid is calculated 
  
  fMaterialPtr = materialPtr;
  fSpectrumPtr = spectrumPtr;
  fDetResponsePtr = detResponsePtr;
  fDim = dim;
  fX = new Double_t[fDim];
  fP = new Double_t[fDim];
  fDx = maxThickness / (dim - 1);
  fS = CalculateS();

  // Fill in the grid
  for (Int_t i = 0; i < dim; i++)
    {
      fX[i] = i * fDx;
      fP[i] = CalculateP(fX[i]);
    }
}

//______________________________________________________________________________
TVpBeamHardeningCorrection::~TVpBeamHardeningCorrection()
{
  // Destructor.
  
  delete[] fX;
  delete[] fP;
}

//______________________________________________________________________________
void TVpBeamHardeningCorrection::SetEffectiveLac(Double_t effLac)
{
  // Set the effective linear attenuation coefficient.
  
  fEffLac = effLac;
}

//______________________________________________________________________________
Double_t TVpBeamHardeningCorrection::GetGridP(Double_t x) const
{
  // Return P for an x.  Use lin-lin interpolation
  //
  // Input parameters:
  // - x - material thickness in cm

  Int_t il = Int_t(x / fDx);
  if (il < 0 || il > fDim-2)
    {
      std::cerr << "Error: TVpBeamHardeningCorrection::GetGridP: x = "
		<< x << " out of range.  Returning 0.0.\n";
      return 0.0;
    }
  Int_t ih = il+1;
  return TVpMath::GetLinLinInterpolation(fX[il], fP[il], fX[ih], fP[ih], x);
}

//______________________________________________________________________________
Double_t TVpBeamHardeningCorrection::GetGridX(Double_t P) const
{
  // Return P for an x.  Use lin-lin interpolation
  //
  // Input parameters:
  // - P - radiological path

  Int_t il = TVpMath::FindIndexByBinarySearch(fDim, fP, P);
  if (il < 0)
    {
      std::cerr << "Error: TVpBeamHardeningCorrection::GetGridX: P = "
		<< P << " out of range.  Returning 0.0.\n";
      return 0.0;
    }
  Int_t ih = il+1;
  return TVpMath::GetLinLinInterpolation(fP[il], fX[il], fP[ih], fX[ih], P);
}

//______________________________________________________________________________
Double_t TVpBeamHardeningCorrection::CalculateP(Double_t x) const
{
  // Calculate the effective radiological path, P, using numerical
  // integration.
  
  Double_t weight, energy, response, rp;
  Double_t sum = 0.0;
  Int_t nc = fSpectrumPtr->GetNumOfChannels();
  
  for (Int_t i = 0; i < nc; i++)
    {
      fSpectrumPtr->GetEnergyAndWeight(i, energy, weight);
      response = fDetResponsePtr->GetResponse(1.0, energy);
      rp = x * fMaterialPtr->GetLac(energy);
      sum += energy * weight * response * exp(-rp);
    }
  sum /= 4*M_PI;
  return -log (sum / GetS());
} 

//______________________________________________________________________________
Double_t TVpBeamHardeningCorrection::CalculateS() const
{
  // Calculate the effective source intensity, S, using numerical integration.
  
  Double_t weight, energy, response;
  Double_t sum = 0.0;
  Int_t nc = fSpectrumPtr->GetNumOfChannels();
  
  for (Int_t i = 0; i < nc; i++)
    {
      fSpectrumPtr->GetEnergyAndWeight(i, energy, weight);
      response = fDetResponsePtr->GetResponse(1.0, energy);
      sum += energy * weight * response;
    }
  return sum / (4*M_PI);
}

//______________________________________________________________________________
Double_t TVpBeamHardeningCorrection::GetI0(Double_t distance) const
{
  // Calculate intensity when only the source and detector is considered.
  // (No bowtie filter and no phantom).
  
  return GetS() / (distance * distance);
}

//______________________________________________________________________________
Double_t TVpBeamHardeningCorrection::Correct(Double_t I1overI0) const
{
  // Correct

  if (I1overI0 > 1.0)  // may happen because of scatter
    I1overI0 = 1.0;    // correct it
  Double_t P = -log(I1overI0);
  Double_t x = GetGridX(P);
  Double_t correction = exp(P - fEffLac * x);
  return I1overI0 * correction;
}

//______________________________________________________________________________
Double_t TVpBeamHardeningCorrection::Correct(Double_t I2overI1, Double_t I1overI0) const
{
  // Correct

  Double_t I2overI0 = I2overI1 * I1overI0;
  return Correct(I2overI0) / Correct(I1overI0);
}

//______________________________________________________________________________
Double_t TVpBeamHardeningCorrection::Correct(Double_t I2, Double_t I1, Double_t I0) const
{
  // Correct

  return Correct(I2/I0) / Correct(I1/I0);
}

//______________________________________________________________________________
TGraph *TVpBeamHardeningCorrection::GetGraphGridP() const
{
  // Return the graph of effective P as a function of x

  return new TGraph(fDim, fX, fP);
}

//______________________________________________________________________________
TGraph *TVpBeamHardeningCorrection::GetGraphMonoP() const
{
  // Return the graph of P for monoenergetic photons as a function of x

  Double_t *y = new Double_t[fDim];
  for (Int_t i = 0; i < fDim; i++)
    y[i] = fEffLac * i * fDx;
  TGraph *g = new TGraph(fDim, fX, y);
  delete[] y;
  return g;
}
