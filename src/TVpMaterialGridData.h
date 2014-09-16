#ifndef TVpMaterialGridData_h
#define TVpMaterialGridData_h

#include "TObject.h"
#include "TGraph.h"
#include <cmath>
#include "TVpSpectrum.h"

class TVpMaterialGridData
{
 protected:
  Int_t       fDimIub;    // Dimension of IUB arrays (the number of intervals)
  
  // Interval Upper Boundary arrays, coefficients of 1st degree polynomials 
  Double_t   *fInIub0;            //[fDimIub] incoherent scattering
  Double_t   *fInIub1;            //[fDimIub] incoherent scattering
  Double_t   *fPhIub0;            //[fDimIub] photoeffect
  Double_t   *fPhIub1;            //[fDimIub] photoeffect
  Double_t   *fCoIub0;            //[fDimIub] coherent scattering
  Double_t   *fCoIub1;            //[fDimIub] coherent scattering

  // Equidistant energy grid
  Double_t   *fInCsg0;            //[fDimIub] incoherent scattering
  Double_t   *fInCsg1;            //[fDimIub] incoherent scattering
  Double_t   *fPhCsg0;            //[fDimIub] photoeffect
  Double_t   *fPhCsg1;            //[fDimIub] photoeffect
  Double_t   *fCoCsg0;            //[fDimIub] coherent scattering
  Double_t   *fCoCsg1;            //[fDimIub] coherent scattering
  Double_t   *fToCsg0;            //[fDimIub] total cross section = In+Ph+Co
  Double_t   *fToCsg1;            //[fDimIub] total cross section = In+Ph+Co

  // Normalization factors for coherent and incoherent scattering (used by PdfTheta)
  Double_t   *fInNtg0;            //[fDimIub] InCs normalization based on Sfg
  Double_t   *fInNtg1;            //[fDimIub] InCs normalization based on Sfg
  Double_t   *fCoNtg0;            //[fDimIub] CoCs normalization based on Ffg
  Double_t   *fCoNtg1;            //[fDimIub] CoCs normalization based on Ffg

  // Normalization factors for coherent and incoherent scattering (used by PdfOmega)
  Double_t   *fInNog0;            //[fDimIub] InCs normalization based on Sfg
  Double_t   *fInNog1;            //[fDimIub] InCs normalization based on Sfg
  Double_t   *fCoNog0;            //[fDimIub] CoCs normalization based on Ffg
  Double_t   *fCoNog1;            //[fDimIub] CoCs normalization based on Ffg

 public:
  enum EGridGraphType {kCoCsg, kInCsg, kPhCsg, kToCsg, kCoIub, kInIub, kPhIub,
		       kSfg, kFfg};

  Double_t    fEnergyMinIUB;     // Min IUB energy
  Double_t    fEnergyMaxIUB;     // Max IUB energy
  Double_t   *fEnergyIUB;        //[fDimIub] IUB energy grid

  Double_t    fEnergyToChannel;  // Energy to channel conversion factor
  Double_t    fChannelToEnergy;  // Channel to energy conversion factor

  // Sfg
  Int_t       fDimLnSfg;     // dimension (the number of intervals)
  Double_t    fLnSfgXMin;    // ln of the first non-zero value of x
  Double_t    fLnSfgXMax;    // ln of the last non-zero value of x
  Double_t    fDLnSfgX;      // = (fLnSfgXMax - fLnSfgXMin)/(fDimLnSfg - 1)
  Double_t   *fLnSfg1;       //[fDimLnSfg]
  Double_t   *fLnSfg0;       //[fDimLnSfg]

  // Ffg
  Int_t       fDimLnFfg;     // dimension (the number of intervals)
  Double_t    fLnFfgXMin;    // ln of the first non-zero value of x
  Double_t    fLnFfgXMax;    // ln of the last non-zero value of x
  Double_t    fDLnFfgX;      // = (fLnFfgXMax-fLnFfgXMin)/(fDimLnFfg-1)
  Double_t   *fLnFfg1;       //[fDimLnFfg]
  Double_t   *fLnFfg0;       //[fDimLnFfg]

  TVpMaterialGridData();
  virtual ~TVpMaterialGridData();

  // Iub
  inline Double_t    GetInIub(Double_t energy);
  inline Double_t    GetPhIub(Double_t energy);
  inline Double_t    GetCoIub(Double_t energy);

  // Csg
  inline Double_t    GetInCsg(Double_t energy);
  inline Double_t    GetPhCsg(Double_t energy);
  inline Double_t    GetCoCsg(Double_t energy);
  inline Double_t    GetToCsg(Double_t energy);
  Double_t           GetToCsg(TVpSpectrum *spectrum);

  // Ntg and Nog
  inline Double_t    GetInNtg(Double_t energy);
  inline Double_t    GetCoNtg(Double_t energy);
  inline Double_t    GetInNog(Double_t energy);
  inline Double_t    GetCoNog(Double_t energy);
  
  inline Double_t    GetLac(Double_t energy);
  Double_t           GetLac(TVpSpectrum *spectrum, Int_t weighting = 0);
  inline Int_t       GetChannel(Double_t energy);
  inline Int_t       GetChannelLnXG(Double_t x);
  inline Double_t    GetEnergy(Int_t channel);
  Double_t           GetSfg(Double_t x) const;
  Double_t           GetFfg(Double_t x) const;
  Double_t           GetFfgValue(Int_t index) const;
  Double_t           GetFfgX(Int_t index) const;

  Double_t           GetGridDataValue(EGridGraphType graphType, Double_t energy);
  
  TH1F *GetHistGridData(EGridGraphType graphType);
  TGraph *GetGraphGridData(EGridGraphType graphType) const;
  void DrawGraphGridData(EGridGraphType graphType, Char_t *opt = "AP") const;

  ClassDef(TVpMaterialGridData,1) // Material grid data, implementation class
};

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetToCsg(Double_t energy)
{
  // Return total macroscopic cross section in 1/cm.  This is the same
  // quantity as the linear attenuation coefficient.
  //
  // Input parameters:
  // - energy - photon energy in keV
  //
  // Method:
  // Linear interpolation of grid data.

  Int_t channel = GetChannel(energy);
  return fToCsg0[channel] + fToCsg1[channel] * energy;
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetLac(Double_t energy)
{
  // Return the linear attenuation coefficient in 1/cm.  This function is an
  // alias for GetToCsg().
  //
  // Input parameters:
  // - energy - photon energy in keV

  return GetToCsg(energy);
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetInIub(Double_t energy)
{
  // Return IUB value for incoherent scattering.
  //
  // Input parameters:
  // - energy - photon energy in keV

  Int_t channel = GetChannel(energy);
  return fInIub0[channel] + fInIub1[channel] * energy;
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetPhIub(Double_t energy)
{
  // Return IUB value for photoelectric effect.
  //
  // Input parameters:
  // - energy - photon energy in keV
  
  Int_t channel = GetChannel(energy);
  return fPhIub0[channel] + fPhIub1[channel] * energy;
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetCoIub(Double_t energy)
{
  // Return IUB value for coherent scattering.
  //
  // Input parameters:
  // - energy - photon energy in keV

  Int_t channel = GetChannel(energy);
  return fCoIub0[channel] + fCoIub1[channel] * energy;
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetInCsg(Double_t energy)
{
  // Return incoherent scattering cross section in cm^2/g.
  //
  // Input parameters:
  // - energy - photon energy in keV

  Int_t channel = GetChannel(energy);
  return fInCsg0[channel] + fInCsg1[channel] * energy;
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetPhCsg(Double_t energy)
{
  // Return photoeffect cross section in cm^2/g.
  //
  // Input parameters:
  // - energy - photon energy in keV

  Int_t channel = GetChannel(energy);
  return fPhCsg0[channel] + fPhCsg1[channel] * energy;
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetCoCsg(Double_t energy)
{
  // Return coherent scattering cross section in cm^2/g.
  //
  // Input parameters:
  // - energy - photon energy in keV

  Int_t channel = GetChannel(energy);
  return fCoCsg0[channel] + fCoCsg1[channel] * energy;
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetInNtg(Double_t energy)
{
  // Return incoherent scattering normalization factor

  Int_t channel = GetChannel(energy);
  return fInNtg0[channel] + fInNtg1[channel] * energy;
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetCoNtg(Double_t energy)
{
  // Return coherent scattering normalization factor

  Int_t channel = GetChannel(energy);
  return fCoNtg0[channel] + fCoNtg1[channel] * energy;
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetInNog(Double_t energy)
{
  // Return incoherent scattering normalization factor

  Int_t channel = GetChannel(energy);
  return fInNog0[channel] + fInNog1[channel] * energy;
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetCoNog(Double_t energy)
{
  // Return coherent scattering normalization factor

  Int_t channel = GetChannel(energy);
  return fCoNog0[channel] + fCoNog1[channel] * energy;
}

//______________________________________________________________________________
inline Int_t TVpMaterialGridData::GetChannel(Double_t energy)
{
  // Return Csg or Iub array channel corresponding to energy

  return (int) (fEnergyToChannel * (energy - fEnergyMinIUB));
}

//______________________________________________________________________________
inline Int_t TVpMaterialGridData::GetChannelLnXG(Double_t x)
{
  Double_t lnX = log(x);
  return (int) ((lnX - fLnSfgXMin) / fDLnSfgX);
}

//______________________________________________________________________________
inline Double_t TVpMaterialGridData::GetEnergy(Int_t channel)
{
  return channel * fChannelToEnergy + fEnergyMinIUB;
}

#endif  // TVpMaterialGridData_h
