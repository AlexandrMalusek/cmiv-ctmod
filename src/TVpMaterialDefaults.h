#ifndef TVpMaterialDefaults_h
#define TVpMaterialDefaults_h

#include <iostream>
#include "TObject.h"

class TVpMaterialDefaults
{
 public:
  Int_t    fDimIub;          // Dimension of IUB arrays (the number of intervals)
  Int_t    fDimA;            // Dimension of the A function (see the theory)
  Int_t    fDimLnSfg;        // dimension (the number of intervals)
  Int_t    fDimLnFfg;        // dimension (the number of intervals)
  Double_t fEnergyMinIUB;    // Min IUB energy [keV]
  Double_t fEnergyMaxIUB;    // Max IUB energy [keV]

  TVpMaterialDefaults();
  virtual ~TVpMaterialDefaults() {};
  inline Int_t    GetDimIub() const;
  inline Int_t    GetDimA() const;
  inline Int_t    GetDimLnSfg() const;
  inline Int_t    GetDimLnFfg() const;
  inline Double_t GetEnergyMinIUB() const;
  inline Double_t GetEnergyMaxIUB() const;
  void            SetDimIub(Int_t dimIub);
  void            SetDimA(Int_t dimA);
  void            SetDimLnSfg(Int_t dimLnSfg);
  void            SetDimLnFfg(Int_t dimLnFfg);
  void            SetEnergyMinIUB(Double_t energyMinIUB);
  void            SetEnergyMaxIUB(Double_t energyMaxIUB);
  void            PrintStatus(std::ostream &out = std::cout) const;

  ClassDef(TVpMaterialDefaults,0) // Physical constants
};

//______________________________________________________________________________
inline Int_t TVpMaterialDefaults::GetDimIub() const
{
  // Return the dimension of the IUB array.

  return fDimIub;
}

//______________________________________________________________________________
inline Int_t TVpMaterialDefaults::GetDimA() const
{
  // Return the dimension of the A array.

  return fDimA;
}

//______________________________________________________________________________
inline Int_t TVpMaterialDefaults::GetDimLnSfg() const
{
  // Return the dimension of the LnSfg arrays

  return fDimLnSfg;
}

//______________________________________________________________________________
inline Int_t TVpMaterialDefaults::GetDimLnFfg() const
{
  // Return the dimension of the LnFfg arrays

  return fDimLnFfg;
}

//______________________________________________________________________________
inline Double_t TVpMaterialDefaults::GetEnergyMinIUB() const
{
  // Return the low-energy cutoff of IUB arrays
 
  return fEnergyMinIUB;
}

//______________________________________________________________________________
inline Double_t TVpMaterialDefaults::GetEnergyMaxIUB() const
{
  // Return the high-energy cutoff of IUB arrays

  return fEnergyMaxIUB;
}

#endif  // TVpMaterialDefaults_h
