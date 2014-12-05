#ifndef TVpDetectorResponse_h
#define TVpDetectorResponse_h

#include "TObject.h"
#include "TH2.h"
#include "TVpAntiscatterGrid.h"

class TVpDetectorResponse : public TObject
{
 public:
  std::string fName;           // Descriptive name
  Int_t       fDimXi;          // Number of xi = cos(theta)
  Int_t       fDimEnergy;      // Number of energies
  Double_t    fMaxEnergy;      // Max energy
  Double_t   *fAer;            // [fDimEnergy*fDimXi] The angle-energy response
  Double_t   *fRer;            // [fDimEnergy*fDimXi] 3sigma abs. err. of aer
  Double_t    fStepXi;         // fStepXi = 1.0 / fDimXi;
  Double_t    fStepEnergy;     // fStepEnergy = fMaxEnergy / fDimEnergy;
  TVpAntiscatterGrid *fAntiscatterGridPtr; //! A pointer to the antiscatter grid

  TVpDetectorResponse();
  virtual ~TVpDetectorResponse();
  inline Int_t     GetIndex(Int_t iXi, Int_t iEnergy) const;
  Double_t         GetResponse(Double_t xi, Double_t energy) const;
  Double_t         GetResponse(Double_t energy, const TVpVector3 &directionL) const;
  inline Double_t  GetGridXi(Int_t iXi) const;
  inline Double_t  GetGridEnergy(Int_t iEnergy) const;
  inline Int_t     GetDimXi() const;
  inline Int_t     GetDimEnergy() const;
  Int_t            ReadAerFile(const Char_t *fileName);
  void             SetAntiscatterGrid(TVpAntiscatterGrid *antiscatterGridPtr);
  virtual void     PrintStatus(std::ostream &out = std::cout) const;
  
  TH2F *GetAer();
  TH2F *GetAerAbsError();
  TH2F *GetAerRelError();
  TH2F *GetAer(Int_t dimXi, Double_t minXi, Double_t maxXi,
	       Int_t dimEnergy, Double_t minEnergy, Double_t maxEnergy);
  
  ClassDef(TVpDetectorResponse,1) // Detector response function
};

//______________________________________________________________________________
inline Int_t TVpDetectorResponse::GetIndex(Int_t iXi, Int_t iEnergy) const
{
  // Return the fAer array index of the grid index pair (iXi, iEnergy).
  //
  // Input parameters:
  // - iXi - Xi index, range: 0 ... fDimXi-1
  // - iEnergy - energy index, range: 0 ... fDimEnergy-1
  
  return iXi * fDimEnergy + iEnergy;
}

//______________________________________________________________________________
inline Double_t TVpDetectorResponse::GetGridXi(Int_t iXi) const
{
  // Return grid value of xi corresponding to an index iXi.
  //
  // Input parameters:
  // - iXi - Xi index, range: 0 ... fDimXi-1

  return fStepXi * (iXi + 1);
}

//______________________________________________________________________________
inline Double_t TVpDetectorResponse::GetGridEnergy(Int_t iEnergy) const
{
  // Return grid value of energy corresponding to an index iEnergy.
  //
  // Input parameters:
  // - iEnergy - energy index, range: 0 ... fDimEnergy-1

  return fStepEnergy * (iEnergy + 1);
}

//______________________________________________________________________________
inline Int_t TVpDetectorResponse::GetDimXi() const
{
  // Return dimension of the xi grid

  return fDimXi;
}

//______________________________________________________________________________
inline Int_t TVpDetectorResponse::GetDimEnergy() const
{
  // Return dimension of the energy grid

  return fDimEnergy;
}

#endif  // TVpDetectorResponse_h
