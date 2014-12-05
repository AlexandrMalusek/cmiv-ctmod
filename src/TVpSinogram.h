#ifndef TVpSinogram_h
#define TVpSinogram_h

#include "TObject.h"
#include "TH2.h"
#include "TVpSetupTomograph.h"

class TVpSinogram
{
 public:
  Int_t              fNumOfProjections;    // number of projection
  Float_t           *fSinogramData;        //! array of values
  TVpSetupTomograph *fSetupTomographPtr;   //! a pointer to an external tomograph setup
  
  TVpSinogram(Int_t numOfProjections, TVpSetupTomograph *setupTomographPtr);
  virtual ~TVpSinogram();

  Int_t    WriteFile(const Char_t *fileName);
  void     CalculatePrimary();
  void     CalculateScatter(Long_t numOfHistories = 10000);
  Int_t    GetNumOfProjections();
  Int_t    GetNumOfDetectors();

  TH2F *GetProbImage(const Char_t *hname = "Prob Sinogram");
  TH2F *GetWopImage(const Char_t *hname = "Wop Sinogram");

  ClassDef(TVpSinogram,1) // Sinogram calculation
};

//______________________________________________________________________________
inline Int_t TVpSinogram::GetNumOfProjections()
{
  // Return the number of projections

  return fNumOfProjections;
}

//______________________________________________________________________________
inline Int_t TVpSinogram::GetNumOfDetectors()
{
  // Return the number of detectors

  return fSetupTomographPtr->GetPointDetectorArrayPtr()->GetNumOfDetectors();
}

#endif  // TVpSinogram_h
