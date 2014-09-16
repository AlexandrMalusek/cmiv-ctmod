//______________________________________________________________________________
//
// TVpSinogram stores a sinogram.

#include <cmath>
#include <iomanip>
#include "TVpSinogram.h"
#include "TVpRunManagerTomograph.h"

ClassImp(TVpSinogram)

//______________________________________________________________________________
TVpSinogram::TVpSinogram(Int_t numOfProjections, TVpSetupTomograph *setupTomographPtr)
{
  // Constructor

  fNumOfProjections = numOfProjections;
  fSetupTomographPtr = setupTomographPtr;
  fSinogramData = new Float_t[GetNumOfDetectors() * numOfProjections];
}

//______________________________________________________________________________
TVpSinogram::~TVpSinogram()
{
  // Destructor
  delete [] fSinogramData;
}

//______________________________________________________________________________
Int_t TVpSinogram::WriteFile(Char_t *fileName)
{
  // Write sinogram to the file fileName

  // empty now
  return 0;
}

//______________________________________________________________________________
void TVpSinogram::CalculatePrimary()
{
  // Calculate the primary sinogram

  const Int_t npro = GetNumOfProjections();
  const Int_t ndet = GetNumOfDetectors();
  Double_t dAngle = 2*M_PI / fNumOfProjections;

  for (Int_t projection = 0; projection < npro; projection++)
    {
      fSetupTomographPtr->SetPosition(projection * dAngle, 0);
      fSetupTomographPtr->GetPointDetectorArrayPtr()->Zero();
      fSetupTomographPtr->AnalyticProjection();
      for (Int_t idet = 0; idet < ndet; idet++)
	fSinogramData[projection*ndet + idet] = 
	  fSetupTomographPtr->GetPointDetectorArrayPtr()->GetMeanEstimate(idet);
    }
}

//______________________________________________________________________________
void TVpSinogram::CalculateScatter(Long_t numOfHistories)
{
  // Calculate the scatter sinogram.
  //
  // Input parameters:
  // - numOfHistories - number of histories per projection (default=10000)

  const Int_t npro = GetNumOfProjections();
  const Int_t ndet = GetNumOfDetectors();
  Double_t dAngle = 2*M_PI / fNumOfProjections;

  fSetupTomographPtr->ActivateSource(1);
  for (Int_t projection = 0; projection < npro; projection++)
    {
      std::cerr << "Info: projection: " << projection << " of " << fNumOfProjections
		<< ", " << std::setprecision(3)
		<< projection/(Double_t)fNumOfProjections*100.0 << "%.\n";
      fSetupTomographPtr->SetPosition(projection * dAngle, 0);
      fSetupTomographPtr->GetPointDetectorArrayPtr()->Zero();
      TVpRunManagerTomograph *trmPtr = new TVpRunManagerTomograph(fSetupTomographPtr);
      trmPtr->Run(numOfHistories);
      for (Int_t idet = 0; idet < ndet; idet++)
	fSinogramData[projection*ndet + idet] = 
	  fSetupTomographPtr->GetPointDetectorArrayPtr()->GetMeanEstimate(idet);
      delete trmPtr;
    }
}

//______________________________________________________________________________
TH2F *TVpSinogram::GetProbImage(Char_t *hname)
{
  // Return 2D image
  
  Int_t nx = GetNumOfDetectors();
  Int_t ny = GetNumOfProjections();

  TH2F *h = new TH2F("sinogram", hname, nx, 0, nx, ny, 0, ny);
  for (Int_t idet = 0; idet < nx; idet++)
    for (Int_t ipro = 0; ipro < ny; ipro++)
      h->SetCellContent(idet+1, ipro+1, fSinogramData[ipro * nx + idet]);
  return h;
}

//______________________________________________________________________________
TH2F *TVpSinogram::GetWopImage(Char_t *hname)
{
  // Return 2D image
  
  Int_t nx = GetNumOfDetectors();
  Int_t ny = GetNumOfProjections();

  TH2F *h = new TH2F("sinogram", hname, nx, 0, nx, ny, 0, ny);
  for (Int_t idet = 0; idet < nx; idet++)
    for (Int_t ipro = 0; ipro < ny; ipro++)
      h->SetCellContent(idet+1, ipro+1, -log(fSinogramData[ipro * nx + idet]));
  return h;
}
