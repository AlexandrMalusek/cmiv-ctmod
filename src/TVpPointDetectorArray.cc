//______________________________________________________________________________
//
// A TVpPointDetectorArray is a fNx x fNy array of point
// detectors.
//
// All point detectors share the same number of channels
// fNumOfChannels and cutoffs fMinEnergy and fMaxEnergy. See
// TVpPointDetector for more information.
//
// Position of each point detector is stored by the point detector
// itself. Derived classes like TVpPointDetectorArrayPlanar and
// TVpPointDetectorArrayCylindrical are provided to simplify
// initialization of positions of individual point detectors.
//
// A DED file example:
// TVpSetupTomograph *setupTomographPtr = getTomograph();
// setupTomographPtr->GetPointDetectorArrayPtr()->Zero();
// setupTomographPtr->GetPointDetectorArrayPtr()->DedEventInitialize(0, 32, "event01.root");
// ... // Run the simulation
// setupTomographPtr->GetPointDetectorArrayPtr()->DedEventClose(0, 32);

#include <cmath>
#include <string>
#include "TVpPointDetectorArray.h"

ClassImp(TVpPointDetectorArray)

//______________________________________________________________________________
TVpPointDetectorArray::TVpPointDetectorArray()
{
  // Default constructor. Set all members to 0.

  fNx = fNy = 0;
  fNumOfChannels = fNumOfDets = 0;
  fPointDetector = 0;
  fMinEnergy = fMaxEnergy = 0.0;
}

//______________________________________________________________________________
TVpPointDetectorArray::TVpPointDetectorArray
(Int_t nx, Int_t ny, Int_t numOfChannels,
 Double_t minEnergy, Double_t maxEnergy, TVpDetectorResponse *detectorResponse)
{
  // Constructor.
  //
  // Input parameters:
  // - nx -
  // - ny - 
  // - numOfChannels -
  // - minEnergy -
  // - maxEnergy -
  // - detectorResponse -

  // Checks
  if (nx <= 0 || ny <= 0)
    std::cerr << "Error: TVpPointDetectorArray::TVpPointDetectorArray: "
	      << "nx = " << nx << " and ny = " << ny << "must be at least 1.\n";

  fNx = nx;
  fNy = ny;
  fNumOfChannels = numOfChannels;
  fMinEnergy = minEnergy;
  fMaxEnergy = maxEnergy;
  
  fNumOfDets = fNx * fNy;
  
  // Allocate the point detector array
  fPointDetector = new TVpPointDetector[fNumOfDets];
  if (fPointDetector == 0)
    {
      std::cerr << "Error: TVpPointDetectorArray::TVpPointDetectorArray: "
		<< "Allocation failed.\n";
      return;
    }
  
  TVpPointDetector *pdPtr;
  for (Int_t iy = 0; iy < fNy; ++iy)
    for (Int_t ix = 0; ix < fNx; ++ix)
      {
	pdPtr = &fPointDetector[GetIndex(iy, ix)];
	pdPtr->SetNumOfChannels(numOfChannels, minEnergy, maxEnergy);
	pdPtr->SetDetectorResponse(detectorResponse);
      }
}

//______________________________________________________________________________
TVpPointDetectorArray::TVpPointDetectorArray(Char_t *fileName)
{
  // Constructor which reads all data from a file name.

  ReadPdaFile(fileName);
}

//______________________________________________________________________________
TVpPointDetectorArray::~TVpPointDetectorArray()
{
  // Destructor

  delete [] fPointDetector;
}

//______________________________________________________________________________
void TVpPointDetectorArray::PrintStatus(std::ostream &out) const
{
  // Print the object status

  out << "<TVpPointDetectorArray>\n"
      << "$Id: TVpPointDetectorArray.cc 62 2009-06-27 10:54:08Z malusek $\n" 
      << "Dimension:\t" << fNx << '\t' << fNy << '\n'
      << "Number of detectors: " << GetNumOfDetectors() << '\n'
      << "Histogram's min energy: " << fMinEnergy << '\n'
      << "Histogram's max energy: " << fMaxEnergy << '\n';
  for (Int_t ix = 0; ix < GetNx(); ++ix)
    for (Int_t iy = 0; iy < GetNy(); ++iy)
      {
	out << "ix = " << ix+1 << ", iy = " << iy+1 << ":\n"; 
	fPointDetector[GetIndex(iy,ix)].PrintStatus(out);
      }
  if (GetNumOfDetectors() > 0 && fPointDetector[0].fDetectorResponse != 0)
    fPointDetector[0].fDetectorResponse->PrintStatus(out);
  out << "</TVpPointDetectorArray>\n";
}


//______________________________________________________________________________
void TVpPointDetectorArray::SetActiveTranslation(TVpVector3 *traVec)
{
  // Set active translation
  
  std::cerr << "Error: TVpPointDetectorArrayPlanar::SetActiveTranslation: "
	    << "Use TVpPointDetectorArrayPlanar::SetActiveTransformation instead.\n";
}

//______________________________________________________________________________
void TVpPointDetectorArray::SetActiveRotation(TVpMatrix3x3 *rotMat)
{
  // Set active rotation
  
  std::cerr << "Error: TVpPointDetectorArrayPlanar::SetActiveRotation: "
	    << "Use TVpPointDetectorArrayPlanar::SetActiveTransformation instead.\n";
}

//______________________________________________________________________________
void TVpPointDetectorArray::SetActiveTransformation(TVpMatrix3x3 *rotMat,
						    TVpVector3 *traVec)
{
  // Set active transformation

  SetDefaultPosition();
  TVpObjectLocation::SetActiveTransformation(rotMat, traVec);

  // Rotation and translation
  TVpPointDetector *pdPtr;
  for (Int_t ix = 0; ix < fNx; ++ix)
    for (Int_t iy = 0; iy < fNy; ++iy)
      {
	pdPtr = &fPointDetector[GetIndex(iy,ix)];
	pdPtr->fTraVecL2u = *rotMat * pdPtr->fTraVecL2u + *traVec;
	pdPtr->fRotMatL2u = *rotMat * pdPtr->fRotMatL2u;
	pdPtr->fRotMatL2uT = transpose(pdPtr->fRotMatL2u);
      }
}

//______________________________________________________________________________
Int_t TVpPointDetectorArray::WritePdaFile(Char_t *fileName, Int_t withErrors)
{
  // Write the file in PDA format

  FILE *fp;

  // Open the file
  if ((fp = fopen(fileName, "w")) == NULL)
    {
      fprintf(stderr, "TVpPointDetectorArray::WritePdaFile: Cannot open the file %s:",
	      fileName);
      perror("");
      return 1;
    }
  WritePdaFile(fp, withErrors);
  
  // Close the file
  fclose(fp);
  
  return 0;
}

//______________________________________________________________________________
Int_t TVpPointDetectorArray::WritePdaFile(FILE *fp, Int_t withErrors)
{
  // Write the file in PDA format

  ULong_t numOfHistories;

  // Write the file's header
  if (withErrors == 0)
    fprintf(fp, "# Format: PDA 1.0\n");
  else
    fprintf(fp, "# Format: PDA 1.0 e\n");
  fprintf(fp, "# Dimension: %d %d\n", fNx, fNy);
  fprintf(fp, "# NumOfChannels: %d\n", fNumOfChannels);
  fprintf(fp, "# Channel range: %e %e\n", fMinEnergy, fMaxEnergy);
  if (withErrors == 1)
    {
      numOfHistories = (fNx > 0 || fNy > 0) ? fPointDetector[0].fNumOfHistories : 0;
      fprintf(fp, "# NumOfHistories: %lu\n", numOfHistories);
    }

  // Write point detector info
  for (Int_t j = 0; j < fNy; j++)
    for (Int_t i = 0; i < fNx; i++)
      fPointDetector[GetIndex(j,i)].WritePDS(fp, withErrors);
   
  return 0;
}

//______________________________________________________________________________
Int_t TVpPointDetectorArray::ReadPdaFile(Char_t *fileName)
{
  // Read the file in PDA format. Allocate space for point detectors.

  FILE *fp;

  // Open the file
  if ((fp = fopen(fileName, "r")) == NULL)
    {
      fprintf(stderr, "TVpPointDetectorArray::ReadPdaFile: ");
      fprintf(stderr, "Cannot open the file %s:", fileName);
      perror("");
      return 1;
    }
  ReadPdaFile(fp);

  // Close the file
  fclose(fp);

  return 0;
}

//______________________________________________________________________________
Int_t TVpPointDetectorArray::ReadPdaFile(FILE *fp)
{
  // Read the file in PDA format. Allocate space for point detectors.

  Char_t line[257];      // one line
  const Char_t *formatSpec = "# Format: PDA 1.0";
  const Char_t *formatSpecE = "# Format: PDA 1.0 e";
  Int_t withErrors;

  // Read the file's header
  fscanf(fp, "%256[^\n]%*[^\n]", line);
  if (strncmp(line, formatSpecE, strlen(formatSpecE)) == 0)
    withErrors = 1;
  else if (strncmp(line, formatSpec, strlen(formatSpec)) == 0)
    withErrors = 0;
  else
    {
      fprintf(stderr, "TVpPointDetectorArray::ReadPdaFile: ");
      fprintf(stderr, "Bad file format.\nIs: %s\nShould be: %s\n", line, formatSpec);
      return 2;
    }
  // Read the dimensions and then skip the rest of the line
  if (fscanf(fp, "\n# Dimension: %d %d%*[^\n]", &fNx, &fNy) != 2)
    {
      fprintf(stderr, "TVpPointDetectorArray::ReadPdaFile: ");
      fprintf(stderr, "Line _Dimension: _ incorrectly scanned.\n");
      return 4;
    } 
  // Read the number of channels and then skip the rest of the line
  if (fscanf(fp, "\n# NumOfChannels: %d%*[^\n]", &fNumOfChannels) != 1)
    {
      fprintf(stderr, "TVpPointDetectorArray::ReadPdaFile: ");
      fprintf(stderr, "Line _NumOfChannels: _ incorrectly scanned.\n");
      return 5;
    } 
  // Read the channel range and then skip the rest of the line
  if (fscanf(fp, "\n# Channel range: %lf %lf%*[^\n]", &fMinEnergy, &fMaxEnergy) != 2)
    {
      fprintf(stderr, "TVpPointDetectorArray::ReadPdaFile: ");
      fprintf(stderr, "Line _Channel range: _ incorrectly scanned.\n");
      return 6;
    } 
  // Read the number of histories if available
  ULong_t numOfHistories;
  if (withErrors == 1)
    {
      if (fscanf(fp, "\n# NumOfHistories: %lu%*[^\n]", &numOfHistories) != 1)
	{
	  fprintf(stderr, "TVpPointDetectorArray::ReadPdaFile: ");
	  fprintf(stderr, "Line _NumOfHistories: _ incorrectly scanned.\n");
	  return 7;
	}
    }
  else
    numOfHistories = 1;  // an arbitrary value
  
  fNumOfDets = fNx * fNy;
  
  // Allocate the point detector array
  fPointDetector = new TVpPointDetector[fNumOfDets];
  if (fPointDetector == 0)
    {
      fprintf(stderr, "TVpPointDetectorArray::ReadPdaFile: ");
      fprintf(stderr, "Allocation of Point Detector Array failed.\n");
      return 7;
    }
  
  // Initialize point detectors and read data
  for (Int_t j = 0; j < fNy; j++)
    for (Int_t i = 0; i < fNx; i++)
      {
	Int_t index = GetIndex(j,i);
	fPointDetector[index].SetNumOfChannels(fNumOfChannels, fMinEnergy, fMaxEnergy);
	fPointDetector[index].ReadPDS(fp, withErrors, numOfHistories);
      }
  
  return 0;
}

//______________________________________________________________________________
void TVpPointDetectorArray::MultiplyChannelContentBy(Double_t factor)
{
  // Multiply the contents of all channels of all point detectors

  for (Int_t index = 0; index < fNumOfDets; index++)
    fPointDetector[index].MultiplyChannelContentBy(factor);
}

//______________________________________________________________________________
void TVpPointDetectorArray::Zero()
{
  // Zero the contents of all channels of all point detectors

  for (Int_t index = 0; index < fNumOfDets; index++)
    fPointDetector[index].Zero();
}

//______________________________________________________________________________
void TVpPointDetectorArray::SetScoredQuantity(TVpPointDetector::EScoredQuantity
					      scoredQuantity)
{
  // Set the scored quantity for each point detector.

  for (Int_t index = 0; index < fNumOfDets; index++)
    fPointDetector[index].SetScoredQuantity(scoredQuantity);
}

//______________________________________________________________________________
void TVpPointDetectorArray::UpdateCountersAtEndOfHistory(ULong_t numOfHistories)
{
  // Update counters of all point detectors at the end of a history

  for (Int_t index = 0; index < fNumOfDets; index++)
    fPointDetector[index].UpdateCountersAtEndOfHistory(numOfHistories);
}

//______________________________________________________________________________
void TVpPointDetectorArray::PrintYProfile(Int_t sliceX, std::ostream &out) const
{
  if (GetNumOfDetectors() == 0)
    return;

  if (sliceX < 0 || sliceX >= fNx)
    {
      std::cerr << "TVpPointDetectorArray::GetYProfile: Slice out of range: "
		<< sliceX << '\n';
      return;
    }
  
  out << std::scientific;
  Double_t sum;
  for (Int_t iy = 0; iy < fNy; iy++)
    {
      sum = GetMeanEstimate(iy, sliceX);
      out << sum << '\n';
    }
}

//______________________________________________________________________________
Int_t TVpPointDetectorArray::WriteValuesToTabFile(Char_t *fileName)
{
  // Write values to a file in Tab format

  FILE *fp;

  // Open the file
  if ((fp = fopen(fileName, "w")) == NULL)
    {
      fprintf(stderr, "TVpPointDetectorArray::WriteValuesToTabFile: ");
      fprintf(stderr, "Cannot open the file %s:", fileName);
      perror("");
      return 1;
    }

  Double_t mean;
  for (Int_t ix = 0; ix < fNy; ix++)
    {
      for (Int_t iy = 0; iy < fNy; iy++)
	{
	  mean = GetMeanEstimate(fNy-iy-1, fNx-ix-1);
	  if (iy != 0)
	    fprintf(fp, " ");
	  fprintf(fp, "%e", mean);
	}
      fprintf(fp, "\n");
    }

  // Close the file
  fclose(fp);

  return 0;
}

//______________________________________________________________________________
Int_t TVpPointDetectorArray::WriteErrorsToTabFile(Char_t *fileName)
{
  // Write errors (standard deviation) to a file in Tab format

  FILE *fp;

  // Open the file
  if ((fp = fopen(fileName, "w")) == NULL)
    {
      fprintf(stderr, "TVpPointDetectorArray::WriteErrorsToTabFile: ");
      fprintf(stderr, "Cannot open the file %s:", fileName);
      perror("");
      return 1;
    }

  Double_t variance;
  for (Int_t ix = 0; ix < fNx; ix++)
    {
      for (Int_t iy = 0; iy < fNy; iy++)
	{
	  variance = GetVarianceEstimate(fNy-iy-1, fNx-ix-1);
	  if (iy != 0)
	    fprintf(fp, " ");
	  fprintf(fp, "%e", sqrt(variance));
	}
      fprintf(fp, "\n");
    }

  // Close the file
  fclose(fp);

  return 0;
}


//______________________________________________________________________________
void TVpPointDetectorArray::FillFomEvent(UInt_t time, TTree *treePtr,
					 FomEvent *fomEvent) const
{
  
}

//______________________________________________________________________________
TH2F *TVpPointDetectorArray::GetImage(Int_t isRotated) const
{
  // Get 2D image.

  if (GetNumOfDetectors() == 0)
    return 0;

  TH2F *h;
  
  switch (isRotated)
    {
      case 0 :
	// Image x-axis is paralel to CT's y-axis. Image y-axis is paralel to CT's x-axis.
	h = new TH2F("image", " PDA Image", fNx, 0, fNx, fNy, 0, fNy);
	for (Int_t ix = 0; ix < fNx; ix++)
	  for (Int_t iy = 0; iy < fNy; iy++)
	    h->SetCellContent(ix+1, iy+1, GetMeanEstimate(iy, ix));
	h->SetXTitle("idx");
	h->SetYTitle("idy");
	break;

      case 1:
	// Image x-axis is paralel to CT's y-axis. Image y-axis is paralel to CT's x-axis.
	h = new TH2F("image", " PDA Image", fNy, 0, fNy, fNx, 0, fNx);
	for (Int_t ix = 0; ix < fNx; ix++)
	  for (Int_t iy = 0; iy < fNy; iy++)
	    h->SetCellContent(iy+1, ix+1, GetMeanEstimate(iy, ix));
	h->SetXTitle("idy");
	h->SetYTitle("idx");
	break;
      default:
	std::cerr << "Error: TVpPointDetectorArray:: GetImage: isRotated is: "
		  << isRotated << " . It should be 0 or 1" << std::endl;
	return 0;
    }

  return h;
}

//______________________________________________________________________________
TH1F *TVpPointDetectorArray::GetYProfile(Int_t sliceX, Int_t isCm, const Char_t *name) const
{
  // Get histogram.  isCm is ignored.

  if (GetNumOfDetectors() == 0)
    return 0;

  if (sliceX < 0 || sliceX >= fNx)
    {
      std::cerr << "TVpPointDetectorArray::GetYProfile: Slice out of range: "
		<< sliceX << '\n';
      return 0;
    }
  
  if (name == 0)
    name = "profile";
  Double_t sum, sum2;
  TH1F *h = new TH1F(name, name, fNy, 0, fNy);
  for (Int_t iy = 0; iy < fNy; iy++)
    {
      sum = GetMeanEstimate(iy, sliceX);
      sum2 = GetVarianceEstimate(iy, sliceX);
      h->SetBinContent(iy+1, sum);
    }
  h->SetXTitle("idy");
  h->SetYTitle("I / (keV cm^{-2})");
  
  return h;
}

//______________________________________________________________________________
TGraph *TVpPointDetectorArray::GetYProfileGraph(Int_t sliceX, Int_t isCm,
						const Char_t *name) const
{
  // Get histogram.  isCm is ignored.

  if (GetNumOfDetectors() == 0)
    return 0;

  if (sliceX < 0 || sliceX >= fNx)
    {
      std::cerr << "TVpPointDetectorArray::GetYProfileGraph: Slice out of range: "
		<< sliceX << '\n';
      return 0;
    }
  
  if (name == 0)
    name = "profile";
  Double_t sum;
  //  TVpVector3 posUni;
  TGraph *g = new TGraph(fNy);
  g->SetTitle(name);
  for (Int_t iy = 0; iy < fNy; iy++)
    {
      sum = GetMeanEstimate(iy, sliceX);
      // sum2 = GetVarianceEstimate(iy, sliceX);  // not used
      // GetPosition(iy, sliceX, &posUni);  // position is not stored in PDA files
      g->SetPoint(iy, iy, sum);
    }

  return g;
}

//______________________________________________________________________________
void TVpPointDetectorArray::Draw(Int_t option) const
{
  // Draw individual point detectors
  
  Int_t elementOption = option & kDrawElementPoints | option & kDrawElementAxes;
  if (option & kDrawPdaAxes)
    DrawLcsAxes(10, 10, 10);
  for (Int_t iy = 0; iy < fNy; ++iy)
    for (Int_t ix = 0; ix < fNx; ++ix)
      fPointDetector[GetIndex(iy,ix)].Draw(elementOption);
}

//______________________________________________________________________________
TH1F *TVpPointDetectorArray::GetHistEnergySpectrum(Int_t ix, Int_t iy) const
{
  // Return histogram containing energy spectrum of the scored quantity for
  // the point detector (ix, iy)
  
  Int_t index = GetIndex(iy, ix);
  return fPointDetector[index].GetHistEnergySpectrum();
}
