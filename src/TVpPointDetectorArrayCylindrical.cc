//______________________________________________________________________________
//
// A TVpPointDetectorArrayCylindrical is a point detector array where
// point dtectors are positioned on a cylindrical surface. The
// cylindrical surface axis coincides with the X direction.


#include <iostream>
#include <string>
#include <cmath>
#include "TVpPointDetectorArrayCylindrical.h"

ClassImp(TVpPointDetectorArrayCylindrical)

//______________________________________________________________________________
TVpPointDetectorArrayCylindrical::TVpPointDetectorArrayCylindrical()
  : TVpPointDetectorArray()
{
  // Default constructor

  fArrayRadius = fArrayWidth = fArrayLength = 0.0;
}

//______________________________________________________________________________
TVpPointDetectorArrayCylindrical::TVpPointDetectorArrayCylindrical
( Double_t arrayWidth, Double_t arrayLength, Double_t arrayRadius,
  Int_t na, Int_t nr, Int_t numOfChannels,
  Double_t minEnergy, Double_t maxEnergy,
  TVpDetectorResponse *detectorResponse)
  : TVpPointDetectorArray(na, nr, numOfChannels, minEnergy, maxEnergy, detectorResponse)
{
  // Constructor with complete initialization. An example:
  // 
  //   TVpPointDetectorArrayCylindrical *pda = new TVpPointDetectorArrayCylindrical(
  //      pdaRadius,      // A vector from the PDA center to the PDA surface
  //      pdaWidth,       // PDA width
  //      pdaLength,       // PDA length
  //      pdaNa,          // Number of point detectors in the axis direction
  //      pdaNr,          // Number of point detectors in tangential direction
  //      1,              // Histogram's number of channels
  //      0.0,            // Histogram's min energy in keV
  //      1e+30,          // Histogram's max energy in keV
  //      aerPtr);        // Pointer to the detector response class

  fArrayWidth = (na == 1) ? 0.0 : arrayWidth;
  fArrayLength = (nr == 1) ? 0.0 : arrayLength;
  fArrayRadius = arrayRadius;
  SetDefaultPosition();
}

//______________________________________________________________________________
TVpPointDetectorArrayCylindrical::TVpPointDetectorArrayCylindrical(const Char_t *fileName)
{
  // Constructor.  Initialize data members from a file.
  //
  // Input parameters:
  // - fileName - name of the file in PDAC format

  ReadPdaFile(fileName);
}

//______________________________________________________________________________
void TVpPointDetectorArrayCylindrical::PrintStatus(std::ostream &out) const
{
  // Print the object status

  out << "<TVpPointDetectorArrayCylindrical>\n"
      << "$Id: TVpPointDetectorArrayCylindrical.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "ArrayWidth: " << fArrayWidth << " cm\n"
      << "ArrayLength: " << fArrayLength << " cm\n"
      << "ArrayRadius: " << fArrayRadius << " cm\n";
  TVpPointDetectorArray::PrintStatus(out);
  out << "</TVpPointDetectorArrayCylindrical>\n";  
}

//______________________________________________________________________________
void TVpPointDetectorArrayCylindrical::SetDefaultPosition()
{
  // Set dafault positions and orientations of all point dectors

  Double_t stepSizeX = (fNx > 1) ? fArrayWidth / (fNx-1) : 0.0;
  Double_t stepAngle = (fNy > 1) ? fArrayLength / ((fNy-1)*fArrayRadius) : 0.0;
  TVpMatrix3x3 rotMat;
  TVpVector3 traVec;
  TVpVector3 axis(1,0,0);
  Double_t x, y, z, angle, angle0;
  for (Int_t ix = 0; ix < fNx; ++ix)
    for (Int_t iy = 0; iy < fNy; ++iy)
      {
	angle0 = 1.5*M_PI - 0.5*fArrayLength / fArrayRadius;
	angle = angle0 + iy*stepAngle;
	x = ix * stepSizeX - 0.5*fArrayWidth;
	y = fArrayRadius * cos(angle);
	z = fArrayRadius * sin(angle);
	traVec.Set(x, y, z);
	rotMat = rotationMatrix(axis, angle - 1.5*M_PI);
	fPointDetector[GetIndex(iy,ix)].SetActiveTransformation(&rotMat, &traVec);
      }
}

//______________________________________________________________________________
Int_t TVpPointDetectorArrayCylindrical::WritePdaFile
(const Char_t *fileName, Int_t withErrors)
{
  // Write the file in PDAC format

  FILE *fp;

  // Open the file
  if ((fp = fopen(fileName, "w")) == NULL)
    {
      fprintf(stderr, "TVpPointDetectorArrayCylindrical::WritePdaFile: ");
      fprintf(stderr, "Cannot open the file %s:", fileName);
      perror("");
      return 1;
    }

  // Write the file's header
  fprintf(fp, "# Format: PDAC 2.0\n");
  fprintf(fp, "# ArrayWidth: %g cm\n", fArrayWidth);
  fprintf(fp, "# ArrayLength: %g cm\n", fArrayLength);
  fprintf(fp, "# ArrayRadius: %g cm\n", fArrayRadius);
  fprintf(fp, "<TVpPointDetectorArray>\n");
  TVpPointDetectorArray::WritePdaFile(fp, withErrors);
  fprintf(fp, "</TVpPointDetectorArray>\n");
  
  // Close the file
  fclose(fp);
  
  return 0;
}

//______________________________________________________________________________
Int_t TVpPointDetectorArrayCylindrical::ReadPdaFile(const Char_t *fileName)
{
  // Read the file in PDAC format

  FILE *fp;
  Char_t line[257];      // one line
  const Char_t *formatSpec = "# Format: PDAC 2.0";

  // Open the file
  if ((fp = fopen(fileName, "r")) == NULL)
    {
      fprintf(stderr, "TVpPointDetectorArrayCylindrical::ReadPdaFile: ");
      fprintf(stderr, "Cannot open the file %s:", fileName);
      perror("");
      return 1;
    }

  // Read the file's header
  fscanf(fp, "%256[^\n]%*[^\n]", line);
  
  if (strncmp(line, formatSpec, strlen(formatSpec)) != 0)
    {
      fprintf(stderr, "TVpPointDetectorArrayCylindrical::ReadPdaFile: ");
      fprintf(stderr, "Bad file format.\nIs: %s\nShould be: %s\n", line, formatSpec);
      return 2;
    }
  if (fscanf(fp, "\n# ArrayWidth: %lf%*[^\n]", &fArrayWidth) != 1)
    {
      fprintf(stderr, "TVpPointDetectorArrayCylindrical::ReadPdaFile: ");
      fprintf(stderr, "Line _ArrayWidth: _ incorrectly scanned.\n");
      return 3;
    } 
  if (fscanf(fp, "\n# ArrayLength: %lf%*[^\n]", &fArrayLength) != 1)
    {
      fprintf(stderr, "TVpPointDetectorArrayCylindrical::ReadPdaFile: ");
      fprintf(stderr, "Line _ArrayLength: _ incorrectly scanned.\n");
      return 4;
    } 
  if (fscanf(fp, "\n# ArrayRadius: %lf%*[^\n]", &fArrayRadius) != 1)
    {
      fprintf(stderr, "TVpPointDetectorArrayCylindrical::ReadPdaFile: ");
      fprintf(stderr, "Line _ArrayRadius: _ incorrectly scanned.\n");
      return 5;
    } 
  fscanf(fp, "%s\n", line);  // Read <TVpPointDetectorArray>
  TVpPointDetectorArray::ReadPdaFile(fp);
  fclose(fp);
  return 0;
}

//______________________________________________________________________________
TVpPointDetectorArray* TVpPointDetectorArrayCylindrical::Clone
(Int_t numDivX, Int_t numDivY) const
{
  // Clone the PDA

  // FIX
#if 0
  TVpPointDetectorArrayCylindrical *pda = new TVpPointDetectorArrayCylindrical(*this);
  return pda;
#endif
  return 0;
}


//______________________________________________________________________________
TH2F *TVpPointDetectorArrayCylindrical::GetImage(Int_t isRotated) const
{
  // Get image

  if (GetNumOfDetectors() == 0)
    return 0;
  TH2F *h;
  
  switch (isRotated)
    {
      case 0: // x-axis is parallel to the cylinder axis.  Use for whole body images.
	h = new TH2F("image", "Cylindrical PDA projection",
		     fNx, -fArrayWidth, fArrayWidth, fNy, -fArrayLength, fArrayLength);
	for (Int_t ix = 0; ix < fNx; ix++)
	  for (Int_t iy = 0; iy < fNy; iy++)
	    h->SetCellContent(ix+1, iy+1, GetMeanEstimate(iy, ix));
	h->SetXTitle("ax / cm");
	h->SetYTitle("cy / cm");
	break;
	
      case 1: // x-axis is perpendicular to the cylinder axis.  Use for short scans.
	h = new TH2F("image", "Cylindrical PDA projection",
		     fNy, -fArrayLength, fArrayLength, fNx, -fArrayWidth, fArrayWidth);
	for (Int_t ix = 0; ix < fNx; ix++)
	  for (Int_t iy = 0; iy < fNy; iy++)
	    h->SetCellContent(iy+1, ix+1, GetMeanEstimate(iy, ix));
	h->SetXTitle("cy / cm");
	h->SetYTitle("ax / cm");
	break;
      default:
	std::cerr << "Error: TVpPointDetectorArrayCylindrical:: GetImage: isRotated is: "
		  << isRotated << " . It should be 0 or 1" << std::endl;
	return 0;
    }
  return h;
}

//______________________________________________________________________________
TH1F *TVpPointDetectorArrayCylindrical::GetYProfile(Int_t sliceX, Int_t isCm) const
{
  // Get histogram
  
  if (GetNumOfDetectors() == 0)
    return 0;

  if (sliceX < 0 || sliceX >= fNx)
    {
      std::cerr << "TVpPointDetectorArrayCylindrical::GetYProfile "
		<< "Slice out of range: " << sliceX << '\n';
      return 0;
    }

  TH1F *h;  
  switch (isCm)
    {
    case 0:
      h = new TH1F("h", "Profile", fNy, 0, fNy);
      h->SetXTitle("iy");
      break;
    case 1:
      h = new TH1F("h", "Profile", fNy, -fArrayLength, fArrayLength);
      h->SetXTitle("cy / cm");
      break;
    }
  Double_t sum, sum2;
  for (Int_t iy = 0; iy < fNy; iy++)
    {
      sum = GetMeanEstimate(iy, sliceX);
      sum2 = GetVarianceEstimate(iy, sliceX);
      h->SetBinContent(iy+1,sum);
      h->SetBinError(iy+1, sqrt(sum2));
    }
  h->SetYTitle("I / (keV cm^{-2})");
  
  return h;
}
