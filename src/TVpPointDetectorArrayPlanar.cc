//______________________________________________________________________________
//
// A TVpPointDetectorArrayPlanar class is a rectangular fNx x fNy grid
// of point detectors.
//
// .  .  .  . fNy-1
//
// .  .  .  .
//
// .  .  .  . 0
// 0        fNx-1
// 
// Position and normal of the point detector with indices (i,j), where
// i=0,...,fNx-1 and  j=0,...,fNy-1 is:
//
// x_i = i * fArraySizeX / (fNx-1) - fArraySizeX/2
// y_j = j * fArraySizeY / (fNy-1) - fArraySizeY/2
//
// Example: 
// TVpPointDetectorArrayPlanar *pda = new TVpPointDetectorArrayPlanar(
//   arraySizeX,   // Array size in X-direction
//   arraySizeY,   // Array size in Y-direction
//   pdaNx,        // Number of point detectors in the X direction (usually CT axis)
//   pdaNy,        // Number of point detectors in the Y direction
//   1,            // Histogram's number of channels
//   0.0,          // Histogram's min energy in keV
//   1e+30,        // Histogram's max energy in keV
//   aerPtr);      // Pointer to a detector response object

#include <string>
#include <iostream>
#include "TVpPointDetectorArrayPlanar.h"

ClassImp(TVpPointDetectorArrayPlanar)

//______________________________________________________________________________
TVpPointDetectorArrayPlanar::TVpPointDetectorArrayPlanar()
  : TVpPointDetectorArray()
{
  // Default constructor

  fArraySizeX = fArraySizeY = 0.0;
}

//______________________________________________________________________________
TVpPointDetectorArrayPlanar::TVpPointDetectorArrayPlanar
(Double_t arraySizeX, Double_t arraySizeY, Int_t nx, Int_t ny, Int_t numOfChannels,
 Double_t minEnergy, Double_t maxEnergy, TVpDetectorResponse *detectorResponse)
  : TVpPointDetectorArray(nx, ny, numOfChannels, minEnergy, maxEnergy,detectorResponse)
{
  // Constructor with complete initialization.
  //
  // Input parameters:
  // - arraySizeX - array size in X-direction
  // - arraySizeY - array size in Y-direction
  // - nx - number of point detectors in the X direction
  // - ny - number of point detectors in the Y direction
  // - numOfChannels - histogram's number of channels
  // - minEnergy - histogram's min energy in keV
  // - maxEnergy - histogram's max energy in keV
  // - detectorResponse - pointer to a detector response object, or 0.

  fArraySizeX = (nx == 1) ? 0.0 : arraySizeX;
  fArraySizeY = (ny == 1) ? 0.0 : arraySizeY;
  SetDefaultPosition();
}

//______________________________________________________________________________
TVpPointDetectorArrayPlanar::TVpPointDetectorArrayPlanar(Char_t *fileName)
{
  // Constructor.  Initialize data members from a file.
  //
  // Input parameters:
  // - fileName - name of the file in PDAP format
  
  ReadPdaFile(fileName);
}

//______________________________________________________________________________
void TVpPointDetectorArrayPlanar::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  // Input parameters:
  // - out - the output stream, cout is the default

  out << "<TVpPointDetectorArrayPlanar>\n"
      << "$Id: TVpPointDetectorArrayPlanar.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Size: " << fArraySizeX << ' ' << fArraySizeY << " cm\n";
  TVpPointDetectorArray::PrintStatus(out);
  out << "</TVpPointDetectorArrayPlanar>\n";

}

//______________________________________________________________________________
void TVpPointDetectorArrayPlanar::SetDefaultPosition()
{
  // Set dafault positions and orientations of all point dectors

  TVpVector3 traVec;
  TVpMatrix3x3 rotMat(1,0,0, 0,1,0, 0,0,1);
  Double_t stepSizeX = (fNx > 1) ? fArraySizeX / (fNx-1) : 0.0;
  Double_t stepSizeY = (fNy > 1) ? fArraySizeY / (fNy-1) : 0.0;
  TVpVector3 ex = TVpVector3(stepSizeX, 0, 0);
  TVpVector3 ey = TVpVector3(0, stepSizeY, 0);
  TVpVector3 shift = TVpVector3(0.5*fArraySizeX, 0.5*fArraySizeY, 0);
  for (Int_t ix = 0; ix < fNx; ++ix)
    for (Int_t iy = 0; iy < fNy; ++iy)
      {
	traVec = ix*ex + iy*ey - shift;
	fPointDetector[GetIndex(iy,ix)].SetActiveTransformation(&rotMat, &traVec);
      }
}

//______________________________________________________________________________
Int_t TVpPointDetectorArrayPlanar::WritePdaFile
(Char_t *fileName, Int_t withErrors)
{
  // Write the file in PDAC format

  FILE *fp;

  // Open the file
  if ((fp = fopen(fileName, "w")) == NULL)
    {
      fprintf(stderr, "TVpPointDetectorArrayPlanar::WritePdaFile: ");
      fprintf(stderr, "Cannot open the file %s:", fileName);
      perror("");
      return 1;
    }

  // Write the file's header
  fprintf(fp, "# Format: PDAP 2.0\n");
  fprintf(fp, "# ArraySizeX: %g cm\n", fArraySizeX);
  fprintf(fp, "# ArraySizeY: %g cm\n", fArraySizeY);
  fprintf(fp, "<TVpPointDetectorArray>\n");
  TVpPointDetectorArray::WritePdaFile(fp, withErrors);
  fprintf(fp, "</TVpPointDetectorArray>\n");
  
  // Close the file
  fclose(fp);
  
  return 0;
}

//______________________________________________________________________________
Int_t TVpPointDetectorArrayPlanar::ReadPdaFile(Char_t *fileName)
{
  // Read the file in PDAP format

  FILE *fp;
  Char_t line[257];      // one line
  const Char_t *formatSpec = "# Format: PDAP 2.0";

  // Open the file
  if ((fp = fopen(fileName, "r")) == NULL)
    {
      fprintf(stderr, "TVpPointDetectorArrayPlanar::ReadPdaFile: ");
      fprintf(stderr, "Cannot open the file %s:", fileName);
      perror("");
      return 1;
    }

  // Read the file's header
  fscanf(fp, "%256[^\n]%*[^\n]", line);
  
  if (strncmp(line, formatSpec, strlen(formatSpec)) != 0)
    {
      fprintf(stderr, "TVpPointDetectorArrayPlanar::ReadPdaFile: ");
      fprintf(stderr, "Bad file format.\nIs: %s\nShould be: %s\n", line, formatSpec);
      return 2;
    }
  if (fscanf(fp, "\n# ArraySizeX: %lf%*[^\n]", &fArraySizeX) != 1)
    {
      fprintf(stderr, "TVpPointDetectorArrayPlanar::ReadPdaFile: ");
      fprintf(stderr, "Line _ArraySizeX: _ incorrectly scanned.\n");
      return 3;
    } 
  if (fscanf(fp, "\n# ArraySizeY: %lf%*[^\n]", &fArraySizeY) != 1)
    {
      fprintf(stderr, "TVpPointDetectorArrayPlanar::ReadPdaFile: ");
      fprintf(stderr, "Line _ArraySizeY: _ incorrectly scanned.\n");
      return 4;
    } 
  fscanf(fp, "%s\n", line);  // Read <TVpPointDetectorArray>
  TVpPointDetectorArray::ReadPdaFile(fp);
  fclose(fp);
  return 0;
}

//______________________________________________________________________________
TVpPointDetectorArray* TVpPointDetectorArrayPlanar::Clone
(Int_t numDivX, Int_t numDivY) const
{
  // Clone the PDA.
  //
  // Input parameters:
  // - numDivX - how many times is the number of point detectors increased in X direction
  // - numDivY - how many times is the number of point detectors increased in Y direction
  //
  // Method:
  // The number of point detectors is multiplied, the distance between them is
  // decreased proportionally.

  // FIX
#if 0
  TVpPointDetectorArrayPlanar *pda = new TVpPointDetectorArrayPlanar
    (
     numDivX * fNx,
     numDivY * fNy,
     fNumOfChannels,
     fMinEnergy,
     fMaxEnergy);
  return pda;
#endif
  return 0;
}

