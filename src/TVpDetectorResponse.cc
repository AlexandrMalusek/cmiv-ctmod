//______________________________________________________________________________
//
// The TVpDetectorResponse class implements the response of a single detector
// element to a photon impact. Angle-energy dependence of the response
// function is considered, f=f(E, xi), where E is the incident photon energy
// in keV and xi = cos(theta) is the cosine of the incident angle.  The value
// theta = 0 (xi = 1) corresponds to a photon trajectory perpendicular to the
// detector element surface, the value theta = Pi/2 (xi = 0) corresponds to a
// a trajectory parallel to the surface.
//
//  normal   incident photon
//       |  / 
//       |-/-- theta      
//       |/
//  ----------------
//  detector element
//  ----------------
//
// The interpretation of the function f(E, xi) is application-dependent.
// Originally, it was the energy absorption efficiency function defined as
// f(E, xi) = epsilon / E where epsilon is the energy imparted to the detector
// by a photon with energy E.  For air kerma calculations, the f(E, xi)
// function specifies mass energy transfer coefficient for air.
//
// Energy and Xi values are equidistant. Their values and corresponding
// indices are:
//
// Xi:      1/fDimXi,                2/fDimXim,               ..., 1
// ix:      0,                       1,                       ..., fDimXi-1
// Energy:  1*fMaxEnergy/fDimEnergy, 2*fMaxEnergy/fDimEnergy, ..., fMaxEnergy
// ie:      0,                       1,                       ..., fDimEnergy-1
//
// Note that values of f(E, xi) for Xi = 0 correspond to photons impinging
// parallel to the surface.  These values do not have a good physical meaning
// since limits for Xi->0 from left and right differ.  Therefore the lowest
// tabulated value is 1/fDimXi and a linear extrapolation is used in the
// interval (0, 1/fDimXi).
//
// Similarly, f(E, xi) for E = 0 do not have physical meaning since photons
// with zero energy are non-existent.  Thus the lowest tabulated value is
// 1*fMaxEnergy/fDimEnergy and a linear extrapolation is used in the interval
// (0, 1*fMaxEnergy/fDimEnergy).
//
// Bilinear interpolation (for a description see Numerical recipes) is used to
// obtain values of f(E, xi).  If fDimXi = 1 then a linear interpolation in E
// is used.  Angle-independent detector response function is used for air
// kerma calculations.  fDimEnergy = 1 does not have a practical meaning and
// is not allowed.
//
// The detector response takes into account an antiscatter grid if
// antiscatterGridPtr != 0. The grid transmission formula then multiplies the
// scored quantity as it reduces the virtual photon's weight.  Note that this
// implementation of the detector response function does not depend on the
// azimuthal angle but the antiscatter grid has its strips oriented.  Thus it
// adds azimuthal dependence in the result.
//
// The detector response function is often calculated by Monte Carlo codes
// which can take electron transport into account like MCNP or Penelope.  It
// is stored in .aer files with given formats:
//
// Format 1.0 .aer file consists of a header which is followed by a
// (fDimEnergy, fDimXi) matrix of values of the detector response function.
// Here, fDimEnergy is the number of rows and fDimXi is the number of columns.
// An example of a typical header is:
//
// # Energy absorption efficiency function, format 1.0
// # Name: Gd2O2S, 133 mg/cm2, MCNP
// # Energy lines: 150
// # Xi columns: 64
// # Max energy: 150 // keV
//
// Errors of the detector response function are initialized to zero.  Current
// implementation does not use these errors for MC simulation.
//
// Format 2.0 .aer file consists of a header, (fDimEnergy, fDimXi) matrix of
// values of the detector response function, and the (fDimEnergy, fDimXi)
// matrix of 3*sigma absolute errors of the detector response function.  The
// header differs only in the version number.  An example:
//
// # Energy absorption efficiency function, format 2.0
// # Name: CsI, 0.6 mm, Penelope
// # Energy lines: 150
// # Xi columns: 64
// # Max energy: 150 // keV
//
// The detector response function and its errors can be displayed by GetAer(),
// GetAerAbsError(), and GetAerRelError().  See mac/detector/x_DrawAer.C.
//
// Streamer: OK
//
// Example:
//
// TVpDetectorResponse *aerPtr = new TVpDetectorResponse();
// aerPtr->ReadAerFile("detector/K_air.aer");
//
// TVpPointDetectorArrayPlanar *pda = new TVpPointDetectorArrayPlanar
//   (TVpVector3(0,0,0),                     // PDA center
//    (pdaSizeX/pdaNx)*TVpVector3(1,0,0),    // StepX
//    (pdaSizeY/pdaNy)*TVpVector3(0,1,0),    // StepY
//    pdaNx,                                 // Dimension in x
//    pdaNy,                                 // Dimension in y
//    1,                                     // Number of energy channels
//    0.0,                                   // Min energy of registered photon[keV]
//    1e+30,                                 // Max energy of registered photon[keV]
//    aerPtr);                               // Pointer to the detector response class
// pda->SetScoredQuantity(kAirKerma);

#include <iostream>
#include <ostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include "TVpDetectorResponse.h"
#include "TVpMath.h"

ClassImp(TVpDetectorResponse)

//______________________________________________________________________________
TVpDetectorResponse::TVpDetectorResponse()
{
  // Constructor
  
  fDimEnergy = fDimXi = 0;
  fMaxEnergy = fStepXi = fStepEnergy = 0.0;
  fAer = fRer = 0;
  fAntiscatterGridPtr = 0;
}

//______________________________________________________________________________
TVpDetectorResponse::~TVpDetectorResponse()
{
  // Destructor.  It does not destroy the antiscatter grid.
  
  delete [] fAer;
  delete [] fRer;
}

//______________________________________________________________________________
Double_t TVpDetectorResponse::GetResponse(Double_t xi, Double_t energy) const
{
  // Returns the detector response.
  //
  // Input:
  // - xi - cosine of the incident angle, (0, 1].
  // - energy - incident photon energy, (0, fMaxEnergy].
  //
  // Output:
  // Returns the interpolated value.
  //
  // Method:
  // If fDimXi = 1 then use linear interpolation for E.  If fDimXi > 1 then
  // use bilinear interpolation for E and Xi.  fDimEnergy must be greater than
  // 1, this should be checked during object's initialization.  The response
  // matrix is given in grid values, i.e. f_{ij}=f(Xi_i, E_j), Xi_i and E_j
  // are not bin edges.
  //
  // References:
  // - Numerical recipes in C.

  Double_t f;
  if (fDimXi == 1)
    // Case 1: Linear interpolation
    f = TVpMath::GetLinLinInterpolation(energy, fDimEnergy, fStepEnergy, fStepEnergy, fAer);
  else
    { // Case 2: Bilinear interpolation
      
      // Find indices ix, ie so that
      // xi1 = xi_{ix} <= xi < xi_{ix+1}
      // energy1 = energy_{ie} <= energy < energy_{ie+1}
      //
      // Note: (Int_t) -0.1 == 0 and therefeore the following algorithm
      // performs extrapolation when ix or energy are bellow minimum values.
      // The corresponding square is just shifted by 1.  The extrapolation may
      // result in a negative value.  This is checked at the end and negative
      // values are replaced by 0.0
      Int_t ix = (Int_t) (xi / fStepXi - 1);
      Int_t ie = (Int_t) (energy / fStepEnergy - 1);

      Double_t xi1 = (ix+1) * fStepXi;
      Double_t energy1 = (ie+1) * fStepEnergy;

      // Checks
      if (ix < 0 || ix >= fDimXi || ie < 0 || ie >= fDimEnergy)
	{
	  std::cerr << "Error: TVpDetectorResponse::GetResponse: index out of range"
		    << ", ix = " << ix << ", " << "ie = " << ie
		    << ", xi = " << xi << ", " << "energy = " << energy
		    << ", Returning 0.0." << std::endl;
	  return 0.0;
	}

      // The corresponding grid square values.  Optimize later.
      Double_t y1 = fAer[GetIndex(ix, ie)];
      Double_t y2 = fAer[GetIndex(ix+1, ie)];
      Double_t y3 = fAer[GetIndex(ix+1, ie+1)];
      Double_t y4 = fAer[GetIndex(ix, ie+1)];      
      // std::cout << y1 << ' ' << y2 << ' ' << y3 << ' ' << y4 << std::endl;

      Double_t t = (xi - xi1) / fStepXi;               // t is in [0,1)
      Double_t u = (energy - energy1) / fStepEnergy;   // u is in [0,1)
      // std::cout << "xi1 = " << xi1 << ", energy1 = " << energy1 << std::endl;
      // std::cout << "t = " << t << ", u = " << u << std::endl;

      f = (1.0 - t)*(1-u)*y1 + t*(1-u)*y2 + t*u*y3 + (1-t)*u*y4;
    }
  if (f < 0.0)
    {

      std::cerr << "Warning: TVpDetectorResponse::GetResponse: bad extrapolation"
		<< ", xi = " << xi << ", energy = " << energy << ", f = " << f
		<< ", Returning 0.0." << std::endl;
      f = 0.0;
    }
  return f;
}

//______________________________________________________________________________
Double_t TVpDetectorResponse::GetResponse(Double_t energy, const TVpVector3 &directionL) const
{
  // Use simple linear interpolation and the antiscatter grid transmission formula
  // to get the response
  
  Double_t xi = std::abs(directionL.fR[2]);   // Works for conventional geometries only
  Double_t response = GetResponse(xi, energy);
  if (fAntiscatterGridPtr != 0)
    response *= fAntiscatterGridPtr->Transmission(energy, directionL);
  return response;
}

//______________________________________________________________________________
void TVpDetectorResponse::PrintStatus(std::ostream &out) const
{
  // Print the object status.
  //
  // Input:
  // - out - the output stream
  //
  // Example output:
  // <TVpDetectorResponse>
  // Name:  CsI, 0.6 mm, Penelope
  // Energy dimension: 150
  // Xi dimension: 64
  // Max energy: 150 keV
  // Antiscater Grid: No
  // </TVpDetectorResponse>

  out << "<TVpDetectorResponse>\n"
      << "$Id: TVpDetectorResponse.cc 62 2009-06-27 10:54:08Z malusek $\n"
      << "Name: " << fName << '\n'
      << "Energy dimension: " << fDimEnergy << '\n'
      << "Xi dimension: " << fDimXi << '\n'
      << "Max energy: " << fMaxEnergy << " keV\n";
  if (fAntiscatterGridPtr != 0)
    fAntiscatterGridPtr->PrintStatus(out);
  else
    out << "Antiscater Grid: No\n";
  out << "</TVpDetectorResponse>\n";
}

//______________________________________________________________________________
void TVpDetectorResponse::SetAntiscatterGrid(TVpAntiscatterGrid *antiscatterGridPtr)
{
  // Set the antiscatter grid

  fAntiscatterGridPtr = antiscatterGridPtr;
}

//______________________________________________________________________________
Int_t TVpDetectorResponse::ReadAerFile(Char_t *fileName)
{
  // Read AER file.
  //
  // Input:
  // - fileName - name of the aer file
  //
  // Format 1.0 header example:
  // # Energy absorption efficiency function, format 1.0
  // # Name: Gd2O2S, 133 mg/cm2, MCNP
  // # Energy lines: 150
  // # Xi columns: 64
  // # Max energy: 150 // keV
  //
  // Format 2.0 header example:
  // # Energy absorption efficiency function, format 2.0
  // # Name: CsI, 0.6 mm, Penelope
  // # Energy lines: 150
  // # Xi columns: 64
  // # Max energy: 150 // keV 
  
  std::string recName;                 // The record name e.g. "# Size"
  std::string recRest;                 // The rest of the record
  Int_t formatVersion;                 // Format version (1.0, 2.0)

  // Read the AER file
  std::ifstream in(fileName);
  if (!in)
    {
      std::cerr << "Error: DetectorResponse::ReadAerFile: Cannot open file: "
		<< fileName << '\n';
      return 1;
    }

  // Read and process the file header
  // # Energy absorption efficiency function, format 1.0
  std::getline(in, recName);
  if (recName == std::string("# Energy absorption efficiency function, format 1.0"))
    formatVersion = 1;
  else if (recName == std::string("# Energy absorption efficiency function, format 2.0"))
    formatVersion = 2;
  else
    {
      std::cerr << "Error: TVpDetectorResponse::ReadAerFile: Record "
		<< "\"# Energy absorption efficiency function, format 1.0\" or \n"
		<< "\"# Energy absorption efficiency function, format 2.0\" is missing.\n";
      return 2;
    }

  // # Name:
  std::getline(in, recName, ':');
  if (recName != std::string("# Name"))
    {
      std::cerr << "Error: TVpDetectorResponse::ReadAerFile: Record "
		<< "\"# Name:\" is missing.\n";
      return 3;
    }
  std::getline(in, fName);

  // # Energy lines:
  std::getline(in, recName, ':');
  if (recName != std::string("# Energy lines"))
    {
      std::cerr << "Error: TVpDetectorResponse::ReadAerFile: Record "
		<< "\"# Energy lines:\" is missing.\n";
      return 4;
    }
  in >> fDimEnergy;
  std::getline(in, recRest);
  
  // # Xi columns:
  std::getline(in, recName, ':');
  if (recName != std::string("# Xi columns"))
    {
      std::cerr << "Error: TVpDetectorResponse::ReadAerFile: Record "
		<< "\"# Xi lines:\" is missing.\n";
      return 5;
    }
  in >> fDimXi;
  std::getline(in, recRest);
 
  // # Max energy:
  std::getline(in, recName, ':');
  if (recName != std::string("# Max energy"))
    {
      std::cerr << "Error: TVpDetectorResponse::ReadAerFile: Record "
		<< "\"# Max energy:\" is missing.\n";
      return 6;
    }
  in >> fMaxEnergy;
  std::getline(in, recRest);
  
  fStepXi = 1.0 / fDimXi;
  fStepEnergy = fMaxEnergy / fDimEnergy;
  
  // Allocate arrays
  fAer = new Double_t[fDimXi*fDimEnergy];
  fRer = new Double_t[fDimXi*fDimEnergy];
  
  // Read Aer table
  for (Int_t ie = 0; ie < fDimEnergy; ie++)
    for (Int_t ix = 0; ix < fDimXi; ix++)
      in >> fAer[GetIndex(ix,ie)];

  if (formatVersion == 2)
    {
      // Read Rer table
      for (Int_t ie = 0; ie < fDimEnergy; ie++)
	for (Int_t ix = 0; ix < fDimXi; ix++)
	  in >> fRer[GetIndex(ix,ie)];
    }

  in.close();
  return 0;
}


//______________________________________________________________________________
TH2F *TVpDetectorResponse::GetAer()
{
  // Return Angle-energy response histogram.  Values from the internal array
  // (fAer) are used, no interpolation routine is called.  The histogram's bin
  // content corresponds to bin's upper edges, see the GetAer(dimXi, ...).
  //
  // Note: The interpolation routine can be checked by comparing this 2D
  // histogram with a histogram created by the GetAer(dimXi, ...) function
  // with arguments.
  //
  // Example:
  // TVpDetectorResponse *aer = new TVpDetectorResponse();
  // aer->ReadAerFile("detector/Gd2O2S.aer");
  // TH2F *h = aer->GetAer();
  // h->Draw("surf0");

  if (fDimEnergy == 0 || fDimXi == 0)
    {
      std::cerr << "Error: GetAer: At least one dimension is 0.\n";
      return 0;
    }
  
  TH2F *h = new TH2F("aer", "aer",  fDimXi, 0, 1, fDimEnergy, 0, fMaxEnergy);
  for (Int_t ie = 0; ie < fDimEnergy; ie++)
    for (Int_t ix = 0; ix < fDimXi; ix++)
      h->SetBinContent(ix+1, ie+1, fAer[GetIndex(ix, ie)]);
  h->SetXTitle("#xi");
  h->GetXaxis()->SetTitleOffset(1.5);
  h->SetYTitle("E / keV");
  h->GetYaxis()->SetTitleOffset(2);
  h->SetZTitle("f");
  return h;
}

//______________________________________________________________________________
TH2F *TVpDetectorResponse::GetAerAbsError()
{
  // Return 3*sigma absolute error of Aer.  Values from the internal array
  // (fRer) are used, no interpolation routine is called.
  //
  // Example:
  // TVpDetectorResponse *aer = new TVpDetectorResponse();
  // aer->ReadAerFile("detector/Gd2O2S.aer");
  // TH2F *h = aer->GetAerAbsError();
  // h->Draw("surf0");

  if (fDimEnergy == 0 || fDimXi == 0)
    {
      std::cerr << "Error: GetAerAbsError: At least one dimension is 0.\n";
      return 0;
    }

  TH2F *h = new TH2F("aerr", "absolute error of aer, 3*sigma",
		     fDimXi, 0, 1, fDimEnergy, 0, fMaxEnergy);
  for (Int_t ie = 0; ie < fDimEnergy; ie++)
    for (Int_t ix = 0; ix < fDimXi; ix++)
      h->SetBinContent(ix+1, ie+1, fRer[GetIndex(ix, ie)]);
  h->SetXTitle("#xi");
  h->GetXaxis()->SetTitleOffset(1.5);
  h->SetYTitle("E / keV");
  h->GetYaxis()->SetTitleOffset(2);
  h->SetZTitle("absolute error");
  return h;
}

//______________________________________________________________________________
TH2F *TVpDetectorResponse::GetAerRelError()
{
  // Return 3*sigma relative error of Aer.  Values from the internal array
  // (fRer) are used, no interpolation routine is called.
  //
  // Example:
  // TVpDetectorResponse *aer = new TVpDetectorResponse();
  // aer->ReadAerFile("detector/Gd2O2S.aer");
  // TH2F *h = aer->GetAerRelError();
  // h->Draw("surf0");

  if (fDimEnergy == 0 || fDimXi == 0)
    {
      std::cerr << "Error: GetAerRelError: At least one dimension is 0.\n";
      return 0;
    }

  Double_t rerror;
  TH2F *h = new TH2F("rerr", "relative error of aer, 3*sigma",
		     fDimXi, 0, 1, fDimEnergy, 0, fMaxEnergy);
  for (Int_t ie = 0; ie < fDimEnergy; ie++)
    for (Int_t ix = 0; ix < fDimXi; ix++)
      {
	if (fAer[GetIndex(ix, ie)] > 0.0)
	  rerror = 100.0 * fRer[GetIndex(ix, ie)] / fAer[GetIndex(ix, ie)];
	else
	  rerror = 0.0;
	h->SetBinContent(ix+1, ie+1, rerror);
      }
  h->SetXTitle("#xi");
  h->GetXaxis()->SetTitleOffset(1.5);
  h->SetYTitle("E / keV");
  h->GetYaxis()->SetTitleOffset(2);
  h->SetZTitle("relative error [%]");
  return h;
}

//______________________________________________________________________________
TH2F *TVpDetectorResponse::GetAer(Int_t dimXi, Double_t minXi, Double_t maxXi,
				  Int_t dimEnergy, Double_t minEnergy, Double_t maxEnergy)
{
  // Return Angle-energy response histogram.  Values are calculated by the
  // interpolation function.  Values included in the grid are:
  //
  // minXi+dXi,         minXi+2*dXi,       ..., maxXi;
  // minEnergy+dEnergy, minEnergy+dEnergy, ..., maxEnergy
  //
  // where 
  // dXi = (maxXi - minXi) / fDimXi
  // dEnergy = (maxEnergy - minEnergy) / fDimEnergy
  //
  // For instance:
  // 0+1/64, 2/64, ..., 1
  // 0+1,    2,    ..., 150
  //
  // Thus the histogram's bin content corresponds to bin's upper edges.
  //
  // Note: The interpolation routine can be checked by comparing this 2D
  // histogram with a histogram created by the GetAer() function without
  // arguments.
  //
  // Example 1:  Compare in grid points.
  // TVpDetectorResponse *aer = new TVpDetectorResponse();
  // aer->ReadAerFile("detector/Gd2O2S.aer");
  // TH2F *h1 = aer->GetAer(64, 0.0, 1.0 , 150, 0.0, 150.0);
  // TH2F *h2 = aer->GetAer();
  // h1->Draw("surf0");

  if (fDimEnergy == 0 || fDimXi == 0)
    {
      std::cerr << "Error: GetAer: At least one dimension is 0.\n";
      return 0;
    }

  TH2F *h = new TH2F("aer", "aer", dimXi, minXi, maxXi,
		     dimEnergy, minEnergy, maxEnergy);
  Double_t dXi = (maxXi - minXi) / fDimXi;
  Double_t dEnergy = (maxEnergy - minEnergy) / fDimEnergy;
  for (Int_t ie = 1; ie <= dimEnergy; ie++)
    for (Int_t ix = 1; ix <= dimXi; ix++)
      {
	Double_t xi = minXi + dXi * ix;
	Double_t energy = minEnergy + dEnergy * ie;
	h->SetBinContent(ix, ie, GetResponse(xi, energy));
      }
  h->SetXTitle("#xi");
  h->GetXaxis()->SetTitleOffset(1.5);
  h->SetYTitle("E / keV");
  h->GetYaxis()->SetTitleOffset(2);
  h->SetZTitle("f");
  return h;
}
