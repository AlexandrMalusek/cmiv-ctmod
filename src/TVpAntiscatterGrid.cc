//______________________________________________________________________________
//
// TVpAntiscatterGrid calculates the transmission formula of a linear
// (unfocussed) antiscatter grid.  The grid consists of upper and lower
// covers, and absorbing strips separated by an interspace material.
//
//    interspace 
//       v
//    ------------------  upper cover
//    ------------------  
//    ||  ||  ||  ||  ||  ^
//    ||  ||  ||  ||  ||  fStripHeight
//    ||  ||  ||  ||  ||  v
//    ------------------
//    ------------------  lower cover
//    ^
//    absorbing strip
//
// In the current implementation, strip normal is parallel to the x-axis
// (head-toe body axis), see TVpTomographSetup.
//
// Typical grid ratio varies between 8:1 and 16:1 (Dendy P.P., Heaton B.:
// Physics for diagnostic radiology).  The former is used for x-ray tube
// potentials less than 85 kV.  For higher potentials the choice is between
// 10:1 and 12:1 grids.

#include <cmath>
#include "TVpAntiscatterGrid.h"

ClassImp(TVpAntiscatterGrid)

//______________________________________________________________________________
TVpAntiscatterGrid::TVpAntiscatterGrid()
{
  // Default constructor

  fStripHeight = fStripWidth = fInterspaceWidth = fUpperCoverWidth = 
    fLowerCoverWidth = 0.0;
  fMaterialCover = fMaterialStrip = fMaterialInterspace = 0;
}

//______________________________________________________________________________
TVpAntiscatterGrid::TVpAntiscatterGrid(Double_t gridHeight, Double_t stripWidth,
				       Double_t interspaceWidth,
				       Double_t upperCoverWidth, Double_t lowerCoverWidth,
				       TVpMaterial *materialCover,
				       TVpMaterial *materialStrip,
				       TVpMaterial *materialInterspace)
{
  // Full initialization constructor.
  //
  // Input:
  // - gridHeight - height of the strips in cm
  // - stripWidth - width of absorbing strips in cm
  // - interspaceWidth - width of the interspace material in cm
  // - upperCoverWidth - width of the upper (closer to the phantom) cover in cm
  // - lowerCoverWidth - width of the lower (closer to the detector) cover in cm
  // - materialCover - material of the cover
  // - materialStrip material of absorbing strips
  // - materialInterspace - interspace material

  fStripHeight = gridHeight;
  fStripWidth = stripWidth;
  fInterspaceWidth = interspaceWidth;
  fUpperCoverWidth = upperCoverWidth;
  fLowerCoverWidth = lowerCoverWidth;
  fMaterialCover = materialCover;
  fMaterialStrip = materialStrip;
  fMaterialInterspace = materialInterspace;
}

//______________________________________________________________________________
TVpAntiscatterGrid::~TVpAntiscatterGrid()
{
  // Destructor. Materials are not deleted.
  
}

//______________________________________________________________________________
void TVpAntiscatterGrid::PrintStatus(std::ostream &out) const
{
  // Print the object status
  //
  // Input:
  // - out - output stream
  //
  // Example output:
  // <TVpAntiscatterGrid>
  // Strip height: 0.16 cm
  // Strip width: 0.004 cm
  // Interspace width: 0.016 cm
  // Upper cover width: 0.05 cm
  // Lower Cover Width: 0.05 cm
  // Cover material: Carbon fibre cover
  // Strip material: Pb, Lead
  // Interspace material: Cotton fibre interspaces
  // Number of strips per cm: 50 1/cm
  // Grid ratio: 10 : 1
  // </TVpAntiscatterGrid>

  out << "<TVpAntiscatterGrid>\n"
      << "Strip height: " << fStripHeight << " cm\n"
      << "Strip width: " << fStripWidth << " cm\n"
      << "Interspace width: " << fInterspaceWidth << " cm\n"
      << "Upper cover width: " << fUpperCoverWidth << " cm\n"
      << "Lower Cover Width: " << fLowerCoverWidth << " cm\n"
      << "Cover material: " << fMaterialCover->GetName() << '\n'
      << "Strip material: " << fMaterialStrip->GetName() << '\n'
      << "Interspace material: " << fMaterialInterspace->GetName() << '\n'
      << "Number of strips per cm: " << 1.0 / (fStripWidth + fInterspaceWidth) 
      << " 1/cm\n"
      << "Grid ratio: " << fStripHeight / fInterspaceWidth << " : 1\n" 
      << "</TVpAntiscatterGrid>\n";
}

//______________________________________________________________________________
Double_t TVpAntiscatterGrid::Transmission(Double_t energy,
					  const TVpVector3 &directionL) const
{
  // Transmission formula for the unfocused antiscatter grid. The same
  // notation as in the G. J. Day's article is used. Grid strips are
  // perpendicular to the x-axis. directionL is expressed in local
  // coordinates.
  //
  // Input:
  // - energy - energy of the impinging photon in keV
  // - directionL - direction of the photon in the coordinate system of the detector

  Double_t B;
  
  Double_t h = fStripHeight;
  Double_t d = fStripWidth;
  Double_t D = fInterspaceWidth;
  Double_t u = std::abs(directionL.fR[0]);
  if (u < 1.0e-10)   // HOT FIX!!! to avoid NaNs
    u = 1.0e-10;
  Double_t w = std::abs(directionL.fR[2]);
  if (w < 1.0e-10)   // HOT FIX!!! to avoid NaNs
    w = 1.0e-10;
  Double_t C = fStripWidth + fInterspaceWidth;
  Double_t muI = fMaterialInterspace->GetLac(energy);
  Double_t muS = fMaterialStrip->GetLac(energy);
  Double_t muC = fMaterialCover->GetLac(energy);
  Double_t dmu = muS - muI;
  Double_t P = h * u / w;
  Int_t n = (P/C > 2000000000) ? 2000000 : (int) (P / C);   // HOT FIX!!!
  Double_t q = P - n*C;
  Double_t a = u / dmu;
  Double_t A = exp(-muI*h/w - n*d/a);
  if (q < d)
    B = D - q + 2*a + (d - q - 2*a)*exp(-q/a);
  else if (q < D)
    B = D - q + 2*a + (q - d - 2*a)*exp(-d/a);
  else
    B = (q - D + 2*a) * exp(-(q - D)/a) + (q - d - 2*a) * exp(-d/a);
  Double_t T = A * B / C;
  
  T *= exp(-(fUpperCoverWidth + fLowerCoverWidth) * muC / w);
  return T;
}

//______________________________________________________________________________
TH2D *TVpAntiscatterGrid::GetHistTransmissionThetaPhi(Double_t energy, Int_t numTheta,
						      Int_t numPhi)
{
  // Return a 2D histogram of the grid transmission function as a function of
  // the polar and azimuthal angle for a photon with a given energy.
  //
  // Input:
  // - energy - energy of the impinging photon in keV
  // - numTheta - number of bins of the polar angle (Pi/2,..,Pi)
  // - numPhi - number of bins ot the azimuthal angle (0, ..,2Pi)
  //
  // Method:
  // Bin (i,j), where 0<=i<numTheta, 0<=j<numPhi contains values for
  // theta = M_PI/2 + i * dTheta, and phi = j * dPhi, where
  // dTheta = M_PI / 2.0 / (numTheta - 1) and
  // dPhi = 2*M_PI / (numPhi - 1).

  TVpVector3 dirL;
  Double_t dTheta = M_PI / 2.0 / (numTheta - 1);
  Double_t dPhi = 2*M_PI / (numPhi - 1);
  TH2D *h2 = new TH2D("TE", "TE", numTheta, M_PI/2, M_PI, numPhi, 0, 2*M_PI);
  for (Int_t it = 0; it < numTheta; it++)
    {
      Double_t theta = M_PI/2 + it * dTheta;
      for (Int_t ip = 0; ip < numPhi; ip++)
	{
	  Double_t phi = ip * dPhi;
	  dirL.SetPolar(1.0, theta, phi);
	  h2->SetBinContent(it+1, ip+1, Transmission(energy, dirL));
	}
    }
  h2->SetXTitle("#theta [rad]");
  h2->SetYTitle("#phi [rad]");
  h2->SetZTitle("T(#theta,#phi)");
  return h2;
}
