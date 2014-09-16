//______________________________________________________________________________
//
// TVpBowTieFilter1dTable represents an ideal (i.e, non-scattering) bowtie
// filter which compensates according to a table specifying a material
// thickness as a function of the fan angle.
//
// The local coordinate system: x-axis points in the direction of the patient
// couch, y-axis in the horizontal direction and z-axis in the vertical
// direction, see the following figure:
//
// *    ^z *
// **  |  **
// ****|****
//     ------> y
//______________________________________________________________________________

#include <iostream>
#include <ostream>
#include <fstream>
#include <cmath>
#include "TVpBowTieFilter1dTable.h"
#include "TVpMath.h"

ClassImp(TVpBowTieFilter1dTable)

//______________________________________________________________________________
TVpBowTieFilter1dTable::TVpBowTieFilter1dTable()
{
  // Default constructor.

  fDim = 0;
  fFanAngleCos = fThickness = 0;
  fFormatVersion = 0;
  fFanAngleCosStep = 0.0;
  fMaterialPtr = 0;
}

//______________________________________________________________________________
TVpBowTieFilter1dTable::TVpBowTieFilter1dTable(const Char_t *fileNameBFT)
{
  // Create the bow-tie filter object.

  fDim = 0;
  fFanAngleCos = fThickness = 0;
  fFormatVersion = 0;
  fFanAngleCosStep = 0.0;
  fMaterialPtr = 0;
  ReadBftFile(fileNameBFT);
}

//______________________________________________________________________________
TVpBowTieFilter1dTable::~TVpBowTieFilter1dTable()
{
  // Destructor.  The material is not deleted.

  delete[] fFanAngleCos;
  delete[] fThickness;
}

//______________________________________________________________________________
Double_t TVpBowTieFilter1dTable::GetThickness(Double_t energy,
					      const TVpVector3 &directionL) const
{
  // Get thickness of the water (by default) filter in a given direction and
  // for a given photon energy.  Currently, the result does not depend on the
  // energy.
  //
  // Input parameters:
  // - energy - photon energy in keV
  // - directionL - direction from source
  //
  // Method:
  // The fThickness array gives the track length projected to the y-z plane.
  // It is supposed that the bowtie filter shape does not change with x
  // coordinate.  It can be shown that the real path length is the projected
  // path lenght devided by cos(beta), where beta is the angle between the the
  // directionL and the y-z plane projected directionL.

  Int_t index;
  Double_t projThickness;
  Double_t u1 = directionL.GetY();
  Double_t u2 = directionL.GetZ();
  Double_t cosBeta = std::sqrt(u1*u1 + u2*u2);  // angle between u and y-z plane projected u
  Double_t angleCos = std::abs(u2 / cosBeta);

  // Get projThick from the table
  if (fFormatVersion == 1)
    {
      index = TVpMath::FindIndexByBinarySearch(fDim, fFanAngleCos, angleCos);
      if (index == -1)
	{
	  std::cerr << "Warning: TVpBowTieFilter1dTable::GetThickness: angleCos = "
		    << angleCos << " is below table minimum.\n";
	  projThickness = fThickness[0];
	}
      else if (index == -2)
	projThickness = fThickness[fDim-1]; // This may happen for fFanAngleCos[fDim-1] == 1.
      else
	projThickness = TVpMath::GetLinLinInterpolation
	  (fFanAngleCos[index], fThickness[index],
	   fFanAngleCos[index+1], fThickness[index+1], angleCos);
    }
  else
    {
      if (angleCos < fFanAngleCos[0])
	{
	  std::cerr << "Warning: TVpBowTieFilter1dTable::GetThickness: angleCos = "
		    << angleCos << " is below table minimum.\n";
	  projThickness = fThickness[0];
	}
      else if (angleCos >= fFanAngleCos[fDim-1])
	projThickness = fThickness[fDim-1]; // This may happen for fFanAngleCos[fDim-1] == 1.
      else 
	{
	  index = (Int_t) ((angleCos - fFanAngleCos[0]) / fFanAngleCosStep);
	  projThickness = TVpMath::GetLinLinInterpolation
	    (fFanAngleCos[index], fThickness[index],
	     fFanAngleCos[index+1], fThickness[index+1], angleCos);
	}
    }
  return projThickness / cosBeta;
}

//______________________________________________________________________________
Double_t TVpBowTieFilter1dTable::GetTransmission(Double_t energy,
						 const TVpVector3 &directionL) const
{
  // Calculate the transmission in a given direction and for a given photon
  // energy.
  //
  // Input parameters:
  // - energy - photon energy in keV
  // - directionL - direction from source
  //
  // Method:
  // The transmission is defined as the probability of passing the bow-tie filter
  // without any interaction, i.e. exp(-radiological_path).
  
  return exp(-fMaterialPtr->GetLac(energy) * GetThickness(energy, directionL));
}

//______________________________________________________________________________
void TVpBowTieFilter1dTable::PrintStatus(std::ostream &out) const
{
  // Print the object's status
  
  out << "<TVpBowTieFilter1dTable>\n"
      <<  "<TVpBowTieFilter1dTable>\n";
}

//______________________________________________________________________________
Int_t TVpBowTieFilter1dTable::ReadBftFile(const Char_t *fileNameBFT)
{
  // Read a BFT file (bowtie filter table).
  //
  // Input parameters:
  // - fileName - name of the BFT file
  //
  // Format description:
  // A format version 1.0 file contains fan angle in radians in the first
  // column and corresponding thickness in the second column.  The grid may be
  // non-equidistant.
  //
  // A format version 2.0 file contains cosine of the fan angle in the first
  // column and the corresponding thickness in the second column.  The grid
  // must be equidistant; othervise the index searching routine will fail.
  //
  // File examples:
  // # Bowtie filter 1D table, format 1.0
  // # Name: GE large, 120 kV, body, QX/i and LS16 (Elly Castellano)
  // # Material: material/C.mat
  // # Density: 1.82 g/cm3
  // # Dimension: 26
  // 0	0
  // 0.018	0.011
  // 0.037	0.016
  // ...
  // 0.433	6.606
  //
  // # Bowtie filter 1D table, format 2.0
  // # Name: Siemens Somatom Definition, formfilter P45, fixed
  // # Material: material/Al.mat
  // # Density: 2.699000e+00 g/cm^3
  // # Dimension: 128
  // 8.716504e-01 3.549542e+00
  // 8.726610e-01 3.544906e+00
  // ...
  // 1.000000e+00 4.996471e-02
  
  std::string recName;        // The record name e.g. "# Material"
  std::string recRest;        // The rest of the record
  std::string materialName;   // Material name
  Double_t materialDensity;   // Material density
    
  // Read the BFT file
  std::ifstream in(fileNameBFT);
  if (!in)
    {
      std::cerr << "Error: TVpBowTieFilter1dTable::ReadBftFile: Cannot open file: "
		<< fileNameBFT << '\n';
      return 1;
    }

  // Read and process the file header
  std::getline(in, recName);
  if (recName == std::string("# Bowtie filter 1D table, format 1.0"))
    fFormatVersion = 1;
  else if (recName == std::string("# Bowtie filter 1D table, format 2.0"))
    fFormatVersion = 2;
  else
    {
      std::cerr << "Error: TVpBowTieFilter1dTable::ReadBftFile: Record "
		<< "\"# Bowtie filter 1D table, format N.N\" is missing.\n";
      return 2;
    }
  
  // # Name:
  std::getline(in, recName, ':');
  if (recName != std::string("# Name"))
    {
      std::cerr << "Error: TVpBowTieFilter1dTable::ReadBftFile: Record "
		<< "\"# Name:\" is missing.\n";
      return 3;
    }
  std::getline(in, fName);
  
  // # Material:
  std::getline(in, recName, ':');
  if (recName != std::string("# Material"))
    {
      std::cerr << "Error: TVpBowTieFilter1dTable::ReadBftFile: Record "
		<< "\"# Material:\" is missing.\n";
      return 4;
    }
  in >> materialName;
  std::getline(in, recRest);

  // # Density:
  std::getline(in, recName, ':');
  if (recName != std::string("# Density"))
    {
      std::cerr << "Error: TVpBowTieFilter1dTable::ReadBftFile: Record "
		<< "\"# Density:\" is missing.\n";
      return 5;
    }
  in >> materialDensity;
  std::getline(in, recRest);

  // # Dimension:
  std::getline(in, recName, ':');
  if (recName != std::string("# Dimension"))
    {
      std::cerr << "Error: TVpBowTieFilter1dTable::ReadBftFile: Record "
		<< "\"# Dimension:\" is missing.\n";
      return 5;
    }
  in >> fDim;
  std::getline(in, recRest);
  if (fDim < 2)
    {
      std::cerr << "Error: TVpBowTieFilter1dTable::ReadBftFile: "
		<< "Dimension must be at least 2.\n";
      return 5;
    }

  // Re-allocate arrays
  delete[] fFanAngleCos;
  delete[] fThickness;
  fFanAngleCos = new Double_t[fDim];
  fThickness = new Double_t[fDim];
  
  // Read Bft table
  Double_t angle;
  if (fFormatVersion == 1)
    {
      for (Int_t i = fDim -1; i >= 0; i--)
	{
	  in >> angle >> fThickness[i];
	  fFanAngleCos[i] = cos(angle);
	}
    }
  else  // fFormatVersion == 2
    {
      for (Int_t i = 0; i < fDim; i++)
	in >> fFanAngleCos[i] >> fThickness[i];
      fFanAngleCosStep = fFanAngleCos[1] - fFanAngleCos[0];
    }
  
  // Read the material
  fMaterialPtr = new TVpMaterial(materialName.c_str());
  if (fMaterialPtr != 0)
    {
      fMaterialPtr->SetDensity(materialDensity);
      fMaterialPtr->Initialize();
    }
  else
    std::cerr << "Error: TVpBowTieFilter1dTable::ReadBftFile: The material is not initialized. "
	      << "The code will crash if the BF is used.\n";
  
  in.close();
  return 0;
}


//______________________________________________________________________________
TGraph *TVpBowTieFilter1dTable::GetBftThickness() const
{
  // Return the bowtie filter thickness (in cm) as a function of the cosine of
  // fan angle

  TGraph *g = new TGraph(fDim, fFanAngleCos, fThickness);
  g->SetTitle(fName.c_str());
  return g;
}
